# by Deus Thindwa, Dan Weinberger, Ginny Pitzer w/ modifications by Katie Kortright
# age-structured mathematical model for pneumococcal transmission
# 09/01/2026

#====================================================================
pacman::p_load(char = c("tidyverse","remotes", "ggplot2", "dplyr","rio","tidyr","plyr","lubridate","reshape2","curl","patchwork", "posterior",
                        "deSolve","adaptivetau","data.table","scales","readr","MASS","rootSolve", "labelVector","PropCIs", "bayesplot", "readr",
                        "binom","coda","Rcpp","gmm","RcppArmadillo","devtools","lattice", "RColorBrewer", "lhs","png", "cmdstanr", "tibble", "GGally",
                        "readxl","viridis","zoo","lattice","latex2exp","ape","cowplot","gridExtra","grid","ggpubr","here","deSolve", "loo"))


#import for posterior analysis
carriage_fit <- rio::import(here::here("results", "carriage_fit_19A.rds"))


#posterior diagnostics
key_pars <- c("log_rho_age", "rho_age", "eps_V", "eps_A", "eps_N", "rho_V", "rho_A", "rho_N", "lambda_V_eq", "lambda_A_eq", "lambda_N_eq")
n_age <- 4
param_names <- c(paste0("rho_age[", 1:n_age, "]"), paste0("eps_V"), paste0("eps_A"), paste0("eps_N"))

summ <-
  carriage_fit$summary(variables = param_names,
                       posterior::default_summary_measures(),
                       extra_quantiles = ~ posterior::quantile2(., probs = c(0.025, 0.975)))
print(summ, n = nrow(summ))


#trace plots
trace_pars <- c("log_rho_age[1]","log_rho_age[2]","log_rho_age[3]","log_rho_age[4]", "eps_V","eps_A","eps_N")
trace_plot <- bayesplot::mcmc_trace(carriage_fit$draws(trace_pars))

print(trace_plot)
ggsave(here::here("results", "trace_plot.png"), trace_plot, width = 9, height = 6, dpi = 330, bg = "white")


#PSIS-LOO and Pareto-k diagnostics
log_lik_draws <- carriage_fit$draws("log_lik", format = "draws_matrix")
log_lik_mat <- as.matrix(log_lik_draws)
loo_out <- loo::loo(log_lik_mat)
print(loo_out)

print(table(cut(loo_out$diagnostics$pareto_k,
                breaks = c(-Inf, 0.5, 0.7, 1, Inf),
                labels = c("k<=0.5 (good)", "0.5<k<=0.7 (ok)", "0.7<k<=1 (bad)", "k>1 (very bad)"))))


#pairs plot (supplementary)
draws_array <- carriage_fit$draws()
p_pairs <-
  bayesplot::mcmc_pairs(draws_array, pars = c("rho_age[1]","rho_age[2]","rho_age[3]","rho_age[4]", "eps_V","eps_A","eps_N"),
                        off_diag_args = list(size = 0.3, alpha = 0.4))

p_pairs
ggsave(here::here('output', "carriage_model", "pairs_plot.png"), p_pairs, width = 8, height = 8, dpi = 150)

######Stopped here
#joint posterior density pairs plot
#use unconstrained log-scale draws for the pairs plot
#space HMC actually explored, where correlations are most informative
post_pairs <-
  as_draws_df(carriage_fit$draws(param_names)) %>%
  dplyr::select(all_of(param_names))

lower_density <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_density_2d_filled(alpha = 0.85, bins = 8, show.legend = FALSE) +
    geom_point(alpha = 0.05, size = 0.4) +
    scale_fill_viridis_d(option = "mako", direction = -1) +
    theme_minimal(base_size = 8)
}

diag_density <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_density(fill = "steelblue", alpha = 0.5, colour = "steelblue4") +
    theme_minimal(base_size = 8)
}

upper_cor <- function(data, mapping, ...) {
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  r <- cor(x, y)
  ggplot(data, mapping) +
    geom_density_2d(colour = "steelblue4", linewidth = 0.4) +
    annotate("text", x = mean(range(x)), y = mean(range(y)), label = sprintf("r = %.2f", r), size = 3, colour = "firebrick", fontface = "bold") +
    theme_minimal(base_size = 8)
}

pairs_density_plot <-
  ggpairs(post_pairs,
          lower = list(continuous = lower_density),
          diag  = list(continuous = diag_density),
          upper = list(continuous = upper_cor),
          axisLabels = "show",
          title = "Joint posterior density of the 7 free log-scale parameters") +
  theme(strip.text = element_text(size = 7))

print(pairs_density_plot)
ggsave(here::here("output", "carriage_model", "pairs_density_plot.png"), pairs_density_plot, width = 14, height = 14, dpi = 330)


#load posterior + reconstruct prior hyperparameters
age_groups <- c("<1y", "1-4y", "5-17y", "18+y")
n_age <- length(age_groups)
sero_grps <- c("V","A","N")
rho_age_logmean <- log(c(0.05, 0.03, 0.015, 0.005)) #MUST match what was passed to stan_data
rho_age_logsd <- 0.5
eps_alpha <- 2
eps_beta <- 2

#observed data (re-declared here so this file is stand-alone)
y_obs <- matrix(c(
  222,  203,  7,  70,   # <1y    (N = 502)
  524,  277,  11, 159,   # 1-4y   (N = 971)
  361,   47,   2,  47,   # 5-17y  (N = 457)
  1789,  79,  3,  70    # 18+y   (N = 1938)
), nrow = n_age, byrow = TRUE)
rownames(y_obs) <- age_groups
colnames(y_obs) <- c("S","V","A","N")
storage.mode(y_obs) <- "integer"
N_total       <- as.integer(rowSums(y_obs))
observed_prev <- sweep(y_obs, 1, N_total, "/")


#prior vs posterior densities
#log-rho_age parameters
post_log_rho_age <- carriage_fit$draws("log_rho_age", format = "draws_matrix")  #iter x 4
prior_panels_age <- lapply(seq_len(n_age), function(a) {
  post <- as.numeric(post_log_rho_age[, a])
  xlim <- range(c(post, rho_age_logmean[a] + c(-1, 1) * 3 * rho_age_logsd))
  xs   <- seq(xlim[1], xlim[2], length.out = 400)

  prior <- dnorm(xs, mean = rho_age_logmean[a], sd = rho_age_logsd)
  ggplot() +
    geom_density(aes(x = post, y = after_stat(density), fill = "Posterior", colour = "Posterior"), alpha = 0.45) +
    geom_line(aes(x = xs, y = prior, colour = "Prior"), linetype = "dashed", linewidth = 0.7) +
    scale_fill_manual(values  = c(Posterior = "steelblue")) +
    scale_colour_manual(values = c(Posterior = "steelblue", Prior = "firebrick")) +
    labs(x = sprintf("log_rho_age[%d]  (%s)", a, age_groups[a]), y = "Density", fill = NULL, colour = NULL) +
    theme(legend.position = "none")
})

#eps_X parameters (on natural scale, prior is Beta(eps_alpha, eps_beta))
eps_draws <- carriage_fit$draws(c("eps_V","eps_A","eps_N"), format = "draws_matrix")
prior_panels_eps <- lapply(seq_len(3), function(s) {
  post <- as.numeric(eps_draws[, s])
  xs    <- seq(0.001, 0.999, length.out = 400)
  prior <- dbeta(xs, eps_alpha, eps_beta)
  ggplot() +
    geom_density(aes(x = post, y = after_stat(density), fill = "Posterior", colour = "Posterior"), alpha = 0.45) +
    geom_line(aes(x = xs, y = prior, colour = "Prior"), linetype = "dashed", linewidth = 0.7) +
    scale_fill_manual(values  = c(Posterior = "steelblue")) +
    scale_colour_manual(values = c(Posterior = "steelblue", Prior     = "firebrick")) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x = sprintf("eps_%s", sero_grps[s]), y = "Density", fill = NULL, colour = NULL) +
    theme(legend.position = "none")
})

prior_post_plot <- (
  prior_panels_age[[1]] | prior_panels_age[[2]] | prior_panels_age[[3]] | prior_panels_age[[4]]) /
  (prior_panels_eps[[1]] | prior_panels_eps[[2]] | prior_panels_eps[[3]] | patchwork::plot_spacer()) +
  patchwork::plot_annotation(title = "Prior (red dashed) vs posterior (blue) density for the 7 free parameters") +
  theme(legend.position = "none")

print(prior_post_plot)
ggsave(here::here("output", "carriage_model", "prior_vs_posterior.png"), prior_post_plot, width = 12, height = 6.5, dpi = 300, bg = "white")


#posterior prevalence and predictive checks (Main Figure)
#credible intervals on equilibrium proportions p_obs[a, k] vs observed
wald_ci <- function(x, n, z = 1.96) {
  p <- x / n
  se <- sqrt(p * (1 - p) / n)
  list(L = pmax(0, p - z * se), U = pmin(1, p + z * se))
}

df_obs <-
  tibble( age_grp = c("<1y", "1-4y", "5-17y", "18+y"),
          S = c(222, 524, 361, 1789),
          V = c(203, 277, 47, 79),
          A = c(7, 11, 2, 3),
          N = c(70, 159, 47, 70)) %>%
  dplyr::mutate(n = S + V + A + N,
                S_prop = S / n,
                V_prop = V / n,
                A_prop = A / n,
                N_prop = N / n,
                S_L = wald_ci(S, n)$L, S_U = wald_ci(S, n)$U,
                V_L = wald_ci(V, n)$L, V_U = wald_ci(V, n)$U,
                A_L = wald_ci(A, n)$L, A_U = wald_ci(A, n)$U,
                N_L = wald_ci(N, n)$L, N_U = wald_ci(N, n)$U) %>%
  dplyr::mutate(n = S + V + A + N) %>%
  tidyr::pivot_longer(cols = c(S, V, A, N), names_to = "state", values_to = "count") %>%
  dplyr::mutate(prop = count / n,
                se = sqrt(prop * (1 - prop) / n),
                L = pmax(0, prop - 1.96 * se),
                U = pmin(1, prop + 1.96 * se)) %>%
  dplyr::select(age_grp, state, prop, L, U)

p_obs_draws <- carriage_fit$draws("p_obs", format = "draws_df") %>% as_tibble()
pp_summary <- tidyr::expand_grid(age = 1:n_age, k = 1:4) %>%
  dplyr::mutate(
    var      = sprintf("p_obs[%d,%d]", age, k),
    mean     = vapply(var, function(v) mean(p_obs_draws[[v]]),     numeric(1)),
    lower    = vapply(var, function(v) quantile(p_obs_draws[[v]], 0.025), numeric(1)),
    upper    = vapply(var, function(v) quantile(p_obs_draws[[v]], 0.975), numeric(1)),
    age_grp  = factor(age_groups[age], levels = c("<1y", "1-4y", "5-17y", "18+y")),
    state    = factor(c("S","V","A","N")[k], levels = c("S","V","A","N")),
    observed = as.numeric(t(observed_prev))[(age - 1) * 4 + k])

plot_df1 <-
  bind_rows(pp_summary %>% dplyr::left_join(df_obs) %>% dplyr::select(age_grp, state, est = mean,  lo = lower, hi = upper) %>% mutate(source = "Model"),
            pp_summary %>% dplyr::left_join(df_obs) %>% dplyr::select(age_grp, state, est = prop,  lo = L, hi = U) %>% mutate(source = "Observed")) %>%
  dplyr::filter(state %in% c('S', "A", "V", "N")) %>%
  dplyr::mutate(state = factor(state,   levels = c('S', "A", "V", "N")),
                age_grp = factor(age_grp, levels = c("<1y", "1-4y", "5-17y", "18+y")),
                source  = factor(source,  levels = c("Model", "Observed")))
dodge <- position_dodge(width = 0.4)   #0ne object shared by both geoms

ppc_prop_plot <-
  ggplot(plot_df1, aes(x = age_grp, y = est, colour = source, group = source)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.8, position = dodge) +
  geom_point(shape = 4, size = 2.8, stroke = 1.2, position = dodge) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_colour_manual(values = c(Model = "steelblue", Observed = "darkred"), name = NULL) +
  facet_wrap(factor(state, levels=c('S', 'A', 'V', 'N')) ~ ., scales = 'free', ncol = 4) +
  labs(x = NULL, y = "Prevalence") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(size = 0), axis.text.y = element_text(size = 14))

ppc_prop_plot

#PPC of y_rep counts vs y_obs
y_rep_draws <- carriage_fit$draws("y_rep", format = "draws_df") %>% as_tibble()
y_rep_summary <- tidyr::expand_grid(age = 1:n_age, k = 1:4) %>%
  dplyr::mutate(
    var      = sprintf("y_rep[%d,%d]", age, k),
    mean     = vapply(var, function(v) mean(y_rep_draws[[v]]),     numeric(1)),
    lower    = vapply(var, function(v) quantile(y_rep_draws[[v]], 0.025), numeric(1)),
    upper    = vapply(var, function(v) quantile(y_rep_draws[[v]], 0.975), numeric(1)),
    age_grp  = factor(age_groups[age], levels = age_groups),
    state    = factor(c("S","V","A","N")[k], levels = c("S","V","A","N")),
    observed = as.numeric(t(y_obs))[(age - 1) * 4 + k])

plot_df2 <-
  bind_rows(y_rep_summary %>% dplyr::select(age_grp, state, est = mean,  lo = lower, hi = upper) %>% mutate(source = "Model"),
            y_rep_summary %>% dplyr::select(age_grp, state, est = observed,  lo = observed, hi = observed) %>% mutate(source = "Observed")) %>%
  dplyr::filter(state %in% c('S', "A", "V", "N")) %>%
  dplyr::mutate(state = factor(state,   levels = c('S', "A", "V", "N")),
                age_grp = factor(age_grp, levels = c("<1y", "1-4y", "5-17y", "18+y")),
                source  = factor(source,  levels = c("Model", "Observed")))
dodge <- position_dodge(width = 0.4)   #0ne object shared by both geoms

ppc_counts_plot <-
  ggplot(plot_df2, aes(x = age_grp, y = est, colour = source, group = source)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.8, position = dodge) +
  geom_point(shape = 4, size = 2.8, stroke = 1.2, position = dodge) +
  scale_colour_manual(values = c(Model = "steelblue", Observed = "darkred"), name = NULL) +
  facet_wrap(factor(state, levels=c('S', 'A', 'V', 'N')) ~ ., scales = 'free', ncol = 4) +
  labs(x = NULL, y = "Counts") +
  theme_bw(base_size = 14, base_family = "Lato") +
  theme(panel.grid.minor  = element_blank(), strip.background  = element_rect(fill = "grey92")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 1))  +
  theme(strip.text.x = element_text(size = 18), strip.background = element_rect(fill = "gray90")) +
  theme(legend.position = 'bottom') +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))

ppc_counts_plot

ppc_combo <- ppc_prop_plot / ppc_counts_plot
print(ppc_combo)
ggsave(here::here("output", "carriage_model", "ppc_prevalence_and_counts.png"), ppc_combo, width = 16, height = 12, dpi = 300, bg = "white")

#posterior of the 4 age + 3 competition parameters
#age baselines
rho_age_draws <- carriage_fit$draws("rho_age", format = "draws_matrix")
age_summary <- tibble::tibble(
  age      = age_groups,
  prior    = exp(rho_age_logmean),
  mean     = colMeans(rho_age_draws),
  lower    = apply(rho_age_draws, 2, quantile, 0.025),
  upper    = apply(rho_age_draws, 2, quantile, 0.975))
print(age_summary)

#competition parameters
eps_draws <- carriage_fit$draws(c("eps_V","eps_A","eps_N"), format = "draws_matrix")
eps_summary <- tibble::tibble(
  serotype = c("V","A","N"),
  prior    = eps_alpha / (eps_alpha + eps_beta),
  mean     = colMeans(eps_draws),
  lower    = apply(eps_draws, 2, quantile, 0.025),
  upper    = apply(eps_draws, 2, quantile, 0.975))
print(eps_summary)


#posterior of cell-level susceptibilities rho_{X}[i]
#rho_V = rho_A = rho_N = rho_age, so the three lines overlap;
rho_long <- bind_rows(carriage_fit$draws("rho_V", format = "draws_df") %>% as_tibble() %>%
                        tidyr::pivot_longer(starts_with("rho_V"), names_to = "var", values_to = "rho") %>%
                        dplyr::mutate(serotype = "V", age = as.integer(sub("rho_V\\[(\\d+)\\]", "\\1", var))),

                      carriage_fit$draws("rho_A", format = "draws_df") %>% as_tibble() %>%
                        tidyr::pivot_longer(starts_with("rho_A"), names_to = "var", values_to = "rho") %>%
                        dplyr::mutate(serotype = "A", age = as.integer(sub("rho_A\\[(\\d+)\\]", "\\1", var))),

                      carriage_fit$draws("rho_N", format = "draws_df") %>% as_tibble() %>%
                        tidyr::pivot_longer(starts_with("rho_N"), names_to = "var", values_to = "rho") %>%
                        dplyr::mutate(serotype = "N", age = as.integer(sub("rho_N\\[(\\d+)\\]", "\\1", var)))) %>%
  dplyr::mutate(age_grp = factor(age_groups[age], levels = age_groups))

rho_summary <-
  rho_long %>%
  dplyr::group_by(age_grp, serotype) %>%
  dplyr::summarise(mean   = mean(rho), lower  = quantile(rho, 0.025), upper  = quantile(rho, 0.975), .groups = "drop")
print(rho_summary)

rho_plot <-
  ggplot(rho_summary, aes(x = age_grp, y = mean, colour = serotype, group = serotype)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3), linetype = "dashed") +
  labs(x = "Age group", y = expression(rho[a*","*s]), title = "Posterior mean (95% CrI) of cell-level susceptibility")

print(rho_plot)
ggsave(here::here("output", "carriage_model", "susceptibility_by_age_sero.png"), rho_plot, width = 7, height = 4, dpi = 300, bg = "white")

#eps posterior plot
eps_plot <-
  ggplot(eps_summary, aes(x = serotype, y = mean)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), colour = "steelblue") +
  geom_hline(yintercept = eps_summary$prior[1], linetype = "dotted") +
  ylim(0, 1) +
  labs(x = "Serotype", y = expression(epsilon), title = "Posterior mean (95% CrI) of competition parameter")

print(eps_plot)
ggsave(here::here("output", "carriage_model", "relative_risk_by_sero.png"), eps_plot, width = 5, height = 4, dpi = 300, bg = "white")


#compartment composition per age group (S, V, A, N, VA, NV, NA) at equilibrium
comp_names <- c("S","V","A","N","VA","NV","NA")

comp_eq_draws <- carriage_fit$draws("comp_eq", format = "draws_df") %>% as_tibble()
comp_summary <- tidyr::expand_grid(age = 1:n_age, k = 1:7) %>%
  dplyr::mutate(
    var       = sprintf("comp_eq[%d,%d]", age, k),
    mean      = vapply(var, function(v) mean(comp_eq_draws[[v]]),     numeric(1)),
    lower     = vapply(var, function(v) quantile(comp_eq_draws[[v]], 0.025), numeric(1)),
    upper     = vapply(var, function(v) quantile(comp_eq_draws[[v]], 0.975), numeric(1)),
    age_grp   = factor(age_groups[age], levels = age_groups),
    compartment = factor(comp_names[k], levels = comp_names))
print(comp_summary %>% dplyr::select(age_grp, compartment, mean, lower, upper), n = nrow(comp_summary))


#stacked-bar layout of posterior means per age group
comp_bar_plot <-
  ggplot(comp_summary, aes(x = age_grp, y = mean, fill = compartment)) +
  geom_col(position = "stack", width = 0.6, colour = "white") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Age group", y = "Equilibrium proportion", fill = "Compartment", title = "Equilibrium compartment composition (posterior mean)", subtitle = "S, V, A, N (single), VA / NV / NA (dual carriage)")

comp_bar_plot


#pointrange layout to show 95% CrI for every compartment
comp_pr_plot <-
  ggplot(comp_summary, aes(x = compartment, y = mean, colour = compartment)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.35) +
  facet_wrap(~ age_grp, nrow = 1, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  scale_colour_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = "Equilibrium proportion", title    = "Posterior mean (95% CrI) of each compartment, by age group", subtitle = "Bars: 95% credible interval")

comp_pr_plot

comp_combo <- comp_bar_plot / comp_pr_plot
print(comp_combo)
ggsave(here::here("output", "carriage_model", "compartments_by_age.png"), comp_combo, width = 16, height = 12, dpi = 300, bg = "white")


#equilibrium force of infection (per year) lambda_V_eq, lambda_A_eq, lambda_N_eq
lambda_draws <-
  carriage_fit$draws(c("lambda_V_eq","lambda_A_eq","lambda_N_eq"), format = "draws_df") %>% as_tibble()

foi_summary <- tidyr::expand_grid(serotype = sero_grps, age = 1:n_age) %>%
  dplyr::mutate(
    var      = sprintf("lambda_%s_eq[%d]", serotype, age),
    mean     = vapply(var, function(v) mean(lambda_draws[[v]]),     numeric(1)),
    lower    = vapply(var, function(v) quantile(lambda_draws[[v]], 0.025), numeric(1)),
    upper    = vapply(var, function(v) quantile(lambda_draws[[v]], 0.975), numeric(1)),
    age_grp  = factor(age_groups[age], levels = age_groups),
    serotype = factor(serotype, levels = sero_grps))
print(foi_summary %>% dplyr::select(age_grp, serotype, mean, lower, upper), n = nrow(foi_summary))

foi_plot <- ggplot(foi_summary, aes(x = age_grp, y = mean, colour = serotype, group = serotype)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.3)) +
  geom_line(position = position_dodge(width = 0.3), linetype = "dashed") +
  scale_y_continuous() +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Age group", y = expression(lambda[s]^{eq}~"(per year)"), colour = "Serotype", title = "Equilibrium force of infection per serotype and age group", subtitle = "Posterior mean and 95% credible interval")

print(foi_plot)
ggsave(here::here("output", "carriage_model", "foi_equilibrium.png"), foi_plot, width = 7.5, height = 4.5, dpi = 300, bg = "white")
