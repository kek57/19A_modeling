# by Deus Thindwa, Dan Weinberger, Ginny Pitzer w/ modifications by Katie Kortright
# age-structured mathematical model for pneumococcal transmission
09/01/2026
#====================================================================

# a deterministic age-structured SIS model of pneumococcal carriage (V, A, N serotype groups across 4 age groups).
#
# free parameters with informative priors:
#
# log_rho_age[a]
# 4 age-specific log-baselines
# susceptibility is shared across serotypes rho_V[i] = rho_A[i] = rho_N[i] = rho_age[i])
# ~ Normal(rho_age_logmean[a], rho_age_logsd)
#
# eps_V, eps_A, eps_N, 3 serotype-specific competition
# relative risk of co-colonisation) in (0,1)
# ~ Beta(eps_alpha, eps_beta)
#
# prior centres for the age-baseline drawn from pneumococcal literature
# e.g. Ojal 2017; Choi 2011/2012, Bottomley 2017
# <1y    ~ 0.05
# 1-4y   ~ 0.03
# 5-17y  ~ 0.015
# 18+y   ~ 0.005
#
# Beta(2, 2) provides a unimodal weakly-informative prior on each
# competition parameter, centred at 0.5 with sd ~0.224.
#
# time and rates in YEARS, explicit aging between groups, independent dual-carriage clearance

# ==============================================================================

#load the require packages
if (!require(pacman)){
  install.packages("pacman")
  install.packages("RcmdrPlugin.KMggplot2")
}

pacman::p_load(
  char = c("tidyverse","remotes", "ggplot2", "dplyr","rio","tidyr","plyr","lubridate","reshape2","curl","patchwork", "posterior",
           "deSolve","adaptivetau","data.table","scales","readr","MASS","rootSolve", "labelVector","PropCIs", "bayesplot", "readr",
           "binom","coda","Rcpp","gmm","RcppArmadillo","devtools","lattice", "RColorBrewer", "lhs","png", "cmdstanr", "tibble", "GGally",
           "readxl","viridis","zoo","lattice","latex2exp","ape","cowplot","gridExtra","grid","ggpubr","here","deSolve", "loo")
)

set.seed(1234)
theme_set(theme_minimal(base_size = 12))
bayesplot::color_scheme_set("brewer-RdYlBu")

# ==============================================================================

# #total contact matrix (run this code only if you want to get the contacts matrix C_day (https://www.epistorm.org/data/epistorm-mix))
# age_groups_5yr <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")
#
# #weighted mean daily in-person contacts (Epistorm-Mix, 2024)
# M <- matrix(c(
#   1.86068,0.70922,0.27294,0.19284,0.38647,0.52152,0.74264,0.65879,0.50745,0.27078,0.13497,0.11624,0.20563,0.22665,0.07370,0.03396,
#   0.21983,4.18650,0.88506,0.37158,0.19047,0.35853,0.66981,0.76950,0.98274,0.32328,0.37390,0.25236,0.15199,0.14732,0.11631,0.08913,
#   0.16316,0.24448,3.97903,0.65597,0.29650,0.11857,0.39405,0.54558,0.78884,0.61383,0.39497,0.11900,0.14986,0.08006,0.07873,0.06910,
#   0.28803,0.66286,0.52931,3.81800,0.63310,0.36918,0.47484,0.51248,0.78651,0.77182,0.60021,0.29266,0.15755,0.13396,0.11359,0.16432,
#   0.24193,0.21859,0.51182,0.80917,1.22247,0.63708,0.48638,0.65223,0.69792,0.49391,0.50489,0.48721,0.35503,0.11931,0.11213,0.19379,
#   0.39399,0.39625,0.16909,0.32176,0.45123,0.95707,0.81098,0.34533,0.35295,0.24431,0.49553,0.34436,0.34734,0.13904,0.08809,0.17117,
#   0.41114,0.40073,0.61230,0.47785,0.54435,0.66417,1.14187,0.95482,0.60610,0.38127,0.50529,0.46032,0.48470,0.27719,0.14004,0.15737,
#   0.60228,0.57665,1.03357,0.39680,0.35656,0.45246,0.74558,0.98167,0.64406,0.59712,0.39807,0.30514,0.35053,0.27646,0.25654,0.17135,
#   0.28746,0.38422,0.48765,0.77529,0.46967,0.51555,0.51592,0.60443,0.93232,0.77716,0.35524,0.46906,0.39576,0.22935,0.29411,0.21664,
#   0.09829,0.36638,0.57683,0.97212,0.37505,0.35279,0.48063,0.60579,0.76518,0.89284,0.77723,0.52809,0.40137,0.21828,0.13325,0.26808,
#   0.19709,0.26375,0.19000,0.39801,0.57346,0.61348,0.86635,0.61171,0.49417,0.62418,0.86114,0.70865,0.51590,0.17011,0.16663,0.18595,
#   0.32640,0.25551,0.48136,0.58423,0.71343,0.70075,0.71431,0.68880,0.70071,0.62075,0.76345,1.05835,0.53343,0.25145,0.17170,0.36944,
#   0.21101,0.17483,0.12477,0.12894,0.27466,0.39818,0.60666,0.33268,0.43297,0.32044,0.38106,0.37968,0.68928,0.33411,0.20188,0.25076,
#   0.17283,0.15241,0.20782,0.17603,0.23082,0.30466,0.39359,0.43447,0.30474,0.27823,0.34361,0.26130,0.46794,0.76394,0.39551,0.36527,
#   0.03191,0.05615,0.12077,0.08700,0.09648,0.15435,0.25678,0.18681,0.32028,0.33532,0.41090,0.24662,0.50743,0.47868,0.65688,0.54496,
#   0.00000,0.02814,0.06460,0.06333,0.11065,0.26669,0.30602,0.17459,0.33583,0.30125,0.33912,0.28585,0.26830,0.41705,0.38301,1.15239
# ), nrow = 16, byrow = TRUE,
# dimnames = list(respondent = age_groups_5yr, contact = age_groups_5yr))
#
# #population sizes (single-year)
# single_yr <- c(
#   `0`=3480117, `1`=3532512, `2`=3672703, `3`=3797741, `4`=3917162,
#   `5`=4001330, `6`=3975522, `7`=4018926, `8`=4059908, `9`=4074737,
#   `10`=4251732,`11`=4268172,`12`=4421435,`13`=4388774,`14`=4297717,
#   `15`=4322527,`16`=4343161,`17`=4281824,`18`=4417941,`19`=4670623)
#
# #population for each 5-yr group used in the matrix
# N <-
#   list(
#     "lt1"  = sum(single_yr[c("0")]),
#     "1to4" = sum(single_yr[c("1","2","3","4")]),
#     "0to4" = sum(single_yr[as.character(0:4)]),
#     "5to9" = sum(single_yr[as.character(5:9)]),
#     "10to14" = sum(single_yr[as.character(10:14)]),
#     "15to17" = sum(single_yr[c("15","16","17")]),
#     "18to19" = sum(single_yr[c("18","19")]),
#     "15to19" = sum(single_yr[as.character(15:19)]),
#     # Approximate 2024 Census estimates for ages 20+
#     "20to24" = 21200000, "25to29" = 22800000, "30to34" = 23500000,
#     "35to39" = 22300000, "40to44" = 21000000, "45to49" = 20100000,
#     "50to54" = 21500000, "55to59" = 22700000, "60to64" = 21100000,
#     "65to69" = 18700000, "70to74" = 16500000, "75plus" = 24300000
#   )
#
# N_18plus <-
#   N[["18to19"]] + N[["20to24"]] + N[["25to29"]] + N[["30to34"]] + N[["35to39"]] + N[["40to44"]] + N[["45to49"]] +
#   N[["50to54"]] + N[["55to59"]] + N[["60to64"]] + N[["65to69"]] + N[["70to74"]] +  N[["75plus"]]
#
# N_5to17 <-
#   N[["5to9"]] + N[["10to14"]] + N[["15to17"]]
#
# #fractional weights for splitting mixed 5-yr groups
# w_lt1_in04    <- N[["lt1"]]    / N[["0to4"]]
# w_1to4_in04   <- N[["1to4"]]  /  N[["0to4"]]
# w_15to17_in19 <- N[["15to17"]] / N[["15to19"]]
# w_18to19_in19 <- N[["18to19"]] / N[["15to19"]]
#
# #aggregate ROWS to broad respondent age groups
# # M_broad[I, ] = Σ_{i∈I} (N_i / N_I) × M[i, ]
#
# #<1y and 1-4y, both mapped from 0-4 row (finest resolution available)
# r_lt1   <- M["0-4", ]
# r_1to4  <- M["0-4", ]
#
# #5-17y: pop-weighted average of 5-9, 10-14, and the 15-17 portion of 15-19
# r_5to17 <- (N[["5to9"]]   * M["5-9",   ] + N[["10to14"]] * M["10-14", ] + N[["15to17"]] * M["15-19", ]) / N_5to17
#
# #18+y: pop-weighted average of 18-19 portion of 15-19, plus 20-24 through 75+
# age20plus_groups <- c("20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75+")
# N_20plus_vec <- c(N[["20to24"]], N[["25to29"]], N[["30to34"]], N[["35to39"]],
#                   N[["40to44"]], N[["45to49"]], N[["50to54"]], N[["55to59"]],
#                   N[["60to64"]], N[["65to69"]], N[["70to74"]], N[["75plus"]])
#
# r_18plus <- (N[["18to19"]] * M["15-19", ] + colSums(N_20plus_vec * M[age20plus_groups, ])) / N_18plus
#
# #aggregate columns to broad contact age groups
# #for contact group J: sum M[i, j] for j∈J
# #(gives total contacts respondent i has with all individuals in group J)
# col_agg <- function(row_vec) {
#   c(`<1y`   = w_lt1_in04   * row_vec["0-4"],
#     `1-4y`  = w_1to4_in04  * row_vec["0-4"],
#     `5-17y` = row_vec["5-9"] + row_vec["10-14"] + w_15to17_in19 * row_vec["15-19"],
#     `18+y`  = w_18to19_in19 * row_vec["15-19"] + sum(row_vec[age20plus_groups]))
# }
#
# M_broad <- rbind( #THIS IS THE CONTACTS WE USE WITH ABCS POP BEFORE SYMMETRISE
#   `<1y`   = col_agg(r_lt1),
#   `1-4y`  = col_agg(r_1to4),
#   `5-17y` = col_agg(r_5to17),
#   `18+y`  = col_agg(r_18plus)
# )
# colnames(M_broad) <- c("<1y","1-4y","5-17y","18+y")
#
# #symmetrise
# #C_ij = (M_ij + M_ji × N_j / N_i) / 2
# pop_broad <- c(`<1y`   = N[["lt1"]],
#                `1-4y`  = N[["1to4"]],
#                `5-17y` = N_5to17,
#                `18+y`  = N_18plus)
#
# broad_labels <- c("<1y","1-4y","5-17y","18+y")
#
# C_sym <- matrix(NA, nrow = 4, ncol = 4,
#                 dimnames = list(respondent = broad_labels,
#                                 contact    = broad_labels))
#
# for (i in broad_labels) {
#   for (j in broad_labels) {
#     C_sym[i, j] <- (M_broad[i, j] + M_broad[j, i] * pop_broad[j] / pop_broad[i]) / 2
#   }
# }

# ==============================================================================

# # observed US ABCs population

# abcs <- readr::read_csv("https://raw.githubusercontent.com/PopHIVE/Ingest/refs/heads/main/data/abcs/raw/abcs_census_age_stratified_pop_full.csv")
#
# # 1. Drop redundant "Total" rows
# abcs <- abcs %>% dplyr::filter(age != "Total")
#
# # 2. Split 0-4 into 0-1 (20%) and 1-4 (80%) within each state-year
# split_04 <- abcs %>%
#   dplyr::filter(age == "0-4") %>%
#   dplyr::mutate(`0-1` = pop * 0.20,
#          `1-4` = pop * 0.80) %>%
#   dplyr::select(state, year, `0-1`, `1-4`) %>%
#   tidyr::pivot_longer(c(`0-1`, `1-4`), names_to = "age", values_to = "pop")
#
# abcs <- abcs %>%
#   dplyr::filter(age != "0-4") %>%
#   dplyr::bind_rows(split_04)
#
# #WITHIN STATE: average the 1998 and 1999 populations per (state, age),
# #relabel the combined year as 1999
# abcs <-
#   abcs %>%
#   dplyr::mutate(year = if_else(year %in% c(1998, 1999), 1999L, as.integer(year))) %>%
#   dplyr::group_by(state, age, year) %>%
#   dplyr::summarise(pop = mean(pop, na.rm = TRUE)) %>%
#   dplyr::ungroup()
#
# #Collapse 18-49, 50-64, 65+ into a single 18+ group, then
# #ACROSS STATES: sum population by year and age group (0-1, 1-4, 5-17, 18+)
# abcs <-
#   abcs %>%
#   dplyr::mutate(age = if_else(age %in% c("18-49", "50-64", "65+"), "18+", age)) %>%
#   dplyr::group_by(year, age) %>%
#   dplyr::summarise(pop_total = round(sum(pop, na.rm = TRUE))) %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(age = factor(age, levels = c("0-1", "1-4", "5-17", "18+"))) %>%
#   dplyr::arrange(year, age)

#data from the abc dataset above
pop_obs <-
  data.frame(
    year = 1999:2009,
    age0_1 = c(258310,331435,335992,366761,370365, 400515,414130,416014,420870,422993,436450),
    age1_4 = c(1033240,1325739,1343966,1467046,1481461, 1602061,1656520,1664058,1683480,1691973,1745798),
    age5_17 = c(3497789,4031650,4065801,4532325,4540649, 4908700,4916630,4956165,4983002,5011344,4968817),
    age18plus = c(14271440,16470672,16694393,18609391,18801641, 20392829,20650525,20984777,21282436,21565745,21395015)
  )

#years to generate
baseline <- pop_obs %>% dplyr::filter(year == 1999)
years_back <- 1980:1998

#empty dataframe
pop_back <- data.frame(year = years_back)

#backward projection function
back_project <- function(pop1999, years, growth = 0.01) {
  # Years before 1999
  t <- 1999 - years
  # Reverse compound growth
  projected <- pop1999 / ((1 + growth)^t)
  return(round(projected))
}

#generate backward populations at 1% annual growth since 1980
pop_back$age0_1 <- back_project(baseline$age0_1, years_back, growth = 0.01)
pop_back$age1_4 <- back_project(baseline$age1_4, years_back, growth = 0.01)
pop_back$age5_17 <- back_project(baseline$age5_17, years_back, growth = 0.01)
pop_back$age18plus <- back_project(baseline$age18plus, years_back, growth = 0.01)
pop_all <- bind_rows(pop_back, pop_obs) %>% arrange(year)

# ==============================================================================

#the M_broad matrix is what the code above gives you after running
#contact matrix in avg number of contacts per person (per capita) per day
M_broad <- matrix(c(
  0.35191855, 1.5087614, 1.095465, 3.958335,
  0.35191855, 1.5087614, 1.095465, 3.958335,
  0.04039284, 0.1731740, 4.581778, 4.792925,
  0.04835139, 0.2072943, 1.006943, 5.701542
), nrow = 4, byrow = TRUE)
n_age      <- 4L
age_groups <- c("<1y", "1-4y", "5-17y", "18+y")
rownames(M_broad) <- colnames(M_broad) <- age_groups

#age groups and ABCs population
pop_total <- c(258310, 1033240, 3497789, 14271441)
names(pop_total) <- age_groups

#symmetrise contact matrix using ABCs pop_total (avg number of contacts per person per day)
C_sym_day <- matrix(NA, n_age, n_age,
                    dimnames = list(age_groups, age_groups))
for (i in seq_len(n_age)) {
  for (j in seq_len(n_age)) {
    C_sym_day[i, j] <- (M_broad[i, j] + M_broad[j, i] * pop_total[j] / pop_total[i]) / 2
  }
}

#avg number of contacts per person per YEAR
C_year <- 365 * C_sym_day

##############################################################KK-Edited##############################################################
#carriage clearance per YEAR (= 365/duration_days)
duration_days <- c(71, 34, 18, 17)
r_per_year <- 365 / duration_days
names(r_per_year) <- age_groups
r_V <- r_per_year
r_A <- r_per_year
r_N <- r_per_year

#aging rates per YEAR
duration_years <- c(1, 4, 13, 59)
l_age          <- 1 / duration_years
names(l_age)   <- age_groups

#birth inflow composition (S only)
p_birth <- c(1, 0, 0, 0, 0, 0, 0)

#relative risk of competition for dual carriage
#eps_V, eps_A, eps_N are estimated in Stan

#observed counts: columns S (non-carriers), V, A, N
# studyPeriod 2001/2

y_obs <- matrix(c(
  222,  203,  7,  70,   # <1y    (N = 502)
  524,  277,  11, 159,   # 1-4y   (N = 971)
  361,   47,   2,  47,   # 5-17y  (N = 457)
  1789,  79,  3,  70    # 18+y   (N = 1938)
), nrow = n_age, byrow = TRUE)
rownames(y_obs) <- age_groups
colnames(y_obs) <- c("S", "V", "A", "N")
storage.mode(y_obs) <- "integer"
#y_obs <- matrix(c(
#  222,  179,  50,  51,   # <1y    (N = 502)
#  524,  294,  50, 103,   # 1-4y   (N = 971)
#  361,   51,   3,  42,   # 5-17y  (N = 457)
#  1789,  76,  24,  49    # 18+y   (N = 1938)
#), nrow = n_age, byrow = TRUE)
#rownames(y_obs) <- age_groups
#colnames(y_obs) <- c("S", "V", "F", "N")
#storage.mode(y_obs) <- "integer"

N_total <- as.integer(rowSums(y_obs))

#Observed prevalence
print(round(sweep(y_obs, 1, N_total, "/") * 100, 1))

#sanity check R model run (eps fixed at 0.5)
sis_rhs <- function(t, y, parms) {
  with(parms, {
    n <- n_age
    S  <- y[(0:(n - 1)) * 7 + 1]
    V  <- y[(0:(n - 1)) * 7 + 2]
    Ac <- y[(0:(n - 1)) * 7 + 3]
    Nc <- y[(0:(n - 1)) * 7 + 4]
    VA <- y[(0:(n - 1)) * 7 + 5]
    NV <- y[(0:(n - 1)) * 7 + 6]
    NAc <- y[(0:(n - 1)) * 7 + 7]

    sumV <- as.numeric(C %*% (V  + NV + VA))
    sumA <- as.numeric(C %*% (Ac + NAc + VA))
    sumN <- as.numeric(C %*% (Nc + NV + NAc))

    lamV <- rho_V * sumV
    lamA <- rho_A * sumA
    lamN <- rho_N * sumN

    dy <- numeric(7 * n)
    for (i in 1:n) {
      rVi <- r_V[i]; rAi <- r_A[i]; rNi <- r_N[i]; li <- l_age[i]

      dS  <- -(lamV[i] + lamA[i] + lamN[i]) * S[i] +
              rVi * V[i] + rAi * Ac[i] + rNi * Nc[i]
      dV  <-  lamV[i] * S[i] - rVi * V[i] -
              eps_V * (lamA[i] + lamN[i]) * V[i] +
              rAi * VA[i] + rNi * NV[i]
      dAc <-  lamA[i] * S[i] - rAi * Ac[i] -
              eps_A * (lamV[i] + lamN[i]) * Ac[i] +
              rVi * VA[i] + rNi * NAc[i]
      dNc <-  lamN[i] * S[i] - rNi * Nc[i] -
              eps_N * (lamV[i] + lamA[i]) * Nc[i] +
              rVi * NV[i] + rAi * NAc[i]
      dVA <-  eps_V * lamA[i] * V[i] +
              eps_A * lamV[i] * Ac[i] - (rVi + rAi) * VA[i]
      dNV <-  eps_V * lamN[i] * V[i] +
              eps_N * lamV[i] * Nc[i] - (rVi + rNi) * NV[i]
      dNAc <-  eps_A * lamN[i] * Ac[i] +
              eps_N * lamA[i] * Nc[i] - (rAi + rNi) * NAc[i]

      if (i == 1) {
        dS  <- dS  + li * (p_birth[1] - S[i])
        dV  <- dV  + li * (p_birth[2] - V[i])
        dAc <- dAc + li * (p_birth[3] - Ac[i])
        dNc <- dNc + li * (p_birth[4] - Nc[i])
        dVA <- dVA + li * (p_birth[5] - VA[i])
        dNV <- dNV + li * (p_birth[6] - NV[i])
        dNAc <- dNAc + li * (p_birth[7] - NAc[i])
      } else {
        dS  <- dS  + li * (S[i - 1]  - S[i])
        dV  <- dV  + li * (V[i - 1]  - V[i])
        dAc <- dAc + li * (Ac[i - 1] - Ac[i])
        dNc <- dNc + li * (Nc[i - 1] - Nc[i])
        dVA <- dVA + li * (VA[i - 1] - VA[i])
        dNV <- dNV + li * (NV[i - 1] - NV[i])
        dNAc <- dNAc + li * (NAc[i - 1] - NAc[i])
      }

      dy[(i - 1) * 7 + 1] <- dS
      dy[(i - 1) * 7 + 2] <- dV
      dy[(i - 1) * 7 + 3] <- dAc
      dy[(i - 1) * 7 + 4] <- dNc
      dy[(i - 1) * 7 + 5] <- dVA
      dy[(i - 1) * 7 + 6] <- dNV
      dy[(i - 1) * 7 + 7] <- dNAc
    }
    list(dy)
  })
}

y0 <- rep(c(0.85, 0.05, 0.05, 0.05, 0, 0, 0), n_age)

#default age-specific rho values (literature-informed prior centres) used for the pre-fit check;
rho_age_default <- c(0.05, 0.03, 0.015, 0.005)

#put all fixed parameters in a list
parms0 <- list(
  n_age   = n_age,
  C       = C_year,
  r_V     = r_V,  r_A = r_A,  r_N = r_N,
  l_age   = l_age,
  p_birth = p_birth,
  rho_V   = rho_age_default,
  rho_A   = rho_age_default,
  rho_N   = rho_age_default,
  eps_V   = 0.5,
  eps_A   = 0.5,
  eps_N   = 0.5
)

#run the model
times <- seq(0, 30, by = 0.25)
sol   <- deSolve::ode(y0, times, sis_rhs, parms = parms0, method = "lsoda")
final <- tail(sol, 1)[-1]

#equilibrium per-age-group totals should still sum to 1
totals <- sapply(1:n_age, function(i) sum(final[(i - 1) * 7 + 1:7]))
print(setNames(round(totals, 6), age_groups))

#predicted prevalence
predicted_prev <- vapply(1:n_age, function(i) {
  Si  <- final[(i - 1) * 7 + 1]
  Vi  <- final[(i - 1) * 7 + 2]
  Ai  <- final[(i - 1) * 7 + 3]
  Ni  <- final[(i - 1) * 7 + 4]
  VAi <- final[(i - 1) * 7 + 5]
  NVi <- final[(i - 1) * 7 + 6]
  NAi <- final[(i - 1) * 7 + 7]
  c(S = Si,
    V = Vi + 0.5 * VAi + 0.5 * NVi,
    A = Ai + 0.5 * VAi + 0.5 * NAi,
    N = Ni + 0.5 * NVi + 0.5 * NAi)
}, numeric(4))

predicted_prev <- t(predicted_prev)
rownames(predicted_prev) <- age_groups
colnames(predicted_prev) <- c("S", "V", "A", "N")

#compare how close observes vs predicted prevalence are
observed_prev <- sweep(y_obs, 1, N_total, "/")
print(round(observed_prev * 100, 1))
print(round(predicted_prev * 100, 1))

# ==============================================================================

#compile and fit Stan model

#age susceptibility to carriage
rho_age_logmean <- log(c(0.05, 0.03, 0.015, 0.005))   #length n_age = 4
rho_age_logsd   <- 0.5

#Beta(2, 2) on eps_X: weak preference for 0.5, sd ~0.224
eps_alpha <- 2
eps_beta  <- 2

#stan list
stan_data <- list(

  #fixed parameters
  n_age           = n_age,
  C               = C_year,
  r_V             = r_V,
  r_A             = r_A,
  r_N             = r_N,
  l_age           = l_age,
  p_birth         = p_birth,
  N_total         = N_total,
  y_obs           = y_obs,
  t0              = 0,
  ts              = array(30, dim = 1),

  #informative priors
  rho_age_logmean = rho_age_logmean,
  rho_age_logsd   = rho_age_logsd,
  eps_alpha       = eps_alpha,
  eps_beta        = eps_beta,

  rel_tol         = 1e-5,
  abs_tol         = 1e-6,
  max_steps       = 20000L
)

#compile stan model
#setwd("/Users/kek57/Desktop/WeinbergerLab/Code/")
#spn_model <- cmdstanr::cmdstan_model(here::here("", "02_pneumocarriage_model.stan"))
spn_model <- cmdstanr::cmdstan_model(here::here("","02_19A_carriage_model.stan"))

#set initial starting values for MCMC
init_fun <- function(chain_id) {
  list(
    #4 age log-baselines, jittered around informative centres
    log_rho_age = rho_age_logmean + rnorm(n_age, 0, 0.1),

    #3 competition parameters near 0.5 with small jitter
    eps_V = runif(1, 0.4, 0.6),
    eps_ = runif(1, 0.4, 0.6),
    eps_N = runif(1, 0.4, 0.6)
  )
}

# #run the stan model under HMC NUTS (takes about 1 hour to run on laptop)
carriage_fit <- spn_model$sample(
  data            = stan_data,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 4000,
  seed            = 1234,
  adapt_delta     = 0.99,
  max_treedepth   = 10,
  init            = init_fun,
  refresh         = 100
)
#
# #persist fit and bundles for post-processing
carriage_fit$save_object(file = "carriage_fit_19A.rds")
