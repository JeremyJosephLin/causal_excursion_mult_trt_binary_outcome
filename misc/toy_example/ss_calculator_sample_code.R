# This is the R_code example on how to use binary sample size calculator 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("wcls_cat_trt_binary_outcome.R")
source("ss_calc.R")
library(dplyr)
library(rootSolve)

# toy dataset with 3 treatment arm with time and SEX covariates
# 50 individual with 15 decision points each
dta <- readRDS("toy_data.RDS")


#================= estimating Marginal Excursion effect =========================
# example 1 where we marginalize the treatment effect for each treatment 
# set moderator and control variable as empty set

id_varname = c("userid")
decision_time_varname = c("time")
treatment_varname = c("A")
avail_varname = c("I")
rand_prob_varname = "prob_A"
trt_level = 3
outcome_varname = c("Y")
moderator_varname = c()
control_varname <- c()

fit1 <- wcls_categorical_treatment(
  dta = dta,
  id_varname = id_varname,
  decision_time_varname = decision_time_varname,
  treatment_varname = treatment_varname,
  outcome_varname = outcome_varname,
  control_varname = control_varname,
  moderator_varname = moderator_varname,
  rand_prob_varname = rand_prob_varname,
  avail_varname = avail_varname,
  trt_level = trt_level,
  estimator_initial_value = NULL
)

# Estimated Marginal Effect 
fit1$beta_hat
# Intercept : trt = 1 Intercept : trt = 2 Intercept : trt = 3 
# -0.1240791          -0.3103041          -0.2582578 

# estimated log ASPN (for sample size estimation)
fit$alpha_hat

#Confidence Interval with small sample size correction
fit1$conf_int_adjusted_t

# example where time is used as moderator variable for treatment effect 
# moderator variable need to be a subset of control variable
# assume both time and sex as the control variable

moderator_varname = c("time")
control_varname <- c("time", "sex")

fit2 <- wcls_categorical_treatment(
  dta = dta,
  id_varname = id_varname,
  decision_time_varname = decision_time_varname,
  treatment_varname = treatment_varname,
  outcome_varname = outcome_varname,
  control_varname = control_varname,
  moderator_varname = moderator_varname,
  rand_prob_varname = rand_prob_varname,
  avail_varname = avail_varname,
  trt_level = trt_level,
  estimator_initial_value = NULL
)

# Estimated Marginal Effect 
fit2$beta_hat

# Intercept : trt = 1      time : trt = 1 Intercept : trt = 2      time : trt = 2 Intercept : trt = 3      time : trt = 3 
# -0.100889148        -0.003976773        -0.101380049        -0.027056126        -0.375833136         0.014001695


# Confidence Interval with small sample size correction
fit2$conf_int_adjusted_t

#============== Sample Size Calculator Example =================================
# it is recommended to use constant spnc and constant mee to obtain more conservative estimate for sample size estimator
# estimated using the first example on the previous section 
ft <- matrix(rep(1, 15), ncol = 1)
gt <-  matrix(rep(1, 15), ncol = 1)
beta <- fit1$beta_hat
alpha <- fit1$alpha_hat


## pt probability assigning treatment are all equal
pt <- matrix(rep(c(0.25, 0.25, 0.25, 0.25),times = 15 ), ncol = 4, byrow = TRUE)

# availability : assume user are always available at any decision point
tau <- rep(1, 15)

#type I error
gamma <- 0.05
# type 2 error
b <-  0.2

# matrix L is used to calculate sample size based on some hypothesis 

# testing to detect statistical difference on treatment arm 1 
L1 <- matrix(c(1,0,0), nrow = 1)

n1 <- mrt_mult_trt_ss(
  f_t = ft,
  g_t = gt,
  beta = beta,
  alpha = alpha,
  p_t = pt,
  gamma = gamma,
  b = b,
  L = L1,
  tau = tau
)

# testing to detect statistical difference between treatment arm1 and treatment arm 2
L2 <- matrix(c(1, -1, 0), nrow = 1)

n2 <- mrt_mult_trt_ss(
  f_t = ft,
  g_t = gt,
  beta = beta,
  alpha = alpha,
  p_t = pt,
  gamma = gamma,
  b = b,
  L = L2,
  tau = tau
)
