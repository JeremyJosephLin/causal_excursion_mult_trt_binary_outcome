rm(list = ls(all = TRUE))
# call functions
#setwd("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/simulations/typeI_error_simulation/1. all working assumptions hold/")
# grab factorial design
SD <- readRDS("factorial_design.RDS")
source("wcls_cat_trt_binary_outcome.R")
source("utillity.R")

# calculate either type1 error(H0) or power(H1) 
truth_hypothesis = "H0"

# get command parameters
args <- commandArgs(trailingOnly = TRUE)
isetting <- as.integer(args[1]) # settings 
nsim <- as.integer(args[2])
nsetting <- as.integer(args[3])

# for debugging purposes isetting any value between 1-6(length of sample size), nsim = 10, nsetting = 4(equal to the number of total_T)
# isetting = 1
# nsim = 1
# nsetting = 1617

setting_start <- (isetting - 1)* nsetting + 1
setting_end <- isetting * nsetting

print(paste0("setting value ", isetting, " nsim ", nsim, " nsetting ", nsetting, " setting start ",setting_start, " setting end ",setting_end))
  
set.seed(isetting)
# 2. Simulation(parallel) -------------------------------------------------------
options(warn = 0)


library(rootSolve) # for solver function multiroot()

dgm <- function(sample_size, total_T, ft, beta_1, beta_2, gt, alpha, tau, pt, mean_var = 1, j_t = NULL, h_At = NULL) {
  # variable A represent the treatment level  (0 , 1 , 2)
  # variable A1 and A2 is indicator function whether A is equal 1 or 2 respectively
  # variable I is a binary indicating availability of the participant at timepoint t 
  # prob_t indicates the probability of someone available at time point t, assume constant for now
  # tau is E(I_t) which we assume is equal 1 (always available at any given time) unless specifiedd otherwise
  # if there is no exp_var provided, assume error is Normal 0,1
  # AA = 0.3
  
  if (is.null(j_t)) {
    j_t <- rep(mean_var, total_T)
  }
  
  if (is.null(h_At)) {
    h_At <- matrix(rep(c(1, 0, 0), total_T), ncol = dim(pt)[2], byrow = TRUE)
  }
  
  
  # treatment name(only for 2 level of treatment for now)
  a_val <- c(0, 1, 2)
  
  df_names <- c("userid", "time", "Y", "A", "A1", "A2", "I", "prob_A","prob_I", "Y_prob1")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  
  for (t in 1:total_T) {
    prob_a <- pt[t,]
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    # generate Availability 
    dta$I[row_index] <- rbinom(sample_size, 1, prob = tau[t])
    dta$prob_I[row_index] <- ifelse(dta$I[row_index] == 0, 1 - tau[t], tau[t])
    # generate treatment 
    dta$A[row_index] <- sample(x = a_val, size = sample_size, replace = TRUE, prob = prob_a)
    # set treatment to 0 when availability is 0 
    dta$A[row_index] <- dta$A[row_index] * dta$I[row_index]
    # record probability for each treatment
    dta$prob_A[row_index] <- ifelse(dta$A[row_index] == 0, prob_a[1], 
                                    ifelse(dta$A[row_index] == 1, prob_a[2], prob_a[3]))
    # create pseudo variable for indicator function
    dta$A1[row_index] <- ifelse(dta$A[row_index] == 1, 1, 0)
    dta$A2[row_index] <- ifelse(dta$A[row_index] == 2, 1, 0)
    # generate random error 
    h_Ai <- h_At[t, 1] + h_At[t, 2] * dta$A1[row_index]+h_At[t, 3] * dta$A2[row_index]
    #dta$eps[row_index] <- rnorm(sample_size, mean = 0, sd = 1) * sqrt(j_t[t]) * sqrt(h_Ai)
    dta$Y_prob1[row_index] <- exp(dta$A1[row_index] * as.numeric(ft[t, ] %*% beta_1)  +
                                    dta$A2[row_index] * as.numeric(ft[t, ] %*% beta_2) +
                                    as.numeric(gt[t, ] %*% alpha))
    if (!(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))) {
      message("min(dta$Y_prob1) = ", min(dta$Y_prob1), 
              "; max(dta$Y_prob1) = ", max(dta$Y_prob1))
    }
    stopifnot(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))
    dta$Y[row_index] <- rbinom(sample_size, 1, prob = dta$Y_prob1[row_index])
  }
  
  return(dta)
}

result_list_collected <- list()

for (i in setting_start:setting_end) {
  result <- list()
  
  start_time <- Sys.time()
  current_time <- Sys.time()
  print(Sys.time())
  print(paste0(round(
    difftime(current_time, start_time, units = "hours"), 2
  ),
  " hours has lapsed."))
  cat("i =", i, "/", nrow(SD), "\n")
  
  #print(SD[i, ])
  
  ## gamma, b, m
  gamma <- SD$gamma[i]
  b <- SD$b[i]
  m <- SD$m[i]
  
  ## pt
  pt <- construct_pt(m, SD$pt_shape[i])
  
  ## taut_t
  taut_t <-
    construct_taut_theta(m , SD$AvgTau_t[i], SD$tau_t_shape[i])
  
  ## gt_t, alpha_t
  gt_t <- construct_ftgt(m, SD$gt_t_shape[i])
  alpha_t <-
    as.numeric(SD[i, c("alpha_t0", "alpha_t1", "alpha_t2")[1:SD$q[i]]])
  
  ## ft_t, beta1_t, beta2_t
  ft_t <- construct_ftgt(m, SD$ft_t_shape[i])
  beta1_t <-
    as.numeric(SD[i, c("beta1_t0", "beta1_t1", "beta1_t2")[1:SD$p[i]]])
  beta2_t <-
    as.numeric(SD[i, c("beta2_t0", "beta2_t1", "beta2_t2")[1:SD$p[i]]])
  beta_t <- c(beta1_t, beta2_t)
  
  L <- construct_L(SD$ft_w_shape[i])
  
  ## taut_w
  # (not used in the following)
  
  ## gt_w, alpha_w
  gt_w <- construct_ftgt(m, SD$gt_w_shape[i])
  alpha_w <-
    as.numeric(SD[i, c("alpha_w0", "alpha_w1", "alpha_w2")[1:SD$q_w[i]]])
  
  ## ft_w
  ft_w <- construct_ftgt(m, SD$ft_w_shape[i])
  beta1_w <-
    as.numeric(SD[i, c("beta1_w0", "beta1_w1", "beta1_w2")[1:SD$p_w[i]]])
  beta2_w <-
    as.numeric(SD[i, c("beta2_w0", "beta2_w1", "beta2_w2")[1:SD$p_w[i]]])
  beta_w <- c(beta1_w, beta2_w)
  
  ## n, p_w, q_w
  n <- SD$n[i]
  p_w <- SD$p_w[i]
  q_w <- SD$q_w[i]
  p_star <- SD$p_star[i] # This is the degree of freedom of T test
  p <- 2 * p_w # this is the dimension of beta
  q <- q_w
  
  ##### simulation begins #####
  beta1_t_H1 <- beta1_t
  beta2_t_H1 <- beta2_t
  
  if (truth_hypothesis == "H0") {
    beta1_t <- rep(0, length(beta1_t_H1))
    beta2_t <- rep(0, length(beta1_t_H1))
  } else if (truth_hypothesis == "H1") {
    beta1_t <- beta1_t_H1
    beta2_t <- beta2_t_H1
  }
  
  for (i_sim in 1:nsim) {
    # Generate data
    dta <- dgm(
      sample_size = SD$n[i],
      total_T = m,
      ft = ft_t,
      beta_1 = beta1_t,
      beta_2 = beta2_t,
      gt = gt_t,
      alpha = alpha_t,
      tau = taut_t,
      pt = pt,
      mean_var = 1,
      h_At = NULL,
      j_t = NULL
    )
    
    dta <- add_analysis_vars(dta, n, ft_w, gt_w)
    control_varname <- get_control_varname(gt_w)
    moderator_varname <- get_moderator_varname(ft_w)
    avail_varname = "I"
    
    # Fit TQ estimator
    fit_TQ  <- wcls_categorical_treatment(
      dta = dta,
      id_varname = "userid",
      decision_time_varname = "time",
      treatment_varname = "A",
      outcome_varname = "Y",
      control_varname = control_varname,
      moderator_varname = moderator_varname,
      rand_prob_varname = "prob_A",
      estimator_initial_value = NULL,
      trt_level = 2,
      pmatrix_tilde = NULL,
      avail_varname = avail_varname
    )
    
    beta_hat <- fit_TQ$beta_hat
    var_beta <- fit_TQ$varcov[(q + 1):(p + q), (q + 1):(p + q)]
    var_beta_adj <-
      fit_TQ$varcov_adjusted[(q + 1):(p + q), (q + 1):(p + q)]
    
    # without small sample correction
    test_stat_unadj <-
      t(L %*% beta_hat) %*% solve(L %*% var_beta %*% t(L)) %*% (L %*% beta_hat)
    
    # with small sample correction
    test_stat_adj <-
      t(L %*% beta_hat) %*% solve(L %*% var_beta_adj %*% t(L)) %*% (L %*% beta_hat)
    
    print(paste0("test statistic adjusted ",test_stat_adj))
    print(paste0("test statistic unadjusted ",test_stat_unadj))
    
    output <- list(
      list(
        fit = fit_TQ,
        test_stat_unadj = test_stat_unadj,
        test_stat_adj = test_stat_adj
      )
    )
    result <- c(result,output)
  }
  result_list_collected <- c(result_list_collected, list(list(
    setting = i,
    result = result
  )))
}


outName <- paste("results_sim_1_1_setting_",isetting,"_",nsim,".RDS",sep = "")
saveRDS(result_list_collected,file=outName)
