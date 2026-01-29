# This simulation is to check consistency of our estimator under control misspecification
# Zt is a moderator variable, here the f(H_t) is correctly specified,
# added a linear time effect on generative model
rm(list = ls(all = TRUE))
# call functions
# source("R_code/wcls_cat_trt_binary_outcome.R")
# source("R_code/functions/utillity.R")
source("wcls_cat_trt_binary_outcome.R")
source("utillity.R")
# grab factorial design

SD <- expand.grid(sample_sizes = c(20, 30, 40, 50, 75, 100), 
                  total_Ts = c( 15, 20, 25, 30))

# get command parameters
args <- commandArgs(trailingOnly = TRUE)
isetting <- as.integer(args[1]) # settings 
nsim <- as.integer(args[2])
nsetting <- as.integer(args[3])

# for debugging purposes isetting any value between 1-6(length of sample size), nsim = 10, nsetting = 4(equal to the number of total_T)
# isetting = 3
# nsim = 10 
# nsetting = 4

setting_start <- (isetting - 1)* nsetting + 1
setting_end <- isetting * nsetting

print(paste0("setting value ", isetting, " nsim ", nsim, " nsetting ", nsetting, " setting start ",setting_start, " setting end ",setting_end))

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


# ------------------------- run simulations ------------------------------------ 

control_vars <- c("time")
moderator_vars <- c("time")
avail_varname = "I"


result_list_collected <- list()

for (i in setting_start:setting_end) {
  result <- list()
  
  start_time <- Sys.time()
  current_time <- Sys.time()
  print(Sys.time())
  print(paste0(round(
    difftime(current_time, start_time, units = "hours"), 2
  ), " hours has lapsed."))
  cat("i =", i, "/", nrow(SD), "\n")
  
  # for each variation of sample sizes and total T do n_sim number of simulation
  total_T = SD[i, "total_Ts"]
  sample_size = SD[i, "sample_sizes"]
  ## pt
  
  pt <- construct_pt(total_T, "constant equal pk")
  
  ## taut_t
  taut_t <-
    construct_taut_theta(total_T , 0.5, "constant")
  
  ## gt_t, alpha_t
  gt_t <- construct_ftgt(total_T, "linear_theta")
  alpha_t <-
    solve_alpha("linear_theta", total_T, 0.2, gt_t, taut_t, theta = 0.3)
  
  ## ft_t, beta1_t, beta2_t
  ft_t <- construct_ftgt(m = total_T, "linear_theta")
  beta1_t <- solve_beta(
    shape = "linear_theta",
    m = total_T,
    ATE = 1.2,
    ft = ft_t,
    taut = taut_t,
    theta = -0.3,
    gt = gt_t,
    alpha = alpha_t
  )
  beta2_t <- solve_beta(
    "linear_theta",
    m = total_T,
    1.6,
    ft = ft_t,
    taut = taut_t,
    theta = -0.3,
    gt = gt_t,
    alpha = alpha_t
  )
  beta_t <- c(beta1_t, beta2_t)
  
  
  print(Sys.time())
  print(paste0("Sample size ", sample_size, " total T ", total_T))
  
  # changing the number of seed for every simulation
  set.seed(i)
  
  for (i_sim in 1:nsim) {
    # Generate data
    dta <- dgm(
      sample_size = sample_size,
      total_T = total_T,
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
    
    # Fit TQ estimator
    fit  <- wcls_categorical_treatment(
      dta = dta,
      id_varname = "userid",
      decision_time_varname = "time",
      treatment_varname = "A",
      outcome_varname = "Y",
      control_varname = c("time"),
      moderator_varname = c("time"),
      rand_prob_varname = "prob_A",
      estimator_initial_value = NULL,
      trt_level = 2,
      pmatrix_tilde = NULL,
      avail_varname = NULL
    )
    output <- list(
      list(
        beta_hat = fit$beta_hat,
        beta_true = beta_t,
        beta_se = fit$beta_se,
        varcov = fit$varcov,
        beta_se_adjusted = fit$beta_se_adjusted,
        varcov_adjusted = fit$varcov_adjusted,
        ci_unadj = fit$conf_int,
        ci_adj_z = fit$conf_int_adjusted_z,
        ci_adj_t = fit$conf_int_adjusted_t,
        p_tilde = fit$p_tilde,
        data = dta
      )
    )
    result <- c(result,output)
  }
      
    # update result list
    result_list_collected <- c(result_list_collected, list(
      list(
      sample_size = sample_size,
      total_T = total_T,
      true_beta = beta_t,
      result = result
    )
  ))
  
}

outName <- paste("results_consistency_sim_setting_",isetting,"_",nsim,".RDS",sep = "")
saveRDS(result_list_collected,file=outName)+
