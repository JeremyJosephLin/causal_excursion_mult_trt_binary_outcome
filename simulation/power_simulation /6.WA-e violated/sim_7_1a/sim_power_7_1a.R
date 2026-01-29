rm(list = ls(all = TRUE))
# call functions
#setwd("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/simulations/power_simulation /6.WA-e violated/sim_7_1a/")
# grab factorial design
SD <- readRDS("factorial_design.RDS")
source("wcls_cat_trt_binary_outcome.R")
source("utillity.R")
source("dgm_endoAvail.R")

# calculate either type1 error(H0) or power(H1) 
truth_hypothesis = "H1"

# get command parameters
args <- commandArgs(trailingOnly = TRUE)
isetting <- as.integer(args[1]) # settings 
nsim <- as.integer(args[2])
nsetting <- as.integer(args[3])

# for debugging purposes isetting any value between 1-6(length of sample size), nsim = 10, nsetting = 4(equal to the number of total_T)
# isetting = 1
# nsim = 1
# nsetting = 729

setting_start <- (isetting - 1)* nsetting + 1
setting_end <- isetting * nsetting

print(paste0("setting value ", isetting, " nsim ", nsim, " nsetting ", nsetting, " setting start ",setting_start, " setting end ",setting_end))
  
set.seed(isetting)
# 2. Simulation(parallel) -------------------------------------------------------
options(warn = 0)


library(rootSolve) # for solver function multiroot()

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
  # coefficients for endoAvail
  gamma3 <- SD$gamma3[i]
  gamma4 <- SD$gamma4[i]
  gamma5 <- SD$gamma5[i]
  
  ## gamma, b, m
  gamma <- SD$gamma[i]
  b <- SD$b[i]
  m <- SD$m[i]
  
  ## pt
  pt <- construct_pt(m, SD$pt_shape[i])
  
  ## gt_t, alpha_t
  gt_t <- construct_ftgt(m, SD$gt_t_shape[i])
  SD$q[i] <- ncol(gt_t)
  alpha_t <- log(SD$ASPNC_t[i]) # only for constant gt
  
  SD[i, c("alpha_t0", "alpha_t1", "alpha_t2")[1:SD$q[i]]] <- alpha_t
  #SD$ASPNC_t_empirical[i] <- compute_ASPNC(gt_t, alpha_t, taut_t)
  #stopifnot(all.equal(SD$ASPNC_t[i], SD$ASPNC_t_empirical[i]))
  
  ## ft_t, beta1_t, beta2_t
  ft_t <- construct_ftgt(m, SD$ft_t_shape[i])
  SD$p[i] <- ncol(ft_t)
  # trt 1 effect
  beta1_t <- log(SD$ATE1_t[i]) #only for constant ft see goodnotes on simulation setting for detail 
  SD[i, c("beta1_t0", "beta1_t1", "beta1_t2")[1:SD$p[i]]] <- beta1_t
  
  
  # trt 2 effect
  beta2_t <- log(SD$ATE2_t[i])
  SD[i, c("beta2_t0", "beta2_t1", "beta2_t2")[1:SD$p[i]]] <- beta2_t
  
  beta_t <- c(beta1_t, beta2_t)
  L <- construct_L(SD$ft_w_shape[i])
  
  ## taut_t
  taut_t <- compute_tau_endoAvail(gt = gt_t, alpha = alpha_t, ft = ft_t, beta = beta_t, pt = pt, gamma3 = gamma3, gamma4 = gamma4, gamma5 = gamma5)
  SD$AvgTau_t[i] <- mean(taut_t)
  
  
  ## gt_w, alpha_w
  gt_w <- gt_t
  
  ## ft_w
  ft_w <- ft_t
  
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
    dta <- dgm_endoAvail(
      sample_size = SD$n[i],
      total_T = m,
      ft = ft_t,
      beta_1 = beta1_t,
      beta_2 = beta2_t,
      gt = gt_t,
      alpha = alpha_t,
      tau = taut_t,
      pt = pt,       
      gamma3  = gamma3,
      gamma4 = gamma4,
      gamma5 = gamma5
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


outName <- paste("results_sim_7_1_setting_",isetting,"_",nsim,".RDS",sep = "")
saveRDS(result_list_collected,file=outName)
