# This simulation is to check consistency of our estimator under control misspecification
# Zt is a moderator variable, here the f(H_t) is correctly specified,
# Simulation for marginal model 

rm(list = ls(all = TRUE))
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# call functions
source("emee_catA.R")
source("dgm.R")

# load library fpr GEE and Lmm
library(geepack) #gee


# grab factorial design
SD <- expand.grid(sample_sizes = c(20, 30, 40, 50, 100),
                  total_Ts = c(30))

# get command parameters
args <- commandArgs(trailingOnly = TRUE)
isetting <- as.integer(args[1]) # settings
nsim <- as.integer(args[2])
nsetting <- as.integer(args[3])

# # for debugging purposes isetting any value between 1-6(length of sample size), nsim = 10, nsetting = 4(equal to the number of total_T)
# isetting = 1
# nsim = 1
# nsetting = dim(SD)[1]

setting_start <- (isetting - 1)* nsetting + 1
setting_end <- isetting * nsetting

print(paste0("setting value ", isetting, " nsim ", nsim, " nsetting ", nsetting, " setting start ",setting_start, " setting end ",setting_end))

library(rootSolve) # for solver function multiroot()


# ------------------------- run simulations ------------------------------------ 

control_vars <- c("S")
moderator_vars <- c("S")
avail_varname = "I"


result_list_collected <- list()

i = 1

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
  
  
  print(Sys.time())
  print(paste0("Sample size ", sample_size, " total T ", total_T))
  
  # changing the number of seed for every simulation
  set.seed(i)
  
  for (i_sim in 1:nsim) {
    # Generate data
    dta <- dgm_1(
      sample_size = sample_size,
      total_T = total_T
    )
    
   
    # Fit TQ estimator
    fit_emee  <- emee_catA(
        dta = dta,
        id_varname = "userid",
        decision_time_varname = "time",
        treatment_varname = "A",
        outcome_varname = "Y",
        control_varname = c("S"),
        moderator_varname = c("S"),
        rand_prob_varname = "prob_A",
        estimator_initial_value = NULL,
        trt_level = 2,
        pmatrix_tilde = NULL,
        avail_varname = NULL
      )
      beta_hat_emee = fit_emee$beta_hat
      beta_se_emee = fit_emee$beta_se
      varcov_emee = fit_emee$varcov
      beta_se_adjusted_emee = fit_emee$beta_se_adjusted
      varcov_adjusted_emee = fit_emee$varcov_adjusted
      ci_unadj_emee = fit_emee$conf_int
      ci_adj_z_emee = fit_emee$conf_int_adjusted_z
      ci_adj_t_emee = fit_emee$conf_int_adjusted_t
      p_tilde = fit_emee$p_tilde
      
      # fit dataset with GEE and Independent Variance covariance matrix
      fit_GEE1  <- geeglm(Y ~ A1*S + A2*S, data = dta, id = userid,
                      corstr = "ind",family = poisson(link = "log"))
      
      fit_summary_GEE1 <- summary(fit_GEE1)
      beta_hat_GEE1 = fit_summary_GEE1$coefficients[c("A1","A1:S", "A2", "S:A2"), "Estimate"]
      beta_se_GEE1 = fit_summary_GEE1$coefficients[c("A1","A1:S", "A2", "S:A2"), "Std.err"]
      varcov_GEE1 = fit_summary_GEE1$cov.scaled[c(2,5,4,6), c(2, 5, 4, 6)]
      ci_lower_GEE1 <- beta_hat_GEE1 - 1.96 * beta_se_GEE1
      ci_upper_GEE1 <- beta_hat_GEE1 + 1.96 * beta_se_GEE1
      ci_GEE1 = cbind(ci_lower_GEE1, ci_upper_GEE1)

    # fit dataset with GEE and Exchangeable Variance covariance matrix
    fit_GEE2  <- geeglm(Y ~ A1*S + A2*S, data = dta, id = userid,
                        corstr = "ind",family = poisson(link = "log"))
    
    fit_summary_GEE2 <- summary(fit_GEE2)
    beta_hat_GEE2 = fit_summary_GEE2$coefficients[c("A1","A1:S", "A2", "S:A2"), "Estimate"]
    beta_se_GEE2 = fit_summary_GEE2$coefficients[c("A1","A1:S", "A2", "S:A2"), "Std.err"]
    varcov_GEE2 = fit_summary_GEE2$cov.scaled[c(2,5,4,6), c(2, 5, 4, 6)]
    ci_lower_GEE2 <- beta_hat_GEE2 - 1.96 * beta_se_GEE2
    ci_upper_GEE2 <- beta_hat_GEE2 + 1.96 * beta_se_GEE2
    ci_GEE2 = cbind(ci_lower_GEE2, ci_upper_GEE2)

    
    output <- list(
      list(
        beta_hat_emee = beta_hat_emee,
        beta_se_emee = beta_se_emee,
        varcov_emee = varcov_emee,
        beta_se_adjusted_emee = beta_se_adjusted_emee,
        varcov_adjusted_emee = varcov_adjusted_emee,
        ci_unadj_emee = ci_unadj_emee,
        ci_adj_z = ci_adj_z_emee,
        ci_adj_t = ci_adj_t_emee,
        p_tilde = p_tilde,
        beta_hat_GEE1 = beta_hat_GEE1, 
        beta_se_GEE1 = beta_se_GEE1, 
        varcov_GEE1 = varcov_GEE1, 
        ci_GEE1 = ci_GEE1, 
        beta_hat_GEE2 = beta_hat_GEE2, 
        beta_se_GEE2 = beta_se_GEE2, 
        varcov_GEE2 = varcov_GEE2, 
        ci_GEE2 = ci_GEE2 
      )
    )
    result <- c(result,output)
  }
      
    # update result list
    result_list_collected <- c(result_list_collected, list(
      list(
      sample_size = sample_size,
      total_T = total_T,
      result = result
    )
  ))
  
}

outName <- paste("results_consistency_sim_setting_",isetting,"_",nsim,".RDS",sep = "")
saveRDS(result_list_collected,file=outName)
