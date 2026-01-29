# This code to fill in SD using sample size calculator no need to run it in HPC since setting is pretty small
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("wcls_cat_trt_binary_outcome.R")
source("utillity.R")
source("ss_calc.R")
source("dgm.R")
library(rootSolve)
SD <- readRDS("factorial_design.RDS")

##### fill the SD #####
# get command parameters
# args <- commandArgs(trailingOnly = TRUE)
# isetting <- as.integer(args[1]) # settings 
# nsim <- as.integer(args[2])
# nsetting <- as.integer(args[3])

isetting = 1
nsim = 1
nsetting = dim(SD)[1]

setting_start <- (isetting - 1)* nsetting + 1
setting_end <- isetting * nsetting

print(paste0("setting value ", isetting, " nsim ", nsim, " nsetting ", nsetting, " setting start ",setting_start, " setting end ",setting_end))


options(warn = 2)
for (i in setting_start:setting_end) {
  
  ## gamma, b, m
  gamma <- SD$gamma[i]
  b <- SD$b[i]
  m <- SD$m[i]
  
  ## pt
  pt <- construct_pt(m, SD$pt_shape[i])
  
  ## taut_t
  taut_t <- construct_taut_theta(m , SD$AvgTau_t[i], SD$tau_t_shape[i], theta = SD$tau_t_theta[i])
  
  ## gt_t, alpha_t
  gt_t <- construct_ftgt(m, SD$gt_t_shape[i])
  SD$q[i] <- ncol(gt_t)
  alpha_t <- solve_alpha(SD$gt_t_shape[i], m, SD$ASPNC_t[i], gt_t, taut_t, 
                         theta = SD$gt_t_theta[i])
  SD[i, c("alpha_t0", "alpha_t1", "alpha_t2")[1:SD$q[i]]] <- alpha_t
  SD$ASPNC_t_empirical[i] <- compute_ASPNC(gt_t, alpha_t, taut_t)
  stopifnot(all.equal(SD$ASPNC_t[i], SD$ASPNC_t_empirical[i]))
  
  ## ft_t, beta1_t, beta2_t
  ft_t <- construct_ftgt(m, SD$ft_t_shape[i])
  SD$p[i] <- ncol(ft_t)
  # trt 1 effect
  beta1_t <- solve_beta(
    shape = SD$ft_t_shape[i],
    m = m,
    ATE = SD$ATE1_t[i],
    ft = ft_t,
    taut = taut_t,
    theta = SD$ft_t_theta[i],
    gt = gt_t,
    alpha = alpha_t
  )
  SD[i, c("beta1_t0", "beta1_t1", "beta1_t2")[1:SD$p[i]]] <- beta1_t
  SD$ATE1_t_empirical[i] <- compute_ATE(gt=gt_t, alpha=alpha_t, ft = ft_t, beta = beta1_t, taut= taut_t)
  stopifnot(all.equal(SD$ATE1_t[i], SD$ATE1_t_empirical[i]))
  
  # trt 2 effect
  beta2_tstar <- solve_beta_star(
    shape = SD$ft_t_shape[i],
    m = m,
    ATE = SD$ATE2_t[i],
    ft = ft_t,
    taut = taut_t,
    theta = SD$ft_t_theta2[i],
    gt = gt_t,
    alpha = alpha_t, 
    beta1 = beta1_t
  )
  
  beta2_t <- beta1_t +beta2_tstar 
  SD[i, c("beta2_t0", "beta2_t1", "beta2_t2")[1:SD$p[i]]] <- beta2_t
  SD$ATE2_t_empirical[i] <- compute_ATE(gt=gt_t, alpha=alpha_t, ft = ft_t, beta = beta2_t, taut= taut_t)
  stopifnot(all.equal(SD$ATE2_t[i], SD$ATE2_t_empirical[i]))
  
  
  ## taut_w
  taut_w <- construct_taut_theta(m, SD$AvgTau_w[i], SD$tau_w_shape[i], theta = SD$tau_w_theta[i])
  
  ## gt_w, alpha_w
  gt_w <- construct_ftgt(m, SD$gt_w_shape[i])
  SD$q_w[i] <- ncol(gt_w)
  alpha_w <- solve_alpha(SD$gt_w_shape[i], m, SD$ASPNC_w[i], gt_w, taut_w,
                         theta = SD$gt_w_theta[i])
  SD[i, c("alpha_w0", "alpha_w1", "alpha_w2")[1:SD$q_w[i]]] <- alpha_w
  SD$ASPNC_w_empirical[i] <- compute_ASPNC(gt_w, alpha_w, taut_w)
  stopifnot(all.equal(SD$ASPNC_w[i], SD$ASPNC_w_empirical[i]))
  
  ## ft_w
  ft_w <- construct_ftgt(m, SD$ft_w_shape[i])
  SD$p_w[i] <- ncol(ft_w)
  # trt 1 eff
  beta1_w <- solve_beta(
    shape = SD$ft_w_shape[i],
    m = m,
    ATE = SD$ATE1_w[i],
    ft = ft_w,
    taut = taut_w,
    theta = SD$ft_w_theta[i],
    gt = gt_w,
    alpha = alpha_w
  )
  SD[i, c("beta1_w0", "beta1_w1", "beta1_w2")[1:SD$p_w[i]]] <- beta1_w
  SD$ATE1_w_empirical[i] <- compute_ATE(gt=gt_w, alpha=alpha_w, ft = ft_w, beta = beta1_w, taut= taut_w)
  stopifnot(all.equal(SD$ATE1_w[i], SD$ATE1_w_empirical[i]))
  
  # trt 2 effect
  beta2_wstar <- solve_beta_star(
    shape = SD$ft_w_shape[i],
    m = m,
    ATE = SD$ATE2_w[i],
    ft = ft_w,
    taut = taut_w,
    theta = SD$ft_w_theta2[i],
    gt = gt_w,
    alpha = alpha_w, 
    beta1 = beta1_w
  )
  
  beta2_w <- beta1_w +beta2_wstar 
  SD[i, c("beta2_w0", "beta2_w1", "beta2_w2")[1:SD$p_w[i]]] <- beta2_w
  SD$ATE2_w_empirical[i] <- compute_ATE(gt=gt_w, alpha=alpha_w, ft = ft_w, beta = beta2_w, taut= taut_w)
  stopifnot(all.equal(SD$ATE2_w[i], SD$ATE2_w_empirical[i]))
  
  ## sample size n
  tryCatch({
    beta_w <- c(beta1_w, beta2_w)
    L <- construct_L(SD$ft_w_shape[i])
    SD$n[i] <- mrt_mult_trt_ss(
      f_t = ft_w,
      g_t = gt_w,
      beta = beta_w,
      alpha = alpha_w,
      p_t = pt,
      gamma = gamma,
      b = b,
      L = L,
      tau = taut_w
    )
    SD$p_star[i] <- dim(L)[1]
  },
  error = function(cond){
    message("iteration ", i, ": ", cond)
    return(NA)
  })
  
  if (!is.na(SD$n[i])) {
    tryCatch({
      ## generate a data set to make sure there are no probabilities outside [0,1]
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
        h_At = NULL,
        j_t = NULL
      )
    },  error = function(cond){
      message("iteration ", i, ": ", cond)
      # need "<<" instead of "<"
      SD$error_indicator[i] <<- 1
    })
  }
}

# extract the SD that being filled up
SD <- SD[setting_start:setting_end,]

# get rid of NA settings
SD <- SD[which(SD$error_indicator == 0),]
SD <- SD[which(!is.na(SD$n)),]


summary(SD)
dim(SD)

outName <- paste("results_sim_5_setting_",isetting,"_",nsim,".RDS",sep = "")
saveRDS(SD,file=outName)

