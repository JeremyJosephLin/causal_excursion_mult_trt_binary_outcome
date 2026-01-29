# This code generalize the simulation for different treatment effect 
# dgm will be able to 
rm(list = ls())
source("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/R_code/wcls_cat_trt_binary_outcome.R")
source("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/R_code/functions/utillity.R")
source("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/R_code/ss_calc.R")
source("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/R_code/functions/dgm.R")

# 1. Simulation design setup ----------------------------------------------

# SD: simulation design
# SD: simulation design
SD <- expand.grid(gamma = 0.05,  # type I error
                  b = 0.2, # type II error
                  # m = c(30, 60),
                  m = c(30),
                  pt_shape = c("constant equal pk"),
                  # pt_shape = c("constant 0.5"),
                  ft_t_shape = c("constant"),
                  ft_t_theta = 0,
                  gt_t_shape = c("quadratic_theta", "linear_theta"),
                  gt_t_theta = 0.3, 
                  # for specifying gt when gt_t_shape is linear or quadratic
                  # ATE_t = c(1.2, 1.3, 1.4, 1.5),
                  ATE1_t = c(1.2, 1.3),
                  ATE1_t_empirical = NA,
                  ATE2_t = c(1.6, 1.8),
                  ATE2_t_empirical = NA,
                  ASPNC_t = c(0.2, 0.3),
                  ASPNC_t_empirical = NA,
                  #AvgTau_t = c(0.5, 0.8, 1),
                  AvgTau_t = c(1),
                  tau_t_shape = c("constant"),
                  beta1_t0 = NA,
                  beta1_t1 = NA,
                  beta1_t2 = NA,
                  beta2_t0 = NA,
                  beta2_t1 = NA,
                  beta2_t2 = NA,
                  alpha_t0 = NA,
                  alpha_t1 = NA,
                  alpha_t2 = NA,
                  p = NA,
                  q = NA,
                  
                  ft_w_shape = c("linear_theta" ),
                  ft_w_theta = seq(from = -1, to = 1, by = 0.1),
                  gt_w_shape = NA,
                  gt_w_theta = NA,
                  beta1_w0 = NA,
                  beta1_w1 = NA,
                  beta1_w2 = NA,
                  beta2_w0 = NA,
                  beta2_w1 = NA,
                  beta2_w2 = NA,
                  alpha_w0 = NA,
                  alpha_w1 = NA,
                  alpha_w2 = NA,
                  ATE1_w = NA,
                  ATE1_w_empirical = NA,
                  ATE2_w = NA,
                  ATE2_w_empirical = NA,
                  ASPNC_w = NA,
                  ASPNC_w_empirical = NA,
                  AvgTau_w = NA,
                  tau_w_shape = NA,
                  p_w = NA,
                  q_w = NA,
                  n = NA,
                  p_star = NA,
                  error_indicator = 0,
                  stringsAsFactors = FALSE)

# remove simulation settings where pt*ft is not in the linear span of gt
SD <- SD[!(SD$ft_t_shape == "linear_theta" &
             SD$gt_t_shape == "constant"), ]
SD <- SD[!(SD$ft_t_shape == "quadratic_theta" &
             SD$gt_t_shape != "quadratic_theta"), ]

SD <- SD[!(SD$ft_t_shape == "constant" &
             SD$ft_t_theta != 0), ]
SD <- SD[!(SD$gt_t_shape == "constant" &
             SD$gt_t_theta != 0), ]

SD <- SD[!(SD$AvgTau_t == 1 &
             SD$tau_t_shape != "constant"), ]

## add stuff that are known by simulation design
# taut
SD$AvgTau_w <- SD$AvgTau_t
SD$tau_w_shape <- SD$tau_t_shape
# gt
SD$gt_w_shape <- SD$gt_t_shape
SD$gt_w_theta <- SD$gt_t_theta
SD$ASPNC_w <- SD$ASPNC_t
# ft
#SD$ft_w_shape <- SD$ft_t_shape
#SD$ft_w_theta <- SD$ft_t_theta
SD$ATE1_w <- SD$ATE1_t
SD$ATE2_w <- SD$ATE2_t


##### fill the SD #####
options(warn = 2)
for (i in 1:nrow(SD)) {
  
  ## gamma, b, m
  gamma <- SD$gamma[i]
  b <- SD$b[i]
  m <- SD$m[i]
  
  ## pt
  pt <- construct_pt(m, SD$pt_shape[i])
  
  ## taut_t
  taut_t <- construct_taut_theta(m , SD$AvgTau_t[i], SD$tau_t_shape[i])
  
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
  beta2_t <- solve_beta(
    shape = SD$ft_t_shape[i],
    m = m,
    ATE = SD$ATE2_t[i],
    ft = ft_t,
    taut = taut_t,
    theta = SD$ft_t_theta[i],
    gt = gt_t,
    alpha = alpha_t
  )
  SD[i, c("beta2_t0", "beta2_t1", "beta2_t2")[1:SD$p[i]]] <- beta2_t
  SD$ATE2_t_empirical[i] <- compute_ATE(gt=gt_t, alpha=alpha_t, ft = ft_t, beta = beta2_t, taut= taut_t)
  stopifnot(all.equal(SD$ATE2_t[i], SD$ATE2_t_empirical[i]))
  
  ## taut_w
  taut_w <- construct_taut_theta(m, SD$AvgTau_w[i], SD$tau_w_shape[i])
  
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
  
  # trt 2 eff
  beta2_w <- solve_beta(
    shape = SD$ft_w_shape[i],
    m = m,
    ATE = SD$ATE2_w[i],
    ft = ft_w,
    taut = taut_w,
    theta = SD$ft_w_theta[i],
    gt = gt_w,
    alpha = alpha_w
  )
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
      tau = taut_t
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

SD$error_indicator
# get rid of SD that has error
SD<- SD[which(SD$error_indicator ==0),]
saveRDS(SD, file = "~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/simulations/power_simulation /2. WA-a violated/sim_2_3/sim_2_3a/factorial_design.RDS")

dim(SD)

