
### dgm with endogenous availability process ###


#source("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/R_code/functions/utillity.R")

dgm_endoAvail <- function(sample_size, total_T, ft, beta_1, beta_2, gt, alpha,
                          gamma3, gamma4, gamma5,
                          tau, pt) {
  # variable A represent the treatment level  (0 , 1 , 2)
  # variable A1 and A2 is indicator function whether A is equal 1 or 2 respectively
  # variable I is a binary indicating availability of the participant at timepoint t 
  # prob_t indicates the probability of someone available at time point t, assume constant for now
  # tau is E(I_t) which we assume is equal 1 (always available at any given time) unless specifiedd otherwise
  # if there is no exp_var provided, assume error is Normal 0,1
  # AA = 0.3
  
  
  
  # treatment name(only for 2 level of treatment for now)
  a_val <- c(0, 1, 2)
  
  df_names <- c("userid", "time", "Y", "A", "A1", "A2", "I", "prob_A","prob_I", "Y_prob1")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  
  for (t in 1:total_T) {
    prob_a <- pt[t,]

    if (t ==1) {
      # generate Availability 
      A1_lag1 <- rep(0, sample_size)
      A2_lag1 <- rep(0, sample_size)
      Y_lag1 <- rep(0, sample_size)
    }
      # row index for the rows corresponding to day t for every subject
      row_index <- seq(from = t, length.out = sample_size, by = total_T)
      
      # generate Availability 
      dta$prob_I[row_index] <- 0.5 + gamma3 * A1_lag1 + gamma4 * A2_lag1 + gamma5 * Y_lag1
      dta$I[row_index] <- rbinom(sample_size, 1, prob = dta$prob_I[row_index])
      
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

      dta$Y_prob1[row_index] <- exp(dta$A1[row_index] * as.numeric(ft[t, ] %*% beta_1)  +
                                      dta$A2[row_index] * as.numeric(ft[t, ] %*% beta_2) +
                                      as.numeric(gt[t, ] %*% alpha))
      
      if (!(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))) {
        message("min(dta$Y_prob1) = ", min(dta$Y_prob1), 
                "; max(dta$Y_prob1) = ", max(dta$Y_prob1))
      }
      stopifnot(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))
      dta$Y[row_index] <- rbinom(sample_size, 1, prob = dta$Y_prob1[row_index])
      
      row_index_lag1 <- row_index
      A1_lag1 <- dta$A1[row_index_lag1]
      A2_lag1 <- dta$A2[row_index_lag1]
      Y_lag1 <- dta$Y[row_index_lag1]
    }
  return(dta)
}


compute_tau_endoAvail <- function(gt, alpha, ft, beta, pt, gamma3, gamma4, gamma5) {
  # function to calculate the mean availability for tau_w 
  m <- nrow(gt)
  tau <- rep(NA, m)
  EYt <- rep(NA, m)
  
  tau[1] <- 0.5
  EYt[1] <- 0 #not going to calculate this, set to 0 for initialization purposes
  for (t in 2:m) {
    sum_pefbeta <- compute_sum_pefbeta(pt = pt[t-1, ], this_f_t = ft_t[t-1, ], beta = beta)
    EYt[t] <- (pt[t-1,1] + sum_pefbeta) * exp(as.numeric(gt[t-1,] %*% alpha))
    tau[t] <- 0.5 + gamma3 * pt[t-1, 2] * tau[t-1] + gamma4 * pt[t-1, 3] * tau[t-1]+
      gamma5 * EYt[t] * tau[t-1]
  }

  return(tau)
}


# # for debugging puposes
# ## m
# m <- 10
# 
# ## pt
# pt <- construct_pt(m, "constant equal pk")
# 
# ## taut_t
# taut_t <-
#   construct_taut_theta(m , 0.5, "constant")
# 
# ## ft_t, beta1_t, beta2_t
# ft_t <- construct_ftgt(m, "constant")
# # trt 1 effect
# beta1_t <- log(1.2) #only for constant ft see goodnotes on simulation setting for detail
# # trt 2 effect
# beta2_t <- log(1.6)
# 
# beta_t <- c(beta1_t, beta2_t)
# 
# ## gt_t, alpha_t
# gt_t <- construct_ftgt(m, "quadratic_theta")
# alpha_t <- solve_alpha("quadratic_theta", m, 0.2, gt_t, taut_t,
#                        theta = 0.3)
# 
# 
# ## taut_w
# # (not used in the following)
# 
# ## n, p_w, q_w
# n <- 15
# sample_size = n
# total_T = m
# ft = ft_t
# beta_1 = beta1_t
# beta_2 = beta2_t
# gt = gt_t
# alpha = alpha_t
# tau = taut_t
# pt = pt
# 
# 
# 
# 
# dta <- dgm_endoAvail(
#   sample_size = 10^6,
#   total_T = m,
#   ft = ft_t,
#   beta_1 = beta1_t,
#   beta_2 = beta2_t,
#   gt = gt_t,
#   alpha = alpha_t,
#   tau = taut_t,
#   pt = pt,
#   gamma3  = 0,
#   gamma4 = 0, 
#   gamma5 = 0.1
# )
# 
# # debug tau 
# A1 <- A2 <- EYt  <- mean_Avail <- mean_tau <- rep(NA, total_T)
# 
# for (time in 1:total_T) {
#   mean_tau[time] <- mean(dta[which(dta$time == time),"prob_I"])
#   mean_Avail[time] <- mean(dta[which(dta$time == time),"I"])
#   A1[time] <- mean(dta[which(dta$time == time),"A1"])
#   A2[time] <- mean(dta[which(dta$time == time),"A2"])
#   EYt[time] <- mean(dta[which(dta$time == time),"Y_prob1"])
# }
# 
# # empirical values from the dataset
# A1 #(E(I_t * A_t), since A_t depends on availability at time t)
# mean_tau
# mean_Avail
# Yt
# 0.5 + 0.1 * A1 + 0.2 * A2
# tau_calc <- compute_tau_endoAvail(gt = gt_t, alpha = alpha_t, ft = ft_t, beta = beta_t, pt = pt, gamma3 = 0, gamma4 = 0, gamma5 = 0.1)
## mean_tau should be similar to tau_calc if the calculation is correct
## in the case of gamma, i.e endo in Y_t, it differ by value of 0.01 ish but the trend is correct. is this just computation error? 
