# dgm function 

source("~/Documents/Research/causal-excursion-multi-treatment-binary-outcome/R_code/functions/utillity.R")

dgm_SerialCorrelationOnY <- function(sample_size, total_T, ft, beta_1, beta_2, gt, alpha, tau, pt, coef_SC) {
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

    if (t == 1) {
      Y_lag1 <- rep(0, sample_size)
    }
    # generate probability Y = 1
    dta$Y_prob1[row_index] <- exp(dta$A1[row_index] * as.numeric(ft[t, ] %*% beta_1)  +
                                    dta$A2[row_index] * as.numeric(ft[t, ] %*% beta_2) +
                                    as.numeric(gt[t, ] %*% alpha)+
                                    Y_lag1 * coef_SC)

    
    #check probability to be in the range of 0,1
    if (!(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))) {
      message("min(dta$Y_prob1) = ", min(dta$Y_prob1), 
              "; max(dta$Y_prob1) = ", max(dta$Y_prob1))
    }
    stopifnot(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))

    # generate binomial Y
    dta$Y[row_index] <- rbinom(sample_size, 1, prob = dta$Y_prob1[row_index])
    # update lag effect 
    Y_lag1 <- dta$Y[row_index]
  }
  
  return(dta)
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
# dgm_SerialCorrelationOnY(
#   sample_size = n,
#   total_T = m,
#   ft = ft_t,
#   beta_1 = beta1_t,
#   beta_2 = beta2_t,
#   gt = gt_t,
#   alpha = alpha_t,
#   tau = taut_t,
#   pt = pt,
#   coef_SC = 0.1
# )
# 
