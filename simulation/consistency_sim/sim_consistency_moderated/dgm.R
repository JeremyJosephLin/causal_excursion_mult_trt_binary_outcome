# Rscript for consistency simulation using Marginal Model 
# E(Y_t) depends on both Treatment assignment and Moderator variable Z_t
# Use for comparing wcls , GEE, and GLMM 
# See GoodNotes for more details

dgm_1 <- function(sample_size, total_T) {
  # same DGM as dgm_binary above, but faster
  # variable A represent the treatment level  (0 , 1 , 2)
  # variable A1 and A2 is indicator function whether A is equal 1 or 2 respectively
  # variable I is a binary indicating availability of the participant at timepoint t 
  # prob_t indicates the probability of someone available at time point t, assume constant for now
  # Z_t is moderator variable, that takes three values with equal probability
  
  # Coefficient on indicator Zt
  baseline_Y_S0 <- 0.2
  baseline_Y_S1 <- 0.5
  baseline_Y_S2 <- 0.4
  
  
  # intercept
  beta_1 <- 0.1 # when At = 1 
  beta_2 <- 0.25  # when At = 2
  
  # Moderator effect 
  beta_3 <- 0.3

  
  
  prob_a <- c(0.2, 0.5, 0.3)
  prob_I <- 0.6
  a_val <- c(0, 1, 2)
  
  
  df_names <- c("userid", "time", "Y", "A", "A1", "A2", "I", "S", "prob_A","prob_I","Y_A0", "Y_prob1")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    # generate treatment assignment
    dta$A[row_index] <- sample(x = a_val, size = sample_size, replace = TRUE, prob = prob_a)
    dta$prob_A[row_index] <- ifelse(dta$A[row_index] == 0, prob_a[1], 
                                    ifelse(dta$A[row_index] == 1, prob_a[2], prob_a[3]))
    # create pseudo variable for indicator function
    dta$A1[row_index] <- ifelse(dta$A[row_index] == 1, 1, 0)
    dta$A2[row_index] <- ifelse(dta$A[row_index] == 2, 1, 0)
    # generate Availability 
    dta$I[row_index] <- rbinom(sample_size, 1, prob = prob_I)
    dta$prob_I[row_index] <- ifelse(dta$I[row_index] == 0, 1 - prob_I, prob_I)
    # generatw moderator variable
    dta$S[row_index] <- sample(c(0,1,2), sample_size, replace = TRUE)
    dta$Y_A0[row_index] <- ifelse(dta$S[row_index] == 0, baseline_Y_S0, 
                                  ifelse(dta$S[row_index] == 1, baseline_Y_S1, baseline_Y_S2))
    
    dta$Y_prob1[row_index] <- dta$Y_A0[row_index] * (exp(dta$A1[row_index] * (beta_1)  +
                                    dta$A2[row_index] * (beta_2) + beta_3 * dta$S[row_index]))
    
    if (!(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))) {
      message("min(dta$Y_prob1) = ", min(dta$Y_prob1), 
              "; max(dta$Y_prob1) = ", max(dta$Y_prob1))
    }
    stopifnot(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))
    dta$Y[row_index] <- rbinom(sample_size, 1, prob = dta$Y_prob1[row_index])
  }
  
  return(dta)
}

# For debugging purposes 
sample_size = 30 
total_T = 10
  
  