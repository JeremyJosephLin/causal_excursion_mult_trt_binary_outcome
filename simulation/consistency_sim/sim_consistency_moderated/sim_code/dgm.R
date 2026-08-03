# Rscript for consistency simulation using Marginal Model 
# E(Y_t) depends on both Treatment assignment and Moderator variable Z_t
# 2026/07/09 updated the interaction effect between two treatment to be different
# Use for comparing wcls , GEE, and GLMM 
# See GoodNotes for more details

dgm_1_nonflat <- function(sample_size, total_T) {
  # same DGM as dgm_1 above, but with p_t(A_t | H_t) made history-dependent
  # (specifically, dependent on the individual's own realized moderator S_t)
  # instead of flat / common-across-individuals at each time point.
  #
  # variable A represent the treatment level  (0 , 1 , 2)
  # variable A1 and A2 is indicator function whether A is equal 1 or 2 respectively
  # variable I is a binary indicating availability of the participant at timepoint t
  # prob_t indicates the probability of someone available at time point t, assume constant for now
  # S is the moderator variable, still takes three values with equal probability
  #
  # CHANGE FROM dgm_1: previously prob_A depended only on t (odd/even), which is
  # common information shared across all individuals at that time point -- this
  # is still a "flat" design in the sense that p_t(A_t | H_t) = p_t(A_t) does not
  # vary with any individual-level history. Here, prob_A depends on each
  # individual's own S_t, so p_t(A_t | H_t) is genuinely history-dependent and
  # not flat. The target/reference distribution used to define the causal
  # estimand (e.g. tilde p_t(A_t | S_t) = 1/3 for WCLS reweighting) should be
  # specified separately downstream -- it is NOT the same as the randomization
  # probability generated here.
  
  # Coefficient on indicator Zt
  baseline_Y_S0 <- 0.2
  baseline_Y_S1 <- 0.5
  baseline_Y_S2 <- 0.4
  
  # intercept
  beta_1 <- 0.1 # when At = 1
  beta_2 <- 0.25  # when At = 2
  
  # Moderator effect
  beta_3 <- 0.3
  beta_4 <- 0.2
  
  # Randomization probabilities p_t(A_t | S_t), one row per value of S = 0, 1, 2
  # columns correspond to a_val = c(0, 1, 2). Rows are deliberately non-flat
  # and distinct across S so that P(A_t | H_t) varies with the individual's
  # own history (here, S_t).
  prob_a_given_S <- rbind(
    c(0.5, 0.3, 0.2),  # S_t = 0
    c(0.2, 0.5, 0.3),  # S_t = 1
    c(0.3, 0.2, 0.5)   # S_t = 2
  )
  
  prob_I <- 1
  a_val <- c(0, 1, 2)
  
  df_names <- c("userid", "time", "Y", "A", "A1", "A2", "I", "S", "prob_A", "prob_I", "Y_A0", "Y_prob1")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  
  for (t in 1:total_T) {
    # row index for the rows corresponding to day t for every subject
    row_index <- seq(from = t, by = total_T, length = sample_size)
    
    # generate moderator variable FIRST: treatment assignment probability now
    # depends on S_t, so S_t must be realized before A_t is drawn
    dta$S[row_index] <- sample(c(0, 1, 2), sample_size, replace = TRUE)
    
    # generate treatment assignment conditional on S_t
    for (s_val in c(0, 1, 2)) {
      s_index <- row_index[dta$S[row_index] == s_val]
      if (length(s_index) == 0) next
      prob_a <- prob_a_given_S[s_val + 1, ]  # +1 for 1-based row indexing
      dta$A[s_index] <- sample(x = a_val, size = length(s_index), replace = TRUE, prob = prob_a)
      dta$prob_A[s_index] <- ifelse(dta$A[s_index] == 0, prob_a[1],
                                    ifelse(dta$A[s_index] == 1, prob_a[2], prob_a[3]))
    }
    
    # create pseudo variable for indicator function
    dta$A1[row_index] <- ifelse(dta$A[row_index] == 1, 1, 0)
    dta$A2[row_index] <- ifelse(dta$A[row_index] == 2, 1, 0)
    # generate Availability
    dta$I[row_index] <- rbinom(sample_size, 1, prob = prob_I)
    dta$prob_I[row_index] <- ifelse(dta$I[row_index] == 0, 1 - prob_I, prob_I)
    
    dta$Y_A0[row_index] <- ifelse(dta$S[row_index] == 0, baseline_Y_S0,
                                  ifelse(dta$S[row_index] == 1, baseline_Y_S1, baseline_Y_S2))
    
    dta$Y_prob1[row_index] <- dta$Y_A0[row_index] * (exp(dta$A1[row_index] * (beta_1 + beta_3 * dta$S[row_index]) +
                                                           dta$A2[row_index] * (beta_2 + beta_4 * dta$S[row_index])))
    
    if (!(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))) {
      message("min(dta$Y_prob1) = ", min(dta$Y_prob1),
              "; max(dta$Y_prob1) = ", max(dta$Y_prob1))
    }
    stopifnot(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))
    dta$Y[row_index] <- rbinom(sample_size, 1, prob = dta$Y_prob1[row_index])
  }
  
  return(dta)
}
# # For debugging purposes
# sample_size = 30
# total_T = 10
# 
# dgm_1(sample_size = sample_size, total_T = total_T)
  