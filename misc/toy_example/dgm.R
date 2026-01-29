# this R script is used to generate synthetic data for toy example
# the true control variable is sex and time
# treatment effect only depended on time

dgm <- function(sample_size, total_T, tau, pt) {
  # variable A represent the treatment level  (0, 1, 2, 3)
  # A = 0 indicate control arm 
  # variable A1, A2, A3 is a binary indicator function whether treatment is equal 1, 2, 3 respectively
  # variable I is a binary indicator for availability of the participant at timepoint t 
  # tau is E(I_t) which we assume is equal 1 (always available at any given time) unless specified otherwise
  # pt is a matrix of total_T by n_A which indicates the probability of assigning each treatment (column) at any given time(row) 
  
  # treatment effect dependend on time as moderator
  beta_1 <- c(-0.1, -0.0005)
  beta_2 <- c(-0.3, -0.001)
  beta_3 <- c(-0.2, -0.0015)
  
  # control effect on both time and sex
  alpha <- c(-0.2, -0.0001, 0.0005)
  
  # treatment name(only for 3 level of treatment for now)
  a_val <- c(0, 1, 2, 3)
  
  df_names <- c("userid", "time", "sex", "Y", "A", "A1", "A2", "A3", "I", "prob_A","prob_I", "Y_prob1")
  
  dta <- data.frame(matrix(NA, nrow = sample_size * total_T, ncol = length(df_names)))
  names(dta) <- df_names
  
  dta$userid <- rep(1:sample_size, each = total_T)
  dta$time <- rep(1:total_T, times = sample_size)
  
  # generate sex for each individual 
  dta$sex <- rep(rbinom(sample_size, 1, prob = 0.5), each = total_T)
  
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
                                    ifelse(dta$A[row_index] == 1, prob_a[2],
                                           ifelse(dta$A[row_index] == 2, prob_a[3], prob_a[4])))
    # create pseudo variable for indicator function
    dta$A1[row_index] <- ifelse(dta$A[row_index] == 1, 1, 0)
    dta$A2[row_index] <- ifelse(dta$A[row_index] == 2, 1, 0)
    dta$A3[row_index] <- ifelse(dta$A[row_index] == 3, 1, 0)

    dta$Y_prob1[row_index] <- exp(dta$A1[row_index] * (beta_1[1] + beta_1[2] * dta$time[row_index])  +
                                    dta$A2[row_index] * (beta_2[1] + beta_2[2] * dta$time[row_index]) +
                                    dta$A3[row_index] * (beta_3[1] + beta_3[2] * dta$time[row_index]) +
                                    alpha[1] + alpha[2] * dta$time[row_index] + alpha[3] * dta$sex[row_index])
    if (!(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))) {
      message("min(dta$Y_prob1) = ", min(dta$Y_prob1), 
              "; max(dta$Y_prob1) = ", max(dta$Y_prob1))
    }
    stopifnot(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))
    dta$Y[row_index] <- rbinom(sample_size, 1, prob = dta$Y_prob1[row_index])
  }
  
  return(dta)
}


# for debugging puposes
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("functions_util.R")
## m
total_T <- 15
sample_size = 50

## pt probability assigning treatment are all equal
pt <- matrix(rep(c(0.25, 0.25, 0.25, 0.25),times = total_T ), ncol = 4, byrow = TRUE)

## taut_t
tau <-
  construct_taut_theta(total_T , 1, "constant")


dta<- dgm(
  sample_size = sample_size,
  total_T = total_T,
  tau = tau,
  pt = pt
)

saveRDS(dta, "toy_data.RDS")
