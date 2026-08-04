# Generative model for the consistency simulation of Supplementary Section C.
#
# One DGM function serves both randomization scenarios, so that the generative model
# displayed in the supplement describes both tables. It consolidates the two divergent
# copies used by the earlier per-scenario pipelines:
#   sim_consistency_moderated/sim_code/dgm.R                       (dgm_1)
#   sim_consistency_marginal_model_with control_non_const_pt/...   (dgm_1_nonflat)
# neither of which matches the generative model printed in the supplement.
#
# E(Y_t | H_t, A_t) = b(Z_t) * exp{ 1(A_t=1)(beta_10 + beta_11 Z_t)
#                                 + 1(A_t=2)(beta_20 + beta_21 Z_t) }
# with b(0) = 0.2, b(1) = 0.5, b(2) = 0.4 and Z_t ~ Unif{0, 1, 2} i.i.d.

dgm_consistency <- function(sample_size,
                            total_T,
                            scenario = c("uniform", "Z_dependent"),
                            beta_10 = 0.1,
                            beta_11 = 0.3,
                            beta_20 = 0.25,
                            beta_21 = 0.3) {
  scenario <- match.arg(scenario)

  baseline_Y <- c(0.2, 0.5, 0.4) # indexed by Z_t = 0, 1, 2

  # Scenario 1: p_t(A_t | H_t) = (1/3, 1/3, 1/3), free of H_t.
  # Scenario 2: p_t(A_t | H_t) depends on the individual's own Z_t.
  prob_a_given_Z <- if (scenario == "uniform") {
    rbind(rep(1 / 3, 3), rep(1 / 3, 3), rep(1 / 3, 3))
  } else {
    rbind(
      c(0.5, 0.3, 0.2), # Z_t = 0
      c(0.2, 0.5, 0.3), # Z_t = 1
      c(0.3, 0.2, 0.5)  # Z_t = 2
    )
  }

  a_val <- c(0, 1, 2)
  n_row <- sample_size * total_T

  dta <- data.frame(
    userid  = rep(1:sample_size, each = total_T),
    time    = rep(1:total_T, times = sample_size),
    Y       = NA_real_,
    A       = NA_real_,
    A1      = NA_real_,
    A2      = NA_real_,
    I       = 1,
    S       = NA_real_,
    prob_A  = NA_real_,
    prob_I  = 1,
    Y_A0    = NA_real_,
    Y_prob1 = NA_real_
  )

  for (t in 1:total_T) {
    row_index <- seq(from = t, by = total_T, length = sample_size)

    # Z_t is drawn before A_t: in Scenario 2 the randomization depends on it.
    dta$S[row_index] <- sample(c(0, 1, 2), sample_size, replace = TRUE)

    for (z_val in c(0, 1, 2)) {
      z_index <- row_index[dta$S[row_index] == z_val]
      if (length(z_index) == 0) next
      prob_a <- prob_a_given_Z[z_val + 1, ]
      dta$A[z_index] <- sample(x = a_val, size = length(z_index), replace = TRUE, prob = prob_a)
      dta$prob_A[z_index] <- prob_a[dta$A[z_index] + 1]
    }

    dta$A1[row_index] <- as.numeric(dta$A[row_index] == 1)
    dta$A2[row_index] <- as.numeric(dta$A[row_index] == 2)

    dta$Y_A0[row_index] <- baseline_Y[dta$S[row_index] + 1]
    dta$Y_prob1[row_index] <- dta$Y_A0[row_index] *
      exp(dta$A1[row_index] * (beta_10 + beta_11 * dta$S[row_index]) +
            dta$A2[row_index] * (beta_20 + beta_21 * dta$S[row_index]))

    stopifnot(all(dta$Y_prob1[row_index] >= 0 & dta$Y_prob1[row_index] <= 1))
    dta$Y[row_index] <- rbinom(sample_size, 1, prob = dta$Y_prob1[row_index])
  }

  dta
}

# True CEE values under the DGM above.
#   moderated (S_t = Z_t):  beta_k = (beta_k0, beta_k1)
#   marginal  (S_t = empty): beta_k = log{ E b(Z) e^{beta_k0 + beta_k1 Z} / E b(Z) }
true_cee <- function(beta_10 = 0.1, beta_11 = 0.3, beta_20 = 0.25, beta_21 = 0.3) {
  b <- c(0.2, 0.5, 0.4)
  z <- c(0, 1, 2)
  list(
    moderated = list(beta1 = c(beta_10, beta_11), beta2 = c(beta_20, beta_21)),
    marginal = c(
      beta1 = log(sum(b * exp(beta_10 + beta_11 * z)) / sum(b)),
      beta2 = log(sum(b * exp(beta_20 + beta_21 * z)) / sum(b))
    )
  )
}
