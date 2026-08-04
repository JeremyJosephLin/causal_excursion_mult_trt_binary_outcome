# reproduce_table2.R
#
# Reproduces the Drink Less data analysis of the paper "Micro-randomized Trials with
# Categorical Treatments and Binary Proximal Outcome": Table 2, the effect statements in
# Section 5.1, and the inputs to the sample size illustration in Section 5.2.
#
# It starts from the public OSF file (see prepare_drinkless_data.R for how to obtain it)
# and needs nothing else beyond this repository.
#
#   Rscript "reproduce_table2.R"
#
# The three analyses, as in equation (16) of the paper:
#   (a) the marginal CEE,
#   (b) the CEE moderated by the decision point index, with t starting at 0,
#   (c) the CEE moderated by the baseline AUDIT risk category.
# All three adjust for the same control variables in the nuisance model g_t: the decision
# point index, age, gender, employment type, AUDIT score, and the lagged proximal outcome.
# Analysis (c) replaces the AUDIT score by the category indicators that define its
# moderator. Every standard error printed is the small-sample adjusted standard error of
# Mancl and DeRouen (2001), and the confidence intervals and p-values are the matching
# t-based quantities, which is what Table 2 reports.

suppressPackageStartupMessages(library(dplyr))

args     <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
HERE     <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else getwd()
source(file.path(HERE, "prepare_drinkless_data.R"))
ROOT <- drinkless_repo_root()
source(file.path(ROOT, "functions", "emee_catA.R"))
source(file.path(ROOT, "functions", "functions_util.R"))
source(file.path(ROOT, "functions", "ss_calc.R"))
source(file.path(ROOT, "functions", "utillity.R"))

options(digits = 10)

dta <- prepare_drinkless_data()

cat("=== DATA ===\n")
cat("n participants          :", length(unique(dta$ID)), "\n")
cat("person-decision points  :", nrow(dta), "\n")
cat("decision index range    :", paste(range(dta$decision_index), collapse = " to "), "\n")
cat("mean(Y)                 :", mean(dta$binary_outcome), "\n")
cat("mean(Y | A = 0)         :", mean(dta$binary_outcome[dta$treatment_cate == 0]), "\n\n")

CONTROL   <- c("decision_index", "age", "gender", "employment_1", "employment_2",
               "AUDIT_score", "binary_outcome.lag")
CONTROL_C <- c("decision_index", "age", "gender", "employment_1", "employment_2",
               "AUDIT_harmful", "AUDIT_atRisk", "binary_outcome.lag")

## tilde p is set to the trial's randomization probability, so that J_t = 1 throughout,
## matching Section 5.1 of the paper. The estimator's own default is 1/3 per level.
## The matrix is decision points by treatment levels, columns ordered (0, 1, 2).
N_TIME <- length(unique(dta$decision_index))
PTILDE <- matrix(rep(c(0.4, 0.3, 0.3), times = N_TIME), ncol = 3, byrow = TRUE)

fitwrap <- function(moderator, control) {
  emee_catA(
    dta = dta, id_varname = "ID", decision_time_varname = "decision_index",
    treatment_varname = "treatment_cate", outcome_varname = "binary_outcome",
    control_varname = control, moderator_varname = moderator,
    rand_prob_varname = "prob_A", avail_varname = "I",
    trt_level = 2, estimator_initial_value = NULL, pmatrix_tilde = PTILDE)
}

beta_rows <- function(fit) {
  s <- emee_catA_summary(fit)
  data.frame(Estimate = fit$beta_hat,
             SE       = fit$beta_se_adjusted,
             CI_lo    = fit$conf_int_adjusted_t[, 1],
             CI_hi    = fit$conf_int_adjusted_t[, 2],
             p        = s$`Pr(>|t|)`)
}
glh_rows <- function(fit, L, labels) {
  g <- emee_catA_glh(fit, L)
  data.frame(row.names = labels,
             Estimate = g$summary$Estimate,
             SE       = g$summary$`Std. Error`,
             CI_lo    = g$conf_int[, 1],
             CI_hi    = g$conf_int[, 2],
             p        = g$summary$`Pr(>|t|)`)
}
# how the paper reports an effect: "X times (= e^y, 95% CI [lo, hi], p)"
show <- function(r, label) {
  cat(sprintf("%-38s est %8.4f  SE %6.4f  CI [%7.4f, %7.4f]  p %.4f\n",
              label, r$Estimate, r$SE, r$CI_lo, r$CI_hi, r$p))
  cat(sprintf("%-38s ratio %.2f (= e^{%.2f}), 95%% CI [%.2f, %.2f]\n", "",
              exp(r$Estimate), r$Estimate, exp(r$CI_lo), exp(r$CI_hi)))
}
tex <- function(r) sprintf("%.3f & %.3f & (%.3f, %.3f) & %s",
                           r$Estimate, r$SE, r$CI_lo, r$CI_hi,
                           ifelse(r$p < 0.001, "$<0.001$", sprintf("%.3f", r$p)))

## ------------------------------------------------------------------- Model (a)
fit_a <- fitwrap(moderator = c(), control = CONTROL)
ra <- beta_rows(fit_a)
da <- glh_rows(fit_a, matrix(c(-1, 1), nrow = 1), "beta2 - beta1")
cat("=== MODEL (a): marginal CEE ===\n")
show(ra[1, ], "beta_1^marginal");  show(ra[2, ], "beta_2^marginal")
show(da[1, ], "beta_2 - beta_1")
cat("\nTable 2 rows, LaTeX:\n")
cat("  b1  &", tex(ra[1, ]), "\n")
cat("  b2  &", tex(ra[2, ]), "\n")
cat("  d   &", tex(da[1, ]), "\n\n")

## ------------------------------------------------------------------- Model (b)
fit_b <- fitwrap(moderator = "decision_index", control = CONTROL)
rb <- beta_rows(fit_b)
db <- glh_rows(fit_b, matrix(c(-1, 0, 1, 0,
                               0, -1, 0, 1), ncol = 4, byrow = TRUE),
               c("int diff", "slp diff"))
cat("=== MODEL (b): CEE moderated by the decision point index (t starts at 0) ===\n")
show(rb[1, ], "beta_1^t-int");  show(rb[2, ], "beta_1^t-slp")
show(rb[3, ], "beta_2^t-int");  show(rb[4, ], "beta_2^t-slp")
show(db[1, ], "int: beta_2 - beta_1"); show(db[2, ], "slp: beta_2 - beta_1")
cat("\nTable 2 rows, LaTeX:\n")
for (i in 1:4) cat("  ", rownames(rb)[i], " &", tex(rb[i, ]), "\n")
cat("   int diff &", tex(db[1, ]), "\n")
cat("   slp diff &", tex(db[2, ]), "\n\n")

## ------------------------------------------------------------------- Model (c)
fit_c <- fitwrap(moderator = c("AUDIT_harmful", "AUDIT_atRisk"), control = CONTROL_C)
rc <- beta_rows(fit_c)
Lc <- rbind(c( 1,  1,  0,  0, 0, 0),    # trt1 harmful
            c( 1,  0,  1,  0, 0, 0),    # trt1 dependent
            c( 0,  0,  0,  1, 1, 0),    # trt2 harmful
            c( 0,  0,  0,  1, 0, 1),    # trt2 dependent
            c(-1,  0,  0,  1, 0, 0),    # trt2 - trt1, hazardous
            c(-1, -1,  0,  1, 1, 0),    # trt2 - trt1, harmful
            c(-1,  0, -1,  1, 0, 1))    # trt2 - trt1, dependent
labc <- c("trt1 harmful", "trt1 dependent", "trt2 harmful", "trt2 dependent",
          "trt2-trt1 hazardous", "trt2-trt1 harmful", "trt2-trt1 dependent")
dc <- do.call(rbind, lapply(1:7, function(i)
  glh_rows(fit_c, matrix(Lc[i, ], ncol = 6), labc[i])))
cat("=== MODEL (c): CEE moderated by AUDIT category ===\n")
show(rc[1, ], "beta_1^int (hazardous)")
show(dc["trt1 harmful", ],   "beta_1^int + beta_1^harmful")
show(dc["trt1 dependent", ], "beta_1^int + beta_1^dependent")
show(rc[4, ], "beta_2^int (hazardous)")
show(dc["trt2 harmful", ],   "beta_2^int + beta_2^harmful")
show(dc["trt2 dependent", ], "beta_2^int + beta_2^dependent")
show(dc["trt2-trt1 hazardous", ], "trt2 - trt1, hazardous")
show(dc["trt2-trt1 harmful", ],   "trt2 - trt1, harmful")
show(dc["trt2-trt1 dependent", ], "trt2 - trt1, dependent")
cat("\nTable 2 rows, LaTeX:\n")
cat("  b1_int              &", tex(rc[1, ]), "\n")
cat("  b1_int + harmful    &", tex(dc["trt1 harmful", ]), "\n")
cat("  b1_int + dependent  &", tex(dc["trt1 dependent", ]), "\n")
cat("  b2_int              &", tex(rc[4, ]), "\n")
cat("  b2_int + harmful    &", tex(dc["trt2 harmful", ]), "\n")
cat("  b2_int + dependent  &", tex(dc["trt2 dependent", ]), "\n")
cat("  b2_int - b1_int     &", tex(dc["trt2-trt1 hazardous", ]), "\n")
cat("  harmful contrast    &", tex(dc["trt2-trt1 harmful", ]), "\n")
cat("  dependent contrast  &", tex(dc["trt2-trt1 dependent", ]), "\n\n")

## ------------------------------------------- Section 5.2: sample size illustration
cat("=== SECTION 5.2 SAMPLE SIZE ILLUSTRATION ===\n")
ATE1 <- exp(fit_a$beta_hat[1]); ATE2 <- exp(fit_a$beta_hat[2])
ASPN_HAT <- mean(dta$binary_outcome[dta$treatment_cate == 0])
cat(sprintf("  ATE_1 = %.4f (prints as %.2f), ATE_2 = %.4f (prints as %.2f)\n",
            ATE1, round(ATE1, 2), ATE2, round(ATE2, 2)))
cat(sprintf("  ASPN  = mean(Y | A = 0) = %.6f (prints as %.3f)\n",
            ASPN_HAT, round(ASPN_HAT, 3)))
cat("  Note that exp of the nuisance intercept is NOT the ASPN: the estimator centers\n")
cat(sprintf("  the treatment indicators, and here exp(alpha_0) = %.4f.\n",
            exp(fit_a$alpha_hat[1])))

m   <- 30
ft  <- matrix(rep(1, m), ncol = 1)
gt  <- matrix(rep(1, m), ncol = 1)
tau <- rep(1, m)
pt  <- construct_pt(m, "constant p0 high")   # (0.4, 0.3, 0.3), the Drink Less design
ss  <- function(a1, a2, aspn) {
  mrt_mult_trt_ss(f_t = ft, g_t = gt, beta = c(log(a1), log(a2)), alpha = log(aspn),
                  p_t = pt, gamma = 0.05, b = 0.2,
                  L = matrix(c(1, -1), nrow = 1), tau = tau)
}
cat(sprintf("\n  required n at the printed inputs: %d\n",
            ss(round(ATE1, 2), round(ATE2, 2), round(ASPN_HAT, 3))))
cat(sprintf("  required n at full precision    : %d\n", ss(ATE1, ATE2, ASPN_HAT)))
cat("\n  (For reference, the Drink Less MRT itself enrolled 349 participants.)\n")

cat("\n=== sessionInfo ===\n")
cat(R.version.string, "\n")
