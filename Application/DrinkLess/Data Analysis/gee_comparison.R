# Drink Less: the proposed estimator (EMEE) against log-link GEE
# ------------------------------------------------------------------
# Written for the Statistics in Biosciences revision (reviewer 2, comment 6): how do the
# findings compare with what a GEE or mixed-model analysis of the same data would give?
#
# Fits, for the marginal model (a) and the AUDIT-moderated model (c):
#   EMEE  the estimator of the paper, with the small-sample adjusted standard error;
#   GEE-ind   Poisson log link, working independence, robust standard errors
#             (the longitudinal form of Zou's modified Poisson regression);
#   GEE-exch  Poisson log link, exchangeable working correlation, robust standard errors;
#   GLMM      Poisson log link with a participant random intercept, model-based standard
#             errors. Its treatment coefficients are comparable with the others because a
#             random intercept factors out under a log link: E(Y | X, b) = exp(x'beta + b)
#             gives E(Y | X) = exp(x'beta) E(e^b), so contrasts of x'beta are unchanged.
# The GEE labels name the working correlation structure rather than numbering the fits:
# GEE1 and GEE2 are taken in the GEE literature, where they mean first- and second-order
# estimating equations, and both fits here are first-order.
# The GEE mean models carry the same covariates as the EMEE nuisance model g_t plus
# treatment indicators (and, for model (c), treatment-by-AUDIT interactions).
#
# It also checks the equivalence claim made in Appendix C: with f_t = g_t = 1 and
# \tilde p_t set equal to the randomization probability, EMEE, working independence GEE
# on treatment indicators alone, and the arm-wise sample mean ratios all coincide.
#
# The analysis follows DrinkLess_DA.Rmd with two corrections adopted in the revision:
#   (i)  employment type enters the nuisance model, as the manuscript states;
#   (ii) the lagged outcome is lagged within participant rather than across the stacked
#        data (the original crossed participant boundaries on 349 of 10470 rows).
# The decision point index starts at 0, which is the manuscript's convention.
#
# DATA must point at the Drink Less MRT dataset shared on OSF.

DATA <- file.path("/Users/tqian/Dropbox (Personal)/data/DrinkLess MRT (Lauren Bell)",
                  "2022.01 - DrinkLess OSF data sharing/FINAL Dataset_A.rds")
REPO <- "~/repos/cat_trt_binary"

suppressPackageStartupMessages({
  library(dplyr)
  library(geepack)
  library(lme4)
})
source(file.path(REPO, "functions/wcls_cat_trt_binary_outcome.R"))
source(file.path(REPO, "functions/functions_util.R"))

options(digits = 10)

## ------------------------------------------------------------------ preprocessing
raw <- readRDS(DATA)
dta <- raw %>%
  mutate(treatment_cate = ifelse(subversion == "Y", 2, ifelse(subversion == "X", 1, 0)),
         decision_index = days_since_download - 1,   # t starts at 0
         binary_outcome = primary_outcome) %>%
  dplyr::select(ID, decision_index, age, AUDIT_score, gender, employment_type,
                binary_outcome, treatment_cate) %>%
  # the raw file keeps one participant's 30 rows in two non-adjacent blocks, so sort
  # before lagging, and give geeglm a contiguous integer cluster id (it does not
  # terminate in reasonable time when handed the character ID)
  arrange(ID, decision_index) %>%
  group_by(ID) %>%
  mutate(binary_outcome.lag = dplyr::lag(binary_outcome, default = 0)) %>%  # within participant
  ungroup() %>%
  mutate(idnum = as.integer(factor(ID)),
         I = 1,
         prob_A = ifelse(treatment_cate == 0, 0.4, 0.3),
         A1 = as.numeric(treatment_cate == 1),
         A2 = as.numeric(treatment_cate == 2),
         employment_1 = as.numeric(employment_type == 1),
         employment_2 = as.numeric(employment_type == 2),
         AUDIT_harmful = as.numeric(AUDIT_score >= 16 & AUDIT_score <= 19),
         AUDIT_atRisk  = as.numeric(AUDIT_score >= 20 & AUDIT_score <= 40)) %>%
  as.data.frame()

cat("n participants          :", length(unique(dta$ID)), "\n")
cat("person-decision points  :", nrow(dta), "\n")
cat("decision index range    :", paste(range(dta$decision_index), collapse = " to "), "\n")
cat("employment_type counts  :", paste(names(table(dta$employment_type)),
                                       table(dta$employment_type), sep = "=", collapse = "  "), "\n\n")

## control sets: decision point index, age, gender, employment type, AUDIT, lagged outcome
CONTROL   <- c("decision_index", "age", "gender", "employment_1", "employment_2",
               "AUDIT_score", "binary_outcome.lag")
CONTROL_C <- c("decision_index", "age", "gender", "employment_1", "employment_2",
               "AUDIT_harmful", "AUDIT_atRisk", "binary_outcome.lag")

fit_emee <- function(moderator, control, d = dta, ptilde = NULL) {
  wcls_categorical_treatment(
    dta = d, id_varname = "ID", decision_time_varname = "decision_index",
    treatment_varname = "treatment_cate", outcome_varname = "binary_outcome",
    control_varname = control, moderator_varname = moderator,
    rand_prob_varname = "prob_A", avail_varname = "I",
    trt_level = 2, estimator_initial_value = NULL, pmatrix_tilde = ptilde)
}

## linear combination of a fitted coefficient vector, with a standard error
lincom <- function(est, vcov, L, labels) {
  e <- as.vector(L %*% est)
  s <- sqrt(diag(L %*% vcov %*% t(L)))
  data.frame(row.names = labels, Estimate = e, SE = s,
             CI_lo = e - 1.96 * s, CI_hi = e + 1.96 * s,
             p = 2 * pnorm(abs(e / s), lower.tail = FALSE))
}

## GLMM fit returning the treatment-related coefficients and their covariance.
## bobyqa is needed: the default optimizer stops with a max|grad| warning on these data.
fit_glmm <- function(formula, terms_wanted) {
  m <- glmer(formula, data = dta, family = poisson(link = "log"),
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  msgs <- m@optinfo$conv$lme4$messages
  if (length(msgs)) warning("glmer: ", paste(msgs, collapse = "; "))
  b <- fixef(m); V <- as.matrix(vcov(m))
  nm <- vapply(terms_wanted, function(tm) {
    if (tm %in% names(b)) tm else paste(rev(strsplit(tm, ":")[[1]]), collapse = ":")
  }, character(1))
  stopifnot(all(nm %in% names(b)))
  list(est = b[nm], vcov = V[nm, nm, drop = FALSE],
       re_sd = as.data.frame(VarCorr(m))$sdcor, fit = m)
}

## GEE fit returning the treatment-related coefficients and their robust covariance
fit_gee <- function(formula, corstr, terms_wanted) {
  f <- geeglm(formula, data = dta, id = idnum, corstr = corstr, family = poisson(link = "log"))
  s <- summary(f)
  V <- f$geese$vbeta                       # robust sandwich covariance
  stopifnot(max(abs(sqrt(diag(V)) - s$coefficients[, "Std.err"])) < 1e-8)
  idx <- match(terms_wanted, rownames(s$coefficients))
  stopifnot(!any(is.na(idx)))
  list(est = s$coefficients[idx, "Estimate"], vcov = V[idx, idx, drop = FALSE], fit = f)
}

## ============================================================ MODEL (a): marginal CEE
cat("================ MODEL (a): marginal CEE ================\n\n")

emee_a <- fit_emee(moderator = c(), control = CONTROL)
La <- rbind(c(1, 0), c(0, 1), c(-1, 1))
lab_a <- c("beta_1", "beta_2", "beta_2 - beta_1")
nb <- names(emee_a$beta_hat)
out_a <- list(
  EMEE = lincom(emee_a$beta_hat, emee_a$varcov_adjusted[nb, nb], La, lab_a))

form_a <- binary_outcome ~ A1 + A2 + decision_index + age + gender +
  employment_1 + employment_2 + AUDIT_score + binary_outcome.lag
for (cs in c("independence", "exchangeable")) {
  g <- fit_gee(form_a, cs, c("A1", "A2"))
  out_a[[if (cs == "independence") "GEE-ind" else "GEE-exch"]] <- lincom(g$est, g$vcov, La, lab_a)
}
m <- fit_glmm(update(form_a, . ~ . + (1 | ID)), c("A1", "A2"))
out_a[["GLMM"]] <- lincom(m$est, m$vcov, La, lab_a)
cat("GLMM model (a): participant random intercept SD =", round(m$re_sd, 4), "\n\n")

for (m in names(out_a)) {
  cat("---", m, "---\n"); print(round(out_a[[m]], 5))
  cat("ratio scale:", sprintf("%.3f", exp(out_a[[m]]$Estimate)), "\n\n")
}

cat("Comparator minus EMEE (log scale):\n")
print(round(data.frame(`GEE-ind`  = out_a[["GEE-ind"]]$Estimate  - out_a$EMEE$Estimate,
                       `GEE-exch` = out_a[["GEE-exch"]]$Estimate - out_a$EMEE$Estimate,
                       GLMM       = out_a$GLMM$Estimate          - out_a$EMEE$Estimate,
                       row.names = lab_a, check.names = FALSE), 5))

## ============================================ MODEL (b): CEE moderated by decision point
## This is the analysis whose mean model is NOT saturated in the moderator: the decision
## point index enters linearly rather than as indicators. The saturation argument of
## Appendix C therefore does not cover it, and the two approaches need not coincide.
cat("\n\n================ MODEL (b): CEE moderated by the decision point index ================\n\n")

emee_b <- fit_emee(moderator = "decision_index", control = CONTROL)
## EMEE coefficient order: (intercept, slope) for trt 1, then for trt 2
Lb <- rbind(c( 1, 0,  0, 0), c( 0, 1,  0, 0),
            c( 0, 0,  1, 0), c( 0, 0,  0, 1),
            c(-1, 0,  1, 0), c( 0,-1,  0, 1))
lab_b <- c("trt1 intercept", "trt1 slope", "trt2 intercept", "trt2 slope",
           "trt2-trt1 intercept", "trt2-trt1 slope")
nbb <- names(emee_b$beta_hat)
out_b <- list(EMEE = lincom(emee_b$beta_hat, emee_b$varcov_adjusted[nbb, nbb], Lb, lab_b))

form_b <- binary_outcome ~ (A1 + A2) * decision_index + age + gender +
  employment_1 + employment_2 + AUDIT_score + binary_outcome.lag
gee_terms_b <- c("A1", "A1:decision_index", "A2", "A2:decision_index")
for (cs in c("independence", "exchangeable")) {
  g <- fit_gee(form_b, cs, gee_terms_b)
  out_b[[if (cs == "independence") "GEE-ind" else "GEE-exch"]] <- lincom(g$est, g$vcov, Lb, lab_b)
}
m <- fit_glmm(update(form_b, . ~ . + (1 | ID)), gee_terms_b)
out_b[["GLMM"]] <- lincom(m$est, m$vcov, Lb, lab_b)
cat("GLMM model (b): participant random intercept SD =", round(m$re_sd, 4), "\n\n")

for (m in names(out_b)) {
  cat("---", m, "---\n"); print(round(out_b[[m]], 5))
  cat("ratio scale:", sprintf("%.3f", exp(out_b[[m]]$Estimate)), "\n\n")
}

cat("Comparator minus EMEE (log scale):\n")
print(round(data.frame(`GEE-ind`  = out_b[["GEE-ind"]]$Estimate  - out_b$EMEE$Estimate,
                       `GEE-exch` = out_b[["GEE-exch"]]$Estimate - out_b$EMEE$Estimate,
                       GLMM       = out_b$GLMM$Estimate          - out_b$EMEE$Estimate,
                       row.names = lab_b, check.names = FALSE), 5))

## ==================================================== MODEL (c): AUDIT-moderated CEE
cat("\n\n================ MODEL (c): CEE moderated by AUDIT category ================\n\n")

emee_c <- fit_emee(moderator = c("AUDIT_harmful", "AUDIT_atRisk"), control = CONTROL_C)
## EMEE coefficient order: (int, harmful, atRisk) for trt 1, then for trt 2
Lc <- rbind(
  c( 1, 0, 0,  0, 0, 0), c( 1, 1, 0,  0, 0, 0), c( 1, 0, 1,  0, 0, 0),
  c( 0, 0, 0,  1, 0, 0), c( 0, 0, 0,  1, 1, 0), c( 0, 0, 0,  1, 0, 1),
  c(-1, 0, 0,  1, 0, 0), c(-1,-1, 0,  1, 1, 0), c(-1, 0,-1,  1, 0, 1))
lab_c <- c("trt1 hazardous", "trt1 harmful", "trt1 dependent",
           "trt2 hazardous", "trt2 harmful", "trt2 dependent",
           "trt2-trt1 hazardous", "trt2-trt1 harmful", "trt2-trt1 dependent")
nbc <- names(emee_c$beta_hat)
out_c <- list(EMEE = lincom(emee_c$beta_hat, emee_c$varcov_adjusted[nbc, nbc], Lc, lab_c))

form_c <- binary_outcome ~ (A1 + A2) * (AUDIT_harmful + AUDIT_atRisk) +
  decision_index + age + gender + employment_1 + employment_2 + binary_outcome.lag
gee_terms_c <- c("A1", "A1:AUDIT_harmful", "A1:AUDIT_atRisk",
                 "A2", "A2:AUDIT_harmful", "A2:AUDIT_atRisk")
for (cs in c("independence", "exchangeable")) {
  g <- fit_gee(form_c, cs, gee_terms_c)
  out_c[[if (cs == "independence") "GEE-ind" else "GEE-exch"]] <- lincom(g$est, g$vcov, Lc, lab_c)
}
m <- fit_glmm(update(form_c, . ~ . + (1 | ID)), gee_terms_c)
out_c[["GLMM"]] <- lincom(m$est, m$vcov, Lc, lab_c)
cat("GLMM model (c): participant random intercept SD =", round(m$re_sd, 4), "\n\n")

for (m in names(out_c)) {
  cat("---", m, "---\n"); print(round(out_c[[m]], 5))
  cat("ratio scale:", sprintf("%.3f", exp(out_c[[m]]$Estimate)), "\n\n")
}

cat("Comparator minus EMEE (log scale):\n")
print(round(data.frame(`GEE-ind`  = out_c[["GEE-ind"]]$Estimate  - out_c$EMEE$Estimate,
                       `GEE-exch` = out_c[["GEE-exch"]]$Estimate - out_c$EMEE$Estimate,
                       GLMM       = out_c$GLMM$Estimate          - out_c$EMEE$Estimate,
                       row.names = lab_c, check.names = FALSE), 5))

## ==================================================== the Appendix C equivalence claim
cat("\n\n================ equivalence check: f_t = g_t = 1, J_t = 1 ================\n\n")

n_time <- length(unique(dta$decision_index))
ptilde_equal_p <- matrix(rep(c(0.4, 0.3, 0.3), times = n_time), ncol = 3, byrow = TRUE)
emee_plain <- fit_emee(moderator = c(), control = c(), ptilde = ptilde_equal_p)
gee_plain <- geeglm(binary_outcome ~ A1 + A2, data = dta, id = idnum,
                    corstr = "independence", family = poisson(link = "log"))
ybar <- tapply(dta$binary_outcome, dta$treatment_cate, mean)

cmp <- data.frame(
  row.names = c("beta_1", "beta_2"),
  EMEE = as.vector(emee_plain$beta_hat),
  `GEE-ind` = as.vector(coef(gee_plain)[c("A1", "A2")]),
  log_ratio_of_arm_means = as.vector(log(ybar[c("1", "2")] / ybar["0"])),
  check.names = FALSE)
print(round(cmp, 10))
cat("\nmax absolute difference:", max(abs(cmp$EMEE - cmp[["GEE-ind"]])),
    max(abs(cmp$EMEE - cmp$log_ratio_of_arm_means)), "\n")
cat("arm-wise means of Y:", sprintf("%.6f", ybar), "\n")

## with tilde p = 1/3 (the default) the weights are not 1 and the equality is only approximate
emee_default_ptilde <- fit_emee(moderator = c(), control = c())
cat("\nEMEE with tilde p = 1/3 (J_t != 1):", sprintf("%.6f", emee_default_ptilde$beta_hat), "\n")
cat("exp():", sprintf("%.4f", exp(emee_default_ptilde$beta_hat)),
    "  alpha_0:", sprintf("%.6f", emee_default_ptilde$alpha_hat[1]), "\n")

cat("\n=== sessionInfo ===\n")
cat(R.version.string, "| geepack", as.character(packageVersion("geepack")),
    "| lme4", as.character(packageVersion("lme4")),
    "| dplyr", as.character(packageVersion("dplyr")), "\n")
