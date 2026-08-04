# Drink Less: analysis (b') for the response to Reviewer 2, Comment 6
# --------------------------------------------------------------------
# Companion to gee_comparison.R, which produces the (a), (b), (c) entries of the
# comparison table in the response letter. This script produces the (b') entries.
#
# Analysis (b') refits analysis (b) (the CEE moderated by the decision point index)
# with the nuisance model deliberately misspecified: the decision point index and the
# lagged outcome are removed from g_t, so that no covariate in g_t tracks the decline
# in app opening over the study. Its purpose is to check that the close EMEE / GEE-ind
# agreement in analysis (b) does not reflect a fortunately specified nuisance model.
# In the GEE and GLMM mean models the main effect of the decision point index remains,
# as the treatment-by-index interaction terms require.
#
# Conventions match gee_comparison.R exactly: EMEE uses the estimator's default
# tilde-p = 1/3 (as in Table 2 of the paper) and the small-sample adjusted standard
# error; GEE-ind and GEE-exch use robust standard errors; the GLMM standard errors are
# model based. For reference the script also reruns analysis (b), whose entries should
# reproduce the response table.
#
# DATA must point at the Drink Less MRT dataset shared on OSF (https://osf.io/w3szp).

DATA <- file.path("/Users/tqian/Dropbox (Personal)/data/DrinkLess MRT (Lauren Bell)",
                  "2022.01 - DrinkLess OSF data sharing/FINAL Dataset_A.rds")
REPO <- "~/repos/cat_trt_binary"

suppressPackageStartupMessages({
  library(dplyr)
  library(geepack)
  library(lme4)
})
source(file.path(REPO, "functions/wcls_cat_trt_binary_outcome.R"))

options(digits = 6, width = 120)

## ------------------------------------------------------------------ preprocessing
## identical to gee_comparison.R
raw <- suppressWarnings(readRDS(DATA))
dta <- raw %>%
  ungroup() %>%
  mutate(treatment_cate = ifelse(subversion == "Y", 2, ifelse(subversion == "X", 1, 0)),
         decision_index = days_since_download - 1,   # t starts at 0
         binary_outcome = primary_outcome) %>%
  dplyr::select(ID, decision_index, age, AUDIT_score, gender, employment_type,
                binary_outcome, treatment_cate) %>%
  arrange(ID, decision_index) %>%
  group_by(ID) %>%
  mutate(binary_outcome.lag = dplyr::lag(binary_outcome, default = 0)) %>%
  ungroup() %>%
  mutate(idnum = as.integer(factor(ID)),
         I = 1,
         prob_A = ifelse(treatment_cate == 0, 0.4, 0.3),
         A1 = as.numeric(treatment_cate == 1),
         A2 = as.numeric(treatment_cate == 2),
         employment_1 = as.numeric(employment_type == 1),
         employment_2 = as.numeric(employment_type == 2)) %>%
  as.data.frame()
stopifnot(nrow(dta) == 10470, length(unique(dta$ID)) == 349)

CONTROL_B  <- c("decision_index", "age", "gender", "employment_1", "employment_2",
                "AUDIT_score", "binary_outcome.lag")           # analysis (b)
CONTROL_BP <- c("age", "gender", "employment_1", "employment_2", "AUDIT_score")  # (b')

## ------------------------------------------------------------------ fitting helpers
## tilde p = the trial's randomization probability, so J_t = 1, as in Section 5.1.
N_TIME <- length(unique(dta$decision_index))
PTILDE <- matrix(rep(c(0.4, 0.3, 0.3), times = N_TIME), ncol = 3, byrow = TRUE)

fit_emee <- function(control) {
  wcls_categorical_treatment(
    dta = dta, id_varname = "ID", decision_time_varname = "decision_index",
    treatment_varname = "treatment_cate", outcome_varname = "binary_outcome",
    control_varname = control, moderator_varname = "decision_index",
    rand_prob_varname = "prob_A", avail_varname = "I",
    trt_level = 2, estimator_initial_value = NULL, pmatrix_tilde = PTILDE)
}
GEE_TERMS <- c("A1", "A1:decision_index", "A2", "A2:decision_index")
pick <- function(b, want) {
  vapply(want, function(tm) if (tm %in% names(b)) tm
         else paste(rev(strsplit(tm, ":")[[1]]), collapse = ":"), character(1))
}
fit_gee <- function(control, corstr) {
  f <- reformulate(c("(A1 + A2) * decision_index", setdiff(control, "decision_index")),
                   response = "binary_outcome")
  g <- geeglm(f, data = dta, id = idnum, corstr = corstr, family = poisson(link = "log"))
  V <- g$geese$vbeta; rownames(V) <- colnames(V) <- names(coef(g))
  idx <- pick(coef(g), GEE_TERMS)
  list(est = coef(g)[idx], vcov = V[idx, idx])
}
fit_glmm <- function(control) {
  f <- reformulate(c("(A1 + A2) * decision_index", setdiff(control, "decision_index"),
                     "(1 | ID)"), response = "binary_outcome")
  m <- glmer(f, data = dta, family = poisson(link = "log"),
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
  msgs <- m@optinfo$conv$lme4$messages
  if (length(msgs)) warning("glmer: ", paste(msgs, collapse = "; "))
  b <- fixef(m); V <- as.matrix(vcov(m)); idx <- pick(b, GEE_TERMS)
  list(est = b[idx], vcov = V[idx, idx], re_sd = as.data.frame(VarCorr(m))$sdcor)
}

## displayed quantities, matching the response table: the two slopes and the
## intercept difference; EMEE coefficient order is (int, slope) for trt 1 then trt 2
L <- rbind(c(0, 1, 0, 0), c(0, 0, 0, 1), c(-1, 0, 1, 0))
labs <- c("beta1 t-slp", "beta2 t-slp", "int diff")
lin <- function(est, vcov) data.frame(row.names = labs,
                                      est = as.vector(L %*% est),
                                      se  = sqrt(diag(L %*% vcov %*% t(L))))
texnum <- function(r) sprintf("%s%.3f (%.3f)", ifelse(r$est < 0, "$-$", "\\phantom{$-$}"),
                              abs(r$est), r$se)

## ------------------------------------------------------------------ run both analyses
for (spec in list(list(ctrl = CONTROL_B,  lab = "(b): g_t = the paper's controls"),
                  list(ctrl = CONTROL_BP,
                       lab = "(b'): g_t without the decision point index and the lagged outcome"))) {
  cat("==============", spec$lab, "==============\n")
  e  <- fit_emee(spec$ctrl)
  nb <- names(e$beta_hat)
  gi <- fit_gee(spec$ctrl, "independence")
  gx <- fit_gee(spec$ctrl, "exchangeable")
  mm <- fit_glmm(spec$ctrl)
  cat("coefficients (trt1 int, trt1 slope, trt2 int, trt2 slope):\n")
  print(round(rbind(EMEE = e$beta_hat, `GEE-ind` = gi$est,
                    gap = gi$est - e$beta_hat,
                    gap_in_adj_SE = (gi$est - e$beta_hat) / e$beta_se_adjusted), 4))
  cat(sprintf("max |GEE-ind - EMEE| over the 4 coefficients: %.4f\n",
              max(abs(gi$est - e$beta_hat))))
  cat(sprintf("max |GEE-exch - EMEE|: %.4f   max |GLMM - EMEE|: %.4f   GLMM RI SD: %.3f\n",
              max(abs(gx$est - e$beta_hat)), max(abs(mm$est - e$beta_hat)), mm$re_sd))
  rE <- lin(e$beta_hat, e$varcov_adjusted[nb, nb])
  rI <- lin(gi$est, gi$vcov); rX <- lin(gx$est, gx$vcov); rM <- lin(mm$est, mm$vcov)
  cat("\nresponse table rows (EMEE & GEE-ind & GEE-exch & GLMM):\n")
  for (i in 1:3) cat(sprintf("  %-12s & %s & %s & %s & %s \\\\\n", labs[i],
                             texnum(rE[i, ]), texnum(rI[i, ]),
                             texnum(rX[i, ]), texnum(rM[i, ])))
  cat("\n")
}

## the descriptive quoted in the response text: app opening under no message, by week
cat("app opening rate under no message, by week of the study:\n")
print(round(with(subset(dta, treatment_cate == 0),
                 tapply(binary_outcome, cut(decision_index, c(-1, 6, 13, 20, 29),
                                            labels = c("week 1", "week 2", "week 3",
                                                       "week 4+")), mean)), 4))

cat("\n=== sessionInfo ===\n")
cat(R.version.string, "| geepack", as.character(packageVersion("geepack")),
    "| lme4", as.character(packageVersion("lme4")),
    "| dplyr", as.character(packageVersion("dplyr")), "\n")
