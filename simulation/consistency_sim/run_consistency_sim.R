# The consistency and asymptotic normality simulation of Supplementary Section C:
# the three tables that compare EMEE-catA with GEE-ind and GEE-exch under a
# misspecified nuisance model, at n = 20, 30, 40, 50, 100.
#
# This driver produces the tables as published. It replaces the two per-scenario
# pipelines in sim_consistency_moderated/ and
# sim_consistency_marginal_model_with control_non_const_pt/, which were written for a
# cluster, use two different generative models, and do not reproduce the published
# numbers; they are kept for the record only. See README.md in this folder.
#
# One generative model, the one displayed in the supplement, serves both scenarios.
#
# Usage:  Rscript run_consistency_sim.R [nsim] [ncores]      (default 1000 and cores - 2)
#         Rscript summarize_consistency.R                    then emits the table bodies

suppressMessages({
  library(rootSolve)
  library(geepack)
  library(parallel)
})

# This script's own folder, whether it is run with Rscript or sourced in RStudio.
HERE <- local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
ROOT <- normalizePath(file.path(HERE, "..", ".."))
source(file.path(ROOT, "functions", "emee_catA.R"))
source(file.path(HERE, "dgm_consistency.R"))

args <- commandArgs(trailingOnly = TRUE)
nsim <- if (length(args) >= 1) as.integer(args[1]) else 1000
ncores <- if (length(args) >= 2) as.integer(args[2]) else max(1, detectCores() - 2)

TOTAL_T <- 30
SAMPLE_SIZES <- c(20, 30, 40, 50, 100)
OUTDIR <- file.path(HERE, "results")
dir.create(OUTDIR, showWarnings = FALSE)

# --- one replicate ----------------------------------------------------------

one_rep <- function(sample_size, scenario, moderated) {
  dta <- dgm_consistency(sample_size, TOTAL_T, scenario = scenario)

  fit <- emee_catA(
    dta = dta, id_varname = "userid", decision_time_varname = "time",
    treatment_varname = "A", outcome_varname = "Y",
    control_varname = c("S"),
    moderator_varname = if (moderated) c("S") else c(),
    rand_prob_varname = "prob_A", estimator_initial_value = NULL,
    trt_level = 2, pmatrix_tilde = NULL, avail_varname = NULL
  )

  gee_terms <- if (moderated) c("A1", "A1:S", "A2", "A2:S") else c("A1", "A2")
  gee_form <- if (moderated) Y ~ A1 * S + A2 * S else Y ~ S + A1 + A2

  fit_gee <- function(corstr) {
    g <- summary(geeglm(gee_form, data = dta, id = userid,
                        corstr = corstr, family = poisson(link = "log")))
    rn <- rownames(g$coefficients)
    # geepack names interaction terms A1:S or S:A1 depending on term order
    idx <- vapply(gee_terms, function(tm) {
      hit <- which(rn == tm)
      if (length(hit) == 0 && grepl(":", tm)) {
        hit <- which(rn == paste(rev(strsplit(tm, ":")[[1]]), collapse = ":"))
      }
      hit[1]
    }, integer(1))
    list(est = g$coefficients[idx, "Estimate"],
         se  = g$coefficients[idx, "Std.err"])
  }

  gee1 <- fit_gee("independence")
  gee2 <- fit_gee("exchangeable")

  list(
    emee_est = fit$beta_hat,
    emee_se = fit$beta_se,
    emee_se_adj = fit$beta_se_adjusted,
    emee_ci_lb = fit$conf_int_adjusted_t[, 1],
    emee_ci_ub = fit$conf_int_adjusted_t[, 2],
    gee1_est = gee1$est, gee1_se = gee1$se,
    gee2_est = gee2$est, gee2_se = gee2$se
  )
}

# --- run --------------------------------------------------------------------

run_setting <- function(scenario, moderated, label, seed) {
  set.seed(seed)   # per setting, so each table is reproducible on its own
  out <- list()
  for (n in SAMPLE_SIZES) {
    t0 <- Sys.time()
    reps <- mclapply(1:nsim, function(i) {
      tryCatch(one_rep(n, scenario, moderated), error = function(e) NULL)
    }, mc.cores = ncores)
    n_fail <- sum(vapply(reps, is.null, logical(1)))
    cat(sprintf("[%s] n = %3d done in %5.1f min, %d failed replicates\n",
                label, n, as.numeric(difftime(Sys.time(), t0, units = "mins")), n_fail))
    flush.console()
    out[[as.character(n)]] <- list(sample_size = n, total_T = TOTAL_T,
                                   n_fail = n_fail, reps = reps[!vapply(reps, is.null, logical(1))])
  }
  saveRDS(out, file.path(OUTDIR, paste0("consistency_", label, "_nsim", nsim, ".RDS")))
  invisible(out)
}

RNGkind("L'Ecuyer-CMRG")
cat(sprintf("nsim = %d, ncores = %d, T = %d\n", nsim, ncores, TOTAL_T))

# Three of the four (estimand, randomization) combinations. The fourth, moderated CEE
# under Z_t-dependent randomization, is omitted: GEE would fail there for both of the
# reasons the other two settings isolate, so it adds nothing.
run_setting("uniform",     TRUE,  "moderated_scenario1", seed = 20260731)
run_setting("Z_dependent", FALSE, "marginal_scenario2",  seed = 20260732)
run_setting("uniform",     FALSE, "marginal_scenario1",  seed = 20260733)

cat("done\n")
