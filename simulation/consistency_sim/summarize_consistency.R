# Summarize the Supplementary Section C consistency simulation and emit the LaTeX
# table bodies. Run run_consistency_sim.R first.
#
# Metrics per coefficient:
#   Bias = mean(estimate) - truth
#   SD   = Monte Carlo standard deviation of the estimates across replicates
#   SE   = mean of the estimated standard errors
#          (EMEE: small-sample adjusted; GEE: robust sandwich)
#   CP   = coverage of the 95% CI
#          (EMEE: t-based with the adjusted SE, as in the paper; GEE: est +- 1.96 SE)
#
# The published tables printed the Monte Carlo SD in the column labelled SE, and
# reported RMSE = sqrt(bias^2 + SD^2); both are recoverable from the columns above.

# This script's own folder, whether it is run with Rscript or sourced in RStudio.
HERE <- local({
  a <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else getwd()
})
ROOT <- normalizePath(file.path(HERE, "..", ".."))
source(file.path(HERE, "dgm_consistency.R"))
nsim <- local({ a <- commandArgs(trailingOnly = TRUE)
                if (length(a)) as.integer(a[1]) else 1000 })

truth <- true_cee()
TRUTH <- list(
  moderated_scenario1 = c(truth$moderated$beta1, truth$moderated$beta2), # b10, b11, b20, b21
  marginal_scenario2 = unname(truth$marginal),
  marginal_scenario1 = unname(truth$marginal)
)
COEF_LABEL <- list(
  moderated_scenario1 = c("trt1 Intercept", "trt1 Z_t", "trt2 Intercept", "trt2 Z_t"),
  marginal_scenario2 = c("trt1", "trt2"),
  marginal_scenario1 = c("trt1", "trt2")
)

metrics <- function(est, se, lb, ub, truth_j) {
  c(bias = mean(est) - truth_j,
    sd = sd(est),
    se = mean(se),
    cp = mean(lb <= truth_j & truth_j <= ub),
    rmse = sqrt((mean(est) - truth_j)^2 + sd(est)^2))
}

collect <- function(label) {
  res <- readRDS(file.path(HERE, "results", paste0("consistency_", label, "_nsim", nsim, ".RDS")))
  tru <- TRUTH[[label]]
  out <- list()
  for (nm in names(res)) {
    reps <- res[[nm]]$reps
    grab <- function(field) do.call(rbind, lapply(reps, function(r) r[[field]]))
    emee_est <- grab("emee_est"); emee_se <- grab("emee_se_adj")
    emee_lb <- grab("emee_ci_lb"); emee_ub <- grab("emee_ci_ub")
    g1e <- grab("gee1_est"); g1s <- grab("gee1_se")
    g2e <- grab("gee2_est"); g2s <- grab("gee2_se")
    for (j in seq_along(tru)) {
      out[[length(out) + 1]] <- data.frame(
        sample_size = as.integer(nm), coef = COEF_LABEL[[label]][j], method = "EMEE-catA",
        t(metrics(emee_est[, j], emee_se[, j], emee_lb[, j], emee_ub[, j], tru[j])))
      out[[length(out) + 1]] <- data.frame(
        sample_size = as.integer(nm), coef = COEF_LABEL[[label]][j], method = "GEE-ind",
        t(metrics(g1e[, j], g1s[, j], g1e[, j] - 1.96 * g1s[, j], g1e[, j] + 1.96 * g1s[, j], tru[j])))
      out[[length(out) + 1]] <- data.frame(
        sample_size = as.integer(nm), coef = COEF_LABEL[[label]][j], method = "GEE-exch",
        t(metrics(g2e[, j], g2s[, j], g2e[, j] - 1.96 * g2s[, j], g2e[, j] + 1.96 * g2s[, j], tru[j])))
    }
  }
  do.call(rbind, out)
}

fmt <- function(x) sprintf("%.5f", x)
fmtcp <- function(x) sprintf("%.3f", x)

## ---------------------------------------------------------------- Scenario 1
s1 <- collect("moderated_scenario1")
cat("### Scenario 1: moderated CEE, equal randomization, truth =",
    fmt(TRUTH$moderated_scenario1), "\n\n")
print(s1, row.names = FALSE, digits = 4)

cat("\n\n%%% LaTeX body for tbl:est_moderator %%%\n")
for (trt in 1:2) {
  ci <- if (trt == 1) c("trt1 Intercept", "trt1 Z_t") else c("trt2 Intercept", "trt2 Z_t")
  for (meth in c("EMEE-catA", "GEE-ind", "GEE-exch")) {
    for (k in seq_along(sort(unique(s1$sample_size)))) {
      n <- sort(unique(s1$sample_size))[k]
      a <- s1[s1$sample_size == n & s1$coef == ci[1] & s1$method == meth, ]
      b <- s1[s1$sample_size == n & s1$coef == ci[2] & s1$method == meth, ]
      lead <- if (k == 1 && meth == "EMEE-catA") sprintf("%d & %s\n ", trt, meth)
              else if (k == 1) sprintf(" & %s\n ", meth) else " & "
      cat(sprintf("%s& %-3d & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n",
                  lead, n, fmt(a$bias), fmt(a$sd), fmt(a$se), fmtcp(a$cp),
                  fmt(b$bias), fmt(b$sd), fmt(b$se), fmtcp(b$cp)))
    }
    if (meth != "GEE-exch") cat("\\cline{2-11}\n")
  }
  cat(if (trt == 1) "\\hline \\hline\n" else "\\hline\n")
}

## ------------------------------------------------- the two marginal-CEE settings
emit_marginal <- function(label, header, tag) {
  d <- collect(label)
  cat("\n\n### ", header, ", truth = ", paste(fmt(TRUTH[[label]]), collapse = " "), "\n\n", sep = "")
  print(d, row.names = FALSE, digits = 4)
  cat("\n\n%%% LaTeX body for ", tag, " %%%\n", sep = "")
  for (trt in 1:2) {
    cf <- paste0("trt", trt)
    for (meth in c("EMEE-catA", "GEE-ind", "GEE-exch")) {
      for (k in seq_along(sort(unique(d$sample_size)))) {
        n <- sort(unique(d$sample_size))[k]
        a <- d[d$sample_size == n & d$coef == cf & d$method == meth, ]
        lead <- if (k == 1 && meth == "EMEE-catA") sprintf("%d & %s\n ", trt, meth)
                else if (k == 1) sprintf(" & %s\n ", meth) else " & "
        cat(sprintf("%s& %-3d & %s & %s & %s & %s \\\\\n",
                    lead, n, fmt(a$bias), fmt(a$sd), fmt(a$se), fmtcp(a$cp)))
      }
      if (meth != "GEE-exch") cat("\\cline{2-7}\n")
    }
    cat(if (trt == 1) "\\hline \\hline\n" else "\\hline\n")
  }
  d
}

s2 <- emit_marginal("marginal_scenario2",
                    "Scenario 2: marginal CEE, Z-dependent randomization", "tbl:est_marginal")
s3 <- emit_marginal("marginal_scenario1",
                    "Scenario 1: marginal CEE, equal randomization", "tbl:est_marginal_scenario1")

saveRDS(list(moderated_scenario1 = s1, marginal_scenario2 = s2, marginal_scenario1 = s3),
        file.path(HERE, "results", "summary_tables.RDS"))
