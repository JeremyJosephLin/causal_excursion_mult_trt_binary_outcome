# prepare_drinkless_data.R
#
# Builds the analysis data set for the Drink Less illustration in the paper
# "Micro-randomized Trials with Categorical Treatments and Binary Proximal Outcome"
# directly from the public Drink Less MRT data, so that no intermediate file has to be
# passed around.
#
# The input is "FINAL Dataset_A.rds" from the Drink Less OSF repository,
# https://osf.io/w3szp. Download it and either
#   * place it at <repo root>/dataset/FINAL Dataset_A.rds, or
#   * point the environment variable DRINKLESS_DATA at it, or
#   * pass its path to prepare_drinkless_data().
#
# Usage:
#   source("prepare_drinkless_data.R")
#   dta <- prepare_drinkless_data()
#
# What the preprocessing does, and why:
#
#   decision_index      days_since_download - 1, so that t runs 0, ..., 29. The paper
#                       starts t at 0 so that the intercept of the time-moderated model
#                       is the effect at the first decision point.
#   binary_outcome.lag  the previous decision point's outcome for the SAME participant,
#                       with 0 at each participant's first decision point. Lagging over
#                       the stacked data instead would give 349 rows the preceding
#                       participant's last outcome.
#   employment_1/2      indicators for the three-level employment type, which enters the
#                       nuisance model g_t along with the decision point index, age,
#                       gender, AUDIT score, and the lagged outcome.
#   AUDIT_harmful,      the AUDIT risk categories used as the moderator in analysis (c):
#   AUDIT_atRisk        hazardous 8-15 (reference), harmful 16-19, dependent 20-40.
#   prob_A              the trial's randomization probabilities: 0.4 for no notification
#                       and 0.3 for each of the two active options.
#   I                   availability, identically 1: every participant was available at
#                       every decision point in this trial.

prepare_drinkless_data <- function(path = NULL) {
  if (is.null(path)) path <- drinkless_data_path()
  if (!file.exists(path)) {
    stop("cannot find the Drink Less data file at\n  ", path,
         "\nDownload 'FINAL Dataset_A.rds' from https://osf.io/w3szp and put it in the ",
         "repository's dataset/ folder, or set DRINKLESS_DATA to its path.")
  }
  raw <- readRDS(path)

  needed <- c("ID", "days_since_download", "subversion", "primary_outcome",
              "age", "AUDIT_score", "gender", "employment_type")
  missing <- setdiff(needed, names(raw))
  if (length(missing)) stop("the data file is missing: ", paste(missing, collapse = ", "))

  dta <- raw
  dta$treatment_cate <- ifelse(dta$subversion == "Y", 2, ifelse(dta$subversion == "X", 1, 0))
  dta$decision_index <- dta$days_since_download - 1
  dta$binary_outcome <- dta$primary_outcome
  dta <- dta[, c("ID", "decision_index", "age", "AUDIT_score", "gender",
                 "employment_type", "binary_outcome", "treatment_cate")]

  # the raw file keeps some participants' rows in two non-adjacent blocks, so sort before
  # lagging
  dta <- dta[order(dta$ID, dta$decision_index), ]
  dta$binary_outcome.lag <- unlist(lapply(split(dta$binary_outcome, dta$ID),
                                          function(y) c(0, head(y, -1))),
                                   use.names = FALSE)

  dta$I             <- 1
  dta$prob_A        <- ifelse(dta$treatment_cate == 0, 0.4, 0.3)
  dta$employment_1  <- as.numeric(dta$employment_type == 1)
  dta$employment_2  <- as.numeric(dta$employment_type == 2)
  dta$AUDIT_harmful <- as.numeric(dta$AUDIT_score >= 16 & dta$AUDIT_score <= 19)
  dta$AUDIT_atRisk  <- as.numeric(dta$AUDIT_score >= 20 & dta$AUDIT_score <= 40)

  dta <- as.data.frame(dta)
  rownames(dta) <- NULL

  stopifnot(nrow(dta) == 10470, length(unique(dta$ID)) == 349,
            range(dta$decision_index) == c(0, 29))
  dta
}

# Where to look for the OSF file, in order of precedence.
drinkless_data_path <- function() {
  from_env <- Sys.getenv("DRINKLESS_DATA", "")
  if (nzchar(from_env)) return(from_env)
  root <- drinkless_repo_root()
  file.path(root, "dataset", "FINAL Dataset_A.rds")
}

drinkless_repo_root <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  here_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else getwd()
  for (up in c("../../..", "../..", "..", ".")) {
    cand <- normalizePath(file.path(here_dir, up), mustWork = FALSE)
    if (dir.exists(file.path(cand, "functions"))) return(cand)
  }
  normalizePath(here_dir, mustWork = FALSE)
}
