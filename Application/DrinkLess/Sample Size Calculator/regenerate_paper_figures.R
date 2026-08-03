# regenerate_paper_figures.R
#
# Regenerates the seven sample size figures of the Drink Less illustration in the paper
# "Micro-randomized Trials with Categorical Treatments and Binary Proximal Outcome".
#
# The seven `* sensitivity to n` folders next to this file hold the original per-figure
# pipelines, which were run on a cluster and whose results are distributed as saved .RDS
# grids. This script instead computes every grid from scratch, so the figures are
# reproducible from a clone of this repository plus the R packages listed below. It takes
# a couple of minutes.
#
# All seven figures use one set of design inputs, matching the Drink Less MRT:
#
#   T = 30 decision points, randomization (0.4, 0.3, 0.3) over
#   {no notification, bank of new messages, standard message},
#   ATE_1 = 3.27, ATE_2 = 3.62 (Model (a) of the paper's Table 2),
#   ASPN = 0.036 (the observed proportion of app openings when no notification was sent),
#   type I error rate 0.05, power 0.8, L = (1, -1), and full availability except where
#   availability is the quantity being varied.
#
# Usage:  Rscript "regenerate_paper_figures.R"
# Output: figures/*.pdf and figures/*.jpg, plus figures/summary.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(latex2exp)
  library(metR)
})

args     <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
here_dir <- if (length(file_arg)) dirname(sub("^--file=", "", file_arg)) else getwd()
ROOT     <- normalizePath(file.path(here_dir, "..", "..", ".."), mustWork = FALSE)
stopifnot(dir.exists(file.path(ROOT, "functions")))

source(file.path(ROOT, "functions", "wcls_cat_trt_binary_outcome.R"))
source(file.path(ROOT, "functions", "utillity.R"))
source(file.path(ROOT, "functions", "ss_calc.R"))

OUT <- file.path(ROOT, "Application", "DrinkLess", "Sample Size Calculator", "figures")
dir.create(OUT, showWarnings = FALSE)

## ------------------------------------------------------------------ design inputs

M        <- 30            # decision points
GAMMA    <- 0.05          # type I error rate
BETA_II  <- 0.2           # 1 - power
PT_SHAPE <- "constant p0 high"   # (0.4, 0.3, 0.3), the Drink Less randomization
ASPN     <- 0.036         # observed P(app opened | no notification)
ATE1     <- 3.27          # Model (a): exp(beta_1)
ATE2     <- 3.62          # Model (a): exp(beta_2)

pt <- construct_pt(M, PT_SHAPE)

## ------------------------------------------------------------------ the calculator
#
# One row of any of the seven grids. `variant` mirrors the two forms the original
# per-figure scripts take: "independent" solves the two treatment effects separately
# (used when both MEEs share a time pattern), "star" solves the second relative to the
# first, with theta_f2 controlling the difference between their slopes.

compute_n <- function(ft_shape, gt_shape,
                      ft_theta = 0, ft_theta2 = 0, gt_theta = 0,
                      ate1 = ATE1, ate2 = ATE2, aspn = ASPN,
                      avg_tau = 1, tau_shape = "constant", tau_theta = 0,
                      variant = c("independent", "star")) {
  variant <- match.arg(variant)
  out <- tryCatch({
    taut  <- construct_taut_theta(M, avg_tau, tau_shape, theta = tau_theta)
    gt    <- construct_ftgt(M, gt_shape)
    alpha <- solve_alpha(gt_shape, M, aspn, gt, taut, theta = gt_theta)
    stopifnot(isTRUE(all.equal(aspn, compute_ASPNC(gt, alpha, taut))))

    ft    <- construct_ftgt(M, ft_shape)
    beta1 <- solve_beta(shape = ft_shape, m = M, ATE = ate1, ft = ft, taut = taut,
                        theta = ft_theta, gt = gt, alpha = alpha)
    beta2 <- if (variant == "independent") {
      solve_beta(shape = ft_shape, m = M, ATE = ate2, ft = ft, taut = taut,
                 theta = ft_theta, gt = gt, alpha = alpha)
    } else {
      beta1 + solve_beta_star(shape = ft_shape, m = M, ATE = ate2, ft = ft, taut = taut,
                              theta = ft_theta2, gt = gt, alpha = alpha, beta1 = beta1)
    }
    stopifnot(isTRUE(all.equal(ate1, compute_ATE(gt, alpha, ft, beta1, taut))),
              isTRUE(all.equal(ate2, compute_ATE(gt, alpha, ft, beta2, taut))))

    # every arm's success probability must stay in (0, 1) for the model to be a model
    pr <- cbind(exp(gt %*% alpha), exp(gt %*% alpha + ft %*% beta1),
                exp(gt %*% alpha + ft %*% beta2))
    stopifnot(all(pr > 0 & pr < 1))

    mrt_mult_trt_ss(f_t = ft, g_t = gt, beta = c(beta1, beta2), alpha = alpha, p_t = pt,
                    gamma = GAMMA, b = BETA_II, L = construct_L(ft_shape), tau = taut)
  }, error = function(e) NA_real_)
  as.numeric(out)
}

# compute_n over a grid, one row at a time
grid_n <- function(grid, variant) {
  grid$n <- vapply(seq_len(nrow(grid)), function(i) {
    do.call(compute_n, c(as.list(grid[i, , drop = FALSE]), list(variant = variant)))
  }, numeric(1))
  grid[!is.na(grid$n), ]
}

## ------------------------------------------------------------------ shared plot theme

myggfont <- function(legend_text_size = 18, legend_title_size = 20, axis_text_size = 22,
                     axis_title_size = 26, plot_title_size = 20, facet_text_size = 16) {
  theme(legend.text  = element_text(size = legend_text_size),
        legend.title = element_text(size = legend_title_size),
        axis.text    = element_text(size = axis_text_size),
        axis.title   = element_text(size = axis_title_size, face = "bold"),
        plot.title   = element_text(size = plot_title_size),
        strip.text.x = element_text(size = facet_text_size),
        strip.text.y = element_text(size = facet_text_size))
}

save_fig <- function(name, height, width) {
  ggsave(file.path(OUT, paste0(name, ".pdf")), height = height, width = width)
  ggsave(file.path(OUT, paste0(name, ".jpg")), height = height, width = width)
}

log_lines <- character(0)
note <- function(...) {
  line <- paste0(...)
  log_lines <<- c(log_lines, line)
  cat(line, "\n")
}

note("design: T = ", M, ", p = (0.4, 0.3, 0.3), ATE = (", ATE1, ", ", ATE2, "), ASPN = ",
     ASPN, ", type I error rate ", GAMMA, ", power ", 1 - BETA_II)
n_reference <- compute_n(ft_shape = "constant", gt_shape = "constant")
note("constant MEE(t) and constant SPNC(t), full availability: n = ", n_reference,
     "   <- the sample size quoted in Section 5.2")

## ------------------------------------------------------------------ DL1: n vs (ATE1, ATE2)

g1 <- expand.grid(ate1 = seq(3.10, 3.50, by = 0.01),
                  ate2 = seq(3.40, 3.80, by = 0.01),
                  ft_shape = "constant", gt_shape = "constant",
                  stringsAsFactors = FALSE)
g1 <- grid_n(g1, "independent")
note("DL1: ", nrow(g1), " grid points, n from ", min(g1$n), " to ", max(g1$n))

brks_line <- seq(250, 8000, by = 250)
brks_lab  <- c(500, 750, 1000, 1250, 1500, 2000, 3000)
# contours and their labels are computed on the displayed window only, so that
# label_placer_flattest() does not park a label outside the plotted region
g1_win <- g1[g1$ate1 >= 3.10 & g1$ate1 <= 3.30 & g1$ate2 >= 3.50 & g1$ate2 <= 3.80, ]
ggplot(g1_win, aes(x = ate1, y = ate2, z = n)) +
  stat_contour(aes(z = n), breaks = brks_line, linewidth = 0.4, na.rm = TRUE,
               colour = "blue") +
  geom_text_contour(aes(z = n, label = after_stat(level)), breaks = brks_lab, skip = 0,
                    stroke = 0.2, size = 4.5, label.placer = label_placer_flattest(),
                    na.rm = TRUE) +
  geom_abline(slope = 1, intercept = ATE2 - ATE1, colour = "black", linetype = "dashed") +
  annotate("point", x = ATE1, y = ATE2, size = 2.5) +
  annotate("text", x = ATE1 - 0.006, y = ATE2 + 0.013, hjust = 1, size = 5,
           label = paste0("n = ", n_reference)) +
  xlab(TeX(r'(ATE$_1$)')) + ylab(TeX(r'(ATE$_2$)')) +
  coord_cartesian(xlim = c(3.10, 3.30), ylim = c(3.50, 3.80)) +
  theme_bw() + myggfont() +
  theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
save_fig("size Drink Less 1 - n vs ASPN, ATE (contour)", height = 6, width = 6)

## ------------------------------------------------------------------ DL2: n vs theta_f1

g2 <- expand.grid(ft_theta = seq(-1, 1, by = 0.01),
                  ft_shape = "linear_theta", gt_shape = "quadratic_theta", gt_theta = 0.2,
                  stringsAsFactors = FALSE)
g2 <- grid_n(g2, "star")
n2_ref <- compute_n(ft_shape = "constant", gt_shape = "quadratic_theta", gt_theta = 0.2)
note("DL2: ", nrow(g2), " points, n from ", min(g2$n), " to ", max(g2$n),
     "; constant-MEE reference n = ", n2_ref)

## ------------------------------------------------------------------ DL3: n vs theta_f2

g3 <- expand.grid(ft_theta2 = seq(-1, 1, by = 0.01),
                  ft_theta = c(-0.3, 0, 0.3),
                  ft_shape = "linear_theta", gt_shape = "quadratic_theta", gt_theta = 0.2,
                  stringsAsFactors = FALSE)
g3 <- grid_n(g3, "star")
n3_ref <- n2_ref
note("DL3: ", nrow(g3), " points, n from ", min(g3$n), " to ", max(g3$n))

## ------------------------------------------------------------------ DL5: n vs theta_g

g5 <- expand.grid(gt_theta = seq(-0.9, 0.9, by = 0.005),
                  gt_shape = c("linear_theta", "quadratic_theta"),
                  ft_shape = "constant",
                  stringsAsFactors = FALSE)
g5 <- grid_n(g5, "star")
n5_ref <- n_reference
note("DL5: ", nrow(g5), " points, n from ", min(g5$n), " to ", max(g5$n),
     "; theta_g feasible over [", min(g5$gt_theta), ", ", max(g5$gt_theta), "]")

# a common vertical scale for the three MEE/SPNC-pattern figures, as in the paper
ylim_pattern <- c(floor(min(g2$n, g3$n, g5$n) / 100) * 100,
                  ceiling(max(g2$n, g3$n, g5$n, n2_ref, n5_ref) / 100) * 100)
note("common y range for DL2, DL3, DL5: ", ylim_pattern[1], " to ", ylim_pattern[2])

ggplot(g2, aes(x = ft_theta, y = n)) +
  geom_line() +
  geom_hline(yintercept = n2_ref, linetype = 2) +
  annotate("text", x = 0, y = n2_ref + diff(ylim_pattern) * 0.035,
           label = paste0("constant MEE(t): n = ", n2_ref), size = 5) +
  xlab(TeX(r'($\theta_{f1}$)')) + ylab("n (sample size)") +
  coord_cartesian(ylim = ylim_pattern) +
  theme_bw() + myggfont() +
  theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
save_fig("size Drinkless n vs ft1_shape", height = 6, width = 8)

g3$mee <- paste0(g3$ft_shape, g3$ft_theta)
ggplot(g3, aes(x = ft_theta2, y = n, colour = mee)) +
  geom_line() +
  geom_hline(yintercept = n3_ref, linetype = 2) +
  scale_colour_discrete(name = "",
                        labels = c(TeX(r'(linear, $\theta_{f1} = -0.3$)'),
                                   TeX(r'(linear, $\theta_{f1} = 0$)'),
                                   TeX(r'(linear, $\theta_{f1} = 0.3$)'))) +
  annotate("text", x = 0, y = n3_ref + diff(ylim_pattern) * 0.035,
           label = paste0("constant MEE(t): n = ", n3_ref), size = 5) +
  xlab(TeX(r'($\theta_{f2}$)')) + ylab("n (sample size)") +
  coord_cartesian(ylim = ylim_pattern) +
  theme_bw() + myggfont(legend_text_size = 14) +
  theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5),
        legend.position = "top", legend.key.width = unit(1.1, "lines"))
save_fig("size Drinkless n vs ft2_shape", height = 6, width = 6)

ggplot(g5, aes(x = gt_theta, y = n, colour = gt_shape)) +
  geom_line() +
  geom_hline(yintercept = n5_ref, linetype = 2) +
  scale_colour_manual(name = "", labels = c("linear g(t)", "quadratic g(t)"),
                      values = c(4, 7)) +
  annotate("text", x = 0, y = n5_ref + diff(ylim_pattern) * 0.035,
           label = paste0("constant SPNC(t): n = ", n5_ref), size = 5) +
  scale_x_continuous(breaks = c(-0.9, -0.5, 0, 0.5, 0.9)) +
  xlab(TeX(r'($\theta_{g}$)')) + ylab("n (sample size)") +
  coord_cartesian(ylim = ylim_pattern) +
  theme_bw() + myggfont() +
  theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5),
        legend.position = "top")
save_fig("size Drinkless n vs gt_shape", height = 8, width = 8)

## ------------------------------------------------------------------ DL4: n vs ASPN

g4 <- expand.grid(aspn = seq(0.03, 0.20, by = 0.001),
                  ft_shape = "constant", gt_shape = "constant",
                  stringsAsFactors = FALSE)
g4 <- grid_n(g4, "independent")
note("DL4: ", nrow(g4), " points, n from ", min(g4$n), " to ", max(g4$n),
     "; at ASPN = ", ASPN, ", n = ", g4$n[which.min(abs(g4$aspn - ASPN))])

ggplot(g4, aes(x = aspn, y = n)) +
  geom_line() +
  xlab(TeX(r'(ASPN)')) + ylab("n (sample size)") +
  theme_bw() + myggfont() +
  theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
save_fig("size Drinkless - n vs ASPNC", height = 6, width = 6)

## ------------------------------------------------------------------ DL6: n vs AA

g6 <- expand.grid(avg_tau = seq(0.5, 1, by = 0.01),
                  ft_shape = "constant", gt_shape = "constant",
                  stringsAsFactors = FALSE)
g6 <- grid_n(g6, "independent")
note("DL6: ", nrow(g6), " points, n from ", min(g6$n), " to ", max(g6$n))

ggplot(g6, aes(x = avg_tau, y = n)) +
  geom_line() +
  xlab(TeX(r'(AA)')) + ylab("n (sample size)") +
  theme_bw() + myggfont() +
  theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
save_fig("size Drinkless 1 - n vs AA", height = 6, width = 6)

## ------------------------------------------------------------------ DL7: n vs theta_tau

g7 <- rbind(
  expand.grid(tau_theta = seq(-0.3, 0.3, by = 0.01),
              tau_shape = c("linear", "sine"),
              ft_theta2 = c(-0.4, 0, 0.4),
              ft_shape = "linear_theta", gt_shape = "quadratic_theta", gt_theta = 0.2,
              avg_tau = 0.7, stringsAsFactors = FALSE),
  expand.grid(tau_theta = seq(-0.3, 0.3, by = 0.01),
              tau_shape = c("linear", "sine"),
              ft_theta2 = 0,
              ft_shape = "constant", gt_shape = "quadratic_theta", gt_theta = 0.2,
              avg_tau = 0.7, stringsAsFactors = FALSE))
g7 <- grid_n(g7, "star")
g7$mee <- ifelse(g7$ft_shape == "constant", "constant",
                 paste0("linear", g7$ft_theta2))
note("DL7: ", nrow(g7), " points, n from ", min(g7$n), " to ", max(g7$n))

ggplot(g7, aes(x = tau_theta, y = n, linetype = tau_shape, colour = mee)) +
  geom_line() +
  scale_colour_discrete(name = "Pattern of MEE(t)",
                        labels = c(TeX(r'(constant)'),
                                   TeX(r'(linear, $\theta_{f2} = -0.4$)'),
                                   TeX(r'(linear, $\theta_{f2} = 0$)'),
                                   TeX(r'(linear, $\theta_{f2} = 0.4$)'))) +
  scale_linetype_discrete(name = TeX(r'(Pattern of $\tau(t)$)'),
                          labels = c("linear", "periodic")) +
  xlab(TeX(r'($\theta_\tau$)')) + ylab("n (sample size)") +
  theme_bw() + myggfont() +
  theme(legend.text.align = 0, plot.title = element_text(hjust = 0.5))
save_fig("size Drinkless 3 - n vs theta_tau", height = 4, width = 7)

## ------------------------------------------------------------------ record

writeLines(log_lines, file.path(OUT, "summary.txt"))
saveRDS(list(DL1 = g1, DL2 = g2, DL3 = g3, DL4 = g4, DL5 = g5, DL6 = g6, DL7 = g7,
             reference = c(constant = n_reference, quadratic_spnc = n2_ref)),
         file.path(OUT, "figure_grids.RDS"))
cat("\nwrote", length(list.files(OUT, pattern = "\\.(pdf|jpg)$")), "figure files to", OUT, "\n")
