# Power / sensitivity check for the disparity-shape null results
#
# The shape analysis concludes "no curvature, no breakpoint, constant capture
# fraction". Those are null results, so the question this script answers is:
# how large a departure from a straight line would these tests have detected?
#
#   1. re-run of the parametric bootstrap null with an independent seed
#      (does the reported p = .35 depend on the seed?)
#   2. power of the maximised-LRT breakpoint test against genuine broken-stick
#      alternatives (slope drop of 25 / 50 / 100 % at 0.40 m)
#   3. minimum detectable effects for the quadratic term and for the slope of
#      the bias ratio on disparity
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/disparity_shape_power_*.csv, disparity_shape_power_log.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

proj_dir    <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
results_dir <- file.path(proj_dir, "results_revision")
log_path    <- file.path(results_dir, "disparity_shape_power_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

cat("Power / sensitivity check for the disparity-shape null results\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds")) %>%
  mutate(participantId = factor(participantId), soundType = factor(soundType)) %>%
  as.data.frame()

breaks <- seq(0.20, 0.60, by = 0.01)

fit_segmented <- function(c_val, data) {
  data$hinge <- pmax(data$stimulusDisparity_m - c_val, 0)
  lmer(signedError_m ~ stimulusDisparity_m + hinge + soundType +
         (1 + stimulusDisparity_m | participantId), data = data, REML = FALSE)
}

lmm_lin  <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                   (1 + stimulusDisparity_m | participantId), data = df, REML = FALSE)
seg_fits <- lapply(breaks, fit_segmented, data = df)
ll_lin   <- as.numeric(logLik(lmm_lin))
lrt_obs  <- max(sapply(seg_fits, function(m) as.numeric(logLik(m)))) - ll_lin
lrt_obs  <- 2 * lrt_obs
cat("observed maximised LRT:", round(lrt_obs, 5), "\n\n")

refit_quiet <- function(model, y) {
  fit <- withCallingHandlers(suppressMessages(refit(model, newresp = y)),
                             warning = function(cond) invokeRestart("muffleWarning"))
  as.numeric(logLik(fit))
}

# maximised LRT over the breakpoint grid for one simulated response vector
max_lrt_one <- function(y) {
  l0  <- refit_quiet(lmm_lin, y)
  lls <- vapply(seg_fits, refit_quiet, numeric(1), y = y)
  c(max_LRT = 2 * (max(lls) - l0), argmax = breaks[which.max(lls)])
}

n_cores <- max(1L, min(4L, parallel::detectCores() - 1L))

# ---------------------------------------------------------------------------
# 1. Null bootstrap with an independent seed
# ---------------------------------------------------------------------------

cat("== 1. Null bootstrap, independent seed (B = 500) ==\n")
B_null <- 500
sim_null <- simulate(lmm_lin, nsim = B_null, seed = 11223344)
t0 <- Sys.time()
null_res <- do.call(rbind, parallel::mclapply(seq_len(B_null),
                                              function(b) max_lrt_one(sim_null[[b]]),
                                              mc.cores = n_cores))
cat("time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
null_max <- null_res[, "max_LRT"]
q_null <- quantile(null_max, c(0.5, 0.9, 0.95, 0.99))
cat("null quantiles 50/90/95/99:", round(q_null, 4), "\n")
cat("bootstrap p:", signif((1 + sum(null_max >= lrt_obs)) / (B_null + 1), 4), "\n")
cat("fraction of argmax at a grid edge:",
    round(mean(null_res[, "argmax"] %in% range(breaks)), 4), "\n\n")
crit95 <- unname(q_null["95%"])

# ---------------------------------------------------------------------------
# 2. Power against real broken-stick alternatives
# ---------------------------------------------------------------------------

cat("== 2. Power of the breakpoint test against true broken-stick data ==\n")
cat("Alternatives have a breakpoint at 0.40 m; the response is simulated with the\n")
cat("random-effect and residual variances estimated from the observed linear LMM.\n")
cat("Rejection uses the 95th percentile of the null bootstrap above (crit =",
    round(crit95, 3), ").\n\n")

seg40 <- fit_segmented(0.40, df)
b_lin <- fixef(lmm_lin)
th    <- getME(lmm_lin, "theta")
sg    <- sigma(lmm_lin)
d_obs <- df$stimulusDisparity_m

# choose the intercept so that the mean simulated signal matches the linear fit
make_beta <- function(s1, s2) {
  hinge_b <- s2 - s1
  target  <- mean(b_lin[["(Intercept)"]] + b_lin[["stimulusDisparity_m"]] * d_obs)
  b0      <- target - mean(s1 * d_obs + hinge_b * pmax(d_obs - 0.40, 0))
  unname(c(b0, s1, hinge_b, b_lin[3:5]))
}

scenarios <- list(
  list(lab = "slope drops 25% at 0.40 m", s1 = 0.40, s2 = 0.30),
  list(lab = "slope drops 50% at 0.40 m", s1 = 0.46, s2 = 0.23),
  list(lab = "slope drops to zero at 0.40 m (full saturation)", s1 = 0.55, s2 = 0.00)
)

B_pow <- 400
pow_tbl <- lapply(seq_along(scenarios), function(k) {
  sc <- scenarios[[k]]
  bb <- make_beta(sc$s1, sc$s2)
  sim <- simulate(seg40, nsim = B_pow, seed = 900000 + k,
                  newparams = list(beta = bb, theta = th, sigma = sg))
  t1 <- Sys.time()
  rr <- do.call(rbind, parallel::mclapply(seq_len(B_pow),
                                          function(b) max_lrt_one(sim[[b]]),
                                          mc.cores = n_cores))
  cat(sc$lab, "| median max LRT", round(median(rr[, "max_LRT"]), 3),
      "| power at crit", round(mean(rr[, "max_LRT"] >= crit95), 3),
      "| time", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 2), "min\n")
  data.frame(scenario = sc$lab, slope_before = sc$s1, slope_after = sc$s2,
             B = B_pow, median_maxLRT = median(rr[, "max_LRT"]),
             power_at_null95 = mean(rr[, "max_LRT"] >= crit95),
             median_argmax_m = median(rr[, "argmax"]),
             mean_signal_drop_cm = 100 * (sc$s1 * 0.70 + (sc$s2 - sc$s1) * 0.30 -
                                            sc$s1 * 0.70) * -1)
}) %>% bind_rows()
cat("\n")
print(pow_tbl, digits = 4, row.names = FALSE)
write.csv(pow_tbl, file.path(results_dir, "disparity_shape_power_breakpoint.csv"),
          row.names = FALSE)

# how big is each alternative in observable terms: signed error at 0.70 m
cat("\nsigned error at 0.70 m implied by each alternative vs the fitted straight line:\n")
lin70 <- 100 * (b_lin[["(Intercept)"]] + b_lin[["stimulusDisparity_m"]] * 0.70)
for (sc in scenarios) {
  bb <- make_beta(sc$s1, sc$s2)
  cat(sprintf("  %-48s %.2f cm (line: %.2f cm, difference %.2f cm)\n", sc$lab,
              100 * (bb[1] + sc$s1 * 0.70 + (sc$s2 - sc$s1) * 0.30), lin70,
              100 * (bb[1] + sc$s1 * 0.70 + (sc$s2 - sc$s1) * 0.30) - lin70))
}
cat("\n")

# ---------------------------------------------------------------------------
# 3. Minimum detectable effects for the other two null tests
# ---------------------------------------------------------------------------

cat("== 3. Minimum detectable effects ==\n")
df$disp2 <- df$stimulusDisparity_m^2
lmm_quad <- lmer(signedError_m ~ stimulusDisparity_m + disp2 + soundType +
                   (1 + stimulusDisparity_m | participantId), data = df, REML = FALSE)
se_q <- summary(lmm_quad)$coefficients["disp2", "Std. Error"]
x1 <- 0.15; x2 <- 0.70; xm <- (x1 + x2) / 2
sag <- function(b2) 100 * b2 * (xm^2 - (x1^2 + x2^2) / 2)
mde_q <- 2.802 * se_q          # 80% power, two-sided alpha = .05
cat("quadratic term SE:", round(se_q, 4), "\n")
cat("MDE (80% power, alpha .05):", round(mde_q, 4),
    " -> mid-range departure from the line:", round(abs(sag(mde_q)), 2), "cm\n")
cat("95% CI bound on the quadratic term implies a departure of up to",
    round(abs(sag(fixef(lmm_quad)["disp2"] - 1.96 * se_q)), 2), "cm\n")
cat("(the fitted linear rise across 0.15-0.70 m is",
    round(100 * fixef(lmm_lin)["stimulusDisparity_m"] * (x2 - x1), 2), "cm)\n\n")

bias_lmm <- lmer(ventriloquistBias ~ stimulusDisparity_m + (1 | participantId), data = df)
se_b <- summary(bias_lmm)$coefficients["stimulusDisparity_m", "Std. Error"]
cat("bias-ratio slope SE:", round(se_b, 4), "\n")
cat("MDE (80% power):", round(2.802 * se_b, 4), "per m ->",
    round(100 * 2.802 * se_b * 0.55, 2), "percentage points across the tested range\n")
ci_b <- confint(bias_lmm, parm = "stimulusDisparity_m", method = "profile")
cat("95% profile CI across the range:", round(100 * ci_b[1, ] * 0.55, 2),
    "percentage points, against a mean capture of",
    round(100 * mean(df$ventriloquistBias), 2), "%\n")
cat("i.e. the capture fraction could fall from",
    round(100 * (mean(df$ventriloquistBias) - 0.5 * ci_b[1, 1] * 0.55), 1), "% to",
    round(100 * (mean(df$ventriloquistBias) + 0.5 * ci_b[1, 1] * 0.55), 1),
    "% across the range and still lie inside the CI\n\n")

mde_tbl <- data.frame(
  test = c("quadratic term on signed error", "slope of bias ratio on disparity"),
  se = c(se_q, se_b),
  mde_80pct = c(2.802 * se_q, 2.802 * se_b),
  mde_in_units = c(paste0(round(abs(sag(mde_q)), 2), " cm mid-range departure"),
                   paste0(round(100 * 2.802 * se_b * 0.55, 2), " pp across the range")))
write.csv(mde_tbl, file.path(results_dir, "disparity_shape_power_mde.csv"), row.names = FALSE)
print(mde_tbl, row.names = FALSE)

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
