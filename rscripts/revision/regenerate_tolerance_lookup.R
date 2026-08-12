# regenerate_tolerance_lookup.R
#
# Regenerates results_revision/tolerance_lookup_nominal.csv under the slope-only
# by-participant covariance. The original table (preserved as
# tolerance_lookup_nominal_fullcov_superseded.csv) computed its percentile-listener
# columns under the full covariance (intercept + slope + correlation); the intercept
# SD (profile CI [0, 0.096]) and the correlation (profile CI [-1, 1]) are not
# identified in these data, which made the 10th/90th-percentile profiles nearly
# coincide at short offsets. The slope-only structure is the one the data support
# (LRT against the full structure reported below).
#
# Primary model: signedError_m ~ stimulusDisparity_m + soundType +
#                (0 + stimulusDisparity_m | participantId), REML.
#
# Columns and CI methods:
#   displacement_avg_*  population-average prediction, marginal over sound type
#                       (design balanced at 186 trials per level, so the marginal
#                       prediction is the equal-weight mean over levels); 95% Wald CI
#                       from the fixed-effect covariance matrix.
#   displacement_p10_*, displacement_p90_*
#                       percentile listener from the fitted slope distribution,
#                       slope = beta +/- qnorm(0.90) * SD_slope, common intercept;
#                       95% CI by parametric bootstrap (bootMer, 2000 replicates).
#   remaining columns   traversed fraction, residual gap, and the empirical
#                       nominal-to-angular crosswalk (trials within +/- 2.5 cm of
#                       each grid point), recomputed from the data as before.
#
# A companion file, tolerance_percentile_definitions.csv, reports the 10th/50th/90th
# percentile listener under both available definitions: (a) normal quantiles of the
# fitted slope distribution and (b) empirical type-7 quantiles of the 31 per-
# participant total slopes (fixed + conditional mode).
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/tolerance_lookup_nominal.csv
#         results_revision/tolerance_percentile_definitions.csv
#         results_revision/regenerate_tolerance_lookup_log.txt
#
# Deterministic: the bootstrap is seeded. Run with R 4.4 (arm64).

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(dplyr)
})

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir  <- file.path(proj_dir, "results_revision")
log_path <- file.path(res_dir, "regenerate_tolerance_lookup_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78))

say("regenerate_tolerance_lookup.R")
say("run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
say("R: ", R.version.string)

N_BOOT <- 2000
Z90    <- qnorm(0.90)   # 1.2816
ctrl   <- lmerControl(optimizer = "bobyqa")

# ---------------------------------------------------------------------------
# 0. Data
# ---------------------------------------------------------------------------
rule("0. DATA")

d <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
say("trials: ", nrow(d), "   participants: ", dplyr::n_distinct(d$participantId))

d <- d %>%
  mutate(disp_cm      = stimulusDisparity_m * 100,
         ang_deg      = angDisp_mean_deg,
         angOnset_deg = angDisp_onset_deg)

# ---------------------------------------------------------------------------
# 1. Primary model: slope-only by-participant covariance
# ---------------------------------------------------------------------------
rule("1. PRIMARY MODEL (SLOPE-ONLY COVARIANCE)")

m <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
            (0 + stimulusDisparity_m | participantId),
          data = d, REML = TRUE, control = ctrl)
print(summary(m))

msgs <- m@optinfo$conv$lme4$messages
say("convergence messages: ", if (length(msgs)) paste(msgs, collapse = " | ") else "none")
say("singular fit: ", isSingular(m))

co   <- summary(m)$coefficients
b1   <- co["stimulusDisparity_m", "Estimate"]
b1se <- co["stimulusDisparity_m", "Std. Error"]
tau1 <- unname(attr(VarCorr(m)$participantId, "stddev")[1])

say("\nfixed slope (per m): ", round(b1, 5), "   SE ", round(b1se, 5),
    "   t(", round(co["stimulusDisparity_m", "df"], 1), ") = ",
    round(co["stimulusDisparity_m", "t value"], 2),
    ", p = ", format.pval(co["stimulusDisparity_m", "Pr(>|t|)"], digits = 3))
say("Wald 95% CI (per m): [", round(b1 - 1.96 * b1se, 4), ", ",
    round(b1 + 1.96 * b1se, 4), "]")

ci_prof <- suppressMessages(confint(m, parm = "stimulusDisparity_m", method = "profile"))
say("profile 95% CI (per m): [", round(ci_prof[1], 4), ", ", round(ci_prof[2], 4), "]")
say("random-slope SD (tau1, per m): ", round(tau1, 5))
say("residual SD: ", round(sigma(m), 5))

# LRT against the full covariance, for the record (ML refits).
say("\nLRT, slope-only vs full covariance (ML):")
m_full_ml <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                    (1 + stimulusDisparity_m | participantId),
                  data = d, REML = FALSE, control = ctrl)
m_slop_ml <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                    (0 + stimulusDisparity_m | participantId),
                  data = d, REML = FALSE, control = ctrl)
print(anova(m_slop_ml, m_full_ml))
say("(2-df test on a boundary; the nominal p-value is conservative)")

# ---------------------------------------------------------------------------
# 2. Predictions
# ---------------------------------------------------------------------------
rule("2. PREDICTIONS")

grid_cm <- seq(15, 70, by = 5)
grid_m  <- grid_cm / 100

# Marginal (equal-weight over sound type) prediction: contrast row for each grid x.
fx     <- fixef(m)
V      <- as.matrix(vcov(m))
st_idx <- grep("^soundType", names(fx))
Cmat <- t(vapply(grid_m, function(x) {
  cc <- setNames(numeric(length(fx)), names(fx))
  cc["(Intercept)"] <- 1
  cc["stimulusDisparity_m"] <- x
  cc[st_idx] <- 1 / (length(st_idx) + 1)   # equal weight incl. reference level
  cc
}, numeric(length(fx))))

avg_m  <- as.vector(Cmat %*% fx)
avg_se <- sqrt(rowSums((Cmat %*% V) * Cmat))
b0_marg <- unname(fx["(Intercept)"]) + mean(c(0, fx[st_idx]))

say("marginal intercept (m): ", round(b0_marg, 5))
say("population-average CI: Wald, +/- 1.96 * SE from the fixed-effect vcov")

# Percentile-listener profiles, definition (a): normal quantiles of the slope
# distribution, intercept common to all listeners under this structure.
pred_pct <- function(b0, b1, tau, z, x) b0 + (b1 + z * tau) * x
p10_m <- pred_pct(b0_marg, b1, tau1, -Z90, grid_m)
p50_m <- pred_pct(b0_marg, b1, tau1,     0, grid_m)   # coincides with the average
p90_m <- pred_pct(b0_marg, b1, tau1,  Z90, grid_m)

# ---------------------------------------------------------------------------
# 3. Parametric bootstrap CIs for the percentile profiles
# ---------------------------------------------------------------------------
rule("3. PARAMETRIC BOOTSTRAP (PERCENTILE PROFILES)")

boot_fn <- function(fit) {
  fxb  <- fixef(fit)
  b0b  <- unname(fxb["(Intercept)"]) + mean(c(0, fxb[grep("^soundType", names(fxb))]))
  b1b  <- unname(fxb["stimulusDisparity_m"])
  taub <- unname(attr(VarCorr(fit)$participantId, "stddev")[1])
  c(pred_pct(b0b, b1b, taub, -Z90, grid_m), pred_pct(b0b, b1b, taub, Z90, grid_m))
}
bb <- bootMer(m, boot_fn, nsim = N_BOOT, seed = 20260805, type = "parametric")
ok <- complete.cases(bb$t)
say(N_BOOT, " replicates, ", sum(!ok), " discarded (failed refit), ", sum(ok), " retained")
if (length(bb$msgs$warning)) {
  say("bootstrap warnings (unique): ", paste(unique(unlist(bb$msgs$warning)), collapse = " | "))
}
tb <- bb$t[ok, , drop = FALSE]
lo <- apply(tb, 2, quantile, 0.025)
hi <- apply(tb, 2, quantile, 0.975)
k  <- length(grid_cm)

# ---------------------------------------------------------------------------
# 4. Lookup table (same column layout as the superseded file)
# ---------------------------------------------------------------------------
rule("4. LOOKUP TABLE")

cross <- bind_rows(lapply(grid_cm, function(g) {
  s <- d[abs(d$disp_cm - g) <= 2.5, ]
  tibble(disparity_cm = g, n_trials = nrow(s),
         ang_mean_deg   = mean(s$ang_deg),
         ang_median_deg = median(s$ang_deg),
         ang_q25_deg    = unname(quantile(s$ang_deg, 0.25)),
         ang_q75_deg    = unname(quantile(s$ang_deg, 0.75)),
         angOnset_median_deg = median(s$angOnset_deg))
}))

tabA <- tibble(
  disparity_cm = grid_cm,
  displacement_avg_cm    = avg_m * 100,
  displacement_avg_lo_cm = (avg_m - 1.96 * avg_se) * 100,
  displacement_avg_hi_cm = (avg_m + 1.96 * avg_se) * 100,
  displacement_p10_cm    = p10_m * 100,
  displacement_p10_lo_cm = lo[1:k] * 100,
  displacement_p10_hi_cm = hi[1:k] * 100,
  displacement_p90_cm    = p90_m * 100,
  displacement_p90_lo_cm = lo[(k + 1):(2 * k)] * 100,
  displacement_p90_hi_cm = hi[(k + 1):(2 * k)] * 100
) %>%
  mutate(
    traversed_avg_pct = 100 * displacement_avg_cm / disparity_cm,
    traversed_p10_pct = 100 * displacement_p10_cm / disparity_cm,
    traversed_p90_pct = 100 * displacement_p90_cm / disparity_cm,
    residual_gap_avg_cm = disparity_cm - displacement_avg_cm,
    residual_gap_p90_cm = disparity_cm - displacement_p90_cm
  ) %>%
  left_join(cross, by = "disparity_cm") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(as.data.frame(tabA))

out_path <- file.path(res_dir, "tolerance_lookup_nominal.csv")
writeLines("# generated by rscripts/revision/regenerate_tolerance_lookup.R (slope-only random-effects structure; supersedes the full-covariance table preserved as tolerance_lookup_nominal_fullcov_superseded.csv)",
           out_path)
suppressWarnings(write.table(tabA, out_path, sep = ",", row.names = FALSE,
                             col.names = TRUE, qmethod = "double", append = TRUE))
say("written: ", out_path)

# ---------------------------------------------------------------------------
# 5. Percentile listeners under both definitions
# ---------------------------------------------------------------------------
rule("5. PERCENTILE DEFINITIONS COMPARED")

# (b) empirical type-7 quantiles of the 31 per-participant total slopes
# (fixed effect + conditional mode).
slopes <- coef(m)$participantId[["stimulusDisparity_m"]]
say("per-participant total slopes, n = ", length(slopes), ":")
say("  min ", round(min(slopes), 4), ", p10 ", round(quantile(slopes, 0.10, type = 7), 4),
    ", median ", round(median(slopes), 4),
    ", p90 ", round(quantile(slopes, 0.90, type = 7), 4), ", max ", round(max(slopes), 4))
say("(conditional modes are shrunk toward the mean, so these quantiles are narrower",
    " than the fitted slope distribution)")

q_emp <- quantile(slopes, c(0.10, 0.50, 0.90), type = 7)
defs <- bind_rows(lapply(seq_along(grid_cm), function(i) {
  x <- grid_m[i]
  tibble(
    disparity_cm = grid_cm[i],
    p10_normal_cm    = pred_pct(b0_marg, b1, tau1, -Z90, x) * 100,
    p50_normal_cm    = pred_pct(b0_marg, b1, tau1,     0, x) * 100,
    p90_normal_cm    = pred_pct(b0_marg, b1, tau1,  Z90, x) * 100,
    p10_empirical_cm = (b0_marg + q_emp[1] * x) * 100,
    p50_empirical_cm = (b0_marg + q_emp[2] * x) * 100,
    p90_empirical_cm = (b0_marg + q_emp[3] * x) * 100
  )
})) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(defs))
write.csv(defs, file.path(res_dir, "tolerance_percentile_definitions.csv"), row.names = FALSE)
say("written: ", file.path(res_dir, "tolerance_percentile_definitions.csv"))

say("\nendpoints, definition (a) normal quantiles: 15 cm ",
    round(defs$p10_normal_cm[1], 1), " / ", round(defs$p90_normal_cm[1], 1),
    " cm; 70 cm ", round(defs$p10_normal_cm[k], 1), " / ",
    round(defs$p90_normal_cm[k], 1), " cm")
say("endpoints, definition (b) empirical quantiles: 15 cm ",
    round(defs$p10_empirical_cm[1], 1), " / ", round(defs$p90_empirical_cm[1], 1),
    " cm; 70 cm ", round(defs$p10_empirical_cm[k], 1), " / ",
    round(defs$p90_empirical_cm[k], 1), " cm")

rule("SESSION INFO")
print(sessionInfo())
