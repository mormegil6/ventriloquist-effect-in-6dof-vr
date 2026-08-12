# ---------------------------------------------------------------------------
# Supporting descriptive statistics and outlier sensitivity for the revision
#
# Inputs (already built; this script never rebuilds them):
#   results_revision/analysis_df_revision.rds        trial-level frame, 744 trials
#   results_revision/trajectory_samples_angular.rds  per-sample head pose + angles
#   results_revision/inter_trial_gaps.csv            run boundaries in session time
#
# Sections:
#   1. Listening distances (near-field question, Reviewer 1)
#   2. Cue eccentricity as a visibility proxy (no gaze tracking, Reviewer 1)
#   3. Inter-trial intervals and trial durations (Reviewer 3)
#   4. Orthogonal error component (bias-free precision proxy)
#   5. Outlier sensitivity of the primary signed-error LMM
#
# All outputs go to results_revision/ as CSV plus a plain-text log.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
})

set.seed(20260805)

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
out_dir  <- file.path(proj_dir, "results_revision")
log_path <- file.path(out_dir, "revision_descriptives_and_sensitivity_log.txt")

log_con <- file(log_path, open = "wt")
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_con)
}
say_df <- function(df, digits = 4, row_names = FALSE) {
  out <- capture.output(print(as.data.frame(df), digits = digits, row.names = row_names))
  for (l in out) say(l)
}

say("Supporting descriptives and sensitivity analyses for the TVCG revision")
say("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  R ", getRversion())
say(strrep("=", 78))

analysis_df <- readRDS(file.path(out_dir, "analysis_df_revision.rds"))
samples     <- readRDS(file.path(out_dir, "trajectory_samples_angular.rds"))
gaps_raw    <- read.csv(file.path(out_dir, "inter_trial_gaps.csv"))

say(sprintf("Trials: %d, participants: %d, samples: %d",
            nrow(analysis_df), n_distinct(analysis_df$participantId), nrow(samples)))

# helper: five-number style summary of a numeric vector
describe <- function(x, label, digits = 4) {
  x <- x[is.finite(x)]
  q <- quantile(x, c(0, 0.25, 0.5, 0.75, 1))
  data.frame(measure = label, n = length(x), mean = mean(x), sd = sd(x),
             min = q[[1]], q25 = q[[2]], median = q[[3]], q75 = q[[4]], max = q[[5]])
}

# ---------------------------------------------------------------------------
# 1. Listening distances
# ---------------------------------------------------------------------------
say("\n", strrep("-", 78))
say("1. LISTENING DISTANCES (head to true source)")
say(strrep("-", 78))

near_field_m <- 0.574  # renderer near-field stage engages below this distance

dist_sample <- describe(samples$dist_source_m, "head-source distance, all samples (m)")
dist_trial_mean <- describe(analysis_df$dist_mean_m, "per-trial time-weighted mean distance (m)")
dist_trial_min  <- describe(analysis_df$dist_min_m, "per-trial minimum distance (m)")
dist_onset      <- describe(analysis_df$dist_at_onset_m, "distance at trial onset (m)")

# time-weighted sample distribution (weights = interval to next sample)
w <- samples$w_time
wq <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  approx(cw, x, xout = p, method = "linear", rule = 2, ties = "ordered")$y
}
dist_tw <- data.frame(
  measure = "head-source distance, time-weighted over samples (m)",
  n = nrow(samples),
  mean = sum(samples$dist_source_m * w) / sum(w),
  sd = NA_real_,
  min = min(samples$dist_source_m),
  q25 = wq(samples$dist_source_m, w, 0.25),
  median = wq(samples$dist_source_m, w, 0.50),
  q75 = wq(samples$dist_source_m, w, 0.75),
  max = max(samples$dist_source_m)
)

dist_tbl <- bind_rows(dist_sample, dist_tw, dist_trial_mean, dist_trial_min, dist_onset)
say_df(dist_tbl)

# fraction of samples inside the two thresholds (unweighted and time-weighted)
frac_05_unw <- mean(samples$dist_source_m < 0.5)
frac_nf_unw <- mean(samples$dist_source_m < near_field_m)
frac_05_tw  <- sum(w * (samples$dist_source_m < 0.5)) / sum(w)
frac_nf_tw  <- sum(w * (samples$dist_source_m < near_field_m)) / sum(w)

say(sprintf("\nSamples within 0.5 m:    %.4f (unweighted), %.4f (time-weighted), n = %d of %d",
            frac_05_unw, frac_05_tw, sum(samples$dist_source_m < 0.5), nrow(samples)))
say(sprintf("Samples within %.3f m: %.4f (unweighted), %.4f (time-weighted), n = %d of %d",
            near_field_m, frac_nf_unw, frac_nf_tw,
            sum(samples$dist_source_m < near_field_m), nrow(samples)))

# per-participant minimum distance
pp_min <- samples %>%
  group_by(participantId) %>%
  summarise(min_dist_m = min(dist_source_m),
            frac_under_0p5 = mean(dist_source_m < 0.5),
            frac_under_nf = mean(dist_source_m < near_field_m),
            .groups = "drop")
say("\nPer-participant minimum head-source distance (m):")
say_df(describe(pp_min$min_dist_m, "per-participant minimum distance (m)"))
say(sprintf("All %d participants came within 0.5 m at some point: %s",
            nrow(pp_min), all(pp_min$min_dist_m < 0.5)))
say(sprintf("All participants came within 0.10 m at some point: %s (n = %d)",
            all(pp_min$min_dist_m < 0.10), sum(pp_min$min_dist_m < 0.10)))
say("\nPer-participant fraction of samples within 0.5 m:")
say_df(describe(pp_min$frac_under_0p5, "per-participant fraction of samples < 0.5 m"))

# trial-level incidence
trial_any_05 <- mean(analysis_df$dist_min_m < 0.5)
trial_any_nf <- mean(analysis_df$dist_min_m < near_field_m)
trial_any_01 <- mean(analysis_df$dist_min_m < 0.10)
say(sprintf("\nTrials with any sample within 0.5 m:    %.4f (%d of %d)",
            trial_any_05, sum(analysis_df$dist_min_m < 0.5), nrow(analysis_df)))
say(sprintf("Trials with any sample within %.3f m: %.4f (%d of %d)",
            near_field_m, trial_any_nf, sum(analysis_df$dist_min_m < near_field_m), nrow(analysis_df)))
say(sprintf("Trials with any sample within 0.10 m:   %.4f (%d of %d)",
            trial_any_01, sum(analysis_df$dist_min_m < 0.10), nrow(analysis_df)))
say(sprintf("Trials that never entered 0.5 m: %d", sum(analysis_df$dist_min_m >= 0.5)))

# upper tail: the attenuation curve mutes sources beyond 3 m
far <- samples$dist_source_m > 3
say(sprintf("Samples beyond the 3 m maximum rendering distance: %d of %d (%.4f%%), from %d participants in %d trials; max distance %.3f m",
            sum(far), nrow(samples), 100 * mean(far),
            n_distinct(samples$participantId[far]),
            nrow(distinct(samples[far, c("participantId", "trialSequenceNum")])),
            max(samples$dist_source_m)))

# per-trial fraction of time inside each threshold (columns from the build script)
say("\nPer-trial fraction of samples inside each threshold:")
say_df(bind_rows(
  describe(analysis_df$frac_samples_under_0p5m, "per-trial fraction of samples < 0.5 m"),
  describe(analysis_df$frac_samples_under_0p574m, "per-trial fraction of samples < 0.574 m")
))

write.csv(dist_tbl, file.path(out_dir, "descriptives_listening_distance.csv"), row.names = FALSE)
write.csv(pp_min, file.path(out_dir, "descriptives_distance_per_participant.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Cue eccentricity as a visibility proxy
# ---------------------------------------------------------------------------
say("\n", strrep("-", 78))
say("2. CUE ECCENTRICITY RELATIVE TO HEAD FORWARD (visibility proxy)")
say(strrep("-", 78))

ecc_tbl <- bind_rows(
  describe(samples$cueEcc_deg, "cue eccentricity, all samples (deg)"),
  describe(analysis_df$cueEcc_mean_deg, "per-trial time-weighted mean cue eccentricity (deg)"),
  describe(analysis_df$cueEcc_median_deg, "per-trial median cue eccentricity (deg)"),
  describe(analysis_df$cueEcc_min_deg, "per-trial minimum cue eccentricity (deg)"),
  describe(analysis_df$cueEcc_onset_deg, "cue eccentricity at trial onset (deg)"),
  describe(analysis_df$frac_time_cue_within_55deg, "per-trial fraction of time cue within 55 deg"),
  describe(samples$srcEcc_deg, "source eccentricity, all samples (deg)"),
  describe(analysis_df$frac_time_src_within_55deg, "per-trial fraction of time source within 55 deg")
)
say_df(ecc_tbl)

frac_ecc_55_samples <- mean(samples$cueEcc_deg < 55)
frac_ecc_55_tw <- sum(w * (samples$cueEcc_deg < 55)) / sum(w)
say(sprintf("\nSamples with cue within 55 deg of head forward: %.4f (unweighted), %.4f (time-weighted)",
            frac_ecc_55_samples, frac_ecc_55_tw))
say(sprintf("Samples with cue within 30 deg: %.4f; within 90 deg: %.4f",
            mean(samples$cueEcc_deg < 30), mean(samples$cueEcc_deg < 90)))
say(sprintf("Trials with cue never inside 55 deg: %d; trials with cue inside 55 deg for >90%% of time: %d",
            sum(analysis_df$frac_time_cue_within_55deg == 0),
            sum(analysis_df$frac_time_cue_within_55deg > 0.9)))

# correlations of cue-looking with bias measures
cue_dat <- analysis_df %>%
  mutate(signedError_cm = signedError_m * 100,
         bias_pct = ventriloquistBias * 100)

cor_rows <- function(x, y, xlab, ylab) {
  p <- cor.test(x, y, method = "pearson")
  s <- suppressWarnings(cor.test(x, y, method = "spearman"))
  data.frame(x = xlab, y = ylab, n = sum(complete.cases(x, y)),
             pearson_r = unname(p$estimate), pearson_ci_lo = p$conf.int[1],
             pearson_ci_hi = p$conf.int[2], pearson_t = unname(p$statistic),
             pearson_df = unname(p$parameter), pearson_p = p$p.value,
             spearman_rho = unname(s$estimate), spearman_p = s$p.value)
}

cor_tbl <- bind_rows(
  cor_rows(cue_dat$frac_time_cue_within_55deg, cue_dat$signedError_cm,
           "frac_time_cue_within_55deg", "signedError_cm"),
  cor_rows(cue_dat$frac_time_cue_within_55deg, cue_dat$bias_pct,
           "frac_time_cue_within_55deg", "ventriloquistBias_pct"),
  cor_rows(cue_dat$cueEcc_mean_deg, cue_dat$signedError_cm,
           "cueEcc_mean_deg", "signedError_cm"),
  cor_rows(cue_dat$cueEcc_mean_deg, cue_dat$bias_pct,
           "cueEcc_mean_deg", "ventriloquistBias_pct"),
  cor_rows(cue_dat$frac_time_cue_within_55deg, cue_dat$participantError_cm,
           "frac_time_cue_within_55deg", "participantError_cm")
)
say("\nTrial-level correlations (pooled across participants):")
say_df(cor_tbl)

# within/between decomposition: person-mean centring separates the two effects
cue_dat <- cue_dat %>%
  group_by(participantId) %>%
  mutate(frac_pm = mean(frac_time_cue_within_55deg)) %>%
  ungroup() %>%
  mutate(frac_within = frac_time_cue_within_55deg - frac_pm,
         frac_between = frac_pm - mean(frac_pm))

m_cue_signed <- lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m +
                       soundType + (1 | participantId),
                     data = cue_dat, REML = TRUE,
                     control = lmerControl(optimizer = "bobyqa"))
m_cue_bias <- lmer(bias_pct ~ frac_within + frac_between + soundType + (1 | participantId),
                   data = cue_dat, REML = TRUE,
                   control = lmerControl(optimizer = "bobyqa"))

say("\nLMM: signed error (cm) ~ within- and between-participant cue-looking + disparity + soundType")
say_df(coef(summary(m_cue_signed)), digits = 5, row_names = TRUE)
say("\nLMM: bias proportion (%) ~ within- and between-participant cue-looking + soundType")
say_df(coef(summary(m_cue_bias)), digits = 5, row_names = TRUE)

# rescale the cue-looking slopes to a 10-percentage-point increase, with Wald CIs
scale_slope <- function(fit, term, unit = 0.1, label) {
  cf <- coef(summary(fit))[term, ]
  ci <- confint(fit, parm = term, method = "Wald")
  data.frame(model = label, term = term,
             est_per_0.1 = cf[["Estimate"]] * unit,
             se_per_0.1 = cf[["Std. Error"]] * unit,
             ci_lo_per_0.1 = ci[1, 1] * unit, ci_hi_per_0.1 = ci[1, 2] * unit,
             df = cf[["df"]], t = cf[["t value"]], p = cf[["Pr(>|t|)"]])
}
cue_slopes <- bind_rows(
  scale_slope(m_cue_signed, "frac_within", label = "signed error (cm)"),
  scale_slope(m_cue_signed, "frac_between", label = "signed error (cm)"),
  scale_slope(m_cue_bias, "frac_within", label = "bias proportion (%)"),
  scale_slope(m_cue_bias, "frac_between", label = "bias proportion (%)")
)
say("\nCue-looking slopes rescaled to a 10-percentage-point increase in time-in-view:")
say_df(cue_slopes, digits = 4)

say(sprintf("\nCue-looking vs stimulus disparity (possible confound): r = %.4f, p = %.4f",
            cor(cue_dat$frac_time_cue_within_55deg, cue_dat$stimulusDisparity_m),
            cor.test(cue_dat$frac_time_cue_within_55deg, cue_dat$stimulusDisparity_m)$p.value))
say(sprintf("Cue-looking vs trial duration: r = %.4f; vs mean head-source distance: r = %.4f",
            cor(cue_dat$frac_time_cue_within_55deg, cue_dat$responseTime_s),
            cor(cue_dat$frac_time_cue_within_55deg, cue_dat$dist_mean_m)))

# per-participant correlation summary (within-participant consistency)
pp_cor <- cue_dat %>%
  group_by(participantId) %>%
  summarise(r_signed = cor(frac_time_cue_within_55deg, signedError_cm),
            r_bias = cor(frac_time_cue_within_55deg, bias_pct),
            .groups = "drop")
say("\nPer-participant correlations of cue-looking with signed error / bias:")
say_df(bind_rows(describe(pp_cor$r_signed, "per-participant r (cue-looking, signed error)"),
                 describe(pp_cor$r_bias, "per-participant r (cue-looking, bias)")))
say(sprintf("Participants with positive r (signed error): %d of %d",
            sum(pp_cor$r_signed > 0), nrow(pp_cor)))

write.csv(ecc_tbl, file.path(out_dir, "descriptives_cue_visibility.csv"), row.names = FALSE)
write.csv(cor_tbl, file.path(out_dir, "descriptives_cue_visibility_correlations.csv"), row.names = FALSE)
write.csv(cue_slopes, file.path(out_dir, "cue_visibility_lmm_slopes.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Inter-trial intervals and trial durations
# ---------------------------------------------------------------------------
say("\n", strrep("-", 78))
say("3. INTER-TRIAL INTERVALS AND TRIAL DURATIONS")
say(strrep("-", 78))

gaps <- gaps_raw$gap_to_next_s[is.finite(gaps_raw$gap_to_next_s)]
gap_tbl <- describe(gaps, "inter-trial gap, run-ordered (s)")
say_df(gap_tbl, digits = 5)
say(sprintf("Gaps: n = %d across %d sessions; IQR width %.4f s",
            length(gaps), n_distinct(gaps_raw$participantId), IQR(gaps)))
say(sprintf("Gaps above 0.2 s: %d; above 1 s: %d; above 5 s: %d",
            sum(gaps > 0.2), sum(gaps > 1), sum(gaps > 5)))
# Tukey outlier rule on the gap distribution
gq <- quantile(gaps, c(0.25, 0.75))
gap_fence <- gq[2] + 1.5 * IQR(gaps)
say(sprintf("Tukey upper fence %.4f s; gaps above it: %d (max %.4f s)",
            gap_fence, sum(gaps > gap_fence), max(gaps)))
say("\nLargest five inter-trial gaps:")
say_df(gaps_raw %>% filter(is.finite(gap_to_next_s)) %>%
         arrange(desc(gap_to_next_s)) %>% head(5) %>%
         select(participantId, level, run_id, run_end, gap_to_next_s), digits = 6)

rt_tbl <- describe(analysis_df$responseTime_s, "trial duration responseTime_s (s)")
say("\nTrial duration:")
say_df(rt_tbl)
say(sprintf("IQR width %.3f s; 5th/95th percentile %.3f / %.3f s",
            IQR(analysis_df$responseTime_s),
            quantile(analysis_df$responseTime_s, 0.05),
            quantile(analysis_df$responseTime_s, 0.95)))
say(sprintf("Trials shorter than 3 s: %d; shorter than 5 s: %d; longer than 120 s: %d",
            sum(analysis_df$responseTime_s < 3), sum(analysis_df$responseTime_s < 5),
            sum(analysis_df$responseTime_s > 120)))

# total session duration implied by trial durations plus gaps
sess <- analysis_df %>%
  group_by(participantId) %>%
  summarise(n_trials = n(), total_trial_time_s = sum(responseTime_s), .groups = "drop") %>%
  left_join(gaps_raw %>% group_by(participantId) %>%
              summarise(session_span_s = max(run_end) - min(run_start), .groups = "drop"),
            by = "participantId")
say("\nPer-participant summed trial time and session span (s):")
say_df(bind_rows(describe(sess$total_trial_time_s, "summed trial time per session (s)"),
                 describe(sess$session_span_s, "session span, first to last sample (s)")))
say(sprintf("Median session span %.1f s (%.1f min); trial time accounts for %.2f%% of the span (median)",
            median(sess$session_span_s), median(sess$session_span_s) / 60,
            100 * median(sess$total_trial_time_s / sess$session_span_s)))

write.csv(bind_rows(gap_tbl, rt_tbl),
          file.path(out_dir, "descriptives_intertrial_timing.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Orthogonal error component
# ---------------------------------------------------------------------------
say("\n", strrep("-", 78))
say("4. ORTHOGONAL ERROR COMPONENT")
say(strrep("-", 78))

orth <- analysis_df %>%
  mutate(orth_cm = orthogonalError_m * 100,
         # the form used in the original pipeline, with its NaN-to-zero fallback
         orth_sqrt_raw = sqrt(participantError_m^2 - signedError_m^2),
         orth_sqrt_imputed = ifelse(is.nan(orth_sqrt_raw), 0, orth_sqrt_raw) * 100,
         ratio_orth_abs_signed = orthogonalError_m / abs(signedError_m),
         ratio_orth_error = orthogonalError_m / participantError_m)

n_nan <- sum(is.nan(orth$orth_sqrt_raw))
say(sprintf("Trials where sqrt(err^2 - signed^2) is NaN (i.e. where the original zero-imputation bit): %d of %d",
            n_nan, nrow(orth)))
say(sprintf("Max absolute difference, vector form vs sqrt form: %.3e cm",
            max(abs(orth$orth_cm - orth$orth_sqrt_imputed))))

orth_tbl <- bind_rows(
  describe(orth$orth_cm, "orthogonal error (cm)"),
  describe(orth$participantError_cm, "unsigned error (cm)"),
  describe(orth$signedError_m * 100, "signed error (cm)"),
  describe(abs(orth$signedError_m) * 100, "absolute signed error (cm)"),
  describe(orth$ratio_orth_error, "orthogonal / unsigned error, per trial"),
  describe(orth$ratio_orth_abs_signed, "orthogonal / |signed| error, per trial")
)
say_df(orth_tbl)

say(sprintf("\nRatio of means: mean(orthogonal) / mean(signed) = %.4f / %.4f = %.4f",
            mean(orth$orth_cm), mean(orth$signedError_m * 100),
            mean(orth$orth_cm) / mean(orth$signedError_m * 100)))
say(sprintf("Ratio of means: mean(orthogonal) / mean(|signed|) = %.4f / %.4f = %.4f",
            mean(orth$orth_cm), mean(abs(orth$signedError_m) * 100),
            mean(orth$orth_cm) / mean(abs(orth$signedError_m) * 100)))
say(sprintf("Ratio of means: mean(orthogonal) / mean(unsigned error) = %.4f",
            mean(orth$orth_cm) / mean(orth$participantError_cm)))
say(sprintf("Trials where orthogonal exceeds |signed|: %d (%.1f%%)",
            sum(orth$orthogonalError_m > abs(orth$signedError_m)),
            100 * mean(orth$orthogonalError_m > abs(orth$signedError_m))))

pp_orth <- orth %>% group_by(participantId) %>%
  summarise(mean_orth_cm = mean(orth_cm), .groups = "drop")
say("\nPer-participant mean orthogonal error (cm):")
say_df(describe(pp_orth$mean_orth_cm, "per-participant mean orthogonal error (cm)"))

write.csv(orth_tbl, file.path(out_dir, "descriptives_orthogonal_error.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Outlier sensitivity of the primary signed-error LMM
# ---------------------------------------------------------------------------
say("\n", strrep("-", 78))
say("5. OUTLIER SENSITIVITY OF THE PRIMARY SIGNED-ERROR LMM")
say(strrep("-", 78))

# extreme cases named explicitly
say("\nFive shortest trials:")
say_df(analysis_df %>% arrange(responseTime_s) %>% head(5) %>%
         select(participantId, trialSequenceNum, soundType, stimulusDisparity_m,
                responseTime_s, participantError_cm, signedError_m, n_samples, total_path_length),
       digits = 5)
say("\nFive largest unsigned errors:")
say_df(analysis_df %>% arrange(desc(participantError_cm)) %>% head(5) %>%
         select(participantId, trialSequenceNum, soundType, stimulusDisparity_m,
                responseTime_s, participantError_cm, signedError_m, ventriloquistBias),
       digits = 5)
say("\nFive most negative signed errors:")
say_df(analysis_df %>% arrange(signedError_m) %>% head(5) %>%
         select(participantId, trialSequenceNum, soundType, stimulusDisparity_m,
                responseTime_s, participantError_cm, signedError_m, ventriloquistBias),
       digits = 5)

err_mean <- mean(analysis_df$participantError_cm)
err_sd   <- sd(analysis_df$participantError_cm)
err_cut  <- err_mean + 3 * err_sd
sgn_mean <- mean(analysis_df$signedError_m)
sgn_sd   <- sd(analysis_df$signedError_m)
say(sprintf("\nUnsigned error: mean %.4f cm, SD %.4f cm, mean + 3 SD = %.4f cm; trials above: %d",
            err_mean, err_sd, err_cut, sum(analysis_df$participantError_cm > err_cut)))
say(sprintf("Signed error: mean %.4f m, SD %.4f m, +/- 3 SD = [%.4f, %.4f] m; trials outside: %d",
            sgn_mean, sgn_sd, sgn_mean - 3 * sgn_sd, sgn_mean + 3 * sgn_sd,
            sum(abs(analysis_df$signedError_m - sgn_mean) > 3 * sgn_sd)))
say(sprintf("Trials with responseTime_s < 3 s: %d", sum(analysis_df$responseTime_s < 3)))

# primary model, refitted on each subset, warnings captured verbatim
fit_primary <- function(dat, label) {
  msgs <- character(0)
  fit <- withCallingHandlers(
    lmer(signedError_m ~ stimulusDisparity_m + soundType +
           (1 + stimulusDisparity_m | participantId),
         data = dat, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      msgs <<- c(msgs, trimws(conditionMessage(m)))
      invokeRestart("muffleMessage")
    })
  cf <- coef(summary(fit))["stimulusDisparity_m", ]
  ci <- confint(fit, parm = "stimulusDisparity_m", method = "Wald")
  vc <- as.data.frame(VarCorr(fit))
  list(
    summary = data.frame(
      subset = label,
      n_trials = nrow(dat),
      n_participants = n_distinct(dat$participantId),
      n_dropped = nrow(analysis_df) - nrow(dat),
      intercept_m = fixef(fit)[["(Intercept)"]],
      slope_m_per_m = cf[["Estimate"]],
      slope_se = cf[["Std. Error"]],
      slope_df = cf[["df"]],
      slope_t = cf[["t value"]],
      slope_p = cf[["Pr(>|t|)"]],
      slope_ci_lo = ci[1, 1],
      slope_ci_hi = ci[1, 2],
      sd_intercept = vc$sdcor[vc$grp == "participantId" & vc$var1 == "(Intercept)" & is.na(vc$var2)],
      sd_slope = vc$sdcor[vc$grp == "participantId" & vc$var1 == "stimulusDisparity_m" & is.na(vc$var2)],
      sigma = sigma(fit),
      is_singular = isSingular(fit),
      warnings = if (length(msgs)) paste(msgs, collapse = " | ") else "none",
      stringsAsFactors = FALSE),
    fit = fit)
}

subsets <- list(
  list(label = "full (published primary model)", dat = analysis_df),
  list(label = "exclude responseTime_s < 3 s",
       dat = filter(analysis_df, responseTime_s >= 3)),
  list(label = "exclude unsigned error > mean + 3 SD",
       dat = filter(analysis_df, participantError_cm <= err_cut)),
  list(label = "exclude signed error beyond +/- 3 SD",
       dat = filter(analysis_df, abs(signedError_m - sgn_mean) <= 3 * sgn_sd)),
  list(label = "exclude both (RT < 3 s and unsigned error > mean + 3 SD)",
       dat = filter(analysis_df, responseTime_s >= 3, participantError_cm <= err_cut))
)

fits <- lapply(subsets, function(s) fit_primary(s$dat, s$label))
sens_tbl <- bind_rows(lapply(fits, `[[`, "summary"))

base_slope <- sens_tbl$slope_m_per_m[1]
sens_tbl <- sens_tbl %>%
  mutate(slope_change_pct = 100 * (slope_m_per_m - base_slope) / base_slope,
         slope_change_in_base_se = (slope_m_per_m - base_slope) / sens_tbl$slope_se[1])

say("\nPrimary model: signedError_m ~ stimulusDisparity_m + soundType + (1 + stimulusDisparity_m | participantId), REML, bobyqa")
say_df(sens_tbl %>% select(subset, n_trials, n_dropped, slope_m_per_m, slope_se, slope_df,
                           slope_t, slope_p, slope_ci_lo, slope_ci_hi,
                           slope_change_pct, slope_change_in_base_se), digits = 5)
say("\nRandom-effect and convergence detail:")
say_df(sens_tbl %>% select(subset, sd_intercept, sd_slope, sigma, is_singular, warnings), digits = 4)

say("\nFull fixed-effect tables per subset:")
for (i in seq_along(fits)) {
  say("\n[", sens_tbl$subset[i], "]  n = ", sens_tbl$n_trials[i])
  say_df(coef(summary(fits[[i]]$fit)), digits = 5, row_names = TRUE)
}

write.csv(sens_tbl, file.path(out_dir, "outlier_sensitivity_signed_lmm.csv"), row.names = FALSE)

say("\n", strrep("=", 78))
say("Files written to ", out_dir, ":")
for (f in c("descriptives_listening_distance.csv",
            "descriptives_distance_per_participant.csv",
            "descriptives_cue_visibility.csv",
            "descriptives_cue_visibility_correlations.csv",
            "cue_visibility_lmm_slopes.csv",
            "descriptives_intertrial_timing.csv",
            "descriptives_orthogonal_error.csv",
            "outlier_sensitivity_signed_lmm.csv",
            basename(log_path))) {
  say("  ", f)
}

close(log_con)
