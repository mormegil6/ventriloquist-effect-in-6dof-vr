#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Did listeners have access to binaural (HRTF) elevation and front/back cues?
#
# The released build ships the atmoky trueSpatial object renderer for Windows
# and macOS but not for Android/arm64, so a standalone Quest build could not
# have loaded it. Without the object renderer the audio-objects bus falls back
# to Wwise built-in panning, which for a stereo endpoint encodes interaural
# level differences only: no spectral elevation cues, and a front/back mirror
# ambiguity. This script asks the data whether listeners behaved as if they
# had elevation information.
#
# Identification. The visual cue was displaced from the true source by a
# randomised vector that is uncorrelated with source position (|r| < .04 on
# every axis), so a response regressed jointly on source position and cue
# position separates the auditory weight (source coefficient) from the visual
# weight (cue coefficient). The horizontal axes, where amplitude panning does
# carry information, act as an internal benchmark for the vertical axis.
#
# Coordinate conventions: Unity world frame, left-handed, y up, metres.
# Quaternions are stored (x, y, z, w); head forward is +Z, up +Y, right +X.
#
# Run with R 4.4 (arm64).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
  library(lme4)
  library(lmerTest)
})

project_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
data_dir    <- file.path(project_dir, "data")
out_dir     <- file.path(project_dir, "results_revision")

log_path <- file.path(out_dir, "hrtf_rendering_evidence_log.txt")
log_con  <- file(log_path, open = "wt")
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_con)
}
say_df <- function(x, digits = 4) {
  out <- capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE))
  cat(out, sep = "\n"); cat("\n")
  cat(out, sep = "\n", file = log_con); cat("\n", file = log_con)
}

say("=== Auditory elevation and front/back evidence for the rendering path ===")
say("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

analysis_df <- readRDS(file.path(out_dir, "analysis_df_revision.rds"))
say("\nTrials: ", nrow(analysis_df), "  participants: ",
    dplyr::n_distinct(analysis_df$participantId))

# ---------------------------------------------------------------------------
# Design check: is the cue displacement orthogonal to source position?
# ---------------------------------------------------------------------------

design_check <- tibble(
  axis = c("x", "y", "z"),
  r_source_vs_cueOffset = c(
    cor(analysis_df$sound_x, analysis_df$soundToFlash_x),
    cor(analysis_df$sound_y, analysis_df$soundToFlash_y),
    cor(analysis_df$sound_z, analysis_df$soundToFlash_z)),
  r_source_vs_cuePosition = c(
    cor(analysis_df$sound_x, analysis_df$flash_x),
    cor(analysis_df$sound_y, analysis_df$flash_y),
    cor(analysis_df$sound_z, analysis_df$flash_z))
)
say("\n--- Design check: source position vs cue displacement ---")
say_df(design_check)
write_csv(design_check, file.path(out_dir, "hrtf_design_orthogonality.csv"))

# ---------------------------------------------------------------------------
# Head-relative geometry from the per-sample trajectories
# ---------------------------------------------------------------------------

json_files <- sort(list.files(data_dir, pattern = "\\.json$", full.names = TRUE))
csv_files  <- sort(list.files(data_dir, pattern = "\\.csv$",  full.names = TRUE))
stopifnot(length(json_files) == length(csv_files))
pids <- sprintf("P%03d", seq_along(json_files))

traj <- map2_dfr(csv_files, pids, function(path, pid) {
  read_csv(path, show_col_types = FALSE, progress = FALSE) %>% mutate(participantId = pid)
})
say("\nTrajectory samples: ", nrow(traj))

trial_static <- analysis_df %>%
  select(participantId, trialSequenceNum,
         sound_x, sound_y, sound_z, flash_x, flash_y, flash_z,
         response_x, response_y, response_z)

# Head-frame basis vectors: columns of the rotation matrix for quaternion
# (x, y, z, w). The forward column reproduces the validated definition used in
# build_analysis_df_revision.R.
quat_basis <- function(x, y, z, w) {
  list(
    right_x = 1 - 2 * (y^2 + z^2), right_y = 2 * (x * y + w * z), right_z = 2 * (x * z - w * y),
    up_x    = 2 * (x * y - w * z), up_y    = 1 - 2 * (x^2 + z^2), up_z    = 2 * (y * z + w * x),
    fwd_x   = 2 * (x * z + w * y), fwd_y   = 2 * (y * z - w * x), fwd_z   = 1 - 2 * (x^2 + y^2)
  )
}

samples <- traj %>%
  select(participantId, trialSequenceNum = level, timestamp,
         px, py, pz, rx, ry, rz, rw, ex, ey, ez) %>%
  inner_join(trial_static, by = c("participantId", "trialSequenceNum")) %>%
  arrange(participantId, timestamp) %>%
  group_by(participantId) %>%
  mutate(run_id = cumsum(trialSequenceNum != lag(trialSequenceNum,
                                                 default = first(trialSequenceNum))) + 1L) %>%
  group_by(participantId, run_id) %>%
  mutate(dt_next = lead(timestamp) - timestamp) %>%
  group_by(participantId, trialSequenceNum) %>%
  mutate(dt_med = median(dt_next, na.rm = TRUE),
         w_time = ifelse(is.na(dt_next), ifelse(is.na(dt_med), 0.111, dt_med), dt_next),
         sample_idx = rank(timestamp, ties.method = "first")) %>%
  ungroup()

basis <- quat_basis(samples$rx, samples$ry, samples$rz, samples$rw)
samples <- bind_cols(samples, as_tibble(basis))

# Sanity check of the basis against the Euler-derived forward vector.
pitch <- samples$ex * pi / 180; yaw <- samples$ey * pi / 180
fwd_check <- acos(pmin(1, pmax(-1,
  samples$fwd_x * (sin(yaw) * cos(pitch)) +
  samples$fwd_y * (-sin(pitch)) +
  samples$fwd_z * (cos(yaw) * cos(pitch))))) * 180 / pi
say(sprintf("Head-frame basis check (quaternion vs Euler forward): median %.3f deg, p99 %.3f deg",
            median(fwd_check), quantile(fwd_check, 0.99)))
stopifnot(median(fwd_check) < 3)

# Unit direction from the head to a world point, expressed in head coordinates.
head_local_unit <- function(df, tx, ty, tz, prefix) {
  vx <- df[[tx]] - df$px; vy <- df[[ty]] - df$py; vz <- df[[tz]] - df$pz
  n  <- sqrt(vx^2 + vy^2 + vz^2)
  lx <- (vx * df$right_x + vy * df$right_y + vz * df$right_z) / n
  ly <- (vx * df$up_x    + vy * df$up_y    + vz * df$up_z)    / n
  lz <- (vx * df$fwd_x   + vy * df$fwd_y   + vz * df$fwd_z)   / n
  out <- tibble(lx, ly, lz,
                elev = asin(pmin(1, pmax(-1, ly))) * 180 / pi,
                azim = atan2(lx, lz) * 180 / pi,
                dist = n)
  setNames(out, paste0(prefix, c("_ux", "_uy", "_uz", "_elev", "_azim", "_dist"))
  )
}

samples <- bind_cols(
  samples,
  head_local_unit(samples, "sound_x", "sound_y", "sound_z", "src"),
  head_local_unit(samples, "flash_x", "flash_y", "flash_z", "cue"),
  head_local_unit(samples, "response_x", "response_y", "response_z", "resp")
)

wmean <- function(x, w) sum(x * w) / sum(w)

# Three head-referenced frames: trial onset, time-weighted over the trial, and
# the final sample (the pose from which the response was confirmed).
geom <- samples %>%
  group_by(participantId, trialSequenceNum) %>%
  arrange(sample_idx, .by_group = TRUE) %>%
  summarise(
    across(c(src_uy, cue_uy, resp_uy, src_ux, cue_ux, resp_ux,
             src_uz, cue_uz, resp_uz, src_elev, cue_elev, resp_elev),
           list(onset = ~ first(.x), final = ~ last(.x),
                tw = ~ wmean(.x, w_time)),
           .names = "{.col}_{.fn}"),
    src_azim_onset = first(src_azim), cue_azim_onset = first(cue_azim),
    resp_azim_onset = first(resp_azim),
    head_y_onset = first(py), head_y_tw = wmean(py, w_time),
    src_dist_onset = first(src_dist), src_dist_tw = wmean(src_dist, w_time),
    # head tilt available for a projection-based elevation strategy
    roll_range_deg  = diff(range(((ez + 180) %% 360) - 180)),
    roll_sd_deg     = sd(((ez + 180) %% 360) - 180),
    pitch_range_deg = diff(range(((ex + 180) %% 360) - 180)),
    pitch_sd_deg    = sd(((ex + 180) %% 360) - 180),
    # front/back exposure of the source over the trial
    frac_time_src_front = wmean(as.numeric(src_uz > 0), w_time),
    frac_time_cue_front = wmean(as.numeric(cue_uz > 0), w_time),
    .groups = "drop"
  )

dat <- analysis_df %>% left_join(geom, by = c("participantId", "trialSequenceNum"))
say("Trials with head-frame geometry: ", sum(!is.na(dat$src_uy_onset)))
saveRDS(dat, file.path(out_dir, "hrtf_geometry_trials.rds"))

# ---------------------------------------------------------------------------
# Model helpers
# ---------------------------------------------------------------------------

# Nakagawa marginal R2 for a Gaussian LMM.
r2_marginal <- function(m) {
  vf <- var(as.vector(model.matrix(m) %*% fixef(m)))
  vr <- sum(unlist(lapply(VarCorr(m), function(v) v[1])))
  vf / (vf + vr + sigma(m)^2)
}

# response ~ source + cue + (1 | participant); the source coefficient is the
# auditory weight, the cue coefficient the visual weight.
cue_weight_model <- function(df, resp, src, cue, label) {
  d <- df %>% select(participantId, y = all_of(resp), s = all_of(src), c = all_of(cue)) %>%
    tidyr::drop_na()
  m  <- lmer(y ~ s + c + (1 | participantId), data = d, REML = TRUE)
  m0 <- lmer(y ~ c + (1 | participantId), data = d, REML = TRUE)
  ci <- suppressMessages(confint(m, parm = c("s", "c"), method = "Wald"))
  co <- summary(m)$coefficients
  lrt <- anova(refitML(m), refitML(m0))
  # standardised (partial) effect size for the source term
  ds <- d %>% mutate(across(c(y, s, c), ~ as.vector(scale(.x))))
  mz <- lmer(y ~ s + c + (1 | participantId), data = ds, REML = TRUE)
  tibble(
    model = label, n = nrow(d),
    b_source = co["s", "Estimate"], se_source = co["s", "Std. Error"],
    ci_lo_source = ci["s", 1], ci_hi_source = ci["s", 2],
    t_source = co["s", "t value"], df_source = co["s", "df"],
    p_source = co["s", "Pr(>|t|)"],
    beta_std_source = fixef(mz)["s"],
    b_cue = co["c", "Estimate"], se_cue = co["c", "Std. Error"],
    ci_lo_cue = ci["c", 1], ci_hi_cue = ci["c", 2],
    p_cue = co["c", "Pr(>|t|)"],
    beta_std_cue = fixef(mz)["c"],
    R2m_full = r2_marginal(m), R2m_no_source = r2_marginal(m0),
    dR2m_source = r2_marginal(m) - r2_marginal(m0),
    chisq_source = lrt$Chisq[2], p_lrt_source = lrt$`Pr(>Chisq)`[2]
  )
}

# ---------------------------------------------------------------------------
# 1. Head-referenced test (primary)
#    Direction cosines in the head frame: uy is sin(elevation), ux the
#    lateral (interaural) axis, uz the front/back axis.
# ---------------------------------------------------------------------------

frames <- c("onset", "tw", "final")
axes <- tribble(
  ~axis,       ~suffix,
  "vertical (uy = sin elevation)",   "uy",
  "lateral (ux = sin lateral angle)", "ux",
  "front/back (uz)",                  "uz"
)

head_frame_res <- expand_grid(frame = frames, axes) %>%
  pmap_dfr(function(frame, axis, suffix) {
    cue_weight_model(dat,
                     sprintf("resp_%s_%s", suffix, frame),
                     sprintf("src_%s_%s", suffix, frame),
                     sprintf("cue_%s_%s", suffix, frame),
                     sprintf("head-frame %s [%s]", axis, frame)) %>%
      mutate(frame = frame, axis = axis, .before = 1)
  })

say("\n--- 1. Head-referenced cue weights (direction cosines) ---")
say("b_source = auditory weight, b_cue = visual weight; 1 = full reliance.")
say_df(head_frame_res %>%
         select(frame, axis, b_source, ci_lo_source, ci_hi_source, p_source,
                b_cue, ci_lo_cue, ci_hi_cue, dR2m_source))
write_csv(head_frame_res, file.path(out_dir, "hrtf_headframe_cue_weights.csv"))

# Same test on elevation in degrees, which is the directly interpretable form.
elev_deg_res <- map_dfr(frames, function(f) {
  cue_weight_model(dat, sprintf("resp_elev_%s", f), sprintf("src_elev_%s", f),
                   sprintf("cue_elev_%s", f),
                   sprintf("head-frame elevation in degrees [%s]", f)) %>%
    mutate(frame = f, .before = 1)
})
say("\n--- 1b. Head-referenced elevation in degrees ---")
say_df(elev_deg_res %>% select(frame, b_source, ci_lo_source, ci_hi_source,
                               p_source, b_cue, dR2m_source))
write_csv(elev_deg_res, file.path(out_dir, "hrtf_headframe_elevation_deg.csv"))

# ---------------------------------------------------------------------------
# 2. World-frame cross-check
# ---------------------------------------------------------------------------

world_res <- pmap_dfr(
  list(c("response_y", "response_x", "response_z"),
       c("sound_y", "sound_x", "sound_z"),
       c("flash_y", "flash_x", "flash_z"),
       c("world vertical y", "world lateral x", "world depth z")),
  function(r, s, c, lab) cue_weight_model(dat, r, s, c, lab))

say("\n--- 2. World-frame cue weights (metres) ---")
say_df(world_res %>% select(model, b_source, ci_lo_source, ci_hi_source, p_source,
                            b_cue, ci_lo_cue, ci_hi_cue, dR2m_source, R2m_full))
write_csv(world_res, file.path(out_dir, "hrtf_worldframe_cue_weights.csv"))

# Random-slope version of the vertical world-frame model, and per-participant
# slopes, to check that a group-level null is not masking a subset of listeners.
m_slope <- lmer(response_y ~ sound_y + flash_y + (1 + sound_y | participantId),
                data = dat, REML = TRUE,
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5)))
say("\nRandom-slope model for world vertical (response_y ~ sound_y + flash_y):")
say(paste(capture.output(print(summary(m_slope)$coefficients, digits = 4)), collapse = "\n"))
say(sprintf("SD of participant-specific sound_y slopes: %.4f",
            attr(VarCorr(m_slope)$participantId, "stddev")["sound_y"]))

part_slopes <- dat %>%
  group_by(participantId) %>%
  group_modify(~ {
    fit <- lm(response_y ~ sound_y + flash_y, data = .x)
    tibble(b_sound_y = coef(fit)["sound_y"], se = summary(fit)$coefficients["sound_y", 2],
           b_flash_y = coef(fit)["flash_y"], n = nrow(.x))
  }) %>% ungroup()
say(sprintf("Per-participant OLS slopes on sound_y: median %.3f, IQR [%.3f, %.3f], %d of %d positive",
            median(part_slopes$b_sound_y), quantile(part_slopes$b_sound_y, .25),
            quantile(part_slopes$b_sound_y, .75),
            sum(part_slopes$b_sound_y > 0), nrow(part_slopes)))
write_csv(part_slopes, file.path(out_dir, "hrtf_participant_vertical_slopes.csv"))

# ---------------------------------------------------------------------------
# 3. Trials where the cue was displaced opposite to the source in elevation
#    If listeners had auditory elevation information, the response should be
#    pulled back towards the source on these trials.
# ---------------------------------------------------------------------------

opp <- dat %>%
  mutate(src_rel_y = sound_y - head_y_tw,
         cue_rel_y = flash_y - head_y_tw,
         cue_offset_y = soundToFlash_y,
         opposite_elev = sign(src_rel_y) != sign(cue_rel_y),
         large_offset = abs(cue_offset_y) > median(abs(cue_offset_y)))

say("\n--- 3a. Subsets by cue displacement in elevation ---")
subset_res <- bind_rows(
  cue_weight_model(opp %>% filter(opposite_elev), "response_y", "sound_y", "flash_y",
                   "cue on the opposite side of the head from the source (vertical)"),
  cue_weight_model(opp %>% filter(!opposite_elev), "response_y", "sound_y", "flash_y",
                   "cue on the same side (vertical)"),
  cue_weight_model(opp %>% filter(large_offset), "response_y", "sound_y", "flash_y",
                   "large vertical cue displacement"),
  cue_weight_model(opp %>% filter(!large_offset), "response_y", "sound_y", "flash_y",
                   "small vertical cue displacement")
)
say_df(subset_res %>% select(model, n, b_source, ci_lo_source, ci_hi_source, p_source, b_cue))
write_csv(subset_res, file.path(out_dir, "hrtf_vertical_subsets.csv"))

# Does the response deviate from the cue in the direction of the source?
# residual pull = (response - cue) projected on the (source - cue) direction.
pull <- dat %>%
  transmute(participantId,
            d_y = soundToFlash_y,                 # cue - source, vertical
            r_y = response_y - flash_y,           # response - cue, vertical
            pull_y = -r_y * sign(d_y),            # positive: response moved towards the source
            d_x = soundToFlash_x, r_x = response_x - flash_x, pull_x = -r_x * sign(d_x),
            d_z = soundToFlash_z, r_z = response_z - flash_z, pull_z = -r_z * sign(d_z)) %>%
  filter(d_y != 0)

pull_res <- map_dfr(c("x", "y", "z"), function(ax) {
  v <- paste0("pull_", ax)
  d <- pull %>% select(participantId, val = all_of(v)) %>% drop_na()
  m <- lmer(val ~ 1 + (1 | participantId), data = d, REML = TRUE)
  ci <- suppressMessages(confint(m, parm = "(Intercept)", method = "Wald"))
  co <- summary(m)$coefficients
  tibble(axis = ax, n = nrow(d),
         mean_pull_towards_source_m = co[1, "Estimate"],
         ci_lo = ci[1, 1], ci_hi = ci[1, 2],
         t = co[1, "t value"], p = co[1, "Pr(>|t|)"])
})
say("\n--- 3b. Displacement of the response from the cue towards the true source (m) ---")
say("Positive = the response moved from the cue in the direction of the source.")
say_df(pull_res)
write_csv(pull_res, file.path(out_dir, "hrtf_response_pull_towards_source.csv"))

# ---------------------------------------------------------------------------
# 4. Front/back
# ---------------------------------------------------------------------------

fb <- dat %>%
  mutate(src_front_onset  = src_uz_onset  > 0,
         cue_front_onset  = cue_uz_onset  > 0,
         resp_front_onset = resp_uz_onset > 0,
         src_front_final  = src_uz_final  > 0,
         cue_front_final  = cue_uz_final  > 0,
         resp_front_final = resp_uz_final > 0,
         reversal_onset   = src_front_onset != resp_front_onset,
         reversal_final   = src_front_final != resp_front_final,
         cue_opposite_onset = src_front_onset != cue_front_onset,
         cue_opposite_final = src_front_final != cue_front_final)

fb_res <- bind_rows(
  tibble(frame = "onset", subset = "all trials", n = nrow(fb),
         reversal_rate = mean(fb$reversal_onset)),
  tibble(frame = "onset", subset = "cue on the opposite front/back side",
         n = sum(fb$cue_opposite_onset),
         reversal_rate = mean(fb$reversal_onset[fb$cue_opposite_onset])),
  tibble(frame = "final", subset = "all trials", n = nrow(fb),
         reversal_rate = mean(fb$reversal_final)),
  tibble(frame = "final", subset = "cue on the opposite front/back side",
         n = sum(fb$cue_opposite_final),
         reversal_rate = mean(fb$reversal_final[fb$cue_opposite_final]))
) %>%
  mutate(ci_lo = map2_dbl(round(reversal_rate * n), n,
                          ~ binom.test(.x, .y)$conf.int[1]),
         ci_hi = map2_dbl(round(reversal_rate * n), n,
                          ~ binom.test(.x, .y)$conf.int[2]))

say("\n--- 4a. Front/back reversal rates (source vs response, head frame) ---")
say_df(fb_res)
write_csv(fb_res, file.path(out_dir, "hrtf_frontback_reversals.csv"))

# On trials where the cue sat on the opposite front/back side, did the response
# follow the source or the cue?
fb_conflict <- fb %>% filter(cue_opposite_onset) %>%
  summarise(n = n(),
            followed_source = mean(resp_front_onset == src_front_onset),
            followed_cue    = mean(resp_front_onset == cue_front_onset))
say("\n--- 4b. Front/back conflict trials (onset frame) ---")
say_df(fb_conflict)

# Time-resolved version: proportion of the trial for which the source lay
# behind the head, which bounds how often a mirror ambiguity could bite.
say(sprintf("\nFraction of trial time the source was in front of the head: mean %.3f (SD %.3f)",
            mean(dat$frac_time_src_front), sd(dat$frac_time_src_front)))
say(sprintf("Trials with the source behind the head for >50%% of the trial: %d of %d",
            sum(dat$frac_time_src_front < 0.5), nrow(dat)))

# ---------------------------------------------------------------------------
# 5. Alternative route to elevation without HRTF
#    (a) distance plus lateral panning constrains the magnitude but not the
#        sign of the elevation offset;
#    (b) head roll/pitch rotates the panning plane, so a listener who tilts
#        the head can in principle read elevation off the lateral panner.
#    If elevation performance were carried by (b) it should scale with head
#    tilt; a genuine HRTF cue should not need it.
# ---------------------------------------------------------------------------

mod_dat <- dat %>%
  mutate(roll_c  = as.vector(scale(roll_sd_deg)),
         pitch_c = as.vector(scale(pitch_sd_deg)),
         sound_y_c = sound_y - mean(sound_y),
         flash_y_c = flash_y - mean(flash_y))

m_roll <- lmer(response_y ~ sound_y_c * roll_c + flash_y_c + (1 | participantId),
               data = mod_dat, REML = TRUE)
m_pitch <- lmer(response_y ~ sound_y_c * pitch_c + flash_y_c + (1 | participantId),
                data = mod_dat, REML = TRUE)
say("\n--- 5. Head-tilt moderation of the auditory vertical weight ---")
say("Roll model:")
say(paste(capture.output(print(summary(m_roll)$coefficients, digits = 4)), collapse = "\n"))
say("Pitch model:")
say(paste(capture.output(print(summary(m_pitch)$coefficients, digits = 4)), collapse = "\n"))

mod_res <- bind_rows(
  tibble(term = "sound_y_c:roll_c",
         estimate = summary(m_roll)$coefficients["sound_y_c:roll_c", "Estimate"],
         se = summary(m_roll)$coefficients["sound_y_c:roll_c", "Std. Error"],
         p = summary(m_roll)$coefficients["sound_y_c:roll_c", "Pr(>|t|)"]),
  tibble(term = "sound_y_c:pitch_c",
         estimate = summary(m_pitch)$coefficients["sound_y_c:pitch_c", "Estimate"],
         se = summary(m_pitch)$coefficients["sound_y_c:pitch_c", "Std. Error"],
         p = summary(m_pitch)$coefficients["sound_y_c:pitch_c", "Pr(>|t|)"])
) %>% mutate(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)
write_csv(mod_res, file.path(out_dir, "hrtf_headtilt_moderation.csv"))
say_df(mod_res)

# Magnitude-only account: does |response elevation offset| track |source
# elevation offset| even when the signed relation is absent?
mag <- dat %>%
  transmute(participantId,
            src_abs = abs(src_elev_tw), cue_abs = abs(cue_elev_tw),
            resp_abs = abs(resp_elev_tw))
m_mag <- lmer(resp_abs ~ src_abs + cue_abs + (1 | participantId), data = mag, REML = TRUE)
say("\nMagnitude-only elevation model (|elevation| in degrees):")
say(paste(capture.output(print(summary(m_mag)$coefficients, digits = 4)), collapse = "\n"))

# ---------------------------------------------------------------------------
# 6. Precision benchmark: residual scatter of the response about the source and
#    about the cue, per axis. Under panning without HRTF the vertical residual
#    should be dominated by the cue.
# ---------------------------------------------------------------------------

resid_res <- tibble(
  axis = c("x", "y", "z"),
  rmse_to_source_m = c(sqrt(mean((dat$response_x - dat$sound_x)^2)),
                       sqrt(mean((dat$response_y - dat$sound_y)^2)),
                       sqrt(mean((dat$response_z - dat$sound_z)^2))),
  rmse_to_cue_m = c(sqrt(mean((dat$response_x - dat$flash_x)^2)),
                    sqrt(mean((dat$response_y - dat$flash_y)^2)),
                    sqrt(mean((dat$response_z - dat$flash_z)^2)))
)
say("\n--- 6. Response scatter about the source and about the cue ---")
say_df(resid_res)
write_csv(resid_res, file.path(out_dir, "hrtf_axis_rmse.csv"))

say("\nOutputs written to ", out_dir)
close(log_con)
