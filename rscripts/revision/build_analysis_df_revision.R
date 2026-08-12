#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Rebuild the trial-level analysis dataset for the 6DoF audiovisual
# localization study, extended with head-trajectory-derived measures that the
# original analysis did not compute.
#
# Steps
#   1. Reproduce the original trial-level frame (JSON responses).
#   2. Per-trial head movement metrics (trajectory CSVs), as originally defined.
#   3. New per-sample measures aggregated per trial:
#        (a) experienced angular disparity source-head-cue,
#        (b) head-to-source distance,
#        (c) cue and source eccentricity relative to head forward,
#        (d) inter-trial gaps in the trajectory timebase.
#   4. Validation against the published descriptives.
#
# Coordinate conventions: Unity world frame, left-handed, y up, metres.
# Quaternions are stored (x, y, z, w); forward is +Z.
#
# Run with R 4.4 (arm64). Deterministic; no random numbers are used.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(readr)
})

set.seed(20260805)  # no stochastic step here, fixed for reproducibility

project_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
data_dir    <- file.path(project_dir, "data")
out_dir     <- file.path(project_dir, "results_revision")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

log_path <- file.path(out_dir, "build_analysis_df_revision_log.txt")
log_con  <- file(log_path, open = "wt")
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_con)
}
say_df <- function(x, digits = 4) {
  out <- capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE))
  cat(out, sep = "\n")
  cat(out, sep = "\n", file = log_con)
  cat("\n"); cat("\n", file = log_con)
}

say("=== Build revision analysis dataset ===")
say("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
say("R version: ", R.version.string)

# ---------------------------------------------------------------------------
# File pairing: participant IDs come from sorted JSON filename order; the
# trajectory CSV exported for the same session carries a slightly different
# filename timestamp, so pairing is by sorted order and validated below.
# ---------------------------------------------------------------------------

json_files <- sort(list.files(data_dir, pattern = "\\.json$", full.names = TRUE))
csv_files  <- sort(list.files(data_dir, pattern = "\\.csv$",  full.names = TRUE))
stopifnot(length(json_files) == length(csv_files))

pids <- sprintf("P%03d", seq_along(json_files))
say("\nSessions found: ", length(json_files),
    " (participants ", pids[1], "..", pids[length(pids)], ")")

file_stamp <- function(paths) sub("^[^_]+_", "", tools::file_path_sans_ext(basename(paths)))
pair_stamps <- tibble(
  participantId = pids,
  json_file  = basename(json_files),
  csv_file   = basename(csv_files),
  json_stamp = as.POSIXct(file_stamp(json_files), format = "%Y%m%d_%H%M%S", tz = "UTC"),
  csv_stamp  = as.POSIXct(file_stamp(csv_files),  format = "%Y%m%d_%H%M%S", tz = "UTC")
) %>%
  mutate(stamp_diff_s = as.numeric(difftime(json_stamp, csv_stamp, units = "secs")))

# --- read JSON responses ---------------------------------------------------
read_experiment <- function(path, pid) {
  entries <- fromJSON(path, flatten = TRUE)$dataEntries
  df <- as.data.frame(entries)
  df$participantId <- pid
  df
}
all_experiments_df <- map2_dfr(json_files, pids, read_experiment)
say("Trials loaded from JSON: ", nrow(all_experiments_df))

# --- read trajectory CSVs --------------------------------------------------
read_trajectory <- function(path, pid) {
  df <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  df$participantId <- pid
  df
}
all_trajectories_df <- map2_dfr(csv_files, pids, read_trajectory)
say("Trajectory samples loaded from CSV: ", nrow(all_trajectories_df))

# ---------------------------------------------------------------------------
# Trajectory runs. A "run" is a maximal block of consecutive samples carrying
# the same level in the session timebase. Normally there is exactly one run per
# level, but a few sessions revisit a level, so runs are the unit that defines
# the sampling timebase (within-run intervals, inter-trial gaps).
# ---------------------------------------------------------------------------

traj_time <- all_trajectories_df %>%
  arrange(participantId, timestamp) %>%
  group_by(participantId) %>%
  mutate(run_id = cumsum(level != lag(level, default = first(level))) + 1L) %>%
  ungroup()

runs <- traj_time %>%
  group_by(participantId, level, run_id) %>%
  summarise(run_start = min(timestamp), run_end = max(timestamp),
            run_samples = n(), .groups = "drop") %>%
  mutate(run_duration = run_end - run_start) %>%
  arrange(participantId, run_start)

# Which run of a level corresponds to the response recorded in the JSON: the
# one whose duration is closest to levelTime.
runs_matched <- runs %>%
  left_join(all_experiments_df %>% select(participantId, level, levelTime),
            by = c("participantId", "level")) %>%
  group_by(participantId, level) %>%
  mutate(n_runs = n(),
         run_dur_diff = run_duration - levelTime,
         is_matched_run = abs(run_dur_diff) == min(abs(run_dur_diff))) %>%
  ungroup()

revisited <- runs_matched %>%
  group_by(participantId, level) %>%
  filter(n_runs > 1) %>%
  ungroup()

# ---------------------------------------------------------------------------
# Pairing validation: identical level sets, and trajectory-derived trial
# duration consistent with the JSON levelTime.
# ---------------------------------------------------------------------------

traj_levels <- all_trajectories_df %>%
  group_by(participantId) %>%
  summarise(csv_levels = n_distinct(level),
            csv_min = min(level), csv_max = max(level),
            csv_samples = n(), .groups = "drop")

json_levels <- all_experiments_df %>%
  group_by(participantId) %>%
  summarise(json_levels = n_distinct(level),
            json_min = min(level), json_max = max(level), .groups = "drop")

traj_dur <- all_trajectories_df %>%
  group_by(participantId, level) %>%
  summarise(traj_dur = max(timestamp) - min(timestamp), .groups = "drop")

dur_check <- all_experiments_df %>%
  select(participantId, level, levelTime) %>%
  left_join(traj_dur, by = c("participantId", "level")) %>%
  left_join(runs_matched %>% filter(is_matched_run) %>%
              select(participantId, level, matched_run_duration = run_duration,
                     n_runs),
            by = c("participantId", "level")) %>%
  mutate(dur_diff = traj_dur - levelTime,
         matched_dur_diff = matched_run_duration - levelTime)

pair_report <- pair_stamps %>%
  left_join(json_levels, by = "participantId") %>%
  left_join(traj_levels, by = "participantId") %>%
  left_join(
    dur_check %>%
      group_by(participantId) %>%
      summarise(n_revisited_trials = sum(n_runs > 1),
                median_dur_diff_s = median(dur_diff),
                max_abs_dur_diff_s = max(abs(dur_diff)),
                dur_cor = cor(traj_dur, levelTime),
                max_abs_matched_dur_diff_s = max(abs(matched_dur_diff)),
                matched_dur_cor = cor(matched_run_duration, levelTime),
                .groups = "drop"),
    by = "participantId"
  ) %>%
  mutate(levels_match = json_levels == 24 & csv_levels == 24 &
           json_min == 1 & json_max == 24 & csv_min == 1 & csv_max == 24,
         duration_ok = max_abs_dur_diff_s < 2 & dur_cor > 0.99,
         matched_duration_ok = max_abs_matched_dur_diff_s < 2 & matched_dur_cor > 0.99)

say("\n--- Session pairing validation ---")
say_df(pair_report %>%
         select(participantId, json_file, csv_file, stamp_diff_s,
                json_levels, csv_levels, csv_samples, n_revisited_trials,
                median_dur_diff_s, max_abs_dur_diff_s, dur_cor,
                max_abs_matched_dur_diff_s, levels_match, duration_ok,
                matched_duration_ok))
say("Pairs failing the level-set check:                         ", sum(!pair_report$levels_match))
say("Pairs failing the duration check (all samples per level):  ", sum(!pair_report$duration_ok))
say("Pairs failing the duration check (matched run only):       ", sum(!pair_report$matched_duration_ok))
say(sprintf("Filename timestamp offset JSON-CSV (s): median %.1f, range %.1f to %.1f",
            median(pair_report$stamp_diff_s), min(pair_report$stamp_diff_s),
            max(pair_report$stamp_diff_s)))
say(sprintf("Trajectory duration minus levelTime (s), all %d trials: median %.3f, IQR %.3f to %.3f, range %.3f to %.3f",
            nrow(dur_check), median(dur_check$dur_diff),
            quantile(dur_check$dur_diff, 0.25), quantile(dur_check$dur_diff, 0.75),
            min(dur_check$dur_diff), max(dur_check$dur_diff)))
say(sprintf("Same, using only the levelTime-matched run: median %.3f, range %.3f to %.3f",
            median(dur_check$matched_dur_diff), min(dur_check$matched_dur_diff),
            max(dur_check$matched_dur_diff)))

say("\n--- Trials whose level is revisited later in the session ---")
say("Trials affected: ", n_distinct(paste(revisited$participantId, revisited$level)),
    " in ", n_distinct(revisited$participantId), " sessions; extra runs: ",
    nrow(revisited) - n_distinct(paste(revisited$participantId, revisited$level)))
say_df(revisited %>%
         select(participantId, level, run_id, run_start, run_end, run_duration,
                run_samples, levelTime, run_dur_diff, is_matched_run))

# ---------------------------------------------------------------------------
# STEP 1: trial-level frame, reproducing the original derived variables
# ---------------------------------------------------------------------------

analysis_df <- all_experiments_df %>%
  mutate(
    participantError_m  = distanceError / 100,
    stimulusDisparity_m = flashDistance,
    # design ranges, classified in cm: Short 15-30, Medium >30-40, Long >40-70
    disparityRange = cut(flashDistance * 100,
                         breaks = c(0, 30, 40, 70),
                         labels = c("Short (15-30cm)", "Medium (30-40cm)", "Long (40-70cm)"),
                         include.lowest = TRUE),
    soundType = factor(sampleId, levels = 1:4,
                       labels = c("Drum", "Flute", "Speech", "Pink Noise")),
    # azimuth wraparound: values above 315 deg map to negative angles
    azimuth_adj = ifelse(flashAzimuth > 315, flashAzimuth - 360, flashAzimuth),
    azimuthSector = cut(azimuth_adj, breaks = c(-45, 45, 135, 225, 315),
                        labels = c("Front", "Right", "Back", "Left"),
                        include.lowest = TRUE),
    elevationCategory = cut(flashElevation, breaks = c(-90, -22.5, 22.5, 90),
                            labels = c("Below", "Level", "Above"),
                            include.lowest = TRUE)
  ) %>%
  transmute(
    participantId,
    trialConditionId = trialId,
    trialSequenceNum = level,
    soundType, disparityRange, azimuthSector, elevationCategory,
    stimulusDisparity_m,
    flashAzimuth_deg   = flashAzimuth,
    flashElevation_deg = flashElevation,
    responseTime_s     = levelTime,
    participantError_m,
    response_x = spherePosition.x, response_y = spherePosition.y, response_z = spherePosition.z,
    sound_x = referencePosition.x, sound_y = referencePosition.y, sound_z = referencePosition.z,
    flash_x = flashSpherePosition.x, flash_y = flashSpherePosition.y, flash_z = flashSpherePosition.z
  )

analysis_df <- analysis_df %>%
  mutate(
    # displacement axis: from the true source towards the visual cue
    soundToFlash_x = flash_x - sound_x,
    soundToFlash_y = flash_y - sound_y,
    soundToFlash_z = flash_z - sound_z,
    soundToFlash_length = sqrt(soundToFlash_x^2 + soundToFlash_y^2 + soundToFlash_z^2),
    soundToFlash_unit_x = soundToFlash_x / soundToFlash_length,
    soundToFlash_unit_y = soundToFlash_y / soundToFlash_length,
    soundToFlash_unit_z = soundToFlash_z / soundToFlash_length,
    # response displacement relative to the true source
    soundToResponse_x = response_x - sound_x,
    soundToResponse_y = response_y - sound_y,
    soundToResponse_z = response_z - sound_z,
    # signed error: positive = shifted towards the visual cue
    signedError_m = soundToResponse_x * soundToFlash_unit_x +
                    soundToResponse_y * soundToFlash_unit_y +
                    soundToResponse_z * soundToFlash_unit_z,
    # orthogonal error: norm of the component perpendicular to the axis.
    # Computed from the residual vector rather than sqrt(err^2 - signed^2),
    # which is undefined whenever rounding pushes the radicand below zero.
    orth_x = soundToResponse_x - signedError_m * soundToFlash_unit_x,
    orth_y = soundToResponse_y - signedError_m * soundToFlash_unit_y,
    orth_z = soundToResponse_z - signedError_m * soundToFlash_unit_z,
    orthogonalError_m = sqrt(orth_x^2 + orth_y^2 + orth_z^2),
    ventriloquistBias = signedError_m / stimulusDisparity_m,
    participantError_cm = participantError_m * 100
  ) %>%
  select(-orth_x, -orth_y, -orth_z)

n_neg_radicand <- sum(analysis_df$participantError_m^2 - analysis_df$signedError_m^2 < 0)
say("\n--- Orthogonal error ---")
say("Trials where sqrt(participantError^2 - signedError^2) would be NaN: ", n_neg_radicand)
say(sprintf("Max |vector orthogonal error - sqrt-form| over the %d valid trials: %.3e m",
            sum(analysis_df$participantError_m^2 - analysis_df$signedError_m^2 >= 0),
            max(abs(analysis_df$orthogonalError_m -
                      suppressWarnings(sqrt(analysis_df$participantError_m^2 -
                                              analysis_df$signedError_m^2))), na.rm = TRUE)))
say(sprintf("Max |soundToFlash_length - stimulusDisparity_m|: %.3e m",
            max(abs(analysis_df$soundToFlash_length - analysis_df$stimulusDisparity_m))))

# ---------------------------------------------------------------------------
# STEP 2: per-trial head movement metrics (original definitions)
# ---------------------------------------------------------------------------

trajectory_clean <- all_trajectories_df %>%
  select(participantId, level, timestamp, px, py, pz, rx, ry, rz, rw, ex, ey, ez) %>%
  rename(trialSequenceNum = level) %>%
  arrange(participantId, trialSequenceNum, timestamp)

trajectory_summary <- trajectory_clean %>%
  group_by(participantId, trialSequenceNum) %>%
  mutate(
    dx = px - lag(px, default = first(px)),
    dy = py - lag(py, default = first(py)),
    dz = pz - lag(pz, default = first(pz)),
    segment_distance = sqrt(dx^2 + dy^2 + dz^2),
    d_yaw   = abs(ey - lag(ey, default = first(ey))),
    d_pitch = abs(ex - lag(ex, default = first(ex))),
    d_roll  = abs(ez - lag(ez, default = first(ez))),
    d_yaw   = pmin(d_yaw, 360 - d_yaw)   # Euler yaw wraparound
  ) %>%
  summarise(
    total_path_length     = sum(segment_distance, na.rm = TRUE),
    horizontal_movement   = sum(sqrt(dx^2 + dz^2), na.rm = TRUE),
    vertical_movement     = sum(abs(dy), na.rm = TRUE),
    total_yaw_rotation    = sum(d_yaw, na.rm = TRUE),
    total_pitch_rotation  = sum(d_pitch, na.rm = TRUE),
    x_range = max(px) - min(px),
    y_range = max(py) - min(py),
    z_range = max(pz) - min(pz),
    exploration_volume = x_range * y_range * z_range,
    trajectory_duration = max(timestamp) - min(timestamp),
    n_samples = n(),
    movement_rate = total_path_length / pmax(trajectory_duration, 0.1),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------
# STEP 3: new per-sample geometry
# ---------------------------------------------------------------------------

deg <- function(rad) rad * 180 / pi

# Angle between two 3-vectors, numerically stable via atan2.
angle_between <- function(ax, ay, az, bx, by, bz) {
  cx <- ay * bz - az * by
  cy <- az * bx - ax * bz
  cz <- ax * by - ay * bx
  deg(atan2(sqrt(cx^2 + cy^2 + cz^2), ax * bx + ay * by + az * bz))
}

# Head forward (+Z) rotated by the stored quaternion (x, y, z, w): third column
# of the rotation matrix. Unity's left-handed convention is already encoded in
# the stored components, so the standard form applies; this is validated below
# against the Euler angles recorded in the same rows.
quat_forward <- function(x, y, z, w) {
  list(fx = 2 * (x * z + w * y),
       fy = 2 * (y * z - w * x),
       fz = 1 - 2 * (x^2 + y^2))
}

# Independent derivation of forward from Unity Euler angles (deg). Roll is
# about the forward axis and therefore does not affect it.
euler_forward <- function(ex, ey, ez) {
  pitch <- ex * pi / 180
  yaw   <- ey * pi / 180
  list(fx = sin(yaw) * cos(pitch),
       fy = -sin(pitch),
       fz = cos(yaw) * cos(pitch))
}

trial_static <- analysis_df %>%
  select(participantId, trialSequenceNum,
         sound_x, sound_y, sound_z, flash_x, flash_y, flash_z)

samples <- traj_time %>%
  select(participantId, trialSequenceNum = level, run_id, timestamp,
         px, py, pz, rx, ry, rz, rw, ex, ey, ez) %>%
  left_join(trial_static, by = c("participantId", "trialSequenceNum")) %>%
  # sampling intervals are defined within a run, so that a level revisited
  # later in the session does not create a spurious multi-second interval
  group_by(participantId, run_id) %>%
  mutate(dt_next = lead(timestamp) - timestamp) %>%
  group_by(participantId, trialSequenceNum) %>%
  mutate(
    dt_med = median(dt_next, na.rm = TRUE),
    # time weight: interval to the next sample; the final sample of each run is
    # given the trial's median interval
    w_time = ifelse(is.na(dt_next), ifelse(is.na(dt_med), 0.111, dt_med), dt_next),
    sample_idx = rank(timestamp, ties.method = "first")
  ) %>%
  ungroup() %>%
  mutate(
    # head -> true source and head -> visual cue
    hs_x = sound_x - px, hs_y = sound_y - py, hs_z = sound_z - pz,
    hc_x = flash_x - px, hc_y = flash_y - py, hc_z = flash_z - pz,
    dist_source_m = sqrt(hs_x^2 + hs_y^2 + hs_z^2),
    dist_cue_m    = sqrt(hc_x^2 + hc_y^2 + hc_z^2),
    angDisp_deg   = angle_between(hs_x, hs_y, hs_z, hc_x, hc_y, hc_z)
  )

fwd_q <- quat_forward(samples$rx, samples$ry, samples$rz, samples$rw)
fwd_e <- euler_forward(samples$ex, samples$ey, samples$ez)
samples <- samples %>%
  mutate(
    fwd_x = fwd_q$fx, fwd_y = fwd_q$fy, fwd_z = fwd_q$fz,
    fwdE_x = fwd_e$fx, fwdE_y = fwd_e$fy, fwdE_z = fwd_e$fz,
    fwd_check_deg = angle_between(fwd_x, fwd_y, fwd_z, fwdE_x, fwdE_y, fwdE_z),
    cueEcc_deg = angle_between(fwd_x, fwd_y, fwd_z, hc_x, hc_y, hc_z),
    srcEcc_deg = angle_between(fwd_x, fwd_y, fwd_z, hs_x, hs_y, hs_z)
  )

say("\n--- Head forward vector validation (quaternion vs Euler) ---")
say(sprintf("Angular difference (deg) over %d samples: median %.4f, mean %.4f, p99 %.4f, max %.4f",
            nrow(samples), median(samples$fwd_check_deg), mean(samples$fwd_check_deg),
            quantile(samples$fwd_check_deg, 0.99), max(samples$fwd_check_deg)))
say(sprintf("Quaternion norm: min %.6f, max %.6f",
            min(sqrt(samples$rx^2 + samples$ry^2 + samples$rz^2 + samples$rw^2)),
            max(sqrt(samples$rx^2 + samples$ry^2 + samples$rz^2 + samples$rw^2))))
if (median(samples$fwd_check_deg) > 3) {
  stop("Quaternion and Euler forward vectors disagree; check handedness/component order.")
}

wmean <- function(x, w) sum(x * w) / sum(w)

summarise_trials <- function(df) {
  df %>%
  group_by(participantId, trialSequenceNum) %>%
  summarise(
    # (a) experienced angular disparity, source-head-cue
    angDisp_onset_deg  = angDisp_deg[which.min(sample_idx)],
    angDisp_mean_deg   = wmean(angDisp_deg, w_time),
    angDisp_median_deg = median(angDisp_deg),
    angDisp_min_deg    = min(angDisp_deg),
    angDisp_max_deg    = max(angDisp_deg),
    # (b) head-to-source distance
    dist_mean_m     = wmean(dist_source_m, w_time),
    dist_min_m      = min(dist_source_m),
    dist_at_onset_m = dist_source_m[which.min(sample_idx)],
    frac_samples_under_0p5m   = mean(dist_source_m < 0.5),
    frac_samples_under_0p574m = mean(dist_source_m < 0.574),
    # (c) eccentricity of cue and source relative to head forward
    cueEcc_mean_deg  = wmean(cueEcc_deg, w_time),
    cueEcc_onset_deg = cueEcc_deg[which.min(sample_idx)],
    cueEcc_median_deg = median(cueEcc_deg),
    cueEcc_min_deg   = min(cueEcc_deg),
    cueEcc_max_deg   = max(cueEcc_deg),
    frac_time_cue_within_55deg = wmean(as.numeric(cueEcc_deg < 55), w_time),
    srcEcc_mean_deg  = wmean(srcEcc_deg, w_time),
    srcEcc_onset_deg = srcEcc_deg[which.min(sample_idx)],
    srcEcc_median_deg = median(srcEcc_deg),
    srcEcc_min_deg   = min(srcEcc_deg),
    srcEcc_max_deg   = max(srcEcc_deg),
    frac_time_src_within_55deg = wmean(as.numeric(srcEcc_deg < 55), w_time),
    n_time_blocks = n_distinct(run_id),
    .groups = "drop"
  ) %>%
  mutate(trial_has_revisit = n_time_blocks > 1)
}

sample_summary <- summarise_trials(samples)

# Sensitivity check for the revisited trials: the same summaries computed on
# the levelTime-matched run alone.
matched_runs <- runs_matched %>% filter(is_matched_run) %>%
  select(participantId, trialSequenceNum = level, run_id)
sample_summary_matched <- summarise_trials(
  samples %>% semi_join(matched_runs, by = c("participantId", "trialSequenceNum", "run_id"))
)

revisit_sensitivity <- sample_summary %>%
  filter(trial_has_revisit) %>%
  select(participantId, trialSequenceNum, angDisp_mean_deg, dist_mean_m,
         cueEcc_mean_deg, frac_time_cue_within_55deg) %>%
  left_join(sample_summary_matched %>%
              select(participantId, trialSequenceNum,
                     angDisp_mean_deg_matched = angDisp_mean_deg,
                     dist_mean_m_matched = dist_mean_m,
                     cueEcc_mean_deg_matched = cueEcc_mean_deg,
                     frac_time_cue_within_55deg_matched = frac_time_cue_within_55deg),
            by = c("participantId", "trialSequenceNum"))

say("\n--- Revisited trials: pooled vs levelTime-matched-run summaries ---")
say_df(revisit_sensitivity)

# --- (d) inter-trial gaps in the trajectory timebase ------------------------
# Primary definition uses runs in session time order, so revisits do not
# produce negative gaps; the level-ordered definition is reported alongside.
inter_trial_gaps <- runs %>%
  group_by(participantId) %>%
  mutate(gap_to_next_s = lead(run_start) - run_end,
         gap_from_prev_s = run_start - lag(run_end)) %>%
  ungroup()

gaps <- inter_trial_gaps$gap_to_next_s[!is.na(inter_trial_gaps$gap_to_next_s)]
say("\n--- Inter-trial gaps (last sample of trial n to first sample of trial n+1) ---")
say(sprintf("Run-based: n = %d gaps across %d sessions", length(gaps),
            n_distinct(inter_trial_gaps$participantId)))
say(sprintf("median %.4f s, IQR %.4f to %.4f s (IQR width %.4f), range %.4f to %.4f s, mean %.4f s",
            median(gaps), quantile(gaps, 0.25), quantile(gaps, 0.75),
            IQR(gaps), min(gaps), max(gaps), mean(gaps)))

level_gaps <- trajectory_clean %>%
  group_by(participantId, trialSequenceNum) %>%
  summarise(t_first = min(timestamp), t_last = max(timestamp), .groups = "drop") %>%
  arrange(participantId, trialSequenceNum) %>%
  group_by(participantId) %>%
  mutate(gap_to_next_s = lead(t_first) - t_last) %>%
  ungroup() %>%
  filter(!is.na(gap_to_next_s))
say(sprintf("Level-ordered (all samples of a level pooled): n = %d, median %.4f s, range %.3f to %.3f s, negative gaps: %d",
            nrow(level_gaps), median(level_gaps$gap_to_next_s),
            min(level_gaps$gap_to_next_s), max(level_gaps$gap_to_next_s),
            sum(level_gaps$gap_to_next_s < 0)))

# --- sampling interval description -----------------------------------------
intervals <- samples$dt_next[!is.na(samples$dt_next)]
say(sprintf("Within-run sampling interval: median %.4f s, IQR %.4f to %.4f, p99 %.4f, max %.4f",
            median(intervals), quantile(intervals, 0.25), quantile(intervals, 0.75),
            quantile(intervals, 0.99), max(intervals)))

# ---------------------------------------------------------------------------
# Assemble the final frame
# ---------------------------------------------------------------------------

analysis_df <- analysis_df %>%
  left_join(trajectory_summary, by = c("participantId", "trialSequenceNum")) %>%
  left_join(sample_summary,     by = c("participantId", "trialSequenceNum")) %>%
  group_by(participantId) %>%
  mutate(participant_mean_rate = mean(movement_rate, na.rm = TRUE),
         trial_rate_deviation  = movement_rate - participant_mean_rate) %>%
  ungroup() %>%
  mutate(obs_id = row_number(),
         stimDisparity_c = stimulusDisparity_m - mean(stimulusDisparity_m, na.rm = TRUE),
         trialSequence_c = trialSequenceNum - mean(trialSequenceNum, na.rm = TRUE))

say("\n--- Merge check ---")
say("Rows: ", nrow(analysis_df),
    " | with trajectory metrics: ", sum(!is.na(analysis_df$total_path_length)),
    " | with angular metrics: ", sum(!is.na(analysis_df$angDisp_mean_deg)))

# ---------------------------------------------------------------------------
# STEP 4: validation against published descriptives
# ---------------------------------------------------------------------------

pp_path <- analysis_df %>%
  group_by(participantId) %>%
  summarise(mean_path = mean(total_path_length), .groups = "drop")

fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
near <- function(a, b, tol) abs(a - b) <= tol

validation <- tibble::tribble(
  ~check, ~published, ~recomputed, ~matches,
  "n trials", "744", as.character(nrow(analysis_df)), nrow(analysis_df) == 744L,
  "n participants", "31", as.character(n_distinct(analysis_df$participantId)),
    n_distinct(analysis_df$participantId) == 31L,
  "mean unsigned error (cm)", "26.2", fmt(mean(analysis_df$participantError_cm), 1),
    near(mean(analysis_df$participantError_cm), 26.2, 0.05),
  "SD unsigned error (cm)", "17.0", fmt(sd(analysis_df$participantError_cm), 1),
    near(sd(analysis_df$participantError_cm), 17.0, 0.05),
  "median unsigned error (cm)", "22.3", fmt(median(analysis_df$participantError_cm), 1),
    near(median(analysis_df$participantError_cm), 22.3, 0.05),
  "mean signed error (cm)", "13.7", fmt(mean(analysis_df$signedError_m) * 100, 1),
    near(mean(analysis_df$signedError_m) * 100, 13.7, 0.05),
  "SD signed error (cm)", "18.5", fmt(sd(analysis_df$signedError_m) * 100, 1),
    near(sd(analysis_df$signedError_m) * 100, 18.5, 0.05),
  "mean bias proportion (%)", "35.9", fmt(mean(analysis_df$ventriloquistBias) * 100, 1),
    near(mean(analysis_df$ventriloquistBias) * 100, 35.9, 0.05),
  "SD bias proportion (%)", "54.6", fmt(sd(analysis_df$ventriloquistBias) * 100, 1),
    near(sd(analysis_df$ventriloquistBias) * 100, 54.6, 0.05),
  "mean path length (m)", "8.36", fmt(mean(analysis_df$total_path_length), 2),
    near(mean(analysis_df$total_path_length), 8.36, 0.005),
  "SD path length (m)", "6.21", fmt(sd(analysis_df$total_path_length), 2),
    near(sd(analysis_df$total_path_length), 6.21, 0.005),
  "per-participant mean path, min (m)", "4.6", fmt(min(pp_path$mean_path), 1),
    near(min(pp_path$mean_path), 4.6, 0.05),
  "per-participant mean path, max (m)", "22.0", fmt(max(pp_path$mean_path), 1),
    near(max(pp_path$mean_path), 22.0, 0.05),
  "distinct disparity values", "56", as.character(n_distinct(analysis_df$stimulusDisparity_m)),
    n_distinct(analysis_df$stimulusDisparity_m) == 56L
)

say("\n--- Validation against published descriptives ---")
say_df(validation)

say("Crosstab soundType x disparityRange:")
crosstab <- as.data.frame.matrix(table(analysis_df$soundType, analysis_df$disparityRange))
crosstab <- cbind(soundType = rownames(crosstab), crosstab)
say_df(crosstab)
write_csv(crosstab, file.path(out_dir, "crosstab_soundtype_disparityrange.csv"))
say(sprintf("responseTime_s: min %.3f, max %.3f, median %.3f, mean %.3f",
            min(analysis_df$responseTime_s), max(analysis_df$responseTime_s),
            median(analysis_df$responseTime_s), mean(analysis_df$responseTime_s)))
say(sprintf("Distinct stimulusDisparity_m values: %d (range %.3f to %.3f m)",
            n_distinct(analysis_df$stimulusDisparity_m),
            min(analysis_df$stimulusDisparity_m), max(analysis_df$stimulusDisparity_m)))

# ---------------------------------------------------------------------------
# Summaries of the new measures
# ---------------------------------------------------------------------------

new_cols <- c("angDisp_onset_deg", "angDisp_mean_deg", "angDisp_median_deg",
              "angDisp_min_deg", "angDisp_max_deg",
              "dist_mean_m", "dist_min_m", "dist_at_onset_m",
              "frac_samples_under_0p5m", "frac_samples_under_0p574m",
              "cueEcc_mean_deg", "cueEcc_onset_deg", "cueEcc_median_deg",
              "cueEcc_min_deg", "cueEcc_max_deg", "frac_time_cue_within_55deg",
              "srcEcc_mean_deg", "srcEcc_onset_deg", "srcEcc_median_deg",
              "srcEcc_min_deg", "srcEcc_max_deg", "frac_time_src_within_55deg")

describe <- function(df, cols) {
  map_dfr(cols, function(v) {
    x <- df[[v]]
    tibble(measure = v, n = sum(!is.na(x)), mean = mean(x, na.rm = TRUE),
           sd = sd(x, na.rm = TRUE), min = min(x, na.rm = TRUE),
           q25 = quantile(x, 0.25, na.rm = TRUE), median = median(x, na.rm = TRUE),
           q75 = quantile(x, 0.75, na.rm = TRUE), max = max(x, na.rm = TRUE))
  })
}

new_summary <- describe(analysis_df, new_cols)
say("\n--- Summary statistics: new trial-level measures ---")
say_df(new_summary, digits = 4)

movement_cols <- c("total_path_length", "horizontal_movement", "vertical_movement",
                   "total_yaw_rotation", "total_pitch_rotation",
                   "x_range", "y_range", "z_range", "exploration_volume",
                   "trajectory_duration", "n_samples", "movement_rate",
                   "participant_mean_rate", "trial_rate_deviation")
movement_summary <- describe(analysis_df, movement_cols)
say("--- Summary statistics: movement measures (Step 2) ---")
say_df(movement_summary, digits = 4)

core_cols <- c("stimulusDisparity_m", "participantError_m", "participantError_cm",
               "signedError_m", "orthogonalError_m", "ventriloquistBias",
               "responseTime_s")
core_summary <- describe(analysis_df, core_cols)
say("--- Summary statistics: core trial measures (Step 1) ---")
say_df(core_summary, digits = 4)

# Relationship between the static design disparity and what was experienced.
say(sprintf("Correlation stimulusDisparity_m with angDisp_mean_deg: r = %.3f",
            cor(analysis_df$stimulusDisparity_m, analysis_df$angDisp_mean_deg)))
say(sprintf("Correlation dist_mean_m with angDisp_mean_deg: r = %.3f",
            cor(analysis_df$dist_mean_m, analysis_df$angDisp_mean_deg)))
say(sprintf("Trials with any sample closer than 0.574 m to the source: %d (%.1f%%)",
            sum(analysis_df$frac_samples_under_0p574m > 0),
            100 * mean(analysis_df$frac_samples_under_0p574m > 0)))
say(sprintf("Trials where the cue was outside 55 deg of head forward for the whole trial: %d",
            sum(analysis_df$frac_time_cue_within_55deg == 0)))

# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

saveRDS(analysis_df, file.path(out_dir, "analysis_df_revision.rds"))
write_csv(analysis_df, file.path(out_dir, "analysis_df_revision.csv"))

# per-sample intermediate for the angular measures (compressed RDS)
samples_out <- samples %>%
  select(participantId, trialSequenceNum, run_id, timestamp, sample_idx, w_time,
         px, py, pz, fwd_x, fwd_y, fwd_z, fwd_check_deg,
         dist_source_m, dist_cue_m, angDisp_deg, cueEcc_deg, srcEcc_deg)
saveRDS(samples_out, file.path(out_dir, "trajectory_samples_angular.rds"), compress = "xz")

write_csv(validation, file.path(out_dir, "validation_published_descriptives.csv"))
write_csv(runs_matched, file.path(out_dir, "trajectory_runs.csv"))
write_csv(revisit_sensitivity, file.path(out_dir, "revisited_trials_sensitivity.csv"))
write_csv(new_summary, file.path(out_dir, "summary_new_measures.csv"))
write_csv(movement_summary, file.path(out_dir, "summary_movement_measures.csv"))
write_csv(core_summary, file.path(out_dir, "summary_core_measures.csv"))
write_csv(pair_report, file.path(out_dir, "session_pairing_report.csv"))
write_csv(inter_trial_gaps, file.path(out_dir, "inter_trial_gaps.csv"))

say("\n--- Files written ---")
for (f in c("analysis_df_revision.rds", "analysis_df_revision.csv",
            "trajectory_samples_angular.rds", "validation_published_descriptives.csv",
            "summary_new_measures.csv", "summary_movement_measures.csv",
            "summary_core_measures.csv", "session_pairing_report.csv",
            "trajectory_runs.csv", "revisited_trials_sensitivity.csv",
            "inter_trial_gaps.csv", "build_analysis_df_revision_log.txt")) {
  p <- file.path(out_dir, f)
  say(sprintf("  %s (%s)", p,
              ifelse(file.exists(p), paste0(round(file.size(p) / 1024), " KB"), "pending")))
}

say("\nFinal frame columns (", ncol(analysis_df), "):")
say(paste(names(analysis_df), collapse = ", "))
say("\nDone.")
close(log_con)
