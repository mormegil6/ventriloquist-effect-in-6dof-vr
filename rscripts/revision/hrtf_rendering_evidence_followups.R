#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Follow-ups to hrtf_rendering_evidence.R.
#
# The primary models show that the response tracks the true source elevation
# well beyond the visual cue. Before that can be read as evidence for binaural
# rendering, three routes to elevation that do not require an HRTF have to be
# examined, all of which exist because the task was free-roaming 6DoF with a
# continuously sounding source and a 3 m distance-attenuation curve:
#
#   (i)   loudness-gradient search: listeners walked to within a median of
#         0.14 m of the source, so proximity alone can fix its position;
#   (ii)  self-motion triangulation: azimuth from amplitude panning plus
#         distance from attenuation, sampled from several vantage points and
#         several ear heights, determines the source in three dimensions;
#   (iii) head tilt: roll and pitch rotate the panning plane, so elevation
#         leaks into the lateral panner even without spectral cues.
#
# Each route predicts that the auditory elevation weight should scale with the
# relevant behaviour (proximity, exploration, tilt). A weight that survives in
# the trials with least movement, least proximity and least tilt cannot be
# produced by any of them. A spectral signature across sound types is the
# positive prediction unique to HRTF rendering.
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
  library(emmeans)
})

project_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
out_dir     <- file.path(project_dir, "results_revision")

log_con <- file(file.path(out_dir, "hrtf_rendering_evidence_followups_log.txt"), open = "wt")
say <- function(...) {
  txt <- paste0(...); cat(txt, "\n", sep = ""); cat(txt, "\n", sep = "", file = log_con)
}
say_df <- function(x, digits = 4) {
  out <- capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE))
  cat(out, sep = "\n"); cat("\n")
  cat(out, sep = "\n", file = log_con); cat("\n", file = log_con)
}

say("=== Follow-ups: can the elevation effect be produced without an HRTF? ===")
say("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

# The primary script caches the head-frame geometry alongside the trial frame.
dat <- readRDS(file.path(out_dir, "hrtf_geometry_trials.rds"))
say("Trials: ", nrow(dat), "  participants: ", n_distinct(dat$participantId))

# ---------------------------------------------------------------------------
# Shared model helper: auditory weight on one axis, optionally in a subset.
# ---------------------------------------------------------------------------

weight_fit <- function(df, resp, src, cue, label) {
  d <- df %>% select(participantId, y = all_of(resp), s = all_of(src), c = all_of(cue)) %>%
    drop_na()
  m <- lmer(y ~ s + c + (1 | participantId), data = d, REML = TRUE)
  co <- summary(m)$coefficients
  ci <- suppressMessages(confint(m, parm = c("s", "c"), method = "Wald"))
  tibble(subset = label, n = nrow(d), n_participants = n_distinct(d$participantId),
         b_source = co["s", "Estimate"], ci_lo = ci["s", 1], ci_hi = ci["s", 2],
         p_source = co["s", "Pr(>|t|)"], b_cue = co["c", "Estimate"])
}

# Vertical axis in two forms: metres in the world frame, and head-referenced
# elevation in degrees (which is invariant to any purely radial correction of
# the response, so it cannot be produced by a distance cue alone).
vertical_forms <- list(
  world_metres = c("response_y", "sound_y", "flash_y"),
  head_elev_deg = c("resp_elev_tw", "src_elev_tw", "cue_elev_tw")
)

subset_scan <- function(df, label) {
  imap_dfr(vertical_forms, function(v, nm) {
    weight_fit(df, v[1], v[2], v[3], label) %>% mutate(measure = nm, .before = 1)
  })
}

# ---------------------------------------------------------------------------
# 1. Behavioural quartiles: does the auditory elevation weight need movement,
#    proximity or head tilt?
# ---------------------------------------------------------------------------

q_label <- function(x) cut(x, breaks = quantile(x, c(0, .25, .5, .75, 1)),
                           labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)

dat <- dat %>%
  mutate(q_path = q_label(total_path_length),
         q_close = q_label(dist_min_m),
         q_yrange = q_label(y_range),
         q_roll = q_label(roll_sd_deg),
         q_pitch = q_label(pitch_sd_deg))

quartile_scan <- bind_rows(
  map_dfr(levels(dat$q_path),   ~ subset_scan(filter(dat, q_path == .x),   paste("path length", .x))   %>% mutate(moderator = "exploration path length")),
  map_dfr(levels(dat$q_close),  ~ subset_scan(filter(dat, q_close == .x),  paste("closest approach", .x)) %>% mutate(moderator = "closest approach to source")),
  map_dfr(levels(dat$q_yrange), ~ subset_scan(filter(dat, q_yrange == .x), paste("head height range", .x)) %>% mutate(moderator = "vertical head excursion")),
  map_dfr(levels(dat$q_roll),   ~ subset_scan(filter(dat, q_roll == .x),   paste("head roll SD", .x))   %>% mutate(moderator = "head roll")),
  map_dfr(levels(dat$q_pitch),  ~ subset_scan(filter(dat, q_pitch == .x),  paste("head pitch SD", .x))  %>% mutate(moderator = "head pitch"))
)

say("\n--- 1. Auditory elevation weight within behavioural quartiles ---")
say("Q1 is the least movement / least proximity / least tilt quartile.")
say_df(quartile_scan %>% select(moderator, subset, measure, n, b_source, ci_lo, ci_hi, p_source))
write_csv(quartile_scan, file.path(out_dir, "hrtf_elevation_weight_quartiles.csv"))

# The conjunction: trials in the lower half on movement, proximity and tilt.
conservative <- dat %>%
  filter(total_path_length < median(total_path_length),
         dist_min_m > median(dist_min_m),
         roll_sd_deg < median(roll_sd_deg))
say("\n--- 1b. Conjunction subset: little walking, never approached closely, little roll ---")
say_df(subset_scan(conservative, "little walking + distant + little roll"))

# ---------------------------------------------------------------------------
# 2. Continuous moderation, so the quartile split is not doing the work
# ---------------------------------------------------------------------------

zs <- function(x) as.vector(scale(x))
mod_dat <- dat %>%
  mutate(sound_y_c = sound_y - mean(sound_y), flash_y_c = flash_y - mean(flash_y),
         log_path = zs(log(total_path_length)), log_close = zs(log(dist_min_m)),
         log_yrange = zs(log(y_range)), roll_z = zs(roll_sd_deg), pitch_z = zs(pitch_sd_deg))

moderators <- c("log_path", "log_close", "log_yrange", "roll_z", "pitch_z")
simple_slopes <- map_dfr(moderators, function(mv) {
  f <- as.formula(sprintf("response_y ~ sound_y_c * %s + flash_y_c + (1 | participantId)", mv))
  m <- lmer(f, data = mod_dat, REML = TRUE)
  co <- summary(m)$coefficients
  ix <- paste0("sound_y_c:", mv)
  # slope of sound_y at the 10th and 90th percentile of the moderator
  qs <- quantile(mod_dat[[mv]], c(.1, .9))
  V <- as.matrix(vcov(m))
  est <- co["sound_y_c", "Estimate"] + qs * co[ix, "Estimate"]
  se  <- sqrt(V["sound_y_c", "sound_y_c"] + qs^2 * V[ix, ix] + 2 * qs * V["sound_y_c", ix])
  tibble(moderator = mv,
         b_interaction = co[ix, "Estimate"], p_interaction = co[ix, "Pr(>|t|)"],
         slope_p10 = est[1], slope_p10_lo = est[1] - 1.96 * se[1], slope_p10_hi = est[1] + 1.96 * se[1],
         slope_p90 = est[2], slope_p90_lo = est[2] - 1.96 * se[2], slope_p90_hi = est[2] + 1.96 * se[2])
})
say("\n--- 2. Continuous moderation of the vertical auditory weight (world metres) ---")
say("slope_p10 / slope_p90: auditory weight at the 10th and 90th percentile of the moderator.")
say_df(simple_slopes)
write_csv(simple_slopes, file.path(out_dir, "hrtf_elevation_weight_moderation.csv"))

# ---------------------------------------------------------------------------
# 3. Spectral signature. HRTF elevation cues live in the high-frequency
#    spectral detail of the pinna transfer function, so broadband pink noise
#    should support elevation better than a narrowband harmonic flute tone.
#    Amplitude panning is frequency independent and predicts no difference.
#    The lateral axis is the control: panning serves it whatever the spectrum.
# ---------------------------------------------------------------------------

sound_scan <- map_dfr(levels(dat$soundType), function(s) {
  d <- filter(dat, soundType == s)
  bind_rows(
    weight_fit(d, "response_y", "sound_y", "flash_y", s) %>% mutate(axis = "vertical"),
    weight_fit(d, "response_x", "sound_x", "flash_x", s) %>% mutate(axis = "lateral x"),
    weight_fit(d, "response_z", "sound_z", "flash_z", s) %>% mutate(axis = "depth z")
  )
})
say("\n--- 3. Auditory weight by sound type ---")
say_df(sound_scan %>% select(axis, subset, n, b_source, ci_lo, ci_hi, p_source))
write_csv(sound_scan, file.path(out_dir, "hrtf_weight_by_soundtype.csv"))

m_sound_y <- lmer(response_y ~ sound_y * soundType + flash_y + (1 | participantId),
                  data = dat, REML = TRUE)
m_sound_x <- lmer(response_x ~ sound_x * soundType + flash_x + (1 | participantId),
                  data = dat, REML = TRUE)
say("Sound type x source interaction, vertical axis:")
say(paste(capture.output(anova(m_sound_y)), collapse = "\n"))
say("Sound type x source interaction, lateral axis:")
say(paste(capture.output(anova(m_sound_x)), collapse = "\n"))

trend_y <- emtrends(m_sound_y, ~ soundType, var = "sound_y")
say("\nVertical auditory weight per sound type (emtrends):")
say(paste(capture.output(print(trend_y)), collapse = "\n"))
write_csv(as_tibble(summary(trend_y)), file.path(out_dir, "hrtf_vertical_trend_by_soundtype.csv"))
write_csv(as_tibble(summary(pairs(trend_y))),
          file.path(out_dir, "hrtf_vertical_trend_soundtype_contrasts.csv"))

# ---------------------------------------------------------------------------
# 4. Median plane vs lateral sources. Spectral elevation cues are strongest
#    near the median plane; the tilt-projection route is weakest there,
#    because a source in the median plane stays in the median plane under
#    pitch and only leaves it under roll.
# ---------------------------------------------------------------------------

plane_scan <- bind_rows(
  subset_scan(filter(dat, abs(src_ux_tw) < quantile(abs(src_ux_tw), .5)),
              "source near the median plane"),
  subset_scan(filter(dat, abs(src_ux_tw) >= quantile(abs(src_ux_tw), .5)),
              "source lateral")
)
say("\n--- 4. Median-plane vs lateral sources ---")
say_df(plane_scan %>% select(measure, subset, n, b_source, ci_lo, ci_hi, p_source))
write_csv(plane_scan, file.path(out_dir, "hrtf_weight_median_plane.csv"))

# ---------------------------------------------------------------------------
# 5. How much response variance does the auditory source position carry once
#    the cue is controlled? Semi-partial R2 per axis, with a bootstrap CI over
#    participants.
# ---------------------------------------------------------------------------

r2m <- function(m) {
  vf <- var(as.vector(model.matrix(m) %*% fixef(m)))
  vr <- sum(unlist(lapply(VarCorr(m), function(v) v[1])))
  vf / (vf + vr + sigma(m)^2)
}
semipartial <- function(d, resp, src, cue) {
  dd <- d %>% select(participantId, y = all_of(resp), s = all_of(src), c = all_of(cue)) %>% drop_na()
  m  <- lmer(y ~ s + c + (1 | participantId), data = dd, REML = TRUE)
  m0 <- lmer(y ~ c + (1 | participantId), data = dd, REML = TRUE)
  r2m(m) - r2m(m0)
}

set.seed(20260810)
boot_participant <- function(d, resp, src, cue, B = 500) {
  ids <- unique(d$participantId)
  replicate(B, {
    take <- sample(ids, length(ids), replace = TRUE)
    db <- bind_rows(lapply(seq_along(take), function(i)
      d %>% filter(participantId == take[i]) %>% mutate(participantId = paste0("b", i))))
    tryCatch(semipartial(db, resp, src, cue), error = function(e) NA_real_)
  })
}

axes_spec <- tribble(
  ~axis, ~resp, ~src, ~cue,
  "vertical y",  "response_y", "sound_y", "flash_y",
  "lateral x",   "response_x", "sound_x", "flash_x",
  "depth z",     "response_z", "sound_z", "flash_z"
)
r2_res <- pmap_dfr(axes_spec, function(axis, resp, src, cue) {
  point <- semipartial(dat, resp, src, cue)
  bs <- boot_participant(dat, resp, src, cue)
  tibble(axis = axis, semipartial_R2_source = point,
         boot_lo = quantile(bs, .025, na.rm = TRUE),
         boot_hi = quantile(bs, .975, na.rm = TRUE))
})
say("\n--- 5. Variance in the response explained by the auditory source ---")
say("Semi-partial marginal R2 for the source term, cue held in the model.")
say("CI from a 500-replicate participant-level bootstrap.")
say_df(r2_res)
write_csv(r2_res, file.path(out_dir, "hrtf_semipartial_r2.csv"))

say("\nOutputs written to ", out_dir)
close(log_con)
