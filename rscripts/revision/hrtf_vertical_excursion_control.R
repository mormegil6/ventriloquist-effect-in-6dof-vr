#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# The elevation weight is strongly moderated by how far the listener moved the
# head vertically within a trial, and not at all by how far they walked. That
# dissociation is the crux of the rendering question, so it is checked here
# against the obvious alternatives:
#
#   - generic engagement: if listeners who barely moved were simply careless,
#     the weight on the visual cue should fall with the weight on the source;
#   - axis specificity: engagement predicts the same moderation on the lateral
#     and depth axes, whereas a geometric elevation-from-ear-height account
#     predicts it on the vertical axis only;
#   - collinearity of vertical excursion with proximity and with walking.
#
# Run with R 4.4 (arm64).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tidyr); library(readr)
  library(lme4); library(lmerTest)
})

project_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
out_dir <- file.path(project_dir, "results_revision")
log_con <- file(file.path(out_dir, "hrtf_vertical_excursion_control_log.txt"), open = "wt")
say <- function(...) { txt <- paste0(...); cat(txt, "\n", sep = ""); cat(txt, "\n", sep = "", file = log_con) }
say_df <- function(x, digits = 4) {
  out <- capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE))
  cat(out, sep = "\n"); cat("\n"); cat(out, sep = "\n", file = log_con); cat("\n", file = log_con)
}

dat <- readRDS(file.path(out_dir, "hrtf_geometry_trials.rds"))
say("=== Vertical head excursion and the auditory elevation weight ===")
say("Run date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

zs <- function(x) as.vector(scale(x))
dat <- dat %>%
  mutate(q_yrange = cut(y_range, quantile(y_range, c(0, .25, .5, .75, 1)),
                        labels = paste0("Q", 1:4), include.lowest = TRUE),
         log_yrange = zs(log(y_range)), log_path = zs(log(total_path_length)),
         log_close = zs(log(dist_min_m)), log_rt = zs(log(responseTime_s)))

say("\nCollinearity of the trial-level movement measures (Pearson r):")
say_df(as.data.frame(round(cor(dat[, c("log_yrange", "log_path", "log_close", "log_rt")]), 3)))

# --- source and cue weights per quartile, all three axes -------------------

axis_spec <- tribble(
  ~axis,        ~resp,        ~src,      ~cue,
  "vertical y", "response_y", "sound_y", "flash_y",
  "lateral x",  "response_x", "sound_x", "flash_x",
  "depth z",    "response_z", "sound_z", "flash_z"
)

fit_weights <- function(df, resp, src, cue) {
  d <- df %>% select(participantId, y = all_of(resp), s = all_of(src), c = all_of(cue)) %>% drop_na()
  m <- lmer(y ~ s + c + (1 | participantId), data = d, REML = TRUE)
  co <- summary(m)$coefficients
  ci <- suppressMessages(confint(m, parm = c("s", "c"), method = "Wald"))
  tibble(n = nrow(d),
         b_source = co["s", "Estimate"], src_lo = ci["s", 1], src_hi = ci["s", 2],
         p_source = co["s", "Pr(>|t|)"],
         b_cue = co["c", "Estimate"], cue_lo = ci["c", 1], cue_hi = ci["c", 2])
}

quart <- expand_grid(q = levels(dat$q_yrange), axis_spec) %>%
  pmap_dfr(function(q, axis, resp, src, cue) {
    fit_weights(filter(dat, q_yrange == q), resp, src, cue) %>%
      mutate(q_yrange = q, axis = axis, .before = 1)
  })
say("\n--- Source and cue weights by quartile of vertical head excursion ---")
say("Q1 = least vertical head movement. A careless-responding account predicts")
say("the cue weight to fall with the source weight; it does not.")
say_df(quart %>% select(axis, q_yrange, n, b_source, src_lo, src_hi, p_source, b_cue, cue_lo, cue_hi))
write_csv(quart, file.path(out_dir, "hrtf_yrange_quartile_weights.csv"))

say("\nMedian vertical head excursion per quartile (m):")
say_df(dat %>% group_by(q_yrange) %>%
         summarise(median_y_range_m = median(y_range), min = min(y_range), max = max(y_range),
                   median_path_m = median(total_path_length),
                   median_min_dist_m = median(dist_min_m), .groups = "drop"))

# --- axis specificity of the moderation, and adjustment for covariates -----

mod_table <- pmap_dfr(axis_spec, function(axis, resp, src, cue) {
  d <- dat %>% mutate(s = .data[[src]] - mean(.data[[src]]),
                      c = .data[[cue]] - mean(.data[[cue]]), y = .data[[resp]])
  m_raw <- lmer(y ~ s * log_yrange + c + (1 | participantId), data = d, REML = TRUE)
  m_adj <- lmer(y ~ s * log_yrange + s * log_path + s * log_close + s * log_rt + c +
                  (1 | participantId), data = d, REML = TRUE)
  co_r <- summary(m_raw)$coefficients; co_a <- summary(m_adj)$coefficients
  tibble(axis = axis,
         b_int_raw = co_r["s:log_yrange", "Estimate"], se_raw = co_r["s:log_yrange", "Std. Error"],
         p_raw = co_r["s:log_yrange", "Pr(>|t|)"],
         b_int_adj = co_a["s:log_yrange", "Estimate"], se_adj = co_a["s:log_yrange", "Std. Error"],
         p_adj = co_a["s:log_yrange", "Pr(>|t|)"])
})
say("\n--- Moderation of the source weight by vertical head excursion, per axis ---")
say("adj = adjusted for walking distance, closest approach and response time,")
say("each also interacting with the source position.")
say_df(mod_table)
write_csv(mod_table, file.path(out_dir, "hrtf_yrange_moderation_by_axis.csv"))

# --- simple slopes on the vertical axis in the adjusted model --------------

d <- dat %>% mutate(s = sound_y - mean(sound_y), c = flash_y - mean(flash_y), y = response_y)
m_adj <- lmer(y ~ s * log_yrange + s * log_path + s * log_close + s * log_rt + c +
                (1 | participantId), data = d, REML = TRUE)
V <- as.matrix(vcov(m_adj)); co <- summary(m_adj)$coefficients
qs <- quantile(d$log_yrange, c(.05, .1, .25, .5, .75, .9))
slopes <- tibble(
  pct = names(qs), log_yrange = as.numeric(qs),
  y_range_m = as.numeric(quantile(d$y_range, c(.05, .1, .25, .5, .75, .9))),
  estimate = co["s", "Estimate"] + as.numeric(qs) * co["s:log_yrange", "Estimate"],
  se = sqrt(V["s", "s"] + as.numeric(qs)^2 * V["s:log_yrange", "s:log_yrange"] +
              2 * as.numeric(qs) * V["s", "s:log_yrange"])
) %>% mutate(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)
say("\n--- Adjusted vertical auditory weight across the range of head excursion ---")
say_df(slopes)
write_csv(slopes, file.path(out_dir, "hrtf_yrange_simple_slopes.csv"))

# --- head-referenced angular version of the same slopes --------------------

da <- dat %>% mutate(s = src_elev_tw - mean(src_elev_tw),
                     c = cue_elev_tw - mean(cue_elev_tw), y = resp_elev_tw)
m_ang <- lmer(y ~ s * log_yrange + s * log_close + c + (1 | participantId), data = da, REML = TRUE)
say("\n--- Head-referenced elevation (degrees), same moderation ---")
say(paste(capture.output(print(summary(m_ang)$coefficients, digits = 4)), collapse = "\n"))

# --- how much elevation information is left with a near-static head? -------
static_head <- dat %>% filter(y_range < quantile(y_range, .2))
say("\n--- Trials in the lowest quintile of vertical head excursion ---")
say(sprintf("n = %d trials, %d participants, y_range < %.3f m",
            nrow(static_head), n_distinct(static_head$participantId),
            quantile(dat$y_range, .2)))
static_res <- pmap_dfr(axis_spec, function(axis, resp, src, cue)
  fit_weights(static_head, resp, src, cue) %>% mutate(axis = axis, .before = 1))
say_df(static_res)
write_csv(static_res, file.path(out_dir, "hrtf_static_head_weights.csv"))

# The same trials in head-referenced angular form, and a check that the source
# elevation range is not restricted in this subset.
static_ang <- bind_rows(
  fit_weights(static_head, "resp_elev_tw", "src_elev_tw", "cue_elev_tw") %>%
    mutate(frame = "time-weighted", .before = 1),
  fit_weights(static_head, "resp_elev_onset", "src_elev_onset", "cue_elev_onset") %>%
    mutate(frame = "onset", .before = 1)
)
say("\nHead-referenced elevation (degrees) in the near-static-head trials:")
say_df(static_ang)
write_csv(static_ang, file.path(out_dir, "hrtf_static_head_elevation_deg.csv"))

say("Source elevation range by quartile of vertical head excursion:")
say_df(dat %>% group_by(q_yrange) %>%
         summarise(sd_sound_y = sd(sound_y), range_sound_y = diff(range(sound_y)),
                   sd_src_elev_deg = sd(src_elev_tw),
                   range_src_elev_deg = diff(range(src_elev_tw)), .groups = "drop"))

# Is the vertical weight in these trials reliably smaller than the lateral one?
d2 <- static_head %>%
  transmute(participantId,
            r_y = response_y, s_y = sound_y, c_y = flash_y,
            r_x = response_x, s_x = sound_x, c_x = flash_x) %>%
  pivot_longer(-participantId, names_to = c(".value", "axis"), names_pattern = "(.)_(x|y)$")
say("\nVertical vs lateral weight contrast within the near-static-head trials")
say("(axis is coded x = lateral as reference, so s:axisy is the difference):")
m_cmp <- lmer(r ~ s * axis + c * axis + (1 | participantId), data = d2, REML = TRUE)
say(paste(capture.output(print(summary(m_cmp)$coefficients, digits = 4)), collapse = "\n"))

# ---------------------------------------------------------------------------
# Is vertical head excursion itself a response to source elevation, i.e. is
# the moderator partly an outcome of knowing where the source was?
# ---------------------------------------------------------------------------

h0 <- median(dat$head_y_tw)
search_dat <- dat %>%
  mutate(src_offset_abs = abs(sound_y - h0), cue_offset_abs = abs(flash_y - h0),
         src_above = sound_y > h0)

say("\n--- Was vertical head excursion driven by the source position? ---")
m_search <- lmer(log(y_range) ~ src_offset_abs + cue_offset_abs + (1 | participantId),
                 data = search_dat, REML = TRUE)
say(paste(capture.output(print(summary(m_search)$coefficients, digits = 4)), collapse = "\n"))
m_height <- lmer(head_y_tw ~ sound_y + flash_y + (1 | participantId), data = dat, REML = TRUE)
say("\nMean head height during the trial as a function of source and cue height:")
say(paste(capture.output(print(summary(m_height)$coefficients, digits = 4)), collapse = "\n"))

say("\n--- Moderation split by whether the source was above or below head height ---")
split_res <- map_dfr(c(TRUE, FALSE), function(ab) {
  d <- search_dat %>% filter(src_above == ab) %>%
    mutate(s = sound_y - mean(sound_y), c = flash_y - mean(flash_y))
  m <- lmer(response_y ~ s * log_yrange + c + (1 | participantId), data = d, REML = TRUE)
  co <- summary(m)$coefficients
  tibble(source = if (ab) "above head height" else "below head height", n = nrow(d),
         b_source_at_mean = co["s", "Estimate"],
         b_interaction = co["s:log_yrange", "Estimate"],
         se = co["s:log_yrange", "Std. Error"], p = co["s:log_yrange", "Pr(>|t|)"])
})
say_df(split_res)
write_csv(split_res, file.path(out_dir, "hrtf_yrange_moderation_by_source_side.csv"))

say("\nOutputs written to ", out_dir)
close(log_con)
