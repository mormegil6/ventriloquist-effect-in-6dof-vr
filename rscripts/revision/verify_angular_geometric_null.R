#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Does the angular re-analysis carry information beyond the metric fit?
#
# The angular predictor and the angular response are both evaluated from the
# SAME onset head pose, so to first order
#
#     theta_D ~ stimulusDisparity / d      theta_R ~ signedError_inline / d
#
# with the same head-to-source distance d in both denominators. Any dispersion
# in d therefore induces a positive theta_R-on-theta_D association even when
# the metric bias does not depend on the stimulus at all. This script measures
# how large that induced association is, so that the observed angular slope and
# the angular-versus-nominal model comparison can be read against the correct
# null rather than against zero.
#
#   G1  reconstructed geometry and the angular pipeline (self-contained)
#   G2  permutation null: the in-line metric displacement is permuted within
#       participant, destroying any dependence on the stimulus while leaving
#       every head pose, source, cue and error magnitude untouched
#   G3  generative placebos with a known metric truth
#   G4  is the one scale-crossing result (time-weighted angle adds to nominal
#       separation on the METRIC response) an angle effect or a distance effect?
#   G5  distance-restricted refits of the primary angular model
#   G6  intercept of the primary angular model at zero angular disparity
#
# Inputs : results_revision/analysis_df_revision.rds
#          results_revision/trajectory_samples_angular.rds
# Outputs: results_revision/verify_angular_null_*.csv, verify_angular_null_log.txt
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir  <- file.path(proj_dir, "results_revision")

log_con <- file(file.path(res_dir, "verify_angular_null_log.txt"), open = "wt")
say <- function(...) { txt <- paste0(...); cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_con) }
say_df <- function(df, digits = 5)
  for (l in capture.output(print(as.data.frame(df), digits = digits, row.names = FALSE))) say(l)
rule <- function(t) { say(""); say(strrep("=", 78)); say(t); say(strrep("=", 78)) }

rule("GEOMETRIC NULL FOR THE ANGULAR RE-ANALYSIS")
say("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", R.version.string)

# ---------------------------------------------------------------------------
# G1  geometry
# ---------------------------------------------------------------------------
rule("G1: GEOMETRY AND THE ANGULAR PIPELINE")

analysis_df <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
samples     <- readRDS(file.path(res_dir, "trajectory_samples_angular.rds"))

onset <- samples %>%
  group_by(participantId, trialSequenceNum) %>%
  slice_min(sample_idx, n = 1, with_ties = FALSE) %>%
  select(participantId, trialSequenceNum, hx = px, hy = py, hz = pz) %>%
  ungroup()

dat <- analysis_df %>%
  left_join(onset, by = c("participantId", "trialSequenceNum")) %>%
  mutate(stimDisp_cm     = stimulusDisparity_m - mean(stimulusDisparity_m),
         angDisp_onset_cm = angDisp_onset_deg - mean(angDisp_onset_deg),
         angDisp_mean_cm  = angDisp_mean_deg  - mean(angDisp_mean_deg),
         invdist_mean_c   = 1 / dist_mean_m   - mean(1 / dist_mean_m),
         dist_mean_c      = dist_mean_m       - mean(dist_mean_m))

unitise <- function(x, y, z) { n <- sqrt(x^2 + y^2 + z^2); list(x = x / n, y = y / n, z = z / n) }
crx <- function(a, b) a$y * b$z - a$z * b$y
cry <- function(a, b) a$z * b$x - a$x * b$z
crz <- function(a, b) a$x * b$y - a$y * b$x
dt  <- function(a, b) a$x * b$x + a$y * b$y + a$z * b$z

# fixed part of the geometry: unit vectors to source and cue at the onset pose
US <- with(dat, unitise(sound_x - hx, sound_y - hy, sound_z - hz))
UC <- with(dat, unitise(flash_x - hx, flash_y - hy, flash_z - hz))
NX <- crx(US, UC); NY <- cry(US, UC); NZ <- crz(US, UC)
NN <- sqrt(NX^2 + NY^2 + NZ^2); NX <- NX / NN; NY <- NY / NN; NZ <- NZ / NN

ang_resp <- function(rx, ry, rz) {
  UR <- unitise(rx - dat$hx, ry - dat$hy, rz - dat$hz)
  atan2(crx(US, UR) * NX + cry(US, UR) * NY + crz(US, UR) * NZ, dt(US, UR)) * 180 / pi
}

# decomposition of the observed response into an in-line and an orthogonal part
u_sf <- with(dat, { dx <- flash_x - sound_x; dy <- flash_y - sound_y; dz <- flash_z - sound_z
  n <- sqrt(dx^2 + dy^2 + dz^2); cbind(dx / n, dy / n, dz / n) })
rel   <- with(dat, cbind(response_x - sound_x, response_y - sound_y, response_z - sound_z))
along <- rowSums(rel * u_sf)
orth  <- rel - along * u_sf
say(sprintf("In-line component reproduces signedError_m to %.3g m",
            max(abs(along - dat$signedError_m))))

build_resp <- function(inline) {
  p <- inline * u_sf + orth
  list(x = dat$sound_x + p[, 1], y = dat$sound_y + p[, 2], z = dat$sound_z + p[, 3])
}
r0 <- build_resp(dat$signedError_m)
dat$angResp_onset_deg <- ang_resp(r0$x, r0$y, r0$z)
say(sprintf("Rebuilt angular response, mean %.4f deg, SD %.4f deg", mean(dat$angResp_onset_deg),
            sd(dat$angResp_onset_deg)))

ctrl  <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
f_ang <- angResp_onset_deg ~ angDisp_onset_cm + soundType + (1 + angDisp_onset_cm || participantId)
f_nom <- angResp_onset_deg ~ stimDisp_cm + soundType + (1 + stimDisp_cm || participantId)

fit_slope <- function(d, f, term) {
  m <- suppressWarnings(lmer(f, data = d, REML = TRUE, control = ctrl))
  s <- summary(m)$coefficients[term, ]
  c(est = unname(s["Estimate"]), se = unname(s["Std. Error"]), t = unname(s["t value"]))
}
aic_ml <- function(d, f) AIC(suppressWarnings(lmer(f, data = d, REML = FALSE, control = ctrl)))

obs_ang  <- fit_slope(dat, f_ang, "angDisp_onset_cm")
obs_dAIC <- aic_ml(dat, f_nom) - aic_ml(dat, f_ang)
say(sprintf("Observed: angular slope %.4f (SE %.4f, t %.2f); dAIC(nominal - angular) = %+.2f",
            obs_ang["est"], obs_ang["se"], obs_ang["t"], obs_dAIC))

# ---------------------------------------------------------------------------
# G2  permutation null
# ---------------------------------------------------------------------------
rule("G2: PERMUTATION NULL FOR THE ANGULAR SLOPE")

say("The in-line metric displacement is permuted among the trials of the same")
say("participant. This destroys every relation between the response and the")
say("stimulus while preserving each head pose, each source and cue position, the")
say("orthogonal scatter and the participant's own distribution of signed errors.")
say("A capture effect of exactly zero therefore holds by construction, on both")
say("scales. The metric model is shown alongside as the control.")

set.seed(20260806)
nperm <- 300
pid <- dat$participantId
perm <- t(vapply(seq_len(nperm), function(i) {
  idx <- seq_len(nrow(dat))
  for (p in unique(pid)) { k <- which(pid == p); idx[k] <- sample(k) }
  d2 <- dat
  d2$signedError_m <- dat$signedError_m[idx]
  rr <- build_resp(d2$signedError_m)
  d2$angResp_onset_deg <- ang_resp(rr$x, rr$y, rr$z)
  a <- fit_slope(d2, f_ang, "angDisp_onset_cm")
  m <- fit_slope(d2, signedError_m ~ stimDisp_cm + soundType + (1 + stimDisp_cm | participantId),
                 "stimDisp_cm")
  c(ang = unname(a["est"]), ang_t = unname(a["t"]),
    dAIC = aic_ml(d2, f_nom) - aic_ml(d2, f_ang),
    metric = unname(m["est"]), metric_t = unname(m["t"]))
}, numeric(5)))

qq <- function(x) sprintf("mean %+.4f, SD %.4f, 2.5%%..97.5%% [%+.4f, %+.4f]",
                          mean(x), sd(x), quantile(x, .025), quantile(x, .975))
say("")
say("Null distribution over ", nperm, " within-participant permutations:")
say("  angular slope (deg/deg)          : ", qq(perm[, "ang"]))
say("  angular t statistic              : ", qq(perm[, "ang_t"]))
say("  dAIC (nominal - angular)         : ", qq(perm[, "dAIC"]))
say("  metric slope, control (m/m)      : ", qq(perm[, "metric"]))
say("  metric t statistic, control      : ", qq(perm[, "metric_t"]))
say("")
say(sprintf("  permutations with an angular slope >= the observed %.4f : %d of %d (p = %.4f)",
            obs_ang["est"], sum(perm[, "ang"] >= obs_ang["est"]), nperm,
            (1 + sum(perm[, "ang"] >= obs_ang["est"])) / (nperm + 1)))
say(sprintf("  permutations with dAIC >= the observed %+.2f            : %d of %d (p = %.4f)",
            obs_dAIC, sum(perm[, "dAIC"] >= obs_dAIC), nperm,
            (1 + sum(perm[, "dAIC"] >= obs_dAIC)) / (nperm + 1)))
say(sprintf("  permutations with a metric slope >= the observed 0.3466  : %d of %d",
            sum(perm[, "metric"] >= 0.34660), nperm))
write.csv(as.data.frame(perm), file.path(res_dir, "verify_angular_null_permutations.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# G3  generative placebos
# ---------------------------------------------------------------------------
rule("G3: GENERATIVE PLACEBOS WITH A KNOWN METRIC TRUTH")

pub_b0 <- 0.0086739; pub_b1 <- 0.34660          # published metric fixed effects
gen <- list(
  "metric truth 0.347 x disparity"        = pub_b0 + pub_b1 * dat$stimulusDisparity_m,
  "constant shift = mean signed error"    = rep(mean(dat$signedError_m), nrow(dat)),
  "constant shift = 0 (no capture)"       = rep(0, nrow(dat)),
  "metric truth with slope 0 (intercept only)" = rep(pub_b0, nrow(dat)),
  "metric truth doubled (0.693 x disparity)"   = pub_b0 + 2 * pub_b1 * dat$stimulusDisparity_m)

gen_tab <- bind_rows(lapply(names(gen), function(nm) {
  d2 <- dat
  d2$signedError_m <- gen[[nm]]
  rr <- build_resp(d2$signedError_m)
  d2$angResp_onset_deg <- ang_resp(rr$x, rr$y, rr$z)
  a <- fit_slope(d2, f_ang, "angDisp_onset_cm")
  dA <- aic_ml(d2, f_nom) - aic_ml(d2, f_ang)
  say(sprintf("  %-44s angular slope %+.4f (SE %.4f, t %5.2f)  dAIC %+7.1f",
              nm, a["est"], a["se"], a["t"], dA))
  data.frame(placebo = nm, angular_slope = a["est"], se = a["se"], t = a["t"], dAIC = dA,
             row.names = NULL)
}))
say(sprintf("  %-44s angular slope %+.4f (SE %.4f, t %5.2f)  dAIC %+7.1f",
            "OBSERVED DATA", obs_ang["est"], obs_ang["se"], obs_ang["t"], obs_dAIC))
write.csv(gen_tab, file.path(res_dir, "verify_angular_null_placebos.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# G4  angle or distance on the metric response?
# ---------------------------------------------------------------------------
rule("G4: THE SCALE-CROSSING RESULT -- EXPERIENCED ANGLE OR LISTENING DISTANCE?")

say("On the METRIC response the response carries no 1/d factor, so the reported")
say("chi2(1) = 12.35 for adding the time-weighted angle to nominal separation is")
say("not vulnerable to the shared-denominator problem. It is however vulnerable to")
say("a simpler reading: the time-weighted angle is largely a function of how close")
say("the listener stood. The tests below add distance instead of, and alongside,")
say("the angle.")

mlfit <- function(f) suppressWarnings(lmer(f, data = dat, REML = FALSE, control = ctrl))
base  <- mlfit(signedError_m ~ stimDisp_cm + soundType + (1 + stimDisp_cm | participantId))
add   <- function(lbl, f) {
  a <- anova(base, mlfit(f))
  say(sprintf("  %-52s chi2(%d) = %7.4f  p = %.4g  dAIC %+7.2f",
              lbl, a$Df[2], a$Chisq[2], a$`Pr(>Chisq)`[2], a$AIC[2] - a$AIC[1]))
  data.frame(added = lbl, chisq = a$Chisq[2], df = a$Df[2], p = a$`Pr(>Chisq)`[2])
}
g4 <- bind_rows(
  add("+ time-weighted angular disparity",
      signedError_m ~ stimDisp_cm + angDisp_mean_cm + soundType + (1 + stimDisp_cm | participantId)),
  add("+ mean head-to-source distance",
      signedError_m ~ stimDisp_cm + dist_mean_c + soundType + (1 + stimDisp_cm | participantId)),
  add("+ inverse mean head-to-source distance",
      signedError_m ~ stimDisp_cm + invdist_mean_c + soundType + (1 + stimDisp_cm | participantId)),
  add("+ inverse distance AND time-weighted angle",
      signedError_m ~ stimDisp_cm + invdist_mean_c + angDisp_mean_cm + soundType +
        (1 + stimDisp_cm | participantId)))
say("")
say("Angle added on top of nominal separation AND inverse distance:")
b2 <- mlfit(signedError_m ~ stimDisp_cm + invdist_mean_c + soundType + (1 + stimDisp_cm | participantId))
f2 <- mlfit(signedError_m ~ stimDisp_cm + invdist_mean_c + angDisp_mean_cm + soundType +
              (1 + stimDisp_cm | participantId))
a2 <- anova(b2, f2)
say(sprintf("  chi2(%d) = %.4f, p = %.4g", a2$Df[2], a2$Chisq[2], a2$`Pr(>Chisq)`[2]))
say("Inverse distance added on top of nominal separation AND the angle:")
b3 <- mlfit(signedError_m ~ stimDisp_cm + angDisp_mean_cm + soundType + (1 + stimDisp_cm | participantId))
a3 <- anova(b3, f2)
say(sprintf("  chi2(%d) = %.4f, p = %.4g", a3$Df[2], a3$Chisq[2], a3$`Pr(>Chisq)`[2]))
write.csv(g4, file.path(res_dir, "verify_angular_null_metric_scale.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# G5  distance-restricted refits
# ---------------------------------------------------------------------------
rule("G5: THE PRIMARY ANGULAR SLOPE RESTRICTED BY ONSET LISTENING DISTANCE")

g5 <- bind_rows(lapply(c(0, 0.3, 0.5, 0.75, 1.0), function(thr) {
  d3 <- dat %>% filter(dist_at_onset_m >= thr)
  m <- suppressWarnings(lmer(f_ang, data = d3, REML = TRUE, control = ctrl))
  s <- summary(m)$coefficients["angDisp_onset_cm", ]
  say(sprintf("  onset distance >= %.2f m : n = %3d, beta = %+.4f, SE = %.4f, t(%.1f) = %5.2f, p = %.4f, SD(angDisp) = %5.1f",
              thr, nobs(m), s["Estimate"], s["Std. Error"], s["df"], s["t value"],
              s["Pr(>|t|)"], sd(d3$angDisp_onset_deg)))
  data.frame(min_distance_m = thr, n = nobs(m), estimate = s["Estimate"], se = s["Std. Error"],
             p = s["Pr(>|t|)"], sd_angDisp = sd(d3$angDisp_onset_deg), row.names = NULL)
}))
write.csv(g5, file.path(res_dir, "verify_angular_null_distance_restricted.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# G6  intercept at zero angular disparity
# ---------------------------------------------------------------------------
rule("G6: PREDICTED RESPONSE AT ZERO ANGULAR DISPARITY")

m_unc <- suppressWarnings(lmer(
  angResp_onset_deg ~ angDisp_onset_deg + soundType + (1 + angDisp_onset_deg || participantId),
  data = dat, REML = TRUE, control = ctrl))
s <- summary(m_unc)$coefficients
say_df(cbind(term = rownames(s), as.data.frame(s)))
say(sprintf("Uncentred intercept = %.4f deg (p = %.4f): the fitted line does not pass",
            s["(Intercept)", "Estimate"], s["(Intercept)", "Pr(>|t|)"]))
say("through the origin, so a proportional reading of the slope is an approximation.")

say(""); say("Done.")
close(log_con)
