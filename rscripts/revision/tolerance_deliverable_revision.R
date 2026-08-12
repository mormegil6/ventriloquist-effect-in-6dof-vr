# tolerance_deliverable_revision.R
#
# Developer-facing disparity-tolerance deliverable for the TVCG revision (Reviewer 2).
#
# (A) Predicted-displacement lookup table across 15-70 cm nominal disparity and the
#     corresponding egocentric angular range, for the population-average listener and
#     for the 10th / 90th percentile listener derived from the random-slope SD.
# (B) Capture-probability model: per-trial binary outcome (response closer to the
#     visual cue than to the true source) as a function of angular and nominal
#     disparity, with 50% / 75% crossing points where they are identified.
#
# Input : results_revision/analysis_df_revision.rds  (built by build_analysis_df_revision.R)
# Output: results_revision/tolerance_*.csv, capture_*.csv, tolerance_deliverable_log.txt
#
# Deterministic: all resampling is seeded. Run with R 4.4 (arm64).

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(DHARMa)
  library(dplyr)
  library(tidyr)
})

set.seed(20260805)

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir  <- file.path(proj_dir, "results_revision")
log_path <- file.path(res_dir, "tolerance_deliverable_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78))

say("tolerance_deliverable_revision.R")
say("run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
say("R: ", R.version.string)

# Number of parametric-bootstrap replicates. 2000 is enough for 95% percentile CIs.
N_BOOT <- 2000
Z90 <- qnorm(0.90)  # 1.2816: multiplier for the 10th / 90th percentile listener

# ---------------------------------------------------------------------------
# 0. Data
# ---------------------------------------------------------------------------
rule("0. DATA")

d <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
say("trials: ", nrow(d), "   participants: ", dplyr::n_distinct(d$participantId))

d <- d %>%
  mutate(
    disp_cm      = stimulusDisparity_m * 100,
    ang_deg      = angDisp_mean_deg,      # time-weighted mean experienced angle (pre-specified primary)
    angOnset_deg = angDisp_onset_deg,     # angle at trial onset (sensitivity)
    ang10        = ang_deg / 10,          # per-10-degree scaling for numerical stability
    angOnset10   = angOnset_deg / 10,
    disp10       = disp_cm / 10           # per-10-cm scaling
  )

mean_disp10 <- mean(d$disp10)
mean_ang10  <- mean(d$ang10)
d$disp10c <- d$disp10 - mean_disp10
d$ang10c  <- d$ang10  - mean_ang10

say("nominal disparity (cm): range ", paste(round(range(d$disp_cm), 1), collapse = "-"),
    ", mean ", round(mean(d$disp_cm), 2))
say("experienced angular disparity, time-weighted mean (deg): range ",
    paste(round(range(d$ang_deg), 1), collapse = "-"), ", mean ", round(mean(d$ang_deg), 2),
    ", median ", round(median(d$ang_deg), 2))
say("angular disparity at onset (deg): mean ", round(mean(d$angOnset_deg), 2),
    ", median ", round(median(d$angOnset_deg), 2))
say("head-to-source distance, per-trial mean (m): mean ", round(mean(d$dist_mean_m), 3),
    ", median ", round(median(d$dist_mean_m), 3))
say("sound type balance: ", paste(names(table(d$soundType)), as.integer(table(d$soundType)),
                                  sep = "=", collapse = ", "))

# ---------------------------------------------------------------------------
# 1. Binary capture outcome, and the identity between the two candidate criteria
# ---------------------------------------------------------------------------
rule("1. CAPTURE OUTCOME DEFINITION")

# Geometric criterion: is the response closer (in 3D) to the visual cue than to the
# true source? This is the nearest-source rule under the perpendicular bisector plane.
dist_to_src   <- sqrt((d$response_x - d$sound_x)^2 + (d$response_y - d$sound_y)^2 +
                        (d$response_z - d$sound_z)^2)
dist_to_flash <- sqrt((d$response_x - d$flash_x)^2 + (d$response_y - d$flash_y)^2 +
                        (d$response_z - d$flash_z)^2)
d$capture_geom <- as.integer(dist_to_flash < dist_to_src)

# Projection criterion: did the response traverse more than half of the sound-to-flash gap?
d$capture_bias <- as.integer(d$ventriloquistBias > 0.5)

# These are algebraically the same event: |r-f| < |r-s|  <=>  (r-s).(f-s) > |f-s|^2 / 2
# <=> signedError > disparity/2 <=> ventriloquistBias > 0.5. Verified numerically here.
say("agreement between geometric and ventriloquistBias > 0.5 criteria: ",
    sum(d$capture_geom == d$capture_bias), "/", nrow(d),
    "  (identical events, as expected algebraically)")
print(table(geometric = d$capture_geom, bias_gt_0.5 = d$capture_bias))

d$capture <- d$capture_geom
say("\nbase rate, primary criterion (capture): ", round(mean(d$capture), 4),
    "  (", sum(d$capture), "/", nrow(d), " trials)")

# Robustness criteria. Because the geometric criterion IS bias > 0.5, a genuinely
# softer/stricter check requires a different traversal fraction.
d$capture_soft   <- as.integer(d$ventriloquistBias > 0.25)
d$capture_strict <- as.integer(d$ventriloquistBias > 0.75)
say("base rate, softer criterion (bias > 0.25): ", round(mean(d$capture_soft), 4))
say("base rate, stricter criterion (bias > 0.75): ", round(mean(d$capture_strict), 4))

pp_cap <- tapply(d$capture, d$participantId, mean)
say("\nper-participant capture rate: min ", round(min(pp_cap), 3),
    ", median ", round(median(pp_cap), 3), ", max ", round(max(pp_cap), 3))

# ---------------------------------------------------------------------------
# 2. Primary signed-error LMMs (metric and angular)
# ---------------------------------------------------------------------------
rule("2. SIGNED-ERROR MIXED MODELS")

# 2a. Published primary model: metric disparity in metres, uncentred (reproduces the
#     manuscript's beta = 0.347, tau_slope = 0.190).
m_metric <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                   (1 + stimulusDisparity_m | participantId),
                 data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
say("\n--- Model M1: signedError_m ~ stimulusDisparity_m + soundType + (1 + stimulusDisparity_m | participantId) ---")
print(summary(m_metric))
msg_metric <- m_metric@optinfo$conv$lme4$messages
say("convergence messages M1: ", if (length(msg_metric)) paste(msg_metric, collapse = " | ") else "none")

# 2b. Angular model. Fitted on a per-10-degree scale: on the raw per-degree scale the
#     identical fit triggers a boundary-gradient warning purely from predictor scaling.
m_ang <- lmer(signedError_m ~ ang10 + soundType + (1 + ang10 | participantId),
              data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
say("\n--- Model M2: signedError_m ~ ang10 + soundType + (1 + ang10 | participantId), ang10 = deg/10 ---")
print(summary(m_ang))
msg_ang <- m_ang@optinfo$conv$lme4$messages
say("convergence messages M2: ", if (length(msg_ang)) paste(msg_ang, collapse = " | ") else "none")

# Documented for the record: the same model on the raw per-degree scale.
m_ang_raw <- lmer(signedError_m ~ ang_deg + soundType + (1 + ang_deg | participantId),
                  data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
msg_ang_raw <- m_ang_raw@optinfo$conv$lme4$messages
say("convergence messages M2 on raw per-degree scale: ",
    if (length(msg_ang_raw)) paste(msg_ang_raw, collapse = " | ") else "none")

# Onset-angle sensitivity model.
m_ang_onset <- lmer(signedError_m ~ angOnset10 + soundType + (1 + angOnset10 | participantId),
                    data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
say("\n--- Model M2b (sensitivity): onset angle instead of time-weighted mean ---")
print(round(summary(m_ang_onset)$coefficients, 5))
msg_ang_onset <- m_ang_onset@optinfo$conv$lme4$messages
say("convergence messages M2b: ",
    if (length(msg_ang_onset)) paste(msg_ang_onset, collapse = " | ") else "none")

# Direct proportional-scaling test: if displacement is a constant fraction of the gap,
# the traversed proportion should not depend on disparity.
m_prop <- lmerTest::lmer(ventriloquistBias ~ disp_cm + (1 | participantId), data = d)
say("\n--- Proportional-scaling check: ventriloquistBias ~ disp_cm + (1 | participantId) ---")
print(round(summary(m_prop)$coefficients, 6))

# ---------------------------------------------------------------------------
# 3. Listener-level prediction machinery
# ---------------------------------------------------------------------------
rule("3. LISTENER PROFILES")

# Predictions marginalise over sound type. The design is perfectly balanced
# (186 trials per level), so the marginal prediction is the equal-weight mean over levels.
sound_offset <- function(fx) {
  st <- fx[grep("^soundType", names(fx))]
  mean(c(0, st))  # reference level contributes 0
}

# Listener profile at percentile p of the random slope. The intercept is set to its
# conditional expectation given that slope, E[b0 | b1] = rho * (tau0 / tau1) * b1, so
# the profile is a coherent draw from the fitted random-effect distribution rather than
# an arbitrary mix of an extreme slope with an average intercept.
listener_pars <- function(b0, b1, tau0, tau1, rho, p) {
  b1_dev <- qnorm(p) * unname(tau1)
  b0_dev <- if (unname(tau1) > 0) unname(rho) * (unname(tau0) / unname(tau1)) * b1_dev else 0
  c(intercept = unname(b0) + b0_dev, slope = unname(b1) + b1_dev)
}

# Extract (intercept, slope, tau0, tau1, rho, soundType offset) from an lmer fit.
lmm_pars <- function(fit) {
  fx <- fixef(fit)
  vc <- VarCorr(fit)$participantId
  sd <- attr(vc, "stddev"); cr <- attr(vc, "correlation")
  c(b0 = unname(fx[1]) + sound_offset(fx), b1 = unname(fx[2]),
    tau0 = unname(sd[1]), tau1 = unname(sd[2]), rho = unname(cr[1, 2]))
}

# Predicted displacement (m) on a grid, for the three listener profiles.
predict_profiles <- function(pars, xgrid) {
  avg <- unname(pars["b0"]) + unname(pars["b1"]) * xgrid
  p10 <- listener_pars(pars["b0"], pars["b1"], pars["tau0"], pars["tau1"], pars["rho"], 0.10)
  p90 <- listener_pars(pars["b0"], pars["b1"], pars["tau0"], pars["tau1"], pars["rho"], 0.90)
  unname(c(avg,
           p10["intercept"] + p10["slope"] * xgrid,
           p90["intercept"] + p90["slope"] * xgrid))
}

pars_metric <- lmm_pars(m_metric)
pars_ang    <- lmm_pars(m_ang)

say("\nMetric model (per metre of disparity):")
say("  population slope        : ", round(pars_metric["b1"], 4),
    "  (", round(pars_metric["b1"] * 100, 2), " cm displacement per 100 cm of offset)")
say("  random-slope SD (tau1)  : ", round(pars_metric["tau1"], 4))
say("  random-intercept SD     : ", round(pars_metric["tau0"], 4),
    "   intercept-slope r: ", round(pars_metric["rho"], 3))
lp10m <- listener_pars(pars_metric["b0"], pars_metric["b1"], pars_metric["tau0"],
                       pars_metric["tau1"], pars_metric["rho"], 0.10)
lp90m <- listener_pars(pars_metric["b0"], pars_metric["b1"], pars_metric["tau0"],
                       pars_metric["tau1"], pars_metric["rho"], 0.90)
say("  10th pct listener slope : ", round(lp10m["slope"], 4), " (resistant)")
say("  90th pct listener slope : ", round(lp90m["slope"], 4), " (susceptible)")

# Empirical check against the 31 fitted random slopes (BLUPs).
blup_slopes <- coef(m_metric)$participantId[["stimulusDisparity_m"]]
say("\nfitted per-listener slopes (BLUPs), n = ", length(blup_slopes), ": ",
    "min ", round(min(blup_slopes), 3), ", 10th pct ", round(quantile(blup_slopes, 0.10), 3),
    ", median ", round(median(blup_slopes), 3),
    ", 90th pct ", round(quantile(blup_slopes, 0.90), 3), ", max ", round(max(blup_slopes), 3))
say("(BLUPs are shrunk toward the mean, so their spread is narrower than the fitted",
    " random-effect distribution used for the table.)")

# ---------------------------------------------------------------------------
# 4. Parametric bootstrap for confidence intervals on all three listener profiles
# ---------------------------------------------------------------------------
rule("4. PARAMETRIC BOOTSTRAP")

grid_cm  <- seq(15, 70, by = 5)          # nominal disparity, cm
grid_deg <- seq(20, 70, by = 5)          # egocentric angular disparity, deg

boot_profiles <- function(fit, xgrid, label) {
  fn <- function(f) predict_profiles(lmm_pars(f), xgrid)
  bb <- bootMer(fit, fn, nsim = N_BOOT, seed = 20260805, type = "parametric")
  n_bad <- sum(!complete.cases(bb$t))
  say(label, ": ", N_BOOT, " parametric bootstrap replicates, ",
      n_bad, " discarded (singular/failed refit), ", N_BOOT - n_bad, " retained")
  if (length(bb$msgs$warning)) {
    say(label, " bootstrap warnings (unique): ",
        paste(unique(unlist(bb$msgs$warning)), collapse = " | "))
  }
  t <- bb$t[complete.cases(bb$t), , drop = FALSE]
  list(lo = apply(t, 2, quantile, 0.025), hi = apply(t, 2, quantile, 0.975), n = nrow(t))
}

bo_metric <- boot_profiles(m_metric, grid_cm / 100, "metric model")
bo_ang    <- boot_profiles(m_ang, grid_deg / 10, "angular model")

# ---------------------------------------------------------------------------
# 5. (A) LOOKUP TABLES
# ---------------------------------------------------------------------------
rule("5. (A) PREDICTED-DISPLACEMENT LOOKUP TABLES")

# Empirical nominal -> angular crosswalk: trials within +/- 2.5 cm of each grid point.
cross <- bind_rows(lapply(grid_cm, function(g) {
  s <- d[abs(d$disp_cm - g) <= 2.5, ]
  tibble(disparity_cm = g, n_trials = nrow(s),
         ang_mean_deg   = mean(s$ang_deg),
         ang_median_deg = median(s$ang_deg),
         ang_q25_deg    = unname(quantile(s$ang_deg, 0.25)),
         ang_q75_deg    = unname(quantile(s$ang_deg, 0.75)),
         angOnset_median_deg = median(s$angOnset_deg))
}))

pm <- predict_profiles(pars_metric, grid_cm / 100)
k  <- length(grid_cm)
tabA <- tibble(
  disparity_cm = grid_cm,
  displacement_avg_cm    = pm[1:k] * 100,
  displacement_avg_lo_cm = bo_metric$lo[1:k] * 100,
  displacement_avg_hi_cm = bo_metric$hi[1:k] * 100,
  displacement_p10_cm    = pm[(k + 1):(2 * k)] * 100,
  displacement_p10_lo_cm = bo_metric$lo[(k + 1):(2 * k)] * 100,
  displacement_p10_hi_cm = bo_metric$hi[(k + 1):(2 * k)] * 100,
  displacement_p90_cm    = pm[(2 * k + 1):(3 * k)] * 100,
  displacement_p90_lo_cm = bo_metric$lo[(2 * k + 1):(3 * k)] * 100,
  displacement_p90_hi_cm = bo_metric$hi[(2 * k + 1):(3 * k)] * 100
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

say("\n--- Table A1: nominal disparity (metric model M1) ---")
print(as.data.frame(tabA %>% select(disparity_cm, displacement_avg_cm, displacement_avg_lo_cm,
                                    displacement_avg_hi_cm, displacement_p10_cm, displacement_p90_cm,
                                    traversed_avg_pct, traversed_p10_pct, traversed_p90_pct,
                                    ang_median_deg)))

pa <- predict_profiles(pars_ang, grid_deg / 10)
ka <- length(grid_deg)
tabB <- tibble(
  angular_disparity_deg  = grid_deg,
  displacement_avg_cm    = pa[1:ka] * 100,
  displacement_avg_lo_cm = bo_ang$lo[1:ka] * 100,
  displacement_avg_hi_cm = bo_ang$hi[1:ka] * 100,
  displacement_p10_cm    = pa[(ka + 1):(2 * ka)] * 100,
  displacement_p10_lo_cm = bo_ang$lo[(ka + 1):(2 * ka)] * 100,
  displacement_p10_hi_cm = bo_ang$hi[(ka + 1):(2 * ka)] * 100,
  displacement_p90_cm    = pa[(2 * ka + 1):(3 * ka)] * 100,
  displacement_p90_lo_cm = bo_ang$lo[(2 * ka + 1):(3 * ka)] * 100,
  displacement_p90_hi_cm = bo_ang$hi[(2 * ka + 1):(3 * ka)] * 100
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Number of observed trials supporting each angular grid row (+/- 2.5 deg).
tabB$n_trials_within_2p5deg <- sapply(grid_deg, function(g) sum(abs(d$ang_deg - g) <= 2.5))

say("\n--- Table A2: egocentric angular disparity (angular model M2) ---")
print(as.data.frame(tabB))

say("\n--- Empirical nominal-to-angular crosswalk ---")
print(as.data.frame(cross %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))

write.csv(tabA, file.path(res_dir, "tolerance_lookup_nominal.csv"), row.names = FALSE)
write.csv(tabB, file.path(res_dir, "tolerance_lookup_angular.csv"), row.names = FALSE)
write.csv(cross, file.path(res_dir, "tolerance_nominal_angular_crosswalk.csv"), row.names = FALSE)

# Sensitivity: percentile listeners with the slope shifted but the intercept left at the
# population value (i.e. ignoring the fitted intercept-slope correlation).
say("\nSensitivity (population intercept retained, slope shifted only), metric model:")
for (x in c(0.15, 0.40, 0.70)) {
  lo <- (pars_metric["b0"] + (pars_metric["b1"] - Z90 * pars_metric["tau1"]) * x) * 100
  hi <- (pars_metric["b0"] + (pars_metric["b1"] + Z90 * pars_metric["tau1"]) * x) * 100
  say("  ", x * 100, " cm: p10 ", round(lo, 2), " cm, p90 ", round(hi, 2), " cm")
}

# ---------------------------------------------------------------------------
# 6. (B) CAPTURE-PROBABILITY MODELS
# ---------------------------------------------------------------------------
rule("6. (B) CAPTURE-PROBABILITY MODELS")

# Note: the fits below name `d` explicitly rather than taking a data argument, so that
# every model records the same data object and anova() can compare them.
fit_capture <- function(yvar, xvar, label) {
  f <- reformulate(c(xvar, sprintf("(1 + %s | participantId)", xvar)), response = yvar)
  fit <- glmmTMB(f, data = d, family = binomial)
  say("\n--- ", label, " ---")
  say("formula: ", deparse(f))
  print(summary(fit)$coefficients$cond)
  print(VarCorr(fit))
  say("convergence code: ", fit$fit$convergence,
      " | message: ", fit$fit$message,
      " | positive-definite Hessian: ", fit$sdr$pdHess)
  fit
}

# Centred predictors: the random intercept then represents a listener's capture
# propensity at the centre of the tested range, which is the axis the heterogeneity
# actually lies on once the slope is near zero.
g_ang <- fit_capture("capture", "ang10c", "B1: capture ~ angular disparity (per 10 deg, centred)")
g_nom <- fit_capture("capture", "disp10c", "B2: capture ~ nominal disparity (per 10 cm, centred)")

# Reduced model: no disparity term at all.
g_null <- glmmTMB(capture ~ 1 + (1 | participantId), data = d, family = binomial)
say("\n--- B0: capture ~ 1 + (1 | participantId) ---")
print(summary(g_null)$coefficients$cond)
print(VarCorr(g_null))
say("convergence code: ", g_null$fit$convergence, " | pdHess: ", g_null$sdr$pdHess)

say("\nLikelihood-ratio tests against the disparity-free model B0:")
print(anova(g_null, g_ang))
print(anova(g_null, g_nom))

# The 3-df test above bundles the fixed slope with the two random-slope variance
# parameters. Decompose it, because the two carry very different messages: a population
# trend versus listener-specific trends around a flat average.
g_ang_ri <- glmmTMB(capture ~ ang10c + (1 | participantId), data = d, family = binomial)
g_nom_ri <- glmmTMB(capture ~ disp10c + (1 | participantId), data = d, family = binomial)
say("\nDecomposition: fixed slope only (1 df), then random slope given the fixed slope (2 df).")
say("Angular:")
print(anova(g_null, g_ang_ri))
print(anova(g_ang_ri, g_ang))
say("Nominal:")
print(anova(g_null, g_nom_ri))
print(anova(g_nom_ri, g_nom))
say("Variance-component tests sit on a parameter-space boundary, so these p-values are conservative.")

# Nonlinearity check on the angular predictor.
g_ang_q <- glmmTMB(capture ~ ang10c + I(ang10c^2) + (1 + ang10c | participantId),
                   data = d, family = binomial)
say("\nQuadratic angular term (checks for a threshold-like rise at large angles):")
print(summary(g_ang_q)$coefficients$cond)
print(anova(g_ang, g_ang_q))

# DHARMa residual checks on the primary capture model.
say("\nDHARMa checks, model B1:")
sim_res <- simulateResiduals(g_ang, n = 1000, seed = 20260805)
print(testUniformity(sim_res, plot = FALSE))
print(testDispersion(sim_res, plot = FALSE))

# ---- Crossing points -------------------------------------------------------
# x at which the linear predictor reaches logit(p). Where the slope is not
# distinguishable from zero the ratio is a Fieller problem and its confidence set is
# unbounded, so the interval is reported as such rather than as a spurious range.
fieller <- function(b, V, target, level = 0.95) {
  b <- unname(b); V <- unname(as.matrix(V))
  z <- qnorm(1 - (1 - level) / 2)
  a2 <- b[2]^2 - z^2 * V[2, 2]
  a1 <- 2 * ((b[1] - target) * b[2] - z^2 * V[1, 2])
  a0 <- (b[1] - target)^2 - z^2 * V[1, 1]
  point <- (target - b[1]) / b[2]
  if (a2 <= 0) return(c(point = point, lo = -Inf, hi = Inf, bounded = 0))
  disc <- a1^2 - 4 * a2 * a0
  if (disc <= 0) return(c(point = point, lo = NA, hi = NA, bounded = -1))  # empty set
  r <- sort((-a1 + c(-1, 1) * sqrt(disc)) / (2 * a2))
  c(point = point, lo = r[1], hi = r[2], bounded = 1)
}

crossings <- function(fit, xname, unit, scale, centre, label) {
  b <- fixef(fit)$cond
  V <- vcov(fit)$cond
  vc <- VarCorr(fit)$cond$participantId
  tau0 <- unname(attr(vc, "stddev")[1]); tau1 <- unname(attr(vc, "stddev")[2])
  rho  <- unname(attr(vc, "correlation")[1, 2])
  out <- list()
  for (p in c(0.50, 0.75)) {
    tgt <- qlogis(p)
    for (who in c("average", "p10_resistant", "p90_susceptible")) {
      if (who == "average") {
        bb <- b
      } else {
        pr <- if (who == "p10_resistant") 0.10 else 0.90
        b0_dev <- qnorm(pr) * tau0
        b1_dev <- if (tau0 > 0) rho * (tau1 / tau0) * b0_dev else 0
        bb <- c(b[1] + b0_dev, b[2] + b1_dev)
      }
      fi <- fieller(bb, V, tgt)
      # back-transform from the centred, scaled predictor to natural units
      bt <- function(v) if (is.finite(v)) (v + centre) * scale else v
      out[[length(out) + 1]] <- tibble(
        model = label, predictor = xname, unit = unit, listener = who,
        probability = p,
        crossing_point = bt(unname(fi["point"])),
        ci_lo = bt(unname(fi["lo"])), ci_hi = bt(unname(fi["hi"])),
        ci_bounded = unname(fi["bounded"])
      )
    }
  }
  bind_rows(out)
}

cross_ang <- crossings(g_ang, "angular disparity", "deg", 10, mean_ang10, "B1 angular")
cross_nom <- crossings(g_nom, "nominal disparity", "cm", 10, mean_disp10, "B2 nominal")
cross_all <- bind_rows(cross_ang, cross_nom) %>%
  mutate(observed_max = ifelse(unit == "deg", max(d$ang_deg), max(d$disp_cm)),
         within_tested_range = crossing_point >= 0 & crossing_point <= observed_max)

say("\n--- 50% and 75% crossing points ---")
say("ci_bounded: 1 = bounded Fieller interval, 0 = unbounded (slope not separable from zero)")
say("Caveat: the percentile-listener intervals treat the fitted random-effect SDs and")
say("correlation as known, propagating only fixed-effect uncertainty, so they are")
say("narrower than a fully propagated interval would be. A parametric bootstrap of the")
say("random-slope logistic models is not usable here: only about half the replicates")
say("return a positive-definite Hessian, and discarding the rest would select on the")
say("very variance components being estimated.")
print(as.data.frame(cross_all %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))

# ---- Listener-level capture probabilities from the disparity-free model ----
# Because no disparity term is supported, the operative quantity is a listener's
# constant capture probability. Parametric bootstrap of B0 (which converges on every
# replicate, unlike the random-slope models) gives the CIs.
say("\n--- Capture probability by listener percentile (model B0, disparity-invariant) ---")
b0_null   <- fixef(g_null)$cond[1]
tau0_null <- attr(VarCorr(g_null)$cond$participantId, "stddev")[1]

set.seed(20260805)
sims_null <- simulate(g_null, nsim = N_BOOT)
boot_null <- vapply(sims_null, function(y) {
  dd <- d; dd$capture <- y
  f <- try(suppressWarnings(update(g_null, data = dd)), silent = TRUE)
  if (inherits(f, "try-error") || !f$sdr$pdHess || f$fit$convergence != 0) return(c(NA, NA))
  c(fixef(f)$cond[1], attr(VarCorr(f)$cond$participantId, "stddev")[1])
}, numeric(2))
keep <- complete.cases(t(boot_null))
say("bootstrap replicates for B0: ", N_BOOT, " requested, ", sum(keep), " retained")

prob_at <- function(b0, tau, p) plogis(b0 + qnorm(p) * tau)
cap_tab <- bind_rows(lapply(c(0.10, 0.50, 0.90), function(p) {
  est <- prob_at(b0_null, tau0_null, p)
  bs  <- prob_at(boot_null[1, keep], boot_null[2, keep], p)
  tibble(listener_percentile = p * 100,
         label = c("10th (resistant)", "50th (median)", "90th (susceptible)")[match(p, c(0.10, 0.50, 0.90))],
         capture_probability = est,
         ci_lo = unname(quantile(bs, 0.025)), ci_hi = unname(quantile(bs, 0.975)))
})) %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
print(as.data.frame(cap_tab))

# The listener percentile at which the constant capture probability itself crosses
# 50% and 75%.
pct_at <- function(target) pnorm((qlogis(target) - b0_null) / tau0_null)
say("\nlistener percentile whose capture probability reaches 50%: ",
    round(100 * pct_at(0.50), 1), "%")
say("listener percentile whose capture probability reaches 75%: ",
    round(100 * pct_at(0.75), 1), "%")
say("(i.e. the 50%/75% thresholds are crossed by moving along the listener axis, not",
    " along the disparity axis.)")

# ---------------------------------------------------------------------------
# 7. Empirical proportions by disparity bin (sanity check)
# ---------------------------------------------------------------------------
rule("7. EMPIRICAL CAPTURE PROPORTIONS BY DISPARITY BIN")

binom_ci <- function(x, n) {
  if (n == 0) return(c(NA, NA))
  as.numeric(binom.test(x, n)$conf.int)
}

emp_nom <- d %>%
  mutate(bin = cut(disp_cm, breaks = c(15, 25, 35, 45, 55, 70), include.lowest = TRUE, right = FALSE)) %>%
  group_by(bin) %>%
  summarise(n = n(), k = sum(capture),
            p_capture = mean(capture),
            p_soft_gt0.25 = mean(capture_soft),
            p_strict_gt0.75 = mean(capture_strict),
            median_ang_deg = median(ang_deg), .groups = "drop") %>%
  rowwise() %>%
  mutate(ci_lo = binom_ci(k, n)[1], ci_hi = binom_ci(k, n)[2]) %>%
  ungroup() %>%
  mutate(metric = "nominal disparity (cm)") %>%
  mutate(bin_label = as.character(bin)) %>%
  select(-bin)

ang_breaks <- c(0, 20, 30, 40, 50, 60, Inf)
emp_ang <- d %>%
  mutate(bin = cut(ang_deg, breaks = ang_breaks, right = FALSE)) %>%
  group_by(bin) %>%
  summarise(n = n(), k = sum(capture),
            p_capture = mean(capture),
            p_soft_gt0.25 = mean(capture_soft),
            p_strict_gt0.75 = mean(capture_strict),
            median_ang_deg = median(ang_deg), .groups = "drop") %>%
  rowwise() %>%
  mutate(ci_lo = binom_ci(k, n)[1], ci_hi = binom_ci(k, n)[2]) %>%
  ungroup() %>%
  mutate(metric = "experienced angular disparity (deg)") %>%
  mutate(bin_label = as.character(bin)) %>%
  select(-bin)

emp <- bind_rows(emp_nom, emp_ang) %>%
  select(metric, bin_label, n, k, p_capture, ci_lo, ci_hi,
         p_soft_gt0.25, p_strict_gt0.75, median_ang_deg) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))
print(as.data.frame(emp))

# ---------------------------------------------------------------------------
# 8. Robustness: softer and stricter criteria
# ---------------------------------------------------------------------------
rule("8. ROBUSTNESS ACROSS CRITERIA")

rob <- bind_rows(lapply(
  list(list("bias > 0.50 (primary; = closer to visual cue)", "capture"),
       list("bias > 0.25 (softer)", "capture_soft"),
       list("bias > 0.75 (stricter)", "capture_strict")),
  function(cr) {
    bind_rows(lapply(list(c("ang10c", "angular, per 10 deg"),
                          c("disp10c", "nominal, per 10 cm")), function(xv) {
      f <- reformulate(c(xv[1], sprintf("(1 + %s | participantId)", xv[1])), response = cr[[2]])
      fit <- suppressWarnings(glmmTMB(f, data = d, family = binomial))
      co <- summary(fit)$coefficients$cond
      tibble(criterion = cr[[1]], base_rate = mean(d[[cr[[2]]]]), predictor = xv[2],
             intercept = co[1, 1], slope = co[2, 1], slope_se = co[2, 2],
             slope_z = co[2, 3], slope_p = co[2, 4],
             slope_ci_lo = co[2, 1] - 1.96 * co[2, 2],
             slope_ci_hi = co[2, 1] + 1.96 * co[2, 2],
             converged = fit$fit$convergence == 0 && isTRUE(fit$sdr$pdHess))
    }))
  })) %>% mutate(across(where(is.numeric), ~ round(.x, 5)))
print(as.data.frame(rob))

# ---------------------------------------------------------------------------
# 9. Write machine-readable outputs
# ---------------------------------------------------------------------------
rule("9. OUTPUT FILES")

coef_tab <- bind_rows(
  tibble(model = "B1 angular (per 10 deg, centred)",
         term = rownames(summary(g_ang)$coefficients$cond),
         estimate = summary(g_ang)$coefficients$cond[, 1],
         se = summary(g_ang)$coefficients$cond[, 2],
         z = summary(g_ang)$coefficients$cond[, 3],
         p = summary(g_ang)$coefficients$cond[, 4]),
  tibble(model = "B2 nominal (per 10 cm, centred)",
         term = rownames(summary(g_nom)$coefficients$cond),
         estimate = summary(g_nom)$coefficients$cond[, 1],
         se = summary(g_nom)$coefficients$cond[, 2],
         z = summary(g_nom)$coefficients$cond[, 3],
         p = summary(g_nom)$coefficients$cond[, 4]),
  tibble(model = "B0 disparity-free",
         term = rownames(summary(g_null)$coefficients$cond),
         estimate = summary(g_null)$coefficients$cond[, 1],
         se = summary(g_null)$coefficients$cond[, 2],
         z = summary(g_null)$coefficients$cond[, 3],
         p = summary(g_null)$coefficients$cond[, 4])
) %>% mutate(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)

lmm_tab <- bind_rows(
  tibble(model = "M1 metric (per m)", term = names(fixef(m_metric)),
         estimate = unname(fixef(m_metric)),
         se = unname(sqrt(diag(as.matrix(vcov(m_metric)))))),
  tibble(model = "M2 angular (per 10 deg)", term = names(fixef(m_ang)),
         estimate = unname(fixef(m_ang)),
         se = unname(sqrt(diag(as.matrix(vcov(m_ang))))))
) %>% mutate(ci_lo = estimate - 1.96 * se, ci_hi = estimate + 1.96 * se)

rand_tab <- bind_rows(
  tibble(model = "M1 metric (per m)", intercept_sd = pars_metric["tau0"],
         slope_sd = pars_metric["tau1"], corr = pars_metric["rho"],
         residual_sd = sigma(m_metric)),
  tibble(model = "M2 angular (per 10 deg)", intercept_sd = pars_ang["tau0"],
         slope_sd = pars_ang["tau1"], corr = pars_ang["rho"],
         residual_sd = sigma(m_ang))
)

files <- c("tolerance_lookup_nominal.csv", "tolerance_lookup_angular.csv",
           "tolerance_nominal_angular_crosswalk.csv")
write.csv(coef_tab,  file.path(res_dir, "capture_model_coefficients.csv"), row.names = FALSE)
write.csv(cross_all, file.path(res_dir, "capture_crossing_points.csv"), row.names = FALSE)
write.csv(cap_tab,   file.path(res_dir, "capture_probability_by_listener.csv"), row.names = FALSE)
write.csv(emp,       file.path(res_dir, "capture_empirical_by_bin.csv"), row.names = FALSE)
write.csv(rob,       file.path(res_dir, "capture_criterion_robustness.csv"), row.names = FALSE)
write.csv(lmm_tab,   file.path(res_dir, "tolerance_lmm_fixed_effects.csv"), row.names = FALSE)
write.csv(rand_tab,  file.path(res_dir, "tolerance_lmm_random_effects.csv"), row.names = FALSE)
files <- c(files, "capture_model_coefficients.csv", "capture_crossing_points.csv",
           "capture_probability_by_listener.csv", "capture_empirical_by_bin.csv",
           "capture_criterion_robustness.csv", "tolerance_lmm_fixed_effects.csv",
           "tolerance_lmm_random_effects.csv")
for (f in files) say("written: ", file.path(res_dir, f))
say("written: ", log_path)

rule("SESSION INFO")
print(sessionInfo())
