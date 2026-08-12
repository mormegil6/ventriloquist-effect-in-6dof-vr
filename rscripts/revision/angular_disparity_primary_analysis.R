#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Primary analysis re-run with EXPERIENCED EGOCENTRIC ANGULAR DISPARITY
#
# Reviewer 1 (point d) observes that the manipulated variable is a world-space
# separation in centimetres, whereas the conflict a freely walking listener
# actually experiences is an angle at the head that changes throughout the
# trial. This script answers that empirically.
#
#   Step 1  Derives an ANGULAR response measure: the in-plane angle, at a given
#           head pose, from the direction of the true source towards the
#           direction of the visual cue, subtended by the response placement.
#           This is the spherical analogue of the published signed error, which
#           is a projection onto the source-to-cue direction.
#   Step 2  Within-trial variation of angular disparity, and the correlation of
#           the nominal predictor with every angular summary.
#   Step 3  Fits the primary model with (a) nominal metric disparity,
#           (b) onset angular disparity, (c) time-weighted mean angular
#           disparity, on both the metric and the angular response.
#   Step 4  Compares nominal against angular by AIC/BIC, marginal and
#           conditional R2, and by likelihood-ratio tests inside models that
#           contain both predictors.
#   Step 5  Reports the angular capture slope with confidence intervals.
#   Step 6  Robustness checks.
#
# PRE-SPECIFICATION, fixed before any model was fitted:
#   Primary predictor summary : angular disparity at trial ONSET.
#   Secondary                 : time-weighted mean angular disparity.
#   Primary response          : angular response displacement (deg).
#   Random structure          : (1 + predictor | participantId); if that fit is
#                               singular or fails to converge, fall back to the
#                               zero-correlation form (1 + predictor ||
#                               participantId), and then to (1 | participantId).
#                               The same structure is used for every model that
#                               shares a response, so that information criteria
#                               compare fixed effects only.
#   Continuous predictors are mean-centred. Centring leaves the slope, the
#   log-likelihood and the correlated random-effects fit unchanged; it makes the
#   intercept the value at mean disparity and makes the zero-correlation
#   fallback a meaningful restriction. The published uncentred model is also
#   refitted verbatim, so its intercept at zero disparity is reproduced.
#
# Inputs : results_revision/analysis_df_revision.rds
#          results_revision/trajectory_samples_angular.rds
# Outputs: results_revision/angular_primary_*.csv and angular_primary_log.txt
#
# R 4.4 arm64 framework build. Deterministic; set.seed is used for DHARMa only.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
  library(DHARMa)
})

set.seed(20260805)

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir  <- file.path(proj_dir, "results_revision")

log_con <- file(file.path(res_dir, "angular_primary_log.txt"), open = "wt")
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_con)
}
say_df <- function(df, digits = 4) {
  for (l in capture.output(print(as.data.frame(df), digits = digits, row.names = FALSE))) say(l)
}
rule <- function(title) {
  say("")
  say(strrep("=", 78))
  say(title)
  say(strrep("=", 78))
}

rule("ANGULAR DISPARITY RE-ANALYSIS OF THE PRIMARY MODEL")
say("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", R.version.string)
say("")
say("PRE-SPECIFICATION (stated before fitting)")
say("  Primary predictor summary : angular disparity at trial ONSET (angDisp_onset_deg)")
say("  Secondary                 : time-weighted mean angular disparity (angDisp_mean_deg)")
say("  Primary response          : angular response displacement at the onset head pose")
say("")
say("  Why onset is primary:")
say("   (i)   Exogeneity. The onset geometry follows from the stimulus placement and")
say("         from whatever pose the listener happens to hold when the trial starts.")
say("         It is not shaped by where the listener subsequently chooses to walk.")
say("         The time-weighted mean is a function of the trajectory, and the")
say("         trajectory may itself be driven by the percept under study, so a slope")
say("         on it mixes a stimulus effect with self-selected sampling.")
say("   (ii)  Frame consistency. An angular response measure must be evaluated from")
say("         some head position. Onset is the only pose available on every trial")
say("         that is independent of the response, so predictor and response are")
say("         defined in one and the same geometric frame.")
say("   (iii) Comparability. Onset angular disparity is the quantity the")
say("         fixed-listener ventriloquism literature manipulates, so it places the")
say("         present slope on the same axis as prior work.")
say("  The time-weighted mean is reported alongside throughout, since it is the")
say("  summary that embodies the 'changes over time' half of the reviewer's point.")

# ---------------------------------------------------------------------------
# STEP 0: data
# ---------------------------------------------------------------------------
rule("STEP 0: DATA")

analysis_df <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
samples     <- readRDS(file.path(res_dir, "trajectory_samples_angular.rds"))

say("Trials: ", nrow(analysis_df), " | participants: ", n_distinct(analysis_df$participantId))
say("Trajectory samples: ", nrow(samples))

# ---------------------------------------------------------------------------
# STEP 1: angular response displacement
# ---------------------------------------------------------------------------
rule("STEP 1: ANGULAR RESPONSE DISPLACEMENT")

say("Let H be a head position and u_s, u_c, u_r the unit vectors from H towards")
say("the true source, the visual cue and the response placement. Let t be the unit")
say("vector in the (u_s, u_c) plane, orthogonal to u_s and pointing towards u_c:")
say("")
say("   angular disparity      theta_D = atan2(|u_s x u_c|, u_s . u_c)")
say("   angular response shift theta_R = atan2(u_r . t,     u_r . u_s)")
say("")
say("theta_R is 0 deg when the response lies in the source direction and equals")
say("theta_D when it lies in the cue direction, so theta_R / theta_D is the angular")
say("analogue of the bias proportion and d(theta_R)/d(theta_D) the angular analogue")
say("of the published displacement slope. The out-of-plane component of the response")
say("direction is discarded, exactly as the published signed error discards the")
say("orthogonal error component.")
say("")
say("Three evaluation frames are computed: the first logged head pose (onset,")
say("primary), the last logged head pose (sensitivity), and the time-weighted")
say("average over all samples of the trial (matched partner of the time-weighted")
say("mean angular disparity).")

deg <- function(rad) rad * 180 / pi

# Signed in-plane angle of vector r away from vector s, positive towards c.
angular_projection <- function(sx, sy, sz, cx, cy, cz, rx, ry, rz) {
  n_s <- sqrt(sx^2 + sy^2 + sz^2)
  us_x <- sx / n_s; us_y <- sy / n_s; us_z <- sz / n_s
  dot_cs <- cx * us_x + cy * us_y + cz * us_z
  tx <- cx - dot_cs * us_x
  ty <- cy - dot_cs * us_y
  tz <- cz - dot_cs * us_z
  n_t <- sqrt(tx^2 + ty^2 + tz^2)
  tx <- tx / n_t; ty <- ty / n_t; tz <- tz / n_t
  deg(atan2(rx * tx + ry * ty + rz * tz, rx * us_x + ry * us_y + rz * us_z))
}

# Per-sample angular response, then per-trial summaries in the three frames.
trial_geom <- analysis_df %>%
  select(participantId, trialSequenceNum,
         sound_x, sound_y, sound_z, flash_x, flash_y, flash_z,
         response_x, response_y, response_z)

samp <- samples %>%
  select(participantId, trialSequenceNum, sample_idx, w_time, px, py, pz, angDisp_deg) %>%
  left_join(trial_geom, by = c("participantId", "trialSequenceNum")) %>%
  mutate(angResp_deg = angular_projection(
    sound_x - px, sound_y - py, sound_z - pz,
    flash_x - px, flash_y - py, flash_z - pz,
    response_x - px, response_y - py, response_z - pz))

wmean <- function(x, w) sum(x * w) / sum(w)

resp_summary <- samp %>%
  group_by(participantId, trialSequenceNum) %>%
  summarise(
    angResp_onset_deg = angResp_deg[which.min(sample_idx)],
    angResp_final_deg = angResp_deg[which.max(sample_idx)],
    angResp_mean_deg  = wmean(angResp_deg, w_time),
    angDisp_final_deg = angDisp_deg[which.max(sample_idx)],
    angDisp_onset_recheck = angDisp_deg[which.min(sample_idx)],
    head_onset_x = px[which.min(sample_idx)],
    head_onset_y = py[which.min(sample_idx)],
    head_onset_z = pz[which.min(sample_idx)],
    .groups = "drop"
  )

dat <- analysis_df %>%
  left_join(resp_summary, by = c("participantId", "trialSequenceNum")) %>%
  mutate(
    angDisp_range_deg = angDisp_max_deg - angDisp_min_deg,
    angBiasProportion_onset = angResp_onset_deg / angDisp_onset_deg,
    angBiasProportion_mean  = angResp_mean_deg / angDisp_mean_deg,
    # mean-centred predictors
    stimDisp_cm      = stimulusDisparity_m - mean(stimulusDisparity_m),
    angDisp_onset_cm = angDisp_onset_deg - mean(angDisp_onset_deg),
    angDisp_mean_cm  = angDisp_mean_deg - mean(angDisp_mean_deg),
    angDisp_final_cm = angDisp_final_deg - mean(angDisp_final_deg)
  )

stopifnot(!anyNA(dat$angResp_onset_deg), !anyNA(dat$angResp_mean_deg))

say("")
say("Reconstruction check, onset angular disparity recomputed here vs the prepared")
say("column: max absolute difference = ",
    signif(max(abs(dat$angDisp_onset_recheck - dat$angDisp_onset_deg)), 3), " deg")
say("Trials with onset angular disparity below 0.5 deg (projection axis ill-conditioned): ",
    sum(dat$angDisp_onset_deg < 0.5))

desc <- data.frame(
  measure = c("angDisp_onset_deg", "angDisp_mean_deg", "angDisp_final_deg",
              "angResp_onset_deg", "angResp_mean_deg", "angResp_final_deg",
              "stimulusDisparity_cm", "signedError_cm"),
  mean = c(mean(dat$angDisp_onset_deg), mean(dat$angDisp_mean_deg), mean(dat$angDisp_final_deg),
           mean(dat$angResp_onset_deg), mean(dat$angResp_mean_deg), mean(dat$angResp_final_deg),
           mean(dat$stimulusDisparity_m) * 100, mean(dat$signedError_m) * 100),
  sd = c(sd(dat$angDisp_onset_deg), sd(dat$angDisp_mean_deg), sd(dat$angDisp_final_deg),
         sd(dat$angResp_onset_deg), sd(dat$angResp_mean_deg), sd(dat$angResp_final_deg),
         sd(dat$stimulusDisparity_m) * 100, sd(dat$signedError_m) * 100),
  median = c(median(dat$angDisp_onset_deg), median(dat$angDisp_mean_deg), median(dat$angDisp_final_deg),
             median(dat$angResp_onset_deg), median(dat$angResp_mean_deg), median(dat$angResp_final_deg),
             median(dat$stimulusDisparity_m) * 100, median(dat$signedError_m) * 100),
  min = c(min(dat$angDisp_onset_deg), min(dat$angDisp_mean_deg), min(dat$angDisp_final_deg),
          min(dat$angResp_onset_deg), min(dat$angResp_mean_deg), min(dat$angResp_final_deg),
          min(dat$stimulusDisparity_m) * 100, min(dat$signedError_m) * 100),
  max = c(max(dat$angDisp_onset_deg), max(dat$angDisp_mean_deg), max(dat$angDisp_final_deg),
          max(dat$angResp_onset_deg), max(dat$angResp_mean_deg), max(dat$angResp_final_deg),
          max(dat$stimulusDisparity_m) * 100, max(dat$signedError_m) * 100)
)
say("")
say("Descriptives (angles in degrees, distances in centimetres):")
say_df(desc)
write.csv(desc, file.path(res_dir, "angular_primary_descriptives.csv"), row.names = FALSE)

say("")
say(sprintf("Angular bias proportion at onset: mean %.3f, median %.3f, SD %.3f",
            mean(dat$angBiasProportion_onset), median(dat$angBiasProportion_onset),
            sd(dat$angBiasProportion_onset)))
say(sprintf("Angular bias proportion, time-weighted: mean %.3f, median %.3f, SD %.3f",
            mean(dat$angBiasProportion_mean), median(dat$angBiasProportion_mean),
            sd(dat$angBiasProportion_mean)))
say("Note: these per-trial ratios are heavy tailed because the denominator can be")
say("small; the model slope, not the ratio, is the quantity to report.")

# ---------------------------------------------------------------------------
# STEP 2: within-trial variation and predictor correlations
# ---------------------------------------------------------------------------
rule("STEP 2: WITHIN-TRIAL VARIATION OF ANGULAR DISPARITY")

qs <- c(0, .05, .25, .5, .75, .95, 1)
rng_q <- quantile(dat$angDisp_range_deg, qs)
say("Within-trial range angDisp_max - angDisp_min, degrees, n = ", nrow(dat), " trials:")
say(sprintf("  mean %.1f  SD %.1f | min %.1f  Q05 %.1f  Q25 %.1f  median %.1f  Q75 %.1f  Q95 %.1f  max %.1f",
            mean(dat$angDisp_range_deg), sd(dat$angDisp_range_deg),
            rng_q[1], rng_q[2], rng_q[3], rng_q[4], rng_q[5], rng_q[6], rng_q[7]))
for (thr in c(30, 60, 90, 120)) {
  say(sprintf("  trials with range > %3d deg: %3d (%.1f%%)",
              thr, sum(dat$angDisp_range_deg > thr), 100 * mean(dat$angDisp_range_deg > thr)))
}
say(sprintf("  ratio of within-trial range to the onset value: median %.2f, IQR %.2f-%.2f",
            median(dat$angDisp_range_deg / dat$angDisp_onset_deg),
            quantile(dat$angDisp_range_deg / dat$angDisp_onset_deg, .25),
            quantile(dat$angDisp_range_deg / dat$angDisp_onset_deg, .75)))
say(sprintf("  within-trial range as a multiple of the between-trial SD of the nominal"))
say(sprintf("  predictor expressed at onset (SD = %.1f deg): median %.2f",
            sd(dat$angDisp_onset_deg),
            median(dat$angDisp_range_deg) / sd(dat$angDisp_onset_deg)))

within_var <- data.frame(
  statistic = c("mean", "sd", names(rng_q), "pct_gt_30deg", "pct_gt_60deg",
                "pct_gt_90deg", "pct_gt_120deg"),
  angDisp_range_deg = c(mean(dat$angDisp_range_deg), sd(dat$angDisp_range_deg),
                        as.numeric(rng_q),
                        100 * mean(dat$angDisp_range_deg > 30),
                        100 * mean(dat$angDisp_range_deg > 60),
                        100 * mean(dat$angDisp_range_deg > 90),
                        100 * mean(dat$angDisp_range_deg > 120))
)
write.csv(within_var, file.path(res_dir, "angular_primary_within_trial_variation.csv"),
          row.names = FALSE)

say("")
say("Correlation of nominal metric disparity with each angular summary.")
say("r_within_participant centres both variables on their participant means.")

ang_cols <- c("angDisp_onset_deg", "angDisp_mean_deg", "angDisp_median_deg",
              "angDisp_min_deg", "angDisp_max_deg", "angDisp_range_deg",
              "angDisp_final_deg")
within_center <- function(x, g) x - ave(x, g, FUN = mean)

cor_tab <- bind_rows(lapply(ang_cols, function(cc) {
  x <- dat$stimulusDisparity_m; y <- dat[[cc]]
  ct <- cor.test(x, y)
  data.frame(angular_summary = cc,
             pearson_r = unname(ct$estimate),
             ci_lo = ct$conf.int[1], ci_hi = ct$conf.int[2],
             p_value = ct$p.value,
             r_squared = unname(ct$estimate)^2,
             spearman_rho = suppressWarnings(cor(x, y, method = "spearman")),
             r_within_participant = cor(within_center(x, dat$participantId),
                                        within_center(y, dat$participantId)))
}))
say_df(cor_tab)
write.csv(cor_tab, file.path(res_dir, "angular_primary_predictor_correlations.csv"),
          row.names = FALSE)

say("")
say("Correlation between the two angular summaries: r = ",
    signif(cor(dat$angDisp_onset_deg, dat$angDisp_mean_deg), 3))
say("Head-to-source distance at onset vs onset angular disparity: r = ",
    signif(cor(dat$dist_at_onset_m, dat$angDisp_onset_deg), 3))
say("Mean head-to-source distance vs mean angular disparity: r = ",
    signif(cor(dat$dist_mean_m, dat$angDisp_mean_deg), 3))

# ---------------------------------------------------------------------------
# STEP 3: model fitting
# ---------------------------------------------------------------------------
rule("STEP 3: MODEL FITS")

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

# Fit with the data argument recorded as the symbol `dat` in every call, so that
# anova() accepts the models as being fitted to the same data object.
fit_lmer <- function(fml, reml) {
  suppressWarnings(do.call("lmer", list(formula = fml, data = quote(dat),
                                        REML = reml, control = quote(ctrl))))
}
fit_ok <- function(fit) {
  is.null(fit@optinfo$conv$lme4$messages) && !isSingular(fit, tol = 1e-4)
}
fit_msgs <- function(fit) paste(fit@optinfo$conv$lme4$messages, collapse = " | ")

# Pre-specified ladder of random-effects structures.
re_forms <- function(resp, pred) {
  c(sprintf("%s ~ %s + soundType + (1 + %s | participantId)", resp, pred, pred),
    sprintf("%s ~ %s + soundType + (1 + %s || participantId)", resp, pred, pred),
    sprintf("%s ~ %s + soundType + (1 | participantId)", resp, pred))
}

# Descend the ladder until both the REML and the ML fit are clean.
select_structure <- function(resp, pred, label) {
  for (fs in re_forms(resp, pred)) {
    f <- as.formula(fs)
    m_reml <- fit_lmer(f, TRUE)
    m_ml   <- fit_lmer(f, FALSE)
    ok <- fit_ok(m_reml) && fit_ok(m_ml)
    say(sprintf("  [%s] %-58s REML ok=%-5s ML ok=%-5s %s", label,
                sub("^.*~ ", "", fs), fit_ok(m_reml), fit_ok(m_ml),
                paste(unique(c(fit_msgs(m_reml), fit_msgs(m_ml))), collapse = " ")))
    if (ok) return(list(reml = m_reml, ml = m_ml, formula = fs, rank = which(re_forms(resp, pred) == fs)))
  }
  stop("no clean random-effects structure for ", label)
}

say("Random-effects structure selection (pre-specified ladder; a structure is")
say("accepted when the REML and the ML fit are both non-singular and warning free):")
say("")

preds <- c(nominal = "stimDisp_cm", onset = "angDisp_onset_cm", timemean = "angDisp_mean_cm")

say("Response = signedError_m (metres, positive towards the cue):")
selA <- lapply(names(preds), function(k) select_structure("signedError_m", preds[[k]], paste0("A_", k)))
names(selA) <- names(preds)

say("")
say("Response = angResp_onset_deg (degrees, positive towards the cue):")
selB <- lapply(names(preds), function(k) select_structure("angResp_onset_deg", preds[[k]], paste0("B_", k)))
names(selB) <- names(preds)

rankA <- max(sapply(selA, `[[`, "rank"))
rankB <- max(sapply(selB, `[[`, "rank"))
say("")
say("Common structure adopted within each response set (the most restrictive one")
say("that any model in the set required, so that information criteria compare")
say("fixed effects only):")
say("  set A (metric response) : ladder level ", rankA)
say("  set B (angular response): ladder level ", rankB)

refit_at <- function(resp, pred, rank) {
  f <- as.formula(re_forms(resp, pred)[rank])
  list(reml = fit_lmer(f, TRUE), ml = fit_lmer(f, FALSE), formula = re_forms(resp, pred)[rank])
}

A0 <- refit_at("signedError_m", "stimDisp_cm", rankA)
A1 <- refit_at("signedError_m", "angDisp_onset_cm", rankA)
A2 <- refit_at("signedError_m", "angDisp_mean_cm", rankA)
B0 <- refit_at("angResp_onset_deg", "stimDisp_cm", rankB)
B1 <- refit_at("angResp_onset_deg", "angDisp_onset_cm", rankB)
B2 <- refit_at("angResp_onset_deg", "angDisp_mean_cm", rankB)

# Matched time-weighted formulation: response and predictor both time-weighted.
B2m <- refit_at("angResp_mean_deg", "angDisp_mean_cm", rankB)
# Published model, verbatim and uncentred, to reproduce the reported intercept.
A0_pub <- list(reml = fit_lmer(signedError_m ~ stimulusDisparity_m + soundType +
                                 (1 + stimulusDisparity_m | participantId), TRUE),
               ml = fit_lmer(signedError_m ~ stimulusDisparity_m + soundType +
                               (1 + stimulusDisparity_m | participantId), FALSE),
               formula = "published, uncentred")

models <- list(
  A0_metric_resp_nominal      = A0,
  A1_metric_resp_ang_onset    = A1,
  A2_metric_resp_ang_timemean = A2,
  B0_ang_resp_nominal         = B0,
  B1_ang_resp_ang_onset       = B1,
  B2_ang_resp_ang_timemean    = B2,
  B2m_angmean_resp_ang_timemean = B2m,
  A0pub_published_uncentred   = A0_pub
)

say("")
say("Final convergence report:")
for (nm in names(models)) {
  say(sprintf("  %-30s REML ok=%-5s ML ok=%-5s %s", nm,
              fit_ok(models[[nm]]$reml), fit_ok(models[[nm]]$ml),
              paste(unique(c(fit_msgs(models[[nm]]$reml), fit_msgs(models[[nm]]$ml))),
                    collapse = " ")))
}

# --- coefficients -----------------------------------------------------------
coef_rows <- bind_rows(lapply(names(models), function(nm) {
  f <- models[[nm]]$reml
  sm <- summary(f)$coefficients
  ci <- confint(f, method = "Wald", parm = rownames(sm))
  data.frame(model = nm, formula = models[[nm]]$formula,
             response = as.character(formula(f))[2], term = rownames(sm),
             estimate = sm[, "Estimate"], se = sm[, "Std. Error"],
             df = sm[, "df"], t = sm[, "t value"], p = sm[, "Pr(>|t|)"],
             ci_lo = ci[rownames(sm), 1], ci_hi = ci[rownames(sm), 2],
             row.names = NULL)
}))
say("")
say("Fixed effects (REML, Satterthwaite df, Wald 95% CI):")
say_df(coef_rows %>% select(-formula) %>% mutate(across(where(is.numeric), ~ signif(.x, 4))))
write.csv(coef_rows, file.path(res_dir, "angular_primary_coefficients.csv"), row.names = FALSE)

say("")
say("Type III F tests (Satterthwaite):")
for (nm in names(models)) {
  a <- anova(models[[nm]]$reml)
  say("  ", nm, ":")
  say_df(cbind(term = rownames(a), as.data.frame(a)))
}

# ---------------------------------------------------------------------------
# STEP 4: nominal versus angular
# ---------------------------------------------------------------------------
rule("STEP 4: NOMINAL METRIC VERSUS EXPERIENCED ANGULAR PREDICTOR")

say("Information criteria come from the ML fits and are comparable only WITHIN a")
say("response scale. Models A0/A1/A2 share signedError_m; B0/B1/B2 share")
say("angResp_onset_deg. Within a set the three models are non-nested, so AIC, BIC")
say("and R2 are the comparison; likelihood-ratio tests inside models holding both")
say("predictors are given afterwards.")

# Nakagawa and Johnson marginal and conditional R2, computed by hand (no MuMIn).
r2_nakagawa <- function(fit) {
  X <- model.matrix(fit)
  var_f <- as.numeric(var(as.vector(X %*% fixef(fit))))
  vc <- VarCorr(fit); mm <- getME(fit, "mmList")
  var_r <- 0
  for (i in seq_along(mm)) {
    M <- as.matrix(mm[[i]]); S <- as.matrix(vc[[i]])
    var_r <- var_r + mean(rowSums((M %*% S) * M))
  }
  var_e <- sigma(fit)^2
  tot <- var_f + var_r + var_e
  tau00 <- as.numeric(vc[[1]][1, 1])
  c(R2m = var_f / tot, R2c = (var_f + var_r) / tot,
    var_fixed = var_f, var_random = var_r, var_resid = var_e,
    ICC_intercept = tau00 / (tau00 + var_e))
}

comp <- bind_rows(lapply(names(models), function(nm) {
  fml <- models[[nm]]$ml; fr <- models[[nm]]$reml
  r2 <- r2_nakagawa(fr)
  data.frame(model = nm, response = as.character(formula(fr))[2],
             predictor = attr(terms(fr), "term.labels")[1],
             n = nobs(fml), df = attr(logLik(fml), "df"),
             logLik = as.numeric(logLik(fml)), AIC = AIC(fml), BIC = BIC(fml),
             R2_marginal = unname(r2["R2m"]), R2_conditional = unname(r2["R2c"]),
             ICC_intercept = unname(r2["ICC_intercept"]), resid_sd = sigma(fr))
}))
comp <- comp %>%
  group_by(response) %>%
  mutate(dAIC_within_response = AIC - min(AIC),
         dBIC_within_response = BIC - min(BIC)) %>%
  ungroup()

say("")
say_df(comp %>% mutate(across(where(is.numeric), ~ signif(.x, 6))))
write.csv(comp, file.path(res_dir, "angular_primary_model_comparison.csv"), row.names = FALSE)

say("")
say("Likelihood-ratio tests inside models containing BOTH predictors (ML fits,")
say("random structure held fixed within each comparison).")

lrt_rows <- list()
add_lrt <- function(label, resp, base_pred, add_pred, rank) {
  f_red  <- as.formula(re_forms(resp, base_pred)[rank])
  f_full <- as.formula(sub(paste0("~ ", base_pred, " \\+"),
                           paste0("~ ", base_pred, " + ", add_pred, " +"),
                           re_forms(resp, base_pred)[rank]))
  a <- anova(fit_lmer(f_red, FALSE), fit_lmer(f_full, FALSE))
  say(""); say("  ", label)
  say("    reduced: ", deparse(f_red))
  say("    full   : ", deparse(f_full))
  say_df(cbind(model = c("reduced", "full"), as.data.frame(a)))
  lrt_rows[[length(lrt_rows) + 1]] <<- data.frame(
    comparison = label, response = resp, base_predictor = base_pred,
    added_predictor = add_pred, chisq = a$Chisq[2], df = a$Df[2],
    p = a$`Pr(>Chisq)`[2], AIC_reduced = a$AIC[1], AIC_full = a$AIC[2])
  invisible(NULL)
}

add_lrt("Set A: onset angular disparity added to nominal metric disparity",
        "signedError_m", "stimDisp_cm", "angDisp_onset_cm", rankA)
add_lrt("Set A: time-weighted mean angular disparity added to nominal metric disparity",
        "signedError_m", "stimDisp_cm", "angDisp_mean_cm", rankA)
add_lrt("Set A: nominal metric disparity added to onset angular disparity",
        "signedError_m", "angDisp_onset_cm", "stimDisp_cm", rankA)
add_lrt("Set B: onset angular disparity added to nominal metric disparity",
        "angResp_onset_deg", "stimDisp_cm", "angDisp_onset_cm", rankB)
add_lrt("Set B: nominal metric disparity added to onset angular disparity",
        "angResp_onset_deg", "angDisp_onset_cm", "stimDisp_cm", rankB)
add_lrt("Set B: time-weighted mean angular disparity added to onset angular disparity",
        "angResp_onset_deg", "angDisp_onset_cm", "angDisp_mean_cm", rankB)

write.csv(bind_rows(lrt_rows), file.path(res_dir, "angular_primary_lrt.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# STEP 5: the developer-facing slope
# ---------------------------------------------------------------------------
rule("STEP 5: CAPTURE SLOPES WITH CONFIDENCE INTERVALS")

slope_report <- function(fitpair, term, label, unit) {
  f <- fitpair$reml
  sm <- summary(f)$coefficients
  est <- sm[term, "Estimate"]; se <- sm[term, "Std. Error"]
  dfv <- sm[term, "df"]; tv <- sm[term, "t value"]; pv <- sm[term, "Pr(>|t|)"]
  wald <- est + c(-1, 1) * qnorm(0.975) * se
  satt <- est + c(-1, 1) * qt(0.975, dfv) * se
  prof <- tryCatch(suppressWarnings(confint(f, parm = term, method = "profile")),
                   error = function(e) matrix(c(NA_real_, NA_real_), 1))
  say(""); say(label)
  say(sprintf("  beta = %.5f %s, SE = %.5f, t(%.1f) = %.2f, p = %s",
              est, unit, se, dfv, tv, format.pval(pv, digits = 3, eps = 1e-16)))
  say(sprintf("  95%% CI Wald          [%.5f, %.5f]", wald[1], wald[2]))
  say(sprintf("  95%% CI Satterthwaite [%.5f, %.5f]", satt[1], satt[2]))
  say(sprintf("  95%% CI profile       [%.5f, %.5f]", prof[1, 1], prof[1, 2]))
  data.frame(model = label, term = term, unit = unit, estimate = est, se = se,
             df = dfv, t = tv, p = pv, wald_lo = wald[1], wald_hi = wald[2],
             satt_lo = satt[1], satt_hi = satt[2],
             prof_lo = prof[1, 1], prof_hi = prof[1, 2])
}

slopes <- bind_rows(
  slope_report(A0_pub, "stimulusDisparity_m",
               "A0pub published: signed error (m) per m of nominal disparity", "m per m"),
  slope_report(A1, "angDisp_onset_cm",
               "A1: signed error (m) per deg of onset angular disparity", "m per deg"),
  slope_report(A2, "angDisp_mean_cm",
               "A2: signed error (m) per deg of time-weighted mean angular disparity", "m per deg"),
  slope_report(B0, "stimDisp_cm",
               "B0: angular response (deg) per m of nominal disparity", "deg per m"),
  slope_report(B1, "angDisp_onset_cm",
               "B1 PRIMARY: angular response (deg) per deg of onset angular disparity", "deg per deg"),
  slope_report(B2, "angDisp_mean_cm",
               "B2: angular response at onset (deg) per deg of time-weighted mean angular disparity", "deg per deg"),
  slope_report(B2m, "angDisp_mean_cm",
               "B2m MATCHED: time-weighted angular response (deg) per deg of time-weighted mean angular disparity", "deg per deg")
)
write.csv(slopes, file.path(res_dir, "angular_primary_slopes.csv"), row.names = FALSE)

say("")
say("Per-participant slopes for the primary angular model B1:")
rs <- coef(B1$reml)$participantId
pp <- data.frame(participantId = rownames(rs),
                 intercept_deg = rs[["(Intercept)"]],
                 slope_deg_per_deg = rs[["angDisp_onset_cm"]])
say(sprintf("  min %.3f, Q25 %.3f, median %.3f, Q75 %.3f, max %.3f, SD %.3f",
            min(pp$slope_deg_per_deg), quantile(pp$slope_deg_per_deg, .25),
            median(pp$slope_deg_per_deg), quantile(pp$slope_deg_per_deg, .75),
            max(pp$slope_deg_per_deg), sd(pp$slope_deg_per_deg)))
vcb <- VarCorr(B1$reml)
say("  random-effect SDs: ",
    paste(sprintf("%s = %.4f", names(vcb), sapply(vcb, function(z) sqrt(z[1, 1]))),
          collapse = ", "))
write.csv(pp, file.path(res_dir, "angular_primary_participant_slopes.csv"), row.names = FALSE)

say("")
say("Random-slope justification for B1 (REML LR test against a random-intercept model):")
B1_ri <- fit_lmer(as.formula(re_forms("angResp_onset_deg", "angDisp_onset_cm")[3]), TRUE)
a_rs <- anova(B1_ri, B1$reml, refit = FALSE)
say_df(cbind(model = c("random intercept", "B1"), as.data.frame(a_rs)))

# ---------------------------------------------------------------------------
# STEP 6: robustness
# ---------------------------------------------------------------------------
rule("STEP 6: ROBUSTNESS")

rob <- function(label, fml, data) {
  m <- suppressWarnings(lmer(fml, data = data, REML = TRUE, control = ctrl))
  tn <- attr(terms(m), "term.labels")[1]
  s <- summary(m)$coefficients[tn, ]
  say(sprintf("  %-46s n = %3d  beta = %+.4f  SE = %.4f  t(%.1f) = %5.2f  p = %s  clean = %s",
              label, nobs(m), s["Estimate"], s["Std. Error"], s["df"], s["t value"],
              format.pval(s["Pr(>|t|)"], digits = 3, eps = 1e-16), fit_ok(m)))
  data.frame(check = label, n = nobs(m), estimate = s["Estimate"], se = s["Std. Error"],
             df = s["df"], t = s["t value"], p = s["Pr(>|t|)"], clean_fit = fit_ok(m),
             row.names = NULL)
}

say("Primary angular slope under alternative specifications:")
f_b1 <- as.formula(re_forms("angResp_onset_deg", "angDisp_onset_cm")[rankB])
rob_tab <- bind_rows(
  rob("full sample (B1)", f_b1, dat),
  rob("onset angular disparity <= 90 deg", f_b1, dat %>% filter(angDisp_onset_deg <= 90)),
  rob("onset angular disparity <= 60 deg", f_b1, dat %>% filter(angDisp_onset_deg <= 60)),
  rob("excluding the 6 revisited trials", f_b1, dat %>% filter(!trial_has_revisit)),
  rob("final-pose frame (response and disparity)",
      as.formula(re_forms("angResp_final_deg", "angDisp_final_cm")[rankB]), dat),
  rob("matched time-weighted frame (B2m)",
      as.formula(re_forms("angResp_mean_deg", "angDisp_mean_cm")[rankB]), dat)
)
write.csv(rob_tab, file.path(res_dir, "angular_primary_robustness.csv"), row.names = FALSE)

say("")
say("Trials retained by the restrictions above: <=90 deg ",
    sum(dat$angDisp_onset_deg <= 90), ", <=60 deg ", sum(dat$angDisp_onset_deg <= 60),
    ", of ", nrow(dat))

say("")
say("Residual diagnostics for the primary angular model B1 (DHARMa, 1000 simulations):")
sim  <- simulateResiduals(B1$reml, n = 1000, seed = 20260805)
ks   <- testUniformity(sim, plot = FALSE)
disp <- testDispersion(sim, plot = FALSE)
outl <- testOutliers(sim, plot = FALSE, type = "bootstrap")
say(sprintf("  KS uniformity : D = %.4f, p = %.4g", ks$statistic, ks$p.value))
say(sprintf("  dispersion    : ratio = %.4f, p = %.4g", disp$statistic, disp$p.value))
say(sprintf("  outliers      : p = %.4g", outl$p.value))
say("Comparison, same diagnostics for the published metric model A0pub:")
simA <- simulateResiduals(A0_pub$reml, n = 1000, seed = 20260805)
ksA  <- testUniformity(simA, plot = FALSE)
dispA <- testDispersion(simA, plot = FALSE)
say(sprintf("  KS uniformity : D = %.4f, p = %.4g", ksA$statistic, ksA$p.value))
say(sprintf("  dispersion    : ratio = %.4f, p = %.4g", dispA$statistic, dispA$p.value))
say("Both responses have heavier tails than Gaussian; the angular response more so,")
say("because a near approach to the source can place the response direction far from")
say("the source direction. The two checks below establish that the slope does not")
say("depend on the Gaussian assumption.")

say("")
say("Heavy-tail check 1: Student-t response via glmmTMB. Three specifications are")
say("attempted. A fit is usable only if the TMB Hessian is positive definite")
say("(sdr$pdHess); a non-positive-definite Hessian yields a degenerate covariance")
say("matrix whose standard errors must not be reported.")

t_specs <- list(
  list(label = "student_t_random_slope",
       call = quote(glmmTMB::glmmTMB(
         angResp_onset_deg ~ angDisp_onset_cm + soundType + (1 | participantId) +
           (0 + angDisp_onset_cm | participantId),
         data = dat, family = glmmTMB::t_family(link = "identity")))),
  list(label = "student_t_random_slope_BFGS",
       call = quote(glmmTMB::glmmTMB(
         angResp_onset_deg ~ angDisp_onset_cm + soundType + (1 | participantId) +
           (0 + angDisp_onset_cm | participantId),
         data = dat, family = glmmTMB::t_family(link = "identity"),
         control = glmmTMB::glmmTMBControl(optimizer = optim,
                                           optArgs = list(method = "BFGS"))))),
  list(label = "student_t_random_intercept_only",
       call = quote(glmmTMB::glmmTMB(
         angResp_onset_deg ~ angDisp_onset_cm + soundType + (1 | participantId),
         data = dat, family = glmmTMB::t_family(link = "identity"))))
)

t_rows <- bind_rows(lapply(t_specs, function(sp) {
  m <- try(suppressWarnings(eval(sp$call)), silent = TRUE)
  if (inherits(m, "try-error")) {
    say(sprintf("  %-32s ERROR: %s", sp$label, attr(m, "condition")$message))
    return(data.frame(check = sp$label, pdHess = NA, estimate = NA_real_, se = NA_real_,
                      z = NA_real_, p = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                      usable = FALSE))
  }
  pd <- isTRUE(m$sdr$pdHess)
  cf <- tryCatch(summary(m)$coefficients$cond["angDisp_onset_cm", ],
                 error = function(e) rep(NA_real_, 4))
  if (!pd) {
    say(sprintf("  %-32s pdHess = FALSE: non-positive-definite Hessian, standard",
                sp$label))
    say(sprintf("  %-32s errors are degenerate and are NOT reported. Point estimate",
                ""))
    say(sprintf("  %-32s only, beta = %.4f deg/deg.", "", cf[1]))
    return(data.frame(check = sp$label, pdHess = FALSE, estimate = cf[1], se = NA_real_,
                      z = NA_real_, p = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
                      usable = FALSE))
  }
  ci <- confint(m, parm = "angDisp_onset_cm")
  say(sprintf("  %-32s pdHess = TRUE | beta = %.4f, SE = %.4f, z = %.2f, p = %s, 95%% CI [%.4f, %.4f]",
              sp$label, cf[1], cf[2], cf[3],
              format.pval(cf[4], digits = 3, eps = 1e-16), ci[1, 1], ci[1, 2]))
  data.frame(check = sp$label, pdHess = TRUE, estimate = cf[1], se = cf[2],
             z = cf[3], p = cf[4], ci_lo = ci[1, 1], ci_hi = ci[1, 2], usable = TRUE)
}))
say("  Caution: the random-intercept-only Student-t model omits between-participant")
say("  slope variation, which is the dominant contributor to the slope standard error")
say("  in this design, so its interval is not comparable with the primary model's and")
say("  its point estimate alone should be read. The cluster bootstrap below is the")
say("  assumption-free check.")
t_row <- t_rows

say("")
say("Heavy-tail check 2: participant-level (cluster) bootstrap of the primary slope,")
say("1000 resamples of participants with replacement, percentile interval.")
set.seed(20260805)
pids <- unique(dat$participantId)
boot_slope <- vapply(seq_len(1000), function(i) {
  take <- sample(pids, length(pids), replace = TRUE)
  bd <- bind_rows(lapply(seq_along(take), function(j) {
    z <- dat[dat$participantId == take[j], , drop = FALSE]
    z$participantId <- paste0("B", j)
    z
  }))
  m <- try(suppressWarnings(lmer(f_b1, data = bd, REML = TRUE, control = ctrl)), silent = TRUE)
  if (inherits(m, "try-error")) NA_real_ else unname(fixef(m)["angDisp_onset_cm"])
}, numeric(1))
n_ok <- sum(!is.na(boot_slope))
bci <- quantile(boot_slope, c(.025, .975), na.rm = TRUE)
say(sprintf("  successful resamples: %d of 1000", n_ok))
say(sprintf("  bootstrap mean %.4f, SD %.4f, 95%% percentile CI [%.4f, %.4f]",
            mean(boot_slope, na.rm = TRUE), sd(boot_slope, na.rm = TRUE), bci[1], bci[2]))
say(sprintf("  resamples with a slope <= 0: %d (%.2f%%)",
            sum(boot_slope <= 0, na.rm = TRUE), 100 * mean(boot_slope <= 0, na.rm = TRUE)))

write.csv(bind_rows(
  t_row,
  data.frame(check = "cluster_bootstrap", pdHess = NA, estimate = mean(boot_slope, na.rm = TRUE),
             se = sd(boot_slope, na.rm = TRUE), z = NA_real_,
             p = NA_real_, ci_lo = bci[1], ci_hi = bci[2], usable = TRUE)
), file.path(res_dir, "angular_primary_heavytail_checks.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# STEP 7: output
# ---------------------------------------------------------------------------
rule("STEP 7: OUTPUT")

out_cols <- dat %>%
  select(participantId, trialSequenceNum, soundType, disparityRange,
         stimulusDisparity_m, signedError_m, participantError_m,
         head_onset_x, head_onset_y, head_onset_z,
         angDisp_onset_deg, angDisp_mean_deg, angDisp_median_deg, angDisp_min_deg,
         angDisp_max_deg, angDisp_range_deg, angDisp_final_deg,
         angResp_onset_deg, angResp_mean_deg, angResp_final_deg,
         angBiasProportion_onset, angBiasProportion_mean,
         dist_at_onset_m, dist_mean_m, dist_min_m, trial_has_revisit)
write.csv(out_cols, file.path(res_dir, "angular_primary_trial_level.csv"), row.names = FALSE)

say("Files written to ", res_dir, ":")
for (f in c("angular_primary_log.txt", "angular_primary_trial_level.csv",
            "angular_primary_descriptives.csv", "angular_primary_coefficients.csv",
            "angular_primary_model_comparison.csv", "angular_primary_lrt.csv",
            "angular_primary_slopes.csv", "angular_primary_participant_slopes.csv",
            "angular_primary_within_trial_variation.csv",
            "angular_primary_predictor_correlations.csv",
            "angular_primary_robustness.csv",
            "angular_primary_heavytail_checks.csv")) say("  ", f)

say("")
say("Done.")
close(log_con)
