# verify_tolerance_deliverable.R
#
# Independent adversarial verification of tolerance_deliverable_revision.R.
# Refits every headline model from the saved trial-level dataset, cross-checks
# optimisers and packages, inspects the conditioning of every variance-covariance
# matrix, and re-derives the lookup-table, capture-model and crossing-point numbers.
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/verify_tolerance_log.txt and verify_tolerance_*.csv
#
# Deterministic. Run with R 4.4 (arm64).

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(DHARMa)
  library(dplyr)
  library(tidyr)
})

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir  <- file.path(proj_dir, "results_revision")
log_path <- file.path(res_dir, "verify_tolerance_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

say  <- function(...) cat(..., "\n", sep = "")
rule <- function(t) say("\n", strrep("=", 78), "\n", t, "\n", strrep("=", 78))

say("verify_tolerance_deliverable.R  |  ", R.version.string)
say("run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))

N_BOOT <- 2000
Z90 <- qnorm(0.90)

# ---------------------------------------------------------------------------
rule("0. DATA AND DERIVED VARIABLES")

d <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
say("trials ", nrow(d), "  participants ", dplyr::n_distinct(d$participantId))
say("complete cases on the modelling variables: ",
    sum(complete.cases(d[, c("signedError_m", "stimulusDisparity_m", "soundType",
                             "angDisp_mean_deg", "participantId")])))

d <- d %>% mutate(
  disp_cm      = stimulusDisparity_m * 100,
  ang_deg      = angDisp_mean_deg,
  angOnset_deg = angDisp_onset_deg,
  ang10        = ang_deg / 10,
  angOnset10   = angOnset_deg / 10,
  disp10       = disp_cm / 10
)
mean_disp10 <- mean(d$disp10); mean_ang10 <- mean(d$ang10)
d$disp10c <- d$disp10 - mean_disp10
d$ang10c  <- d$ang10  - mean_ang10

say("sound type counts: ", paste(names(table(d$soundType)), as.integer(table(d$soundType)),
                                 sep = "=", collapse = ", "))
say("disparity cm: range ", paste(round(range(d$disp_cm), 2), collapse = "-"),
    "  mean ", round(mean(d$disp_cm), 3))
say("ang_deg: range ", paste(round(range(d$ang_deg), 2), collapse = "-"),
    "  mean ", round(mean(d$ang_deg), 3), "  median ", round(median(d$ang_deg), 3))

# ---------------------------------------------------------------------------
rule("1. CAPTURE OUTCOME: ARE THE TWO CRITERIA THE SAME EVENT?")

dist_to_src   <- sqrt((d$response_x - d$sound_x)^2 + (d$response_y - d$sound_y)^2 +
                        (d$response_z - d$sound_z)^2)
dist_to_flash <- sqrt((d$response_x - d$flash_x)^2 + (d$response_y - d$flash_y)^2 +
                        (d$response_z - d$flash_z)^2)
d$capture_geom <- as.integer(dist_to_flash < dist_to_src)
d$capture_bias <- as.integer(d$ventriloquistBias > 0.5)
say("agreement geometric vs bias>0.5: ", sum(d$capture_geom == d$capture_bias), "/", nrow(d))
say("closest trial to the boundary: min |bias - 0.5| = ",
    signif(min(abs(d$ventriloquistBias - 0.5)), 4))
d$capture <- d$capture_geom
d$capture_soft   <- as.integer(d$ventriloquistBias > 0.25)
d$capture_strict <- as.integer(d$ventriloquistBias > 0.75)
say("base rates: >0.5 ", round(mean(d$capture), 4), " (", sum(d$capture), "/", nrow(d), ")",
    "   >0.25 ", round(mean(d$capture_soft), 4), "   >0.75 ", round(mean(d$capture_strict), 4))
pp <- tapply(d$capture, d$participantId, mean)
say("per-participant capture rate: min ", round(min(pp), 4), " median ", round(median(pp), 4),
    " max ", round(max(pp), 4))

# ---------------------------------------------------------------------------
rule("2. SIGNED-ERROR LMMs, REFITTED ACROSS OPTIMISERS")

fit_lmm <- function(form, opt) {
  eval(bquote(lmer(.(form), data = d, REML = TRUE,
                   control = lmerControl(optimizer = .(opt)))))
}
report_lmm <- function(fit, label) {
  fx <- fixef(fit); se <- sqrt(diag(as.matrix(vcov(fit))))
  vc <- VarCorr(fit)$participantId
  sd <- attr(vc, "stddev"); cr <- attr(vc, "correlation")
  msg <- fit@optinfo$conv$lme4$messages
  say(sprintf("%-34s b1=%.6f se=%.6f  tau0=%.6f tau1=%.6f rho=%+.4f sigma=%.6f  logLik=%.4f",
              label, fx[2], se[2], sd[1], sd[2], cr[1, 2], sigma(fit), as.numeric(logLik(fit))))
  say("   convergence messages: ", if (length(msg)) paste(msg, collapse = " | ") else "none",
      "   singular: ", isSingular(fit))
  invisible(NULL)
}

f_metric <- signedError_m ~ stimulusDisparity_m + soundType + (1 + stimulusDisparity_m | participantId)
f_ang10  <- signedError_m ~ ang10 + soundType + (1 + ang10 | participantId)
f_angraw <- signedError_m ~ ang_deg + soundType + (1 + ang_deg | participantId)

say("\n-- M1 metric --")
m_metric <- fit_lmm(f_metric, "bobyqa")
report_lmm(m_metric, "bobyqa")
m1_nm <- fit_lmm(f_metric, "Nelder_Mead"); report_lmm(m1_nm, "Nelder_Mead")
m1_nl <- fit_lmm(f_metric, "nloptwrap"); report_lmm(m1_nl, "nloptwrap")
print(round(summary(m_metric)$coefficients, 6))

# Fully independent implementation: the same Gaussian LMM in glmmTMB (REML, TMB
# Laplace/AD machinery rather than lme4's PLS), as a check on the SEs.
tmb_metric <- glmmTMB(f_metric, data = d, REML = TRUE)
say("glmmTMB REML refit of M1: b1 = ", signif(fixef(tmb_metric)$cond[2], 8),
    "  se = ", signif(summary(tmb_metric)$coefficients$cond[2, 2], 8),
    "  tau0 = ", signif(attr(VarCorr(tmb_metric)$cond$participantId, "stddev")[1], 6),
    "  tau1 = ", signif(attr(VarCorr(tmb_metric)$cond$participantId, "stddev")[2], 6),
    "  rho = ", signif(attr(VarCorr(tmb_metric)$cond$participantId, "correlation")[1, 2], 5),
    "  pdHess ", tmb_metric$sdr$pdHess)

say("\n-- M2 angular (per 10 deg) --")
m_ang <- fit_lmm(f_ang10, "bobyqa")
report_lmm(m_ang, "bobyqa")
m2_nm <- fit_lmm(f_ang10, "Nelder_Mead"); report_lmm(m2_nm, "Nelder_Mead")
print(round(summary(m_ang)$coefficients, 6))

say("\n-- M2 angular on the raw per-degree scale (scaling-artifact claim) --")
m_ang_raw <- fit_lmm(f_angraw, "bobyqa")
report_lmm(m_ang_raw, "bobyqa raw deg")
say("raw-scale slope x 10 = ", signif(fixef(m_ang_raw)[2] * 10, 8),
    "   per-10-deg slope = ", signif(fixef(m_ang)[2], 8),
    "   relative difference = ", signif(abs(fixef(m_ang_raw)[2] * 10 / fixef(m_ang)[2] - 1), 3))

# ML refits, to make sure the REML/ML choice does not move the headline values.
say("\n-- ML refits (REML = FALSE) --")
m1_ml <- update(m_metric, REML = FALSE); report_lmm(m1_ml, "M1 ML")
m2_ml <- update(m_ang, REML = FALSE); report_lmm(m2_ml, "M2 ML")

# vcov conditioning of the fixed effects.
for (nm in c("m_metric", "m_ang")) {
  V <- as.matrix(vcov(get(nm)))
  say(nm, " fixed-effect vcov: condition number ", signif(kappa(V), 5),
      "  min eigenvalue ", signif(min(eigen(V, only.values = TRUE)$values), 5))
}

# ---------------------------------------------------------------------------
rule("3. PROFILE LIKELIHOOD ON THE RANDOM-EFFECT SDs (IS tau1 WELL IDENTIFIED?)")

pp_metric <- try(profile(update(m_metric, REML = TRUE), which = "theta_", signames = FALSE),
                 silent = TRUE)
if (!inherits(pp_metric, "try-error")) {
  print(round(confint(pp_metric, level = 0.95), 5))
} else {
  say("profile failed: ", as.character(pp_metric))
}

# ---------------------------------------------------------------------------
rule("4. LISTENER PROFILES AND LOOKUP TABLE, RECOMPUTED")

sound_offset <- function(fx) mean(c(0, fx[grep("^soundType", names(fx))]))
lmm_pars <- function(fit) {
  fx <- fixef(fit); vc <- VarCorr(fit)$participantId
  sd <- attr(vc, "stddev"); cr <- attr(vc, "correlation")
  c(b0 = unname(fx[1]) + sound_offset(fx), b1 = unname(fx[2]),
    tau0 = unname(sd[1]), tau1 = unname(sd[2]), rho = unname(cr[1, 2]))
}
listener_pars <- function(p, pars) {
  b1d <- qnorm(p) * pars["tau1"]
  b0d <- if (pars["tau1"] > 0) pars["rho"] * (pars["tau0"] / pars["tau1"]) * b1d else 0
  c(intercept = unname(pars["b0"] + b0d), slope = unname(pars["b1"] + b1d))
}
predict_profiles <- function(pars, xg) {
  p10 <- listener_pars(0.10, pars); p90 <- listener_pars(0.90, pars)
  unname(c(pars["b0"] + pars["b1"] * xg,
           p10["intercept"] + p10["slope"] * xg,
           p90["intercept"] + p90["slope"] * xg))
}

pars_metric <- lmm_pars(m_metric); pars_ang <- lmm_pars(m_ang)
say("metric pars: ", paste(names(pars_metric), signif(pars_metric, 6), sep = "=", collapse = "  "))
say("angular pars: ", paste(names(pars_ang), signif(pars_ang, 6), sep = "=", collapse = "  "))
say("10th/90th pct listener slopes (metric): ",
    signif(listener_pars(0.10, pars_metric)["slope"], 5), " / ",
    signif(listener_pars(0.90, pars_metric)["slope"], 5))
say("10th/90th pct listener intercepts (metric, cm): ",
    signif(listener_pars(0.10, pars_metric)["intercept"] * 100, 5), " / ",
    signif(listener_pars(0.90, pars_metric)["intercept"] * 100, 5))
say("population intercept incl. sound offset (cm): ", signif(pars_metric["b0"] * 100, 5))

blup <- coef(m_metric)$participantId[["stimulusDisparity_m"]]
say("BLUP slopes: min ", round(min(blup), 4), " q10 ", round(quantile(blup, .10), 4),
    " median ", round(median(blup), 4), " q90 ", round(quantile(blup, .90), 4),
    " max ", round(max(blup), 4), "  SD ", round(sd(blup), 4))

grid_cm  <- seq(15, 70, by = 5)
grid_deg <- seq(20, 70, by = 5)

# ---------------------------------------------------------------------------
rule("5. PARAMETRIC BOOTSTRAP, REPRODUCED, PLUS THE EFFECT OF DISCARDING REPLICATES")

boot_profiles <- function(fit, xg, label) {
  fn <- function(f) predict_profiles(lmm_pars(f), xg)
  bb <- bootMer(fit, fn, nsim = N_BOOT, seed = 20260805, type = "parametric")
  n_bad <- sum(!complete.cases(bb$t))
  say(label, ": ", N_BOOT, " replicates, ", n_bad, " incomplete, ", N_BOOT - n_bad, " retained")
  tt <- bb$t
  list(all = tt,
       kept = tt[complete.cases(tt), , drop = FALSE],
       n_bad = n_bad)
}
bo_metric <- boot_profiles(m_metric, grid_cm / 100, "metric")
bo_ang    <- boot_profiles(m_ang, grid_deg / 10, "angular")

qlo <- function(m) apply(m, 2, quantile, 0.025, na.rm = TRUE)
qhi <- function(m) apply(m, 2, quantile, 0.975, na.rm = TRUE)

k <- length(grid_cm)
tabA <- tibble(
  disparity_cm = grid_cm,
  avg_cm  = predict_profiles(pars_metric, grid_cm / 100)[1:k] * 100,
  avg_lo  = qlo(bo_metric$kept)[1:k] * 100,
  avg_hi  = qhi(bo_metric$kept)[1:k] * 100,
  p10_cm  = predict_profiles(pars_metric, grid_cm / 100)[(k + 1):(2 * k)] * 100,
  p10_lo  = qlo(bo_metric$kept)[(k + 1):(2 * k)] * 100,
  p10_hi  = qhi(bo_metric$kept)[(k + 1):(2 * k)] * 100,
  p90_cm  = predict_profiles(pars_metric, grid_cm / 100)[(2 * k + 1):(3 * k)] * 100,
  p90_lo  = qlo(bo_metric$kept)[(2 * k + 1):(3 * k)] * 100,
  p90_hi  = qhi(bo_metric$kept)[(2 * k + 1):(3 * k)] * 100
) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
say("\n-- Table A1 recomputed (complete-case bootstrap, as in the original) --")
print(as.data.frame(tabA))

# Average-listener CI is unaffected by NA in the percentile columns; the original
# nevertheless dropped whole rows. Quantify what that costs.
say("\n-- Average-listener CI: all 2000 replicates vs the 1896 complete-case subset --")
avg_all  <- bo_metric$all[, 1:k, drop = FALSE]
cmp_avg <- tibble(disparity_cm = grid_cm,
                  lo_all  = round(qlo(avg_all) * 100, 3),  hi_all  = round(qhi(avg_all) * 100, 3),
                  lo_kept = round(qlo(bo_metric$kept)[1:k] * 100, 3),
                  hi_kept = round(qhi(bo_metric$kept)[1:k] * 100, 3))
print(as.data.frame(cmp_avg))
say("max absolute difference in the average-listener CI limits (cm): ",
    round(max(abs(c(cmp_avg$lo_all - cmp_avg$lo_kept, cmp_avg$hi_all - cmp_avg$hi_kept))), 4))

# Which replicates were dropped, and what were their tau1 values?
fn_tau <- function(f) { p <- lmm_pars(f); c(p["tau1"], p["tau0"], p["rho"]) }
set.seed(1)
bb_tau <- bootMer(m_metric, fn_tau, nsim = 400, seed = 20260805, type = "parametric")
tau_t <- bb_tau$t
say("\nbootstrap tau1 distribution (400 replicates): ",
    "n with tau1 < 1e-6 = ", sum(tau_t[, 1] < 1e-6),
    "; n with NA rho = ", sum(is.na(tau_t[, 3])),
    "; tau1 quantiles ", paste(round(quantile(tau_t[, 1], c(.025, .5, .975), na.rm = TRUE), 4),
                               collapse = " "))

ka <- length(grid_deg)
tabB <- tibble(
  angle_deg = grid_deg,
  avg_cm = predict_profiles(pars_ang, grid_deg / 10)[1:ka] * 100,
  avg_lo = qlo(bo_ang$kept)[1:ka] * 100, avg_hi = qhi(bo_ang$kept)[1:ka] * 100,
  p90_cm = predict_profiles(pars_ang, grid_deg / 10)[(2 * ka + 1):(3 * ka)] * 100,
  p90_lo = qlo(bo_ang$kept)[(2 * ka + 1):(3 * ka)] * 100,
  p90_hi = qhi(bo_ang$kept)[(2 * ka + 1):(3 * ka)] * 100,
  n_within_2p5deg = sapply(grid_deg, function(g) sum(abs(d$ang_deg - g) <= 2.5))
) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
say("\n-- Table A2 recomputed --")
print(as.data.frame(tabB))

# Crosswalk medians quoted in the report's Table A1 last column.
cw <- bind_rows(lapply(grid_cm, function(g) {
  s <- d[abs(d$disp_cm - g) <= 2.5, ]
  tibble(disparity_cm = g, n = nrow(s), ang_median_deg = median(s$ang_deg),
         ang_mean_deg = mean(s$ang_deg))
}))
say("\n-- Nominal-to-angular crosswalk --")
print(as.data.frame(cw %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))

# ---------------------------------------------------------------------------
rule("6. CAPTURE MODELS: glmmTMB REFIT, DEGENERACY AND CROSS-PACKAGE CHECK")

show_tmb <- function(fit, label) {
  co <- summary(fit)$coefficients$cond
  vc <- VarCorr(fit)$cond$participantId
  sd <- attr(vc, "stddev"); cr <- attr(vc, "correlation")
  say("\n-- ", label, " --")
  print(round(co, 6))
  say("random SDs: ", paste(round(sd, 6), collapse = ", "),
      if (length(sd) > 1) paste0("   corr = ", round(cr[1, 2], 5)) else "")
  say("convergence ", fit$fit$convergence, " | ", fit$fit$message,
      " | pdHess ", fit$sdr$pdHess, " | logLik ", round(as.numeric(logLik(fit)), 5))
  # conditioning of the FULL joint Hessian, including the variance parameters
  h <- fit$sdr$cov.fixed
  if (!is.null(h)) {
    ev <- eigen(h, only.values = TRUE)$values
    say("joint cov of all estimated parameters: dim ", nrow(h),
        " | min eigenvalue ", signif(min(ev), 4),
        " | condition number ", signif(max(ev) / min(ev), 5))
  }
  invisible(co)
}

g_ang  <- glmmTMB(capture ~ ang10c  + (1 + ang10c  | participantId), data = d, family = binomial)
g_nom  <- glmmTMB(capture ~ disp10c + (1 + disp10c | participantId), data = d, family = binomial)
g_null <- glmmTMB(capture ~ 1 + (1 | participantId), data = d, family = binomial)
show_tmb(g_ang, "B1 angular, glmmTMB")
show_tmb(g_nom, "B2 nominal, glmmTMB")
show_tmb(g_null, "B0 disparity-free, glmmTMB")

# Independent cross-check with lme4::glmer (Laplace and 15-point AGQ where possible).
say("\n-- lme4::glmer cross-check --")
gl_ang <- glmer(capture ~ ang10c + (1 + ang10c | participantId), data = d, family = binomial,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
gl_nom <- glmer(capture ~ disp10c + (1 + disp10c | participantId), data = d, family = binomial,
                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
gl_null <- glmer(capture ~ 1 + (1 | participantId), data = d, family = binomial,
                 control = glmerControl(optimizer = "bobyqa"))
for (nm in c("gl_ang", "gl_nom", "gl_null")) {
  f <- get(nm); co <- summary(f)$coefficients
  vcc <- VarCorr(f)$participantId
  say(nm, ": ", paste(rownames(co), " b=", round(co[, 1], 5), " se=", round(co[, 2], 5),
                      " p=", signif(co[, 4], 4), sep = "", collapse = " | "))
  say("   SDs ", paste(round(attr(vcc, "stddev"), 5), collapse = ", "),
      if (nrow(vcc) > 1) paste0("  corr ", round(attr(vcc, "correlation")[1, 2], 5)) else "",
      "  singular ", isSingular(f), "  logLik ", round(as.numeric(logLik(f)), 5))
  msg <- f@optinfo$conv$lme4$messages
  say("   messages: ", if (length(msg)) paste(msg, collapse = " | ") else "none")
}
gl_null_agq <- glmer(capture ~ 1 + (1 | participantId), data = d, family = binomial, nAGQ = 15)
say("B0 with 15-point AGQ: b0 = ", round(fixef(gl_null_agq)[1], 5),
    "  se = ", round(sqrt(vcov(gl_null_agq)[1, 1]), 5),
    "  tau0 = ", round(attr(VarCorr(gl_null_agq)$participantId, "stddev")[1], 5))

# ---------------------------------------------------------------------------
rule("7. LIKELIHOOD-RATIO TESTS AND DECOMPOSITION")

g_ang_ri <- glmmTMB(capture ~ ang10c  + (1 | participantId), data = d, family = binomial)
g_nom_ri <- glmmTMB(capture ~ disp10c + (1 | participantId), data = d, family = binomial)
say("angular, 3 df:"); print(anova(g_null, g_ang))
say("nominal, 3 df:"); print(anova(g_null, g_nom))
say("angular fixed slope, 1 df:"); print(anova(g_null, g_ang_ri))
say("angular random slope given fixed, 2 df:"); print(anova(g_ang_ri, g_ang))
say("nominal fixed slope, 1 df:"); print(anova(g_null, g_nom_ri))
say("nominal random slope given fixed, 2 df:"); print(anova(g_nom_ri, g_nom))

g_ang_q <- glmmTMB(capture ~ ang10c + I(ang10c^2) + (1 + ang10c | participantId),
                   data = d, family = binomial)
say("\nquadratic angular:"); print(summary(g_ang_q)$coefficients$cond)
print(anova(g_ang, g_ang_q))

# Within-between decomposition: is the angular effect within participants or purely
# a between-participant contrast? This is the check the report's narrative assumes.
d <- d %>% group_by(participantId) %>%
  mutate(ang10_pm = mean(ang10), ang10_wi = ang10 - ang10_pm) %>% ungroup()
d$ang10_pmc <- d$ang10_pm - mean(d$ang10_pm)
g_wb <- glmmTMB(capture ~ ang10_wi + ang10_pmc + (1 | participantId), data = d, family = binomial)
say("\nwithin/between decomposition of the angular effect (random intercept only):")
print(round(summary(g_wb)$coefficients$cond, 6))
g_wb_rs <- glmmTMB(capture ~ ang10_wi + ang10_pmc + (1 + ang10_wi | participantId),
                   data = d, family = binomial)
say("same with a random within-slope:")
print(round(summary(g_wb_rs)$coefficients$cond, 6))
say("convergence ", g_wb_rs$fit$convergence, " pdHess ", g_wb_rs$sdr$pdHess)

# ---------------------------------------------------------------------------
rule("8. DHARMa ON THE ANGULAR CAPTURE MODEL, SEED STABILITY")

for (s in c(20260805, 1, 2, 3)) {
  sr <- simulateResiduals(g_ang, n = 1000, seed = s)
  u <- testUniformity(sr, plot = FALSE); di <- testDispersion(sr, plot = FALSE)
  say("seed ", s, ": KS D = ", signif(unname(u$statistic), 5), " p = ", signif(u$p.value, 4),
      " | dispersion ", signif(unname(di$statistic), 5), " p = ", signif(di$p.value, 4))
}

# ---------------------------------------------------------------------------
rule("9. CROSSING POINTS AND FIELLER SETS")

fieller <- function(b, V, target, level = 0.95) {
  b <- unname(b); V <- unname(as.matrix(V)); z <- qnorm(1 - (1 - level) / 2)
  a2 <- b[2]^2 - z^2 * V[2, 2]
  a1 <- 2 * ((b[1] - target) * b[2] - z^2 * V[1, 2])
  a0 <- (b[1] - target)^2 - z^2 * V[1, 1]
  point <- (target - b[1]) / b[2]
  if (a2 <= 0) return(c(point = point, lo = -Inf, hi = Inf, bounded = 0))
  disc <- a1^2 - 4 * a2 * a0
  if (disc <= 0) return(c(point = point, lo = NA, hi = NA, bounded = -1))
  r <- sort((-a1 + c(-1, 1) * sqrt(disc)) / (2 * a2))
  c(point = point, lo = r[1], hi = r[2], bounded = 1)
}
crossings <- function(fit, scale, centre, label) {
  b <- fixef(fit)$cond; V <- vcov(fit)$cond
  vc <- VarCorr(fit)$cond$participantId
  tau0 <- unname(attr(vc, "stddev")[1]); tau1 <- unname(attr(vc, "stddev")[2])
  rho <- unname(attr(vc, "correlation")[1, 2])
  out <- list()
  for (p in c(0.50, 0.75)) for (who in c("average", "p10", "p90")) {
    if (who == "average") bb <- b else {
      pr <- if (who == "p10") 0.10 else 0.90
      b0d <- qnorm(pr) * tau0
      b1d <- rho * (tau1 / tau0) * b0d
      bb <- c(b[1] + b0d, b[2] + b1d)
    }
    fi <- fieller(bb, V, qlogis(p))
    bt <- function(v) if (is.finite(v)) (v + centre) * scale else v
    out[[length(out) + 1]] <- tibble(model = label, listener = who, prob = p,
                                     point = bt(unname(fi["point"])),
                                     lo = bt(unname(fi["lo"])), hi = bt(unname(fi["hi"])),
                                     bounded = unname(fi["bounded"]))
  }
  bind_rows(out)
}
cr_all <- bind_rows(crossings(g_ang, 10, mean_ang10, "B1 angular deg"),
                    crossings(g_nom, 10, mean_disp10, "B2 nominal cm"))
print(as.data.frame(cr_all %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
say("observed max angle ", round(max(d$ang_deg), 2), " deg; max disparity ",
    round(max(d$disp_cm), 2), " cm")

# How sensitive is the 90th-percentile angular crossing to propagating tau0/tau1/rho?
# Delta-method-free check: bootstrap the whole fit and see how many replicates are usable.
say("\nprobe: 300-replicate parametric bootstrap of B1, retention of usable refits")
set.seed(20260805)
sims_b1 <- simulate(g_ang, nsim = 300)
res_b1 <- vapply(sims_b1, function(y) {
  dd <- d; dd$capture <- y
  f <- try(suppressWarnings(update(g_ang, data = dd)), silent = TRUE)
  if (inherits(f, "try-error") || is.null(f$sdr) || !isTRUE(f$sdr$pdHess) ||
      f$fit$convergence != 0) return(rep(NA_real_, 6))
  vc <- VarCorr(f)$cond$participantId
  b <- fixef(f)$cond
  c(b[1], b[2], attr(vc, "stddev")[1], attr(vc, "stddev")[2],
    attr(vc, "correlation")[1, 2], 1)
}, numeric(6))
ok <- !is.na(res_b1[6, ])
say("usable replicates: ", sum(ok), "/300 = ", round(mean(ok), 4))
if (sum(ok) > 30) {
  b0b <- res_b1[1, ok]; b1b <- res_b1[2, ok]
  t0b <- res_b1[3, ok]; t1b <- res_b1[4, ok]; rob <- res_b1[5, ok]
  x50 <- ((qlogis(0.5) - (b0b + Z90 * t0b)) / (b1b + rob * (t1b / t0b) * Z90 * t0b) + mean_ang10) * 10
  say("fully propagated 90th-pct 50% crossing (deg): median ", round(median(x50, na.rm = TRUE), 2),
      "  2.5% ", round(quantile(x50, .025, na.rm = TRUE), 2),
      "  97.5% ", round(quantile(x50, .975, na.rm = TRUE), 2),
      "  fraction of replicates outside 0-115.5 deg: ",
      round(mean(x50 < 0 | x50 > 115.5, na.rm = TRUE), 3))
  say("bootstrap correlation rho: median ", round(median(rob, na.rm = TRUE), 3),
      "  fraction |rho| > 0.99: ", round(mean(abs(rob) > 0.99, na.rm = TRUE), 3))
}

# ---------------------------------------------------------------------------
rule("10. DISPARITY-FREE LISTENER PROBABILITIES, FULL 2000-REPLICATE BOOTSTRAP")

b0_null <- fixef(g_null)$cond[1]
tau0_null <- attr(VarCorr(g_null)$cond$participantId, "stddev")[1]
say("b0 = ", signif(b0_null, 6), "  se = ", signif(summary(g_null)$coefficients$cond[1, 2], 6),
    "  tau0 = ", signif(tau0_null, 6))
say("profile CI on tau0:")
print(try(confint(g_null, parm = "theta_", method = "profile"), silent = TRUE))
print(try(confint(g_null), silent = TRUE))

set.seed(20260805)
sims_null <- simulate(g_null, nsim = N_BOOT)
boot_null <- vapply(sims_null, function(y) {
  dd <- d; dd$capture <- y
  f <- try(suppressWarnings(update(g_null, data = dd)), silent = TRUE)
  if (inherits(f, "try-error") || !isTRUE(f$sdr$pdHess) || f$fit$convergence != 0)
    return(c(NA, NA))
  c(fixef(f)$cond[1], attr(VarCorr(f)$cond$participantId, "stddev")[1])
}, numeric(2))
keep <- complete.cases(t(boot_null))
say("B0 bootstrap: ", N_BOOT, " requested, ", sum(keep), " retained")
prob_at <- function(b0, tau, p) plogis(b0 + qnorm(p) * tau)
cap_tab <- bind_rows(lapply(c(0.10, 0.50, 0.90), function(p) {
  bs <- prob_at(boot_null[1, keep], boot_null[2, keep], p)
  tibble(pct = p * 100, est = prob_at(b0_null, tau0_null, p),
         lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975)))
})) %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
print(as.data.frame(cap_tab))
pct_at <- function(t) 100 * pnorm((qlogis(t) - b0_null) / tau0_null)
say("listener percentile reaching 50%: ", round(pct_at(0.50), 2),
    "   reaching 75%: ", round(pct_at(0.75), 2))
# bootstrap uncertainty on those percentiles
p50b <- 100 * pnorm((qlogis(0.5) - boot_null[1, keep]) / boot_null[2, keep])
p75b <- 100 * pnorm((qlogis(0.75) - boot_null[1, keep]) / boot_null[2, keep])
say("bootstrap 95% CI on the 50% percentile: [", round(quantile(p50b, .025), 1), ", ",
    round(quantile(p50b, .975), 1), "]   on the 75% percentile: [",
    round(quantile(p75b, .025), 1), ", ", round(quantile(p75b, .975), 1), "]")
# observed per-participant rates for comparison
say("observed per-participant capture rates, quantiles 10/50/90: ",
    paste(round(quantile(pp, c(.1, .5, .9)), 4), collapse = " "))
say("observed fraction of participants with rate > 0.5: ", round(mean(pp > 0.5), 4),
    "  > 0.75: ", round(mean(pp > 0.75), 4))

# ---------------------------------------------------------------------------
rule("11. EMPIRICAL BINS AND CRITERION ROBUSTNESS")

binci <- function(x, n) as.numeric(binom.test(x, n)$conf.int)
emp_nom <- d %>%
  mutate(bin = cut(disp_cm, c(15, 25, 35, 45, 55, 70), include.lowest = TRUE, right = FALSE)) %>%
  group_by(bin) %>% summarise(n = n(), k = sum(capture), p = mean(capture), .groups = "drop") %>%
  rowwise() %>% mutate(lo = binci(k, n)[1], hi = binci(k, n)[2]) %>% ungroup()
print(as.data.frame(emp_nom %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))
emp_ang <- d %>%
  mutate(bin = cut(ang_deg, c(0, 20, 30, 40, 50, 60, Inf), right = FALSE)) %>%
  group_by(bin) %>% summarise(n = n(), k = sum(capture), p = mean(capture), .groups = "drop") %>%
  rowwise() %>% mutate(lo = binci(k, n)[1], hi = binci(k, n)[2]) %>% ungroup()
print(as.data.frame(emp_ang %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))

rob <- bind_rows(lapply(list(c("capture", "bias>0.50"), c("capture_soft", "bias>0.25"),
                             c("capture_strict", "bias>0.75")), function(cr) {
  bind_rows(lapply(c("ang10c", "disp10c"), function(xv) {
    f <- reformulate(c(xv, sprintf("(1 + %s | participantId)", xv)), response = cr[1])
    fit <- suppressWarnings(glmmTMB(f, data = d, family = binomial))
    co <- summary(fit)$coefficients$cond
    vc <- VarCorr(fit)$cond$participantId
    tibble(criterion = cr[2], predictor = xv, base = mean(d[[cr[1]]]),
           slope = co[2, 1], se = co[2, 2], z = co[2, 3], p = co[2, 4],
           tau0 = attr(vc, "stddev")[1], tau1 = attr(vc, "stddev")[2],
           rho = attr(vc, "correlation")[1, 2],
           conv = fit$fit$convergence, pdHess = isTRUE(fit$sdr$pdHess))
  }))
})) %>% mutate(across(where(is.numeric), ~ signif(.x, 5)))
print(as.data.frame(rob))

# Proportional-scaling test.
m_prop <- lmerTest::lmer(ventriloquistBias ~ disp_cm + (1 | participantId), data = d)
say("\nproportional-scaling test:"); print(round(summary(m_prop)$coefficients, 6))
# variance of ventriloquistBias by disparity bin, the mechanism the report invokes
vb <- d %>% mutate(bin = cut(disp_cm, c(15, 25, 35, 45, 55, 70), include.lowest = TRUE, right = FALSE)) %>%
  group_by(bin) %>% summarise(n = n(), sd_bias = sd(ventriloquistBias),
                              sd_signed_cm = sd(signedError_m) * 100, .groups = "drop")
print(as.data.frame(vb %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))

# ---------------------------------------------------------------------------
rule("12. EQUIVALENCE / POWER: WHAT SLOPE COULD THE NOMINAL CAPTURE MODEL RULE OUT?")

co_nom <- summary(g_nom)$coefficients$cond
say("nominal slope per 10 cm: ", signif(co_nom[2, 1], 5), " se ", signif(co_nom[2, 2], 5),
    "  95% CI [", signif(co_nom[2, 1] - 1.96 * co_nom[2, 2], 4), ", ",
    signif(co_nom[2, 1] + 1.96 * co_nom[2, 2], 4), "]")
p_at_centre <- plogis(co_nom[1, 1])
for (b in c(co_nom[2, 1] - 1.96 * co_nom[2, 2], co_nom[2, 1] + 1.96 * co_nom[2, 2])) {
  say("  at the CI limit b = ", signif(b, 4), ": predicted capture at 15 cm = ",
      round(plogis(co_nom[1, 1] + b * (1.5 - mean_disp10)), 4),
      ", at 70 cm = ", round(plogis(co_nom[1, 1] + b * (7.0 - mean_disp10)), 4))
}
co_a <- summary(g_ang)$coefficients$cond
for (b in c(co_a[2, 1] - 1.96 * co_a[2, 2], co_a[2, 1] + 1.96 * co_a[2, 2])) {
  say("angular CI limit b = ", signif(b, 4), ": capture at 20 deg = ",
      round(plogis(co_a[1, 1] + b * (2.0 - mean_ang10)), 4), ", at 70 deg = ",
      round(plogis(co_a[1, 1] + b * (7.0 - mean_ang10)), 4))
}

# ---------------------------------------------------------------------------
rule("13. TRAJECTORY CLAIMS: DISTANCE AND EXPERIENCED ANGLE")

say("head-to-source distance, per-trial time-weighted mean (m): mean ",
    round(mean(d$dist_mean_m), 4), " median ", round(median(d$dist_mean_m), 4),
    " range ", paste(round(range(d$dist_mean_m), 3), collapse = "-"))
say("distance at onset (m): mean ", round(mean(d$dist_at_onset_m), 4),
    " median ", round(median(d$dist_at_onset_m), 4))
say("experienced angle (time-weighted mean, deg): median ", round(median(d$ang_deg), 3),
    " mean ", round(mean(d$ang_deg), 3), " max ", round(max(d$ang_deg), 3))
say("fraction of trials with mean angle > 30 deg: ", round(mean(d$ang_deg > 30), 4))
say("fraction of trials with onset angle > 30 deg: ", round(mean(d$angOnset_deg > 30), 4))
say("angle implied by 70 cm at the observed mean distance 0.50 m: ",
    round(2 * atan(0.35 / 0.5005) * 180 / pi, 1), " deg (two-point subtense), ",
    round(atan(0.70 / 0.5005) * 180 / pi, 1), " deg (one-sided)")
say("distance that would make 70 cm subtend 24 deg: ",
    round(0.70 / tan(24 * pi / 180), 3), " m (one-sided), ",
    round(0.35 / tan(12 * pi / 180), 3), " m (two-point subtense)")
say("median experienced angle among trials with disparity >= 60 cm: ",
    round(median(d$ang_deg[d$disp_cm >= 60]), 2))

# ---------------------------------------------------------------------------
rule("14. SOUND TYPE IN THE CAPTURE MODEL (OMITTED IN THE ORIGINAL)")
g_st <- glmmTMB(capture ~ ang10c + soundType + (1 + ang10c | participantId),
                data = d, family = binomial)
print(round(summary(g_st)$coefficients$cond, 5))
print(anova(g_ang, g_st))

# ---------------------------------------------------------------------------
rule("15. WRITE OUTPUTS")
write.csv(tabA, file.path(res_dir, "verify_tolerance_lookup_nominal.csv"), row.names = FALSE)
write.csv(tabB, file.path(res_dir, "verify_tolerance_lookup_angular.csv"), row.names = FALSE)
write.csv(cr_all, file.path(res_dir, "verify_capture_crossings.csv"), row.names = FALSE)
write.csv(cap_tab, file.path(res_dir, "verify_capture_by_listener.csv"), row.names = FALSE)
write.csv(rob, file.path(res_dir, "verify_capture_criterion_robustness.csv"), row.names = FALSE)
say("written to ", res_dir)

rule("SESSION INFO")
print(sessionInfo())
