#!/usr/bin/env Rscript
# Independent verification of the Reviewer 1 primary-model extensions.
#
# Refits every model reported in rscripts/revision/reviewer1_primary_model_extensions.R
# from the prepared trial-level dataset, without reusing that script's code, and
# additionally probes:
#   - optimiser sensitivity and conditioning of the fixed-effect vcov,
#   - the Johnson (2014) random-slope R2 via an independent matrix-trace route,
#   - whether the reported learning effect on unsigned error survives adjustment
#     for stimulus disparity (which the supplementary model omits),
#   - whether the ML LRT cited as "confirming" the centred parameterisation is
#     in fact parameterisation-invariant.
#
# Deterministic; set.seed only matters for the parametric bootstrap section.

library(lme4)
library(lmerTest)
library(emmeans)

set.seed(20260805)

root <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res  <- file.path(root, "results_revision")
d    <- readRDS(file.path(res, "analysis_df_revision.rds"))

log_path <- file.path(res, "verify_reviewer1_primary_model_extensions_log.txt")
con <- file(log_path, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

ctrl <- lmerControl(optimizer = "bobyqa")

rule("0. DATA INTEGRITY")
cat("n trials:", nrow(d), " n participants:", length(unique(d$participantId)), "\n")
cat("trials per participant:", paste(unique(table(d$participantId)), collapse = ","), "\n")
cat("disparity range:", paste(range(d$stimulusDisparity_m), collapse = " - "), " mean:",
    mean(d$stimulusDisparity_m), "\n")
cat("trialSequence_c is grand-mean centred:",
    isTRUE(all.equal(d$trialSequence_c, d$trialSequenceNum - mean(d$trialSequenceNum))), "\n")
cat("stimDisparity_c is grand-mean centred:",
    isTRUE(all.equal(d$stimDisparity_c, d$stimulusDisparity_m - mean(d$stimulusDisparity_m))), "\n")
cat("cor(trialSequenceNum, stimulusDisparity_m):",
    cor(d$trialSequenceNum, d$stimulusDisparity_m), "\n")
print(table(d$soundType, d$trialSequenceNum <= 12))

# ---------------------------------------------------------------- helpers ----

# Nakagawa/Johnson R2, implemented two independent ways.
r2_johnson <- function(m) {
  X    <- getME(m, "X")
  beta <- fixef(m)
  vf   <- as.numeric(var(as.vector(X %*% beta)))
  vres <- sigma(m)^2

  # route A: mean over observations of z_i' Sigma z_i, per grouping term
  vc   <- VarCorr(m)
  mmL  <- getME(m, "mmList")
  vrA  <- 0
  for (k in seq_along(vc)) {
    Sig <- as.matrix(vc[[k]])
    Zg  <- as.matrix(mmL[[k]])
    vrA <- vrA + mean(rowSums((Zg %*% Sig) * Zg))
  }

  # route B: trace form, (1/n) * sigma^2 * tr(Z Lambda Lambda' Z'), using lme4's
  # full sparse Z and relative covariance factor. Independent of route A.
  Z   <- getME(m, "Z")
  Lam <- Matrix::t(getME(m, "Lambdat"))   # Lambda, q x q
  M   <- Z %*% Lam                        # n x q
  vrB <- sum(M^2) / nrow(X) * vres

  c(var_fixed = vf, var_rand_A = vrA, var_rand_B = vrB, var_resid = vres,
    R2m_A = vf / (vf + vrA + vres), R2c_A = (vf + vrA) / (vf + vrA + vres),
    R2m_B = vf / (vf + vrB + vres), R2c_B = (vf + vrB) / (vf + vrB + vres))
}

conv_info <- function(m, label) {
  g <- m@optinfo$derivs
  mg <- tryCatch(max(abs(solve(g$Hessian, g$gradient))), error = function(e) NA_real_)
  V <- as.matrix(vcov(m))
  ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  th <- as.matrix(VarCorr(m)[[1]])
  thev <- eigen(th, symmetric = TRUE, only.values = TRUE)$values
  data.frame(model = label,
             singular = isSingular(m),
             max_abs_rel_grad = mg,
             vcov_min_eig = min(ev),
             vcov_cond = max(ev) / min(ev),
             RE_min_eig = min(thev),
             RE_cond = max(thev) / min(thev),
             n_warn = length(m@optinfo$conv$lme4$messages),
             warns = paste(m@optinfo$conv$lme4$messages, collapse = " | "))
}

fx <- function(m, term) {
  s <- summary(m)$coefficients
  s[term, , drop = FALSE]
}

# ------------------------------------------------------- M0 reproduction -----

rule("1. M0 BASELINE REPRODUCTION")
m0 <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
             (1 + stimulusDisparity_m | participantId),
           data = d, REML = TRUE, control = ctrl)
print(summary(m0)$coefficients)
cat("\nWald 95% CI for disparity:\n")
b <- fixef(m0)["stimulusDisparity_m"]; se <- sqrt(vcov(m0)["stimulusDisparity_m", "stimulusDisparity_m"])
cat(sprintf("  b = %.6f  SE = %.6f  CI = [%.4f, %.4f]\n", b, se, b - 1.96 * se, b + 1.96 * se))
print(VarCorr(m0))
cat("residual SD:", sigma(m0), "\n")
cat("AIC(REML):", AIC(m0), "\n")

rule("1b. OPTIMISER SENSITIVITY / VCOV CONDITIONING FOR M0")
opts <- list(
  bobyqa      = lmerControl(optimizer = "bobyqa"),
  NelderMead  = lmerControl(optimizer = "Nelder_Mead"),
  nloptwrap   = lmerControl(optimizer = "nloptwrap"),
  nlminbwrap  = lmerControl(optimizer = "nlminbwrap")
)
opt_tab <- do.call(rbind, lapply(names(opts), function(nm) {
  fit <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                (1 + stimulusDisparity_m | participantId),
              data = d, REML = TRUE, control = opts[[nm]])
  s <- summary(fit)$coefficients
  vcs <- as.data.frame(VarCorr(fit))
  data.frame(optimizer = nm,
             logLik = as.numeric(logLik(fit)),
             b_disp = s["stimulusDisparity_m", "Estimate"],
             se_disp = s["stimulusDisparity_m", "Std. Error"],
             t_disp = s["stimulusDisparity_m", "t value"],
             df_disp = s["stimulusDisparity_m", "df"],
             sd_int = vcs$sdcor[1], sd_slope = vcs$sdcor[2],
             corr = vcs$sdcor[3], sigma = sigma(fit),
             singular = isSingular(fit))
}))
print(opt_tab, row.names = FALSE)
cat("\nmax spread of b_disp across optimisers:", diff(range(opt_tab$b_disp)), "\n")
cat("max spread of se_disp across optimisers:", diff(range(opt_tab$se_disp)), "\n")

cat("\nM0 fixed-effect vcov eigenvalues:\n")
print(eigen(as.matrix(vcov(m0)), symmetric = TRUE, only.values = TRUE)$values)
cat("M0 random-effect covariance eigenvalues:\n")
print(eigen(as.matrix(VarCorr(m0)$participantId), symmetric = TRUE, only.values = TRUE)$values)

cat("\nProfile CI for disparity and for the RE parameters (cross-check on Wald):\n")
prof_ci <- tryCatch(confint(m0, method = "profile", oldNames = FALSE),
                    error = function(e) paste("profile failed:", conditionMessage(e)))
print(prof_ci)

rule("1c. M0 RANDOM SLOPE JUSTIFICATION + R2")
m0_ri <- lmer(signedError_m ~ stimulusDisparity_m + soundType + (1 | participantId),
              data = d, REML = TRUE, control = ctrl)
print(anova(m0_ri, m0, refit = FALSE))
cat("\nR2 (two independent routes):\n")
print(round(r2_johnson(m0), 6))
cat("\nR2 helper sanity check on the random-intercept model (var_rand should equal tau00 =",
    as.numeric(VarCorr(m0_ri)$participantId), "):\n")
print(round(r2_johnson(m0_ri), 6))

# ------------------------------------------------------------ TEST 1 ---------

rule("2. TEST 1: disparity x soundType")
m1 <- lmer(signedError_m ~ stimulusDisparity_m * soundType +
             (1 + stimulusDisparity_m | participantId),
           data = d, REML = TRUE, control = ctrl)
cat("Type III (lmerTest anova, Satterthwaite):\n")
print(anova(m1, type = 3))

d_sum <- d
contrasts(d_sum$soundType) <- contr.sum(4)
m1_sum <- lmer(signedError_m ~ stimulusDisparity_m * soundType +
                 (1 + stimulusDisparity_m | participantId),
               data = d_sum, REML = TRUE, control = ctrl)
cat("\nType III under contr.sum:\n")
print(anova(m1_sum, type = 3))

m0_ml <- update(m0, REML = FALSE)
m1_ml <- update(m1, REML = FALSE)
cat("\nML LRT M0 vs M1:\n")
print(anova(m0_ml, m1_ml))

cat("\nemmeans default df class for this object:\n")
cat("pbkrtest installed:", requireNamespace("pbkrtest", quietly = TRUE), "\n")
et_default <- emtrends(m1, ~ soundType, var = "stimulusDisparity_m")
cat("emtrends WITHOUT setting lmer.df (i.e. emmeans default):\n")
print(summary(et_default, infer = c(TRUE, TRUE)))

emm_options(lmer.df = "satterthwaite")
et <- emtrends(m1, ~ soundType, var = "stimulusDisparity_m")
cat("\nemtrends with Satterthwaite:\n")
et_s <- summary(et, infer = c(TRUE, TRUE))
print(et_s)
cat("\nCI widths:", paste(round(et_s$upper.CL - et_s$lower.CL, 4), collapse = ", "), "\n")
cat("mean half-width:", mean((et_s$upper.CL - et_s$lower.CL) / 2), "\n")

cat("\nPairwise slope contrasts (Tukey):\n")
pw <- summary(pairs(et), infer = c(TRUE, TRUE))
print(pw)
cat("widest Tukey CI:", max(pw$upper.CL - pw$lower.CL), "\n")
cat("\nPairwise unadjusted:\n")
print(summary(pairs(et, adjust = "none"), infer = c(TRUE, TRUE)))

cat("\nM1 R2:\n"); print(round(r2_johnson(m1), 6))
print(conv_info(m1, "M1"))

# ------------------------------------------------------------ TEST 2 ---------

rule("3. TEST 2: trial and trial x disparity")
m2_raw <- lmer(signedError_m ~ stimulusDisparity_m * trialSequence_c + soundType +
                 (1 + stimulusDisparity_m | participantId),
               data = d, REML = TRUE, control = ctrl)
cat("RAW-disparity parameterisation:\n")
print(summary(m2_raw)$coefficients)
print(anova(m2_raw, type = 3))

m2_c <- lmer(signedError_m ~ stimDisparity_c * trialSequence_c + soundType +
               (1 + stimDisparity_c | participantId),
             data = d, REML = TRUE, control = ctrl)
cat("\nCENTRED-disparity parameterisation:\n")
print(summary(m2_c)$coefficients)
print(anova(m2_c, type = 3))
cat("\nlogLik raw vs centred (must be identical - same model):",
    as.numeric(logLik(m2_raw)), as.numeric(logLik(m2_c)), "\n")

cat("\nTrial main effect at mean disparity, converted to cm over the 23-trial span:\n")
sc <- summary(m2_c)$coefficients
bt <- sc["trialSequence_c", "Estimate"]; st <- sc["trialSequence_c", "Std. Error"]
dfl <- sc["trialSequence_c", "df"]; tq <- qt(0.975, dfl)
cat(sprintf("  b = %.7f m/trial, SE = %.7f, t(%.2f) = %.3f, p = %.4f\n",
            bt, st, dfl, sc["trialSequence_c", "t value"], sc["trialSequence_c", "Pr(>|t|)"]))
cat(sprintf("  per-trial CI (cm): [%.4f, %.4f]\n", 100 * (bt - tq * st), 100 * (bt + tq * st)))
cat(sprintf("  over 23 trials (cm): %.3f, CI [%.3f, %.3f]\n",
            2300 * bt, 2300 * (bt - tq * st), 2300 * (bt + tq * st)))

cat("\nTrial slope at the disparity extremes:\n")
print(summary(emtrends(m2_raw, ~ stimulusDisparity_m, var = "trialSequence_c",
                       at = list(stimulusDisparity_m = c(0.15, 0.70))),
              infer = c(TRUE, TRUE)))

cat("\nCapture slope with trial in the model (raw and centred must agree):\n")
print(fx(m2_raw, "stimulusDisparity_m"))
print(fx(m2_c, "stimDisparity_c"))
b2 <- fixef(m2_raw)["stimulusDisparity_m"]
se2 <- sqrt(vcov(m2_raw)["stimulusDisparity_m", "stimulusDisparity_m"])
cat(sprintf("  CI [%.4f, %.4f]; change from M0: %.7f (%.4f%%)\n",
            b2 - 1.96 * se2, b2 + 1.96 * se2, b2 - b, 100 * (b2 - b) / b))

rule("3b. IS THE ML LRT PARAMETERISATION-INVARIANT?")
# The report claims the ML LRT "confirms the centred parameterisation is correct".
# Test that claim: run the trial-main-effect LRT under BOTH parameterisations,
# both against M0 (no trial at all) and against the interaction model with the
# trial main effect dropped.
m2raw_ml  <- update(m2_raw, REML = FALSE)
m2c_ml    <- update(m2_c,   REML = FALSE)
m_trial_raw_ml <- lmer(signedError_m ~ stimulusDisparity_m + trialSequence_c + soundType +
                         (1 + stimulusDisparity_m | participantId),
                       data = d, REML = FALSE, control = ctrl)
cat("M0(ML) vs main-effects+trial(ML):\n"); print(anova(m0_ml, m_trial_raw_ml))
cat("\nlogLik of the two interaction models under ML (identical?):",
    as.numeric(logLik(m2raw_ml)), as.numeric(logLik(m2c_ml)), "\n")
# drop only the trial MAIN effect from the interaction model, under each coding
m2raw_notrial <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                        stimulusDisparity_m:trialSequence_c +
                        (1 + stimulusDisparity_m | participantId),
                      data = d, REML = FALSE, control = ctrl)
m2c_notrial   <- lmer(signedError_m ~ stimDisparity_c + soundType +
                        stimDisparity_c:trialSequence_c +
                        (1 + stimDisparity_c | participantId),
                      data = d, REML = FALSE, control = ctrl)
cat("\nLRT dropping ONLY the trial main effect, RAW coding:\n")
print(anova(m2raw_notrial, m2raw_ml))
cat("\nLRT dropping ONLY the trial main effect, CENTRED coding:\n")
print(anova(m2c_notrial, m2c_ml))

rule("3c. RANDOM SLOPE FOR TRIAL")
m2_rs <- lmer(signedError_m ~ stimulusDisparity_m * trialSequence_c + soundType +
                (1 + stimulusDisparity_m + trialSequence_c | participantId),
              data = d, REML = TRUE, control = ctrl)
print(VarCorr(m2_rs))
print(conv_info(m2_rs, "M2rs"))
cat("\nREML LRT M2 vs M2rs:\n"); print(anova(m2_raw, m2_rs, refit = FALSE))
cat("\nM2rs capture slope:\n"); print(fx(m2_rs, "stimulusDisparity_m"))
cat("M2rs interaction:\n"); print(fx(m2_rs, "stimulusDisparity_m:trialSequence_c"))

cat("\nM2 R2:\n"); print(round(r2_johnson(m2_raw), 6))
cat("M2rs R2:\n"); print(round(r2_johnson(m2_rs), 6))

m2_ri <- lmer(signedError_m ~ stimulusDisparity_m * trialSequence_c + soundType +
                (1 | participantId), data = d, REML = TRUE, control = ctrl)
cat("\nREML LRT: disparity random slope still justified in M2:\n")
print(anova(m2_ri, m2_raw, refit = FALSE))

# ------------------------------------------------ supplementary learning -----

rule("4. SUPPLEMENTARY: LEARNING ON UNSIGNED ERROR")
sup <- lmer(participantError_m ~ trialSequence_c + soundType + (1 | participantId),
            data = d, REML = TRUE, control = ctrl)
print(summary(sup)$coefficients)
ss <- summary(sup)$coefficients
bs <- ss["trialSequence_c", "Estimate"]; ses <- ss["trialSequence_c", "Std. Error"]
cat(sprintf("Wald CI: [%.6f, %.6f]\n", bs - 1.96 * ses, bs + 1.96 * ses))

cat("\nPredicted unsigned error at trial 1 and 24, AVERAGED over sound type (emmeans):\n")
pr <- summary(emmeans(sup, ~ trialSequence_c, at = list(trialSequence_c = c(-11.5, 11.5))))
print(pr)
cat(sprintf("  trial1 = %.2f cm, trial24 = %.2f cm, change = %.2f%%\n",
            100 * pr$emmean[1], 100 * pr$emmean[2],
            100 * (pr$emmean[2] - pr$emmean[1]) / pr$emmean[1]))
cat("\nPrediction using the raw intercept (i.e. the Drum reference level only):\n")
ic <- fixef(sup)["(Intercept)"]
cat(sprintf("  trial1 = %.2f cm, trial24 = %.2f cm, change = %.2f%%\n",
            100 * (ic - 11.5 * bs), 100 * (ic + 11.5 * bs),
            100 * (23 * bs) / (ic - 11.5 * bs)))

rule("4b. IS THE LEARNING EFFECT CONFOUNDED WITH STIMULUS DISPARITY?")
# Unsigned error depends strongly on disparity; the supplementary model omits it.
sup_adj <- lmer(participantError_m ~ trialSequence_c + stimulusDisparity_m + soundType +
                  (1 | participantId), data = d, REML = TRUE, control = ctrl)
print(summary(sup_adj)$coefficients)
cat("\nAlso with the by-participant disparity slope, matching the primary RE structure:\n")
sup_adj2 <- lmer(participantError_m ~ trialSequence_c + stimulusDisparity_m + soundType +
                   (1 + stimulusDisparity_m | participantId), data = d, REML = TRUE, control = ctrl)
print(summary(sup_adj2)$coefficients)
cat("\nBy-participant mean disparity across the session (trend?):\n")
print(summary(lm(stimulusDisparity_m ~ trialSequenceNum, data = d))$coefficients)

cat("\nDistributional check on the supplementary Gaussian LMM:\n")
cat("  skewness of unsigned error:",
    mean((d$participantError_m - mean(d$participantError_m))^3) /
      sd(d$participantError_m)^3, "\n")
cat("  min unsigned error (m):", min(d$participantError_m), "\n")
r <- residuals(sup)
cat("  residual skewness:", mean((r - mean(r))^3) / sd(r)^3, "\n")
cat("  Shapiro-Wilk on residuals:\n"); print(shapiro.test(r))
cat("\nLog-transformed refit (robustness of the learning claim):\n")
sup_log <- lmer(log(participantError_m) ~ trialSequence_c + soundType + (1 | participantId),
                data = d, REML = TRUE, control = ctrl)
print(summary(sup_log)$coefficients)

rule("5. CONVERGENCE SUMMARY ACROSS ALL FITS")
all_fits <- list(M0 = m0, M0_ri = m0_ri, M1 = m1, M1_sum = m1_sum,
                 M2_raw = m2_raw, M2_c = m2_c, M2_rs = m2_rs, M2_ri = m2_ri,
                 Sup = sup, Sup_adj = sup_adj, Sup_adj_rs = sup_adj2, Sup_log = sup_log)
conv <- do.call(rbind, lapply(names(all_fits), function(n) conv_info(all_fits[[n]], n)))
print(conv, row.names = FALSE)
write.csv(conv, file.path(res, "verify_reviewer1_convergence.csv"), row.names = FALSE)

rule("6. HEADLINE NUMBER COMPARISON TABLE")
cmp <- rbind(
  data.frame(quantity = "M0 b_disparity",       reported = 0.34660,  recomputed = unname(fixef(m0)["stimulusDisparity_m"])),
  data.frame(quantity = "M0 SE_disparity",      reported = 0.05351,  recomputed = unname(se)),
  data.frame(quantity = "M0 t",                 reported = 6.4767,   recomputed = unname(summary(m0)$coefficients["stimulusDisparity_m", "t value"])),
  data.frame(quantity = "M0 df",                reported = 30.505,   recomputed = unname(summary(m0)$coefficients["stimulusDisparity_m", "df"])),
  data.frame(quantity = "M0 SD_int",            reported = 0.04449,  recomputed = as.data.frame(VarCorr(m0))$sdcor[1]),
  data.frame(quantity = "M0 SD_slope",          reported = 0.19015,  recomputed = as.data.frame(VarCorr(m0))$sdcor[2]),
  data.frame(quantity = "M0 corr",              reported = -0.556,   recomputed = as.data.frame(VarCorr(m0))$sdcor[3]),
  data.frame(quantity = "M0 sigma",             reported = 0.16412,  recomputed = sigma(m0)),
  data.frame(quantity = "M0 R2m",               reported = 0.08173,  recomputed = unname(r2_johnson(m0)["R2m_A"])),
  data.frame(quantity = "M0 R2c",               reported = 0.21080,  recomputed = unname(r2_johnson(m0)["R2c_A"])),
  data.frame(quantity = "T1 F",                 reported = 0.8059,   recomputed = anova(m1, type = 3)["stimulusDisparity_m:soundType", "F value"]),
  data.frame(quantity = "T1 LRT chisq",         reported = 2.437,    recomputed = anova(m0_ml, m1_ml)$Chisq[2]),
  data.frame(quantity = "T2c b_trial",          reported = 0.0001369, recomputed = unname(bt)),
  data.frame(quantity = "T2 interaction b",     reported = 0.008868, recomputed = unname(fixef(m2_raw)["stimulusDisparity_m:trialSequence_c"])),
  data.frame(quantity = "T2 capture slope",     reported = 0.34451,  recomputed = unname(b2)),
  data.frame(quantity = "Sup b_trial",          reported = -0.002094, recomputed = unname(bs))
)
cmp$abs_diff <- abs(cmp$reported - cmp$recomputed)
cmp$agrees <- cmp$abs_diff < pmax(1e-4 * abs(cmp$reported), 5e-5)
print(cmp, row.names = FALSE)
write.csv(cmp, file.path(res, "verify_reviewer1_headline_comparison.csv"), row.names = FALSE)

cat("\nR session:\n"); print(R.version.string)
cat("lme4", as.character(packageVersion("lme4")),
    "lmerTest", as.character(packageVersion("lmerTest")),
    "emmeans", as.character(packageVersion("emmeans")), "\n")

sink(type = "message"); sink(); close(con)
