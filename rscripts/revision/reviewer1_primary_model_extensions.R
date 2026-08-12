#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Reviewer 1 extensions to the primary Gaussian LMM on signed error.
#
# Baseline (published primary model, structurally unchanged):
#   signedError_m ~ stimulusDisparity_m + soundType + (1 + stimulusDisparity_m | participantId)
#
# Test 1: does sound type MODULATE the capture slope?
#         stimulusDisparity_m * soundType, Type III F test + ML likelihood-ratio test,
#         per-sound-type slopes (emtrends) with 95% CIs and pairwise slope contrasts.
#
# Test 2: does trial number matter, and does it interact with disparity?
#         stimulusDisparity_m * trialSequence_c + soundType, with an attempt at a
#         random slope for trial.
#
# Marginal and conditional R2 are computed by hand using the Nakagawa/Johnson
# formulation (MuMIn is not installed in this environment).
#
# Run with the framework R 4.4-arm64 build:
#   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/bin/Rscript \
#     rscripts/revision/reviewer1_primary_model_extensions.R
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(dplyr)
})

set.seed(20260805)

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir  <- file.path(proj_dir, "results_revision")
data_rds <- file.path(res_dir, "analysis_df_revision.rds")
log_path <- file.path(res_dir, "reviewer1_primary_model_extensions_log.txt")

dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

# Mirror all console output into the log file.
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

# Satterthwaite df throughout, matching the df method used for the published model.
emm_options(lmerTest.limit = 1000, pbkrtest.limit = 1000, lmer.df = "satterthwaite")

ctrl <- lmerControl(optimizer = "bobyqa")  # optimizer used in the original analysis

rule <- function(title) cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
sub  <- function(title) cat("\n--- ", title, " ---\n", sep = "")

# Fit a model and record any convergence warnings or messages verbatim.
fit_records <- list()
fit_lmm <- function(label, ...) {
  msgs <- character(0)
  m <- withCallingHandlers(
    tryCatch(lmer(...), error = function(e) {
      msgs <<- c(msgs, paste0("ERROR: ", conditionMessage(e)))
      NULL
    }),
    warning = function(w) {
      msgs <<- c(msgs, paste0("WARNING: ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    },
    message = function(m2) {
      msgs <<- c(msgs, paste0("MESSAGE: ", trimws(conditionMessage(m2))))
      invokeRestart("muffleMessage")
    }
  )
  if (length(msgs) == 0) msgs <- "none"
  fit_records[[label]] <<- msgs
  cat("[fit] ", label, ": ", paste(msgs, collapse = " | "), "\n", sep = "")
  m
}

# --- Nakagawa (2013) / Johnson (2014) R2 for LMMs, computed by hand -------
# var_fixed  = variance of the fixed-effect linear predictor across observations
# var_random = mean over observations of z_i' Sigma z_i, summed over grouping terms
#              (Johnson's extension, required because the model has a random slope)
# var_resid  = residual variance
r2_nakagawa <- function(model) {
  X <- model.matrix(model)
  var_fixed <- as.numeric(var(as.vector(X %*% fixef(model))))
  vc <- VarCorr(model)
  # mmList() is named by random-effects term, VarCorr() by grouping factor;
  # both follow the cnms ordering, so pair them positionally.
  zmats <- getME(model, "mmList")
  stopifnot(length(zmats) == length(vc))
  var_random <- 0
  for (i in seq_along(vc)) {
    Sigma <- as.matrix(vc[[i]])
    Z <- as.matrix(zmats[[i]])
    stopifnot(identical(colnames(Z), colnames(Sigma)))
    var_random <- var_random + mean(rowSums((Z %*% Sigma) * Z))
  }
  var_resid <- sigma(model)^2
  total <- var_fixed + var_random + var_resid
  c(var_fixed = var_fixed, var_random = var_random, var_resid = var_resid,
    R2m = var_fixed / total, R2c = (var_fixed + var_random) / total)
}

# Fixed-effect table with Satterthwaite df and Wald 95% CIs.
coef_table <- function(model, model_label) {
  cf <- as.data.frame(coef(summary(model)))
  names(cf) <- c("estimate", "SE", "df", "t", "p")
  ci <- confint(model, method = "Wald", parm = rownames(cf))
  data.frame(model = model_label, term = rownames(cf), cf,
             ci_low = ci[, 1], ci_high = ci[, 2],
             ci_width = ci[, 2] - ci[, 1], row.names = NULL)
}

re_table <- function(model, model_label) {
  vc <- as.data.frame(VarCorr(model))
  data.frame(model = model_label, group = vc$grp, term1 = vc$var1,
             term2 = ifelse(is.na(vc$var2), "", vc$var2),
             variance = vc$vcov, sd_or_cor = vc$sdcor, row.names = NULL)
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
rule("DATA")
dat <- readRDS(data_rds)
cat("rows:", nrow(dat), " participants:", dplyr::n_distinct(dat$participantId), "\n")
cat("soundType levels:", paste(levels(dat$soundType), collapse = ", "), "\n")
cat("trials per sound type:\n"); print(table(dat$soundType))
cat("trialSequenceNum range:", range(dat$trialSequenceNum), "\n")
cat("trialSequence_c = trialSequenceNum - ", mean(dat$trialSequenceNum), " (grand-mean centred)\n", sep = "")
cat("stimulusDisparity_m range:", range(dat$stimulusDisparity_m), "\n")

# ---------------------------------------------------------------------------
# Baseline: published primary model
# ---------------------------------------------------------------------------
rule("BASELINE: published primary model (M0)")
m0 <- fit_lmm("M0 baseline",
              signedError_m ~ stimulusDisparity_m + soundType +
                (1 + stimulusDisparity_m | participantId),
              data = dat, REML = TRUE, control = ctrl)
print(summary(m0))
sub("Type III F tests (Satterthwaite)")
print(anova(m0))
tab_m0 <- coef_table(m0, "M0 baseline")
sub("Fixed effects with Wald 95% CIs")
print(tab_m0, digits = 4)
r2_m0 <- r2_nakagawa(m0)
sub("Nakagawa/Johnson R2 (by hand)")
print(round(r2_m0, 5))

sub("Baseline check against published values")
b0 <- tab_m0[tab_m0$term == "stimulusDisparity_m", ]
check <- data.frame(
  quantity = c("beta_disparity", "SE", "df (Satterthwaite)", "t", "p", "CI low", "CI high"),
  published = c(0.347, 0.054, 30.5, 6.48, NA, 0.24, 0.45),
  recomputed = c(b0$estimate, b0$SE, b0$df, b0$t, b0$p, b0$ci_low, b0$ci_high))
check$abs_diff <- abs(check$recomputed - check$published)
print(check, digits = 6)

# Is the random slope justified?
rule("RANDOM-EFFECTS STRUCTURE: is the random slope justified?")
m0_ri <- fit_lmm("M0 random intercept only",
                 signedError_m ~ stimulusDisparity_m + soundType + (1 | participantId),
                 data = dat, REML = TRUE, control = ctrl)
sub("REML LRT: (1 | pid) vs (1 + disparity | pid)")
lrt_slope <- anova(m0_ri, m0, refit = FALSE)
print(lrt_slope)
sub("Random effects, M0")
print(re_table(m0, "M0 baseline"), digits = 4)

# ---------------------------------------------------------------------------
# TEST 1: disparity x sound type
# ---------------------------------------------------------------------------
rule("TEST 1: does sound type MODULATE the capture slope?")
m1 <- fit_lmm("M1 disparity x soundType",
              signedError_m ~ stimulusDisparity_m * soundType +
                (1 + stimulusDisparity_m | participantId),
              data = dat, REML = TRUE, control = ctrl)
print(summary(m1))

sub("Type III F tests, treatment contrasts (Satterthwaite)")
print(anova(m1, type = 3))

# Type III main effects are only marginal-interpretable under sum-to-zero coding;
# the interaction F is identical either way. Refit with contr.sum as a check.
dat_sum <- dat
contrasts(dat_sum$soundType) <- contr.sum(nlevels(dat_sum$soundType))
m1_sum <- fit_lmm("M1 with contr.sum",
                  signedError_m ~ stimulusDisparity_m * soundType +
                    (1 + stimulusDisparity_m | participantId),
                  data = dat_sum, REML = TRUE, control = ctrl)
sub("Type III F tests, sum-to-zero contrasts (Satterthwaite)")
print(anova(m1_sum, type = 3))

sub("Likelihood-ratio test, both models refitted with ML")
m0_ml <- fit_lmm("M0 (ML)",
                 signedError_m ~ stimulusDisparity_m + soundType +
                   (1 + stimulusDisparity_m | participantId),
                 data = dat, REML = FALSE, control = ctrl)
m1_ml <- fit_lmm("M1 (ML)",
                 signedError_m ~ stimulusDisparity_m * soundType +
                   (1 + stimulusDisparity_m | participantId),
                 data = dat, REML = FALSE, control = ctrl)
lrt_t1 <- anova(m0_ml, m1_ml)
print(lrt_t1)

sub("Per-sound-type capture slopes (emtrends, Satterthwaite, 95% CI)")
et <- emtrends(m1, ~ soundType, var = "stimulusDisparity_m")
et_df <- as.data.frame(summary(et, infer = c(TRUE, TRUE)))
names(et_df)[names(et_df) == "stimulusDisparity_m.trend"] <- "slope"
et_df$ci_width <- et_df$upper.CL - et_df$lower.CL
et_df$n_trials <- as.integer(table(dat$soundType)[as.character(et_df$soundType)])
print(et_df, digits = 4)

sub("Pairwise differences between slopes (Tukey-adjusted)")
pw <- as.data.frame(summary(pairs(et), infer = c(TRUE, TRUE)))
pw$ci_width <- pw$upper.CL - pw$lower.CL
print(pw, digits = 4)

sub("Pairwise differences between slopes (unadjusted, for CI width reporting)")
pw_none <- as.data.frame(summary(pairs(et, adjust = "none"), infer = c(TRUE, TRUE)))
pw_none$ci_width <- pw_none$upper.CL - pw_none$lower.CL
print(pw_none, digits = 4)

sub("Precision summary for the null")
cat("slope CI widths (m/m): ", paste(sprintf("%.3f", et_df$ci_width), collapse = ", "), "\n")
cat("mean slope CI half-width: ", sprintf("%.3f", mean(et_df$ci_width) / 2), "\n")
cat("largest observed slope difference: ",
    sprintf("%.3f", max(abs(pw_none$estimate))), "\n")
cat("widest pairwise-difference CI: ", sprintf("%.3f", max(pw_none$ci_width)),
    " (half-width ", sprintf("%.3f", max(pw_none$ci_width) / 2), ")\n", sep = "")
r2_m1 <- r2_nakagawa(m1)
sub("Nakagawa/Johnson R2, M1")
print(round(r2_m1, 5))
sub("Random effects, M1")
print(re_table(m1, "M1 disparity x soundType"), digits = 4)

# ---------------------------------------------------------------------------
# TEST 2: trial number and trial x disparity
# ---------------------------------------------------------------------------
rule("TEST 2: trial number and trial x disparity")

# First attempt: random slope for trial as well as disparity.
m2_rs <- fit_lmm("M2 with random slope for trial",
                 signedError_m ~ stimulusDisparity_m * trialSequence_c + soundType +
                   (1 + stimulusDisparity_m + trialSequence_c | participantId),
                 data = dat, REML = TRUE, control = ctrl)
if (!is.null(m2_rs)) {
  sub("M2 with random trial slope: convergence diagnostics")
  print(VarCorr(m2_rs))
  cat("singular fit: ", isSingular(m2_rs), "\n")
  cat("max |relative gradient|: ",
      max(abs(with(m2_rs@optinfo$derivs, solve(Hessian, gradient)))), "\n")
  sub("M2 with random trial slope: fixed effects")
  print(coef_table(m2_rs, "M2 random trial slope"), digits = 4)
  sub("M2 with random trial slope: Type III F tests")
  print(anova(m2_rs, type = 3))
  sub("Nakagawa/Johnson R2, M2 with random trial slope")
  print(round(r2_nakagawa(m2_rs), 5))
}

# Fallback / comparison: fixed effect for trial only.
m2 <- fit_lmm("M2 trial as fixed effect",
              signedError_m ~ stimulusDisparity_m * trialSequence_c + soundType +
                (1 + stimulusDisparity_m | participantId),
              data = dat, REML = TRUE, control = ctrl)
print(summary(m2))
sub("Type III F tests (Satterthwaite)")
print(anova(m2, type = 3))
tab_m2 <- coef_table(m2, "M2 trial")
sub("Fixed effects with Wald 95% CIs")
print(tab_m2, digits = 4)

sub("REML LRT: is the random slope for trial justified?")
lrt_trial_rs <- anova(m2, m2_rs, refit = FALSE)
print(lrt_trial_rs)

# stimulusDisparity_m is on its raw scale (0.15-0.70 m), so in the interaction model
# the trialSequence_c coefficient is the trial slope extrapolated to zero disparity,
# which is outside the stimulus range. Refit on centred disparity so that the trial
# main effect is evaluated at the mean disparity. The disparity slope, the
# interaction and the model fit are unchanged by this reparameterisation.
sub("Same model on centred disparity (trial main effect at mean disparity)")
cat("mean stimulusDisparity_m = ", sprintf("%.4f", mean(dat$stimulusDisparity_m)), " m\n", sep = "")
m2_c <- fit_lmm("M2 on centred disparity",
                signedError_m ~ stimDisparity_c * trialSequence_c + soundType +
                  (1 + stimDisparity_c | participantId),
                data = dat, REML = TRUE, control = ctrl)
print(coef_table(m2_c, "M2 centred disparity"), digits = 4)
sub("Type III F tests on centred disparity")
print(anova(m2_c, type = 3))

sub("Trial slope at representative disparities (emtrends)")
disp_levels <- c(min(dat$stimulusDisparity_m), mean(dat$stimulusDisparity_m),
                 max(dat$stimulusDisparity_m))
et_trial <- emtrends(m2, ~ stimulusDisparity_m, var = "trialSequence_c",
                     at = list(stimulusDisparity_m = disp_levels))
et_trial_df <- as.data.frame(summary(et_trial, infer = c(TRUE, TRUE)))
names(et_trial_df)[names(et_trial_df) == "trialSequence_c.trend"] <- "trial_slope_m_per_trial"
et_trial_df$change_over_23_trials_cm <- et_trial_df$trial_slope_m_per_trial * 23 * 100
print(et_trial_df, digits = 4)

sub("Likelihood-ratio tests, ML refits")
m2_ml <- fit_lmm("M2 (ML)",
                 signedError_m ~ stimulusDisparity_m * trialSequence_c + soundType +
                   (1 + stimulusDisparity_m | participantId),
                 data = dat, REML = FALSE, control = ctrl)
m2_main_ml <- fit_lmm("M2 main effects only (ML)",
                      signedError_m ~ stimulusDisparity_m + trialSequence_c + soundType +
                        (1 + stimulusDisparity_m | participantId),
                      data = dat, REML = FALSE, control = ctrl)
cat("\nM0 (no trial) vs M2main (trial main effect added):\n")
lrt_trial_main <- anova(m0_ml, m2_main_ml)
print(lrt_trial_main)
cat("\nM2main vs M2 (disparity x trial interaction added):\n")
lrt_trial_int <- anova(m2_main_ml, m2_ml)
print(lrt_trial_int)

sub("Does adding trial change the disparity slope?")
b2 <- tab_m2[tab_m2$term == "stimulusDisparity_m", ]
slope_cmp <- data.frame(
  model = c("M0 (no trial)", "M2 (trial + trial x disparity)"),
  beta_disparity = c(b0$estimate, b2$estimate),
  SE = c(b0$SE, b2$SE), df = c(b0$df, b2$df), t = c(b0$t, b2$t),
  p = c(b0$p, b2$p), ci_low = c(b0$ci_low, b2$ci_low), ci_high = c(b0$ci_high, b2$ci_high))
slope_cmp$change_vs_M0 <- slope_cmp$beta_disparity - b0$estimate
slope_cmp$pct_change <- 100 * slope_cmp$change_vs_M0 / b0$estimate
print(slope_cmp, digits = 5)

sub("Trial effect on the response scale (at mean disparity)")
tab_m2c <- coef_table(m2_c, "M2 centred disparity")
bt <- tab_m2c[tab_m2c$term == "trialSequence_c", ]
span <- diff(range(dat$trialSequenceNum))
cat("trial slope at mean disparity: ", sprintf("%.6f", bt$estimate), " m/trial (",
    sprintf("%.4f", bt$estimate * 100), " cm/trial), 95% CI [",
    sprintf("%.6f", bt$ci_low), ", ", sprintf("%.6f", bt$ci_high), "]\n", sep = "")
cat("change across the ", span, "-trial span: ",
    sprintf("%.5f", bt$estimate * span), " m (",
    sprintf("%.3f", bt$estimate * span * 100), " cm), 95% CI [",
    sprintf("%.3f", bt$ci_low * span * 100), ", ",
    sprintf("%.3f", bt$ci_high * span * 100), "] cm\n", sep = "")
cat("mean signed error: ", sprintf("%.4f", mean(dat$signedError_m)), " m (",
    sprintf("%.2f", mean(dat$signedError_m) * 100), " cm)\n", sep = "")
cat("change as a percentage of the mean signed error: ",
    sprintf("%.2f", 100 * bt$estimate * span / mean(dat$signedError_m)), "%\n", sep = "")

# Supplementary: the reviewer's "19-20% learning effect" refers to unsigned error,
# so quantify it on that measure for completeness. Not part of the primary model.
sub("Supplementary: trial effect on UNSIGNED error (context for the reviewer's 19-20% claim)")
m_unsigned <- fit_lmm("Supplementary unsigned-error trial model",
                      participantError_m ~ trialSequence_c + soundType +
                        (1 | participantId),
                      data = dat, REML = TRUE, control = ctrl)
print(coef_table(m_unsigned, "unsigned error ~ trial"), digits = 4)
bu <- coef_table(m_unsigned, "unsigned error ~ trial")
bu <- bu[bu$term == "trialSequence_c", ]
pred_first <- fixef(m_unsigned)[["(Intercept)"]] + bu$estimate * (1 - mean(dat$trialSequenceNum))
pred_last  <- fixef(m_unsigned)[["(Intercept)"]] + bu$estimate * (24 - mean(dat$trialSequenceNum))
cat("model-predicted unsigned error, trial 1: ", sprintf("%.2f", pred_first * 100), " cm\n", sep = "")
cat("model-predicted unsigned error, trial 24: ", sprintf("%.2f", pred_last * 100), " cm\n", sep = "")
cat("relative change over 24 trials: ",
    sprintf("%.1f", 100 * (pred_last - pred_first) / pred_first), "%\n", sep = "")
obs_tr <- dat %>% group_by(block = ifelse(trialSequenceNum <= 12, "trials 1-12", "trials 13-24")) %>%
  summarise(mean_unsigned_cm = 100 * mean(participantError_m),
            mean_signed_cm = 100 * mean(signedError_m), n = dplyr::n(), .groups = "drop")
print(as.data.frame(obs_tr), digits = 4)

# The size of the descriptive "learning effect" depends strongly on how early and late
# trials are windowed, so report the sweep rather than a single figure.
sub("Learning effect on unsigned error under different windowings")
u_cm <- dat$participantError_m * 100
tn <- dat$trialSequenceNum
windows <- list("trials 1-12 vs 13-24" = list(tn <= 12, tn > 12),
                "trials 1-6 vs 19-24"  = list(tn <= 6,  tn >= 19),
                "trials 1-4 vs 21-24"  = list(tn <= 4,  tn >= 21),
                "trial 1 vs trial 24"  = list(tn == 1,  tn == 24))
learn <- do.call(rbind, lapply(names(windows), function(nm) {
  w <- windows[[nm]]
  a <- mean(u_cm[w[[1]]]); b <- mean(u_cm[w[[2]]])
  data.frame(window = nm, early_cm = a, late_cm = b, pct_change = 100 * (b - a) / a)
}))
learn <- rbind(learn, data.frame(window = "LMM-predicted trial 1 vs 24",
                                 early_cm = pred_first * 100, late_cm = pred_last * 100,
                                 pct_change = 100 * (pred_last - pred_first) / pred_first))
print(learn, digits = 4)
write.csv(learn, file.path(res_dir, "reviewer1_learning_effect_windows.csv"), row.names = FALSE)
r2_m2 <- r2_nakagawa(m2)
sub("Nakagawa/Johnson R2, M2")
print(round(r2_m2, 5))
sub("Random effects, M2")
print(re_table(m2, "M2 trial"), digits = 4)

# Random slope still justified in M2?
m2_ri <- fit_lmm("M2 random intercept only",
                 signedError_m ~ stimulusDisparity_m * trialSequence_c + soundType +
                   (1 | participantId),
                 data = dat, REML = TRUE, control = ctrl)
sub("REML LRT in M2: (1 | pid) vs (1 + disparity | pid)")
lrt_slope_m2 <- anova(m2_ri, m2, refit = FALSE)
print(lrt_slope_m2)

# ---------------------------------------------------------------------------
# Machine-readable output
# ---------------------------------------------------------------------------
rule("WRITING CSV OUTPUT")

fixed_all <- rbind(tab_m0, coef_table(m1, "M1 disparity x soundType"), tab_m2, tab_m2c,
                   coef_table(m2_rs, "M2 random trial slope"))
write.csv(fixed_all, file.path(res_dir, "reviewer1_fixed_effects.csv"), row.names = FALSE)

anova_to_df <- function(a, label) {
  d <- as.data.frame(a)
  data.frame(model = label, term = rownames(d), d, row.names = NULL)
}
anova_all <- rbind(anova_to_df(anova(m0), "M0 baseline"),
                   anova_to_df(anova(m1, type = 3), "M1 disparity x soundType"),
                   anova_to_df(anova(m2, type = 3), "M2 trial"),
                   anova_to_df(anova(m2_c, type = 3), "M2 centred disparity"),
                   anova_to_df(anova(m2_rs, type = 3), "M2 random trial slope"))
write.csv(anova_all, file.path(res_dir, "reviewer1_type3_ftests.csv"), row.names = FALSE)

lrt_to_df <- function(a, label) {
  d <- as.data.frame(a)
  data.frame(comparison = label, model = rownames(d), d, row.names = NULL)
}
lrt_all <- rbind(
  lrt_to_df(lrt_slope, "random slope for disparity (REML, M0)"),
  lrt_to_df(lrt_t1, "disparity x soundType (ML)"),
  lrt_to_df(lrt_trial_main, "trial main effect (ML)"),
  lrt_to_df(lrt_trial_int, "disparity x trial (ML)"),
  lrt_to_df(lrt_slope_m2, "random slope for disparity (REML, M2)"),
  lrt_to_df(lrt_trial_rs, "random slope for trial (REML, M2)"))
write.csv(lrt_all, file.path(res_dir, "reviewer1_likelihood_ratio_tests.csv"), row.names = FALSE)

write.csv(et_df, file.path(res_dir, "reviewer1_slopes_by_soundtype.csv"), row.names = FALSE)
write.csv(et_trial_df, file.path(res_dir, "reviewer1_trial_slopes_by_disparity.csv"),
          row.names = FALSE)
write.csv(rbind(cbind(adjust = "tukey", pw), cbind(adjust = "none", pw_none)),
          file.path(res_dir, "reviewer1_slope_pairwise_contrasts.csv"), row.names = FALSE)

r2_all <- rbind(data.frame(model = "M0 baseline", t(r2_m0)),
                data.frame(model = "M1 disparity x soundType", t(r2_m1)),
                data.frame(model = "M2 trial", t(r2_m2)),
                data.frame(model = "M2 random trial slope", t(r2_nakagawa(m2_rs))))
write.csv(r2_all, file.path(res_dir, "reviewer1_r2_nakagawa.csv"), row.names = FALSE)

re_all <- rbind(re_table(m0, "M0 baseline"), re_table(m1, "M1 disparity x soundType"),
                re_table(m2, "M2 trial"))
if (!is.null(m2_rs)) re_all <- rbind(re_all, re_table(m2_rs, "M2 with random trial slope"))
write.csv(re_all, file.path(res_dir, "reviewer1_random_effects.csv"), row.names = FALSE)

write.csv(data.frame(model = names(fit_records),
                     messages = vapply(fit_records, paste, character(1), collapse = " | ")),
          file.path(res_dir, "reviewer1_convergence.csv"), row.names = FALSE)

write.csv(check, file.path(res_dir, "reviewer1_baseline_check.csv"), row.names = FALSE)

cat("\nCSV files written to ", res_dir, "\n", sep = "")
rule("SESSION INFO")
print(sessionInfo())
cat("\nDone.\n")
