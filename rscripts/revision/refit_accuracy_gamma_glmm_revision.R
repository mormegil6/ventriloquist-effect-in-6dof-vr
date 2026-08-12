# Re-estimation of the Analysis 2 accuracy models (localisation error, cm).
#
# The originally reported Gamma GLMM did not converge and returned a degenerate
# variance-covariance matrix (every fixed effect received an SE of about 0.0035).
# This script refits the same model specification with glmmTMB, reproduces the
# lme4::glmer pathology for the record, runs Tukey-adjusted sound-type contrasts,
# checks residuals with DHARMa, and repeats the fit under a lognormal alternative.
#
# Model: participantError_cm ~ stimDisparity_c + soundType + trialSequence_c +
#                              azimuthSector + elevationCategory + (1 | participantId)
#
# Input : results_revision/analysis_df_revision.rds  (built by build_analysis_df_revision.R)
# Output: results_revision/accuracy_gamma_glmm_*.csv and a plain-text log.

# NOTE (2026-08-07): glmmTMB::lognormal() is BROKEN in this R installation. It reports a
# log-likelihood it cannot attain and returns a nonsensical dispersion (sigma = 16.16 against
# sd(log y) = 0.63). Any lognormal() output below is retained only to document the defect and
# must not be reported. The correct log-normal model is a Gaussian fit on log(y); see
# refit_signed_error_lmm_revision.R and the audit in results_revision/.

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(emmeans)
  library(DHARMa)
})

set.seed(20260805)

proj_dir    <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
results_dir <- file.path(proj_dir, "results_revision")
log_path    <- file.path(results_dir, "accuracy_gamma_glmm_refit_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

cat("Re-estimation of Analysis 2 accuracy GLMMs\n")
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "\n")
cat("glmmTMB:", as.character(packageVersion("glmmTMB")),
    "| lme4:", as.character(packageVersion("lme4")),
    "| emmeans:", as.character(packageVersion("emmeans")),
    "| DHARMa:", as.character(packageVersion("DHARMa")), "\n\n")

# ---------------------------------------------------------------- data -------

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds"))

model_vars <- c("participantError_cm", "stimDisparity_c", "soundType",
                "trialSequence_c", "azimuthSector", "elevationCategory",
                "participantId")
dat <- df[stats::complete.cases(df[, model_vars]), model_vars]
dat$participantId <- factor(dat$participantId)

cat("Trials:", nrow(dat), "| participants:", nlevels(dat$participantId), "\n")
cat("Reference levels: soundType =", levels(dat$soundType)[1],
    ", azimuthSector =", levels(dat$azimuthSector)[1],
    ", elevationCategory =", levels(dat$elevationCategory)[1], "\n")
cat("Error (cm): min", round(min(dat$participantError_cm), 4),
    "mean", round(mean(dat$participantError_cm), 4),
    "max", round(max(dat$participantError_cm), 4), "\n\n")

model_formula <- participantError_cm ~ stimDisparity_c + soundType +
  trialSequence_c + azimuthSector + elevationCategory + (1 | participantId)

# --------------------------------------------------- (a) lme4 pathology ------

cat("=== (a) lme4::glmer refits, to document the original pathology ===\n")

# Helper: fit with glmer under a given optimiser, capturing warnings verbatim.
fit_glmer_capture <- function(optimizer) {
  warns <- character(0)
  ctrl <- if (is.null(optimizer)) glmerControl() else glmerControl(optimizer = optimizer)
  fit <- withCallingHandlers(
    glmer(model_formula, data = dat, family = Gamma(link = "log"), control = ctrl),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warns)
}

glmer_summary <- function(fit, label) {
  tab <- as.data.frame(coef(summary(fit)))
  names(tab) <- c("estimate", "se", "z", "p")
  tab <- tibble::rownames_to_column(tab, "term")
  tab$fit <- label
  tab[, c("fit", "term", "estimate", "se", "z", "p")]
}

# The published analysis used the bobyqa optimiser; the current lme4 default
# (nloptwrap) reaches a different, well-behaved solution. Both are reported so
# the discrepancy can be described honestly.
glmer_bobyqa  <- fit_glmer_capture("bobyqa")
glmer_default <- fit_glmer_capture(NULL)

for (nm in c("bobyqa", "default (nloptwrap)")) {
  res <- if (nm == "bobyqa") glmer_bobyqa else glmer_default
  cat("\n-- glmer, optimiser:", nm, "--\n")
  cat("Warnings (verbatim):\n")
  if (length(res$warnings) == 0) {
    cat("  <none>\n")
  } else {
    for (w in res$warnings) cat("  ", w, "\n", sep = "")
  }
  tab <- glmer_summary(res$fit, nm)
  print(as.data.frame(tab), digits = 6)
  cat("SE range:", format(range(tab$se)),
      "| ratio max/min SE:", round(max(tab$se) / min(tab$se), 2), "\n")
  cat("Random intercept SD:", as.data.frame(VarCorr(res$fit))$sdcor[1],
      "| sigma:", sigma(res$fit), "\n")
  cat("logLik:", as.numeric(logLik(res$fit)), "| AIC:", AIC(res$fit), "\n")
}

glmer_tab <- rbind(glmer_summary(glmer_bobyqa$fit, "bobyqa"),
                   glmer_summary(glmer_default$fit, "default (nloptwrap)"))

# Tukey contrasts from the pathological fit, i.e. the numbers the paper reports.
con_bobyqa <- as.data.frame(summary(
  contrast(emmeans(glmer_bobyqa$fit, ~ soundType), method = "pairwise", adjust = "tukey")))
cat("\nTukey sound-type contrasts from the pathological bobyqa fit (link scale):\n")
print(con_bobyqa, digits = 6)
cat("\n")

# ------------------------------------------------ glmmTMB Gamma GLMM ---------

cat("=== glmmTMB Gamma(log) refit ===\n")

tmb_warnings <- character(0)
fit_tmb <- withCallingHandlers(
  glmmTMB(model_formula, data = dat, family = Gamma(link = "log")),
  warning = function(w) {
    tmb_warnings <<- c(tmb_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

cat("Convergence code:", fit_tmb$fit$convergence,
    "| message:", fit_tmb$fit$message, "\n")
cat("pdHess (positive-definite Hessian):", fit_tmb$sdr$pdHess, "\n")
cat("Warnings during fit:",
    if (length(tmb_warnings) == 0) "<none>" else paste(tmb_warnings, collapse = " | "), "\n")
cat("logLik:", as.numeric(logLik(fit_tmb)), "| AIC:", AIC(fit_tmb),
    "| df:", attr(logLik(fit_tmb), "df"), "\n\n")

tmb_coefs <- as.data.frame(summary(fit_tmb)$coefficients$cond)
names(tmb_coefs) <- c("estimate", "se", "z", "p")
tmb_ci <- confint(fit_tmb, parm = "beta_", method = "wald")

fixed_tab <- tibble::tibble(
  term       = rownames(tmb_coefs),
  estimate   = tmb_coefs$estimate,
  se         = tmb_coefs$se,
  z          = tmb_coefs$z,
  p          = tmb_coefs$p,
  ci_low     = tmb_ci[rownames(tmb_coefs), 1],
  ci_high    = tmb_ci[rownames(tmb_coefs), 2],
  ratio      = exp(tmb_coefs$estimate),
  ratio_low  = exp(tmb_ci[rownames(tmb_coefs), 1]),
  ratio_high = exp(tmb_ci[rownames(tmb_coefs), 2])
)

cat("glmmTMB fixed effects (log link, with Wald 95% CI and multiplicative ratios):\n")
print(as.data.frame(fixed_tab), digits = 6)

vc <- VarCorr(fit_tmb)$cond$participantId
sd_id <- sqrt(as.numeric(vc))
shape <- sigma(fit_tmb)
cat("\nRandom intercept SD (participantId):", sd_id, "| variance:", as.numeric(vc), "\n")
cat("Gamma shape parameter:", shape, "\n")

# Marginal / conditional R2 for the Gamma-log model (Nakagawa et al. 2017,
# delta method: distribution-specific variance = trigamma(shape) on the log scale).
var_fixed <- var(as.vector(model.matrix(fit_tmb)[, names(fixef(fit_tmb)$cond), drop = FALSE] %*%
                             fixef(fit_tmb)$cond))
var_rand  <- as.numeric(vc)
var_dist  <- trigamma(shape)
r2m <- var_fixed / (var_fixed + var_rand + var_dist)
r2c <- (var_fixed + var_rand) / (var_fixed + var_rand + var_dist)
icc <- var_rand / (var_rand + var_dist)
cat("Marginal R2:", r2m, "| conditional R2:", r2c, "| ICC:", icc, "\n\n")

# ---------------------------------------------- (b) learning effect ----------

cat("=== (b) Learning effect (trialSequence_c) ===\n")
lrn <- fixed_tab[fixed_tab$term == "trialSequence_c", ]
pct_per_trial <- (exp(lrn$estimate) - 1) * 100
pct_per_trial_ci <- (exp(c(lrn$ci_low, lrn$ci_high)) - 1) * 100
# Trials 1 to 24 span 23 unit increments of trialSequence_c.
span <- 23
pct_cum <- (exp(lrn$estimate * span) - 1) * 100
pct_cum_ci <- (exp(c(lrn$ci_low, lrn$ci_high) * span) - 1) * 100
cat("beta =", lrn$estimate, "| SE =", lrn$se, "| z =", lrn$z, "| p =", lrn$p,
    "| 95% CI [", lrn$ci_low, ",", lrn$ci_high, "]\n")
cat("Percent change per trial:", pct_per_trial,
    "% [", pct_per_trial_ci[1], ",", pct_per_trial_ci[2], "]\n")
cat("Cumulative change trial 1 -> 24 (23 increments):", pct_cum,
    "% [", pct_cum_ci[1], ",", pct_cum_ci[2], "]\n")
cat("For reference, over 24 unit increments:", (exp(lrn$estimate * 24) - 1) * 100, "%\n\n")

# ------------------------------------- (c) spatial factors, omnibus tests ----

cat("=== (c) Omnibus Wald tests for the categorical predictors (type II) ===\n")
anova_tab <- as.data.frame(emmeans::joint_tests(fit_tmb))
print(anova_tab, digits = 6)
cat("\n")

# ------------------------------------------- Tukey sound-type contrasts ------

cat("=== Sound-type estimated marginal means and Tukey contrasts ===\n")

emm_link <- emmeans(fit_tmb, ~ soundType)                    # log scale
emm_resp <- emmeans(fit_tmb, ~ soundType, type = "response") # cm

emm_link_df <- as.data.frame(summary(emm_link, infer = c(TRUE, TRUE)))
emm_resp_df <- as.data.frame(summary(emm_resp, infer = c(TRUE, TRUE)))
cat("EMMs on the link (log) scale:\n"); print(emm_link_df, digits = 6)
cat("\nEMMs back-transformed to cm:\n"); print(emm_resp_df, digits = 6)

con_link <- as.data.frame(summary(contrast(emm_link, method = "pairwise", adjust = "tukey"),
                                  infer = c(TRUE, TRUE)))
con_resp <- as.data.frame(summary(contrast(emm_link, method = "pairwise", adjust = "tukey"),
                                  type = "response", infer = c(TRUE, TRUE)))
cat("\nTukey-adjusted pairwise contrasts, link scale (log differences):\n")
print(con_link, digits = 6)
cat("\nTukey-adjusted pairwise contrasts, response scale (ratios):\n")
print(con_resp, digits = 6)
cat("\n")

# ------------------------------------------------ (d) DHARMa diagnostics -----

cat("=== (d) DHARMa residual diagnostics on the glmmTMB Gamma fit ===\n")
sim_res <- simulateResiduals(fit_tmb, n = 1000, seed = 20260805)

ks_test   <- testUniformity(sim_res, plot = FALSE)
disp_test <- testDispersion(sim_res, plot = FALSE)
out_test  <- testOutliers(sim_res, type = "bootstrap", nBoot = 500, plot = FALSE)
zi_test   <- NULL
qq_slope  <- testQuantiles(sim_res, plot = FALSE)

cat("Uniformity (one-sample KS): D =", unname(ks_test$statistic),
    ", p =", ks_test$p.value, "\n")
cat("Dispersion (simulated): ratio =", unname(disp_test$statistic),
    ", p =", disp_test$p.value,
    ", alternative =", disp_test$alternative, "\n")
cat("Outliers (bootstrap): observed outliers =", unname(out_test$statistic),
    "of", unname(out_test$parameter), "| observed frequency =",
    unname(out_test$estimate), ", p =", out_test$p.value, "\n")
cat("Quantile deviations (testQuantiles): p =", qq_slope$p.value, "\n\n")

png(file.path(results_dir, "accuracy_gamma_glmm_dharma_gamma.png"),
    width = 1400, height = 700, res = 130)
plot(sim_res)
dev.off()

# DHARMa p-values depend on the simulation draw, so the uniformity test is
# repeated over several seeds for both the glmmTMB fit and the pathological
# glmer fit. The original manuscript reported a failing KS test (p = .002);
# this shows where that came from.
ks_seeds <- c(1L, 7L, 42L, 999L, 20260805L)
ks_stability <- do.call(rbind, lapply(ks_seeds, function(s) {
  k_tmb <- testUniformity(simulateResiduals(fit_tmb, n = 1000, seed = s), plot = FALSE)
  k_bad <- testUniformity(simulateResiduals(glmer_bobyqa$fit, n = 1000, seed = s), plot = FALSE)
  data.frame(seed = s,
             glmmTMB_D = unname(k_tmb$statistic), glmmTMB_p = k_tmb$p.value,
             glmer_bobyqa_D = unname(k_bad$statistic), glmer_bobyqa_p = k_bad$p.value)
}))
cat("KS uniformity across simulation seeds (glmmTMB vs pathological glmer):\n")
print(ks_stability, digits = 5, row.names = FALSE)
write.csv(ks_stability, file.path(results_dir, "accuracy_glmm_dharma_ks_seed_stability.csv"),
          row.names = FALSE)
cat("\n")

# ------------------------------------- (e) lognormal robustness check --------

cat("=== (e) Lognormal robustness check ===\n")
fit_ln <- glmmTMB(model_formula, data = dat, family = lognormal(link = "log"))
cat("Convergence code:", fit_ln$fit$convergence, "| pdHess:", fit_ln$sdr$pdHess, "\n")
cat("AIC Gamma:", AIC(fit_tmb), "| AIC lognormal:", AIC(fit_ln), "\n")

ln_coefs <- as.data.frame(summary(fit_ln)$coefficients$cond)
names(ln_coefs) <- c("estimate", "se", "z", "p")
ln_ci <- confint(fit_ln, parm = "beta_", method = "wald")
ln_tab <- tibble::tibble(
  term       = rownames(ln_coefs),
  estimate   = ln_coefs$estimate,
  se         = ln_coefs$se,
  z          = ln_coefs$z,
  p          = ln_coefs$p,
  ci_low     = ln_ci[rownames(ln_coefs), 1],
  ci_high    = ln_ci[rownames(ln_coefs), 2],
  ratio      = exp(ln_coefs$estimate),
  ratio_low  = exp(ln_ci[rownames(ln_coefs), 1]),
  ratio_high = exp(ln_ci[rownames(ln_coefs), 2])
)
cat("\nLognormal fixed effects:\n"); print(as.data.frame(ln_tab), digits = 6)

emm_ln <- emmeans(fit_ln, ~ soundType)
con_ln <- as.data.frame(summary(contrast(emm_ln, method = "pairwise", adjust = "tukey"),
                                type = "response", infer = c(TRUE, TRUE)))
emm_ln_df <- as.data.frame(summary(emmeans(fit_ln, ~ soundType, type = "response"),
                                   infer = c(TRUE, TRUE)))
cat("\nLognormal Tukey contrasts (ratios):\n"); print(con_ln, digits = 6)
cat("\nLognormal EMMs (cm, median scale):\n"); print(emm_ln_df, digits = 6)

sim_ln <- simulateResiduals(fit_ln, n = 1000, seed = 20260805)
ks_ln   <- testUniformity(sim_ln, plot = FALSE)
disp_ln <- testDispersion(sim_ln, plot = FALSE)
out_ln  <- testOutliers(sim_ln, type = "bootstrap", nBoot = 500, plot = FALSE)
cat("\nLognormal DHARMa: KS D =", unname(ks_ln$statistic), ", p =", ks_ln$p.value,
    "| dispersion ratio =", unname(disp_ln$statistic), ", p =", disp_ln$p.value,
    "| outlier p =", out_ln$p.value, "\n")

png(file.path(results_dir, "accuracy_gamma_glmm_dharma_lognormal.png"),
    width = 1400, height = 700, res = 130)
plot(sim_ln)
dev.off()

# Gaussian LMM on log error, as a second robustness variant.
fit_loglmer <- lmerTest::lmer(log(participantError_cm) ~ stimDisparity_c + soundType +
                      trialSequence_c + azimuthSector + elevationCategory +
                      (1 | participantId), data = dat, REML = TRUE)
loglmer_tab <- as.data.frame(coef(summary(fit_loglmer)))
names(loglmer_tab) <- c("estimate", "se", "df", "t", "p")
loglmer_tab <- tibble::rownames_to_column(loglmer_tab, "term")
cat("\nGaussian LMM on log(error) fixed effects:\n"); print(loglmer_tab, digits = 6)

con_loglmer <- as.data.frame(summary(
  contrast(emmeans(fit_loglmer, ~ soundType), method = "pairwise", adjust = "tukey"),
  type = "response", infer = c(TRUE, TRUE)))
cat("\nGaussian LMM on log(error), Tukey contrasts (ratios):\n")
print(con_loglmer, digits = 6)

# --------------------------------------------------------------- outputs -----

write.csv(fixed_tab, file.path(results_dir, "accuracy_gamma_glmm_fixed_effects.csv"),
          row.names = FALSE)
write.csv(glmer_tab, file.path(results_dir, "accuracy_gamma_glmer_degenerate_fixed_effects.csv"),
          row.names = FALSE)
write.csv(con_bobyqa, file.path(results_dir, "accuracy_gamma_glmer_bobyqa_soundtype_contrasts.csv"),
          row.names = FALSE)
write.csv(con_link, file.path(results_dir, "accuracy_gamma_glmm_soundtype_contrasts_link.csv"),
          row.names = FALSE)
write.csv(con_resp, file.path(results_dir, "accuracy_gamma_glmm_soundtype_contrasts_ratio.csv"),
          row.names = FALSE)
write.csv(emm_resp_df, file.path(results_dir, "accuracy_gamma_glmm_soundtype_emmeans_cm.csv"),
          row.names = FALSE)
write.csv(emm_link_df, file.path(results_dir, "accuracy_gamma_glmm_soundtype_emmeans_link.csv"),
          row.names = FALSE)
write.csv(anova_tab, file.path(results_dir, "accuracy_gamma_glmm_joint_tests.csv"),
          row.names = FALSE)
write.csv(ln_tab, file.path(results_dir, "accuracy_lognormal_glmm_fixed_effects.csv"),
          row.names = FALSE)
write.csv(con_ln, file.path(results_dir, "accuracy_lognormal_glmm_soundtype_contrasts_ratio.csv"),
          row.names = FALSE)
write.csv(loglmer_tab, file.path(results_dir, "accuracy_loglmer_fixed_effects.csv"),
          row.names = FALSE)
write.csv(con_loglmer, file.path(results_dir, "accuracy_loglmer_soundtype_contrasts_ratio.csv"),
          row.names = FALSE)

diag_tab <- tibble::tibble(
  model = c(rep("Gamma(log) glmmTMB", 4), rep("lognormal(log) glmmTMB", 3)),
  test  = c("KS uniformity", "dispersion", "outliers (bootstrap)", "quantile deviations",
            "KS uniformity", "dispersion", "outliers (bootstrap)"),
  statistic = c(unname(ks_test$statistic), unname(disp_test$statistic),
                unname(out_test$statistic), NA_real_,
                unname(ks_ln$statistic), unname(disp_ln$statistic),
                unname(out_ln$statistic)),
  p_value = c(ks_test$p.value, disp_test$p.value, out_test$p.value, qq_slope$p.value,
              ks_ln$p.value, disp_ln$p.value, out_ln$p.value)
)
write.csv(diag_tab, file.path(results_dir, "accuracy_glmm_dharma_diagnostics.csv"),
          row.names = FALSE)

fit_tab <- tibble::tibble(
  model = c("Gamma(log) glmmTMB", "lognormal(log) glmmTMB", "Gaussian LMM on log(error)"),
  logLik = c(as.numeric(logLik(fit_tmb)), as.numeric(logLik(fit_ln)),
             as.numeric(logLik(fit_loglmer))),
  AIC = c(AIC(fit_tmb), AIC(fit_ln), AIC(fit_loglmer)),
  random_intercept_sd = c(sd_id, sqrt(as.numeric(VarCorr(fit_ln)$cond$participantId)),
                          as.data.frame(VarCorr(fit_loglmer))$sdcor[1])
)
write.csv(fit_tab, file.path(results_dir, "accuracy_glmm_model_comparison.csv"),
          row.names = FALSE)

cat("\nSession info:\n")
print(sessionInfo())
cat("\nDone.\n")
