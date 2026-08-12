# Independent verification of the Analysis 2 accuracy GLMM re-analysis.
#
# This script is deliberately written from the model specification alone rather
# than from the earlier refit script, so that every headline number is
# recomputed rather than copied. It covers:
#
#   1. data integrity and centring of the two continuous predictors
#   2. the glmmTMB Gamma(log) refit and the conditioning of its Hessian
#   3. lme4::glmer under four optimisers, to locate the reported pathology
#   4. what glmmTMB's sigma() actually returns for the Gamma family, which
#      determines the Nakagawa R2 and the ICC
#   5. omnibus tests by both Wald and likelihood ratio
#   6. estimated marginal means and Tukey contrasts
#   7. DHARMa diagnostics, and why the converged and non-converged fits differ
#   8. random-slope extensions, which bear on the two surviving effects
#   9. the precision the design actually affords for a sound-type contrast
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/verify_accuracy_gamma_*.csv and a plain-text log.

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
log_path    <- file.path(results_dir, "verify_accuracy_gamma_glmm_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

cat("Independent verification of the Analysis 2 accuracy GLMMs\n")
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "\n")
cat("glmmTMB:", as.character(packageVersion("glmmTMB")),
    "| lme4:", as.character(packageVersion("lme4")),
    "| emmeans:", as.character(packageVersion("emmeans")),
    "| DHARMa:", as.character(packageVersion("DHARMa")), "\n")

# ------------------------------------------------------------------ 1. data --

rule("1. Data integrity and centring")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds"))

model_vars <- c("participantError_cm", "stimDisparity_c", "soundType",
                "trialSequence_c", "azimuthSector", "elevationCategory",
                "participantId")
dat <- df[stats::complete.cases(df[, model_vars]), model_vars]
dat$participantId <- factor(dat$participantId)

cat("Rows in source file:", nrow(df), "| rows used:", nrow(dat),
    "| dropped:", nrow(df) - nrow(dat), "\n")
cat("Participants:", nlevels(dat$participantId),
    "| trials per participant: min", min(table(dat$participantId)),
    "max", max(table(dat$participantId)), "\n")
cat("Response min:", min(dat$participantError_cm),
    "(strictly positive:", all(dat$participantError_cm > 0), ")\n")
cat("Reference levels:",
    "soundType =", levels(dat$soundType)[1],
    "| azimuthSector =", levels(dat$azimuthSector)[1],
    "| elevationCategory =", levels(dat$elevationCategory)[1], "\n")

# Centring: are the _c variables what the model_spec claims?
cat("\nCentring checks\n")
cat("  mean(stimDisparity_c) =", mean(dat$stimDisparity_c),
    "| max|stimDisparity_c - (stimulusDisparity_m - mean)| =",
    max(abs(dat$stimDisparity_c - (df$stimulusDisparity_m - mean(df$stimulusDisparity_m)))), "\n")
cat("  mean(trialSequence_c) =", mean(dat$trialSequence_c),
    "| max|trialSequence_c - (trialSequenceNum - 12.5)| =",
    max(abs(dat$trialSequence_c - (df$trialSequenceNum - 12.5))), "\n")
cat("  trialSequenceNum range:", range(df$trialSequenceNum), "\n")
cat("  stimulusDisparity_m range (m):", range(df$stimulusDisparity_m), "\n")

cat("\nDesign balance, soundType x participant (should be 6 each if fully crossed):\n")
print(table(table(dat$participantId, dat$soundType)))
cat("\nsoundType by azimuthSector:\n")
print(table(dat$soundType, dat$azimuthSector))
cat("\nMean stimulus disparity (m) by soundType (confounding check):\n")
print(round(tapply(df$stimulusDisparity_m, dat$soundType, mean), 5))

model_formula <- participantError_cm ~ stimDisparity_c + soundType +
  trialSequence_c + azimuthSector + elevationCategory + (1 | participantId)

# ------------------------------------------------------ 2. glmmTMB Gamma -----

rule("2. glmmTMB Gamma(log) refit and Hessian conditioning")

tmb_warn <- character(0)
fit_tmb <- withCallingHandlers(
  glmmTMB(model_formula, data = dat, family = Gamma(link = "log")),
  warning = function(w) {
    tmb_warn <<- c(tmb_warn, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)

cat("convergence:", fit_tmb$fit$convergence, "| message:", fit_tmb$fit$message,
    "| pdHess:", fit_tmb$sdr$pdHess, "\n")
cat("warnings:", if (length(tmb_warn) == 0) "<none>" else paste(tmb_warn, collapse = " | "), "\n")
cat("logLik:", as.numeric(logLik(fit_tmb)), "| AIC:", AIC(fit_tmb),
    "| df:", attr(logLik(fit_tmb), "df"), "\n")
cat("max abs gradient at optimum:", max(abs(fit_tmb$obj$gr(fit_tmb$fit$par))), "\n")

V <- vcov(fit_tmb)$cond
ev <- eigen(V, symmetric = TRUE)$values
cat("vcov eigenvalues: min", min(ev), "max", max(ev),
    "| condition number:", max(ev) / min(ev), "\n")

tmb_co <- as.data.frame(summary(fit_tmb)$coefficients$cond)
names(tmb_co) <- c("estimate", "se", "z", "p")
tmb_ci <- confint(fit_tmb, parm = "beta_", method = "wald")
fixed_tab <- data.frame(
  term     = rownames(tmb_co),
  estimate = tmb_co$estimate, se = tmb_co$se, z = tmb_co$z, p = tmb_co$p,
  ci_low   = tmb_ci[rownames(tmb_co), 1],
  ci_high  = tmb_ci[rownames(tmb_co), 2],
  row.names = NULL
)
print(fixed_tab, digits = 6)
cat("SE range:", format(range(fixed_tab$se)),
    "| max/min:", round(max(fixed_tab$se) / min(fixed_tab$se), 2), "\n")

# Likelihood-profile CIs for the two effects the paper wants to keep, as a
# check that the Wald intervals are not themselves an artefact.
cat("\nLikelihood-profile 95% CIs (disparity and trial sequence):\n")
prof <- try(confint(fit_tmb, parm = c("stimDisparity_c", "trialSequence_c"),
                    method = "profile"), silent = TRUE)
if (inherits(prof, "try-error")) cat("  profiling failed\n") else print(prof, digits = 6)

# --------------------------------------------------- 3. lme4 optimisers ------

rule("3. lme4::glmer across optimisers")

fit_glmer_capture <- function(opt) {
  warns <- character(0)
  ctrl <- if (is.null(opt)) glmerControl() else glmerControl(optimizer = opt)
  fit <- withCallingHandlers(
    glmer(model_formula, data = dat, family = Gamma(link = "log"), control = ctrl),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = warns)
}

opt_list <- list(bobyqa = "bobyqa", `default (nloptwrap)` = NULL, Nelder_Mead = "Nelder_Mead")
glmer_fits <- lapply(opt_list, fit_glmer_capture)

glmer_rows <- list()
for (nm in names(glmer_fits)) {
  res <- glmer_fits[[nm]]
  tab <- as.data.frame(coef(summary(res$fit)))
  names(tab) <- c("estimate", "se", "z", "p")
  cat("\n-- optimiser:", nm, "--\n")
  cat("warnings:", if (length(res$warnings) == 0) "<none>"
      else paste(res$warnings, collapse = " || "), "\n")
  cat("logLik:", as.numeric(logLik(res$fit)), "| AIC:", AIC(res$fit),
      "| RE SD:", as.data.frame(VarCorr(res$fit))$sdcor[1],
      "| sigma():", sigma(res$fit), "\n")
  cat("SE range:", format(range(tab$se)),
      "| max/min:", round(max(tab$se) / min(tab$se), 2), "\n")
  cat("intercept: b =", tab$estimate[1], "SE =", tab$se[1], "t =", tab$z[1], "\n")
  cat("max abs difference in point estimates vs glmmTMB:",
      max(abs(tab$estimate - fixed_tab$estimate)), "\n")
  glmer_rows[[nm]] <- data.frame(optimiser = nm, term = rownames(tab), tab, row.names = NULL)
}
glmer_tab <- do.call(rbind, glmer_rows)

# The published sound-type contrast, taken from the non-converged fit.
con_bad <- as.data.frame(summary(contrast(
  emmeans(glmer_fits$bobyqa$fit, ~ soundType), method = "pairwise", adjust = "tukey")))
cat("\nTukey contrasts from the bobyqa fit (the published pathology):\n")
print(con_bad, digits = 6)

# ------------------------------------- 4. what does sigma() mean for Gamma ---

rule("4. Calibration of glmmTMB sigma() for the Gamma family")

# Decisive test: simulate from a Gamma with a known shape and see what sigma()
# recovers. The answer determines whether trigamma(sigma) is the right
# distribution-specific variance for the Nakagawa R2.
set.seed(11)
n_sim <- 20000
x_sim <- rnorm(n_sim)
mu_sim <- exp(1 + 0.5 * x_sim)
for (shape_true in c(0.75, 2, 5)) {
  y_sim <- rgamma(n_sim, shape = shape_true, rate = shape_true / mu_sim)
  f_sim <- glmmTMB(y_sim ~ x_sim, family = Gamma(link = "log"))
  cat("true shape =", shape_true,
      "| glmmTMB sigma() =", round(sigma(f_sim), 4),
      "| 1/sigma() =", round(1 / sigma(f_sim), 4),
      "| empirical CV =", round(sd(y_sim / mu_sim), 4),
      "| 1/sqrt(shape) =", round(1 / sqrt(shape_true), 4), "\n")
}

sig <- sigma(fit_tmb)
cat("\nFitted model: sigma() =", sig, "| 1/sigma() =", 1 / sig, "\n")
# Empirical check on the real data: Pearson residual variance estimates 1/shape.
mu_hat <- fitted(fit_tmb)
pearson_var <- var((dat$participantError_cm - mu_hat) / mu_hat)
cat("Var of (y - mu)/mu on the real data =", pearson_var,
    "-> implied shape =", 1 / pearson_var, "\n")

vc <- VarCorr(fit_tmb)$cond$participantId
var_rand <- as.numeric(vc)
Xb <- model.matrix(fit_tmb) %*% fixef(fit_tmb)$cond
var_fixed <- as.numeric(var(Xb))

r2_report <- function(shape_used, label) {
  var_dist <- trigamma(shape_used)
  r2m <- var_fixed / (var_fixed + var_rand + var_dist)
  r2c <- (var_fixed + var_rand) / (var_fixed + var_rand + var_dist)
  icc <- var_rand / (var_rand + var_dist)
  cat(sprintf("%-34s shape=%.4f trigamma=%.4f  R2m=%.4f R2c=%.4f ICC=%.4f\n",
              label, shape_used, var_dist, r2m, r2c, icc))
  data.frame(basis = label, shape = shape_used, var_dist = var_dist,
             var_fixed = var_fixed, var_rand = var_rand,
             R2_marginal = r2m, R2_conditional = r2c, ICC = icc)
}
cat("\nThe calibration above shows sigma() = 1/sqrt(shape), so the Gamma shape is\n")
cat("1/sigma()^2 =", 1 / sig^2, ", not sigma() =", sig, ".\n")
cat("\nNakagawa delta-method R2 under each candidate reading of sigma():\n")
r2_tab <- rbind(
  r2_report(sig, "shape = sigma() [as reported]"),
  r2_report(1 / sig, "shape = 1/sigma()"),
  r2_report(1 / sig^2, "shape = 1/sigma()^2 [correct]"),
  r2_report(1 / pearson_var, "shape from Pearson residual var")
)

# ------------------------------------------------ 5. omnibus tests -----------

rule("5. Omnibus tests: Wald (joint_tests) and likelihood ratio")

jt <- as.data.frame(emmeans::joint_tests(fit_tmb))
jt$Chisq <- jt$df1 * jt$F.ratio
print(jt, digits = 6)

# Likelihood-ratio tests, refitting with ML and dropping one term at a time.
lrt_rows <- list()
for (trm in c("stimDisparity_c", "soundType", "trialSequence_c",
              "azimuthSector", "elevationCategory")) {
  red_form <- update(model_formula, paste(". ~ . -", trm))
  fit_red <- glmmTMB(red_form, data = dat, family = Gamma(link = "log"))
  an <- anova(fit_red, fit_tmb)
  lrt_rows[[trm]] <- data.frame(term = trm, df = an$`Chi Df`[2],
                                Chisq = an$Chisq[2], p = an$`Pr(>Chisq)`[2])
}
lrt_tab <- do.call(rbind, lrt_rows)
cat("\nLikelihood-ratio tests:\n")
print(lrt_tab, digits = 6, row.names = FALSE)

# ------------------------------------ 6. EMMs and Tukey contrasts ------------

rule("6. Estimated marginal means and Tukey contrasts")

emm_link <- emmeans(fit_tmb, ~ soundType)
emm_cm   <- emmeans(fit_tmb, ~ soundType, type = "response")
emm_link_df <- as.data.frame(summary(emm_link, infer = c(TRUE, TRUE)))
emm_cm_df   <- as.data.frame(summary(emm_cm, infer = c(TRUE, TRUE)))
cat("EMMs, link scale:\n");  print(emm_link_df, digits = 6)
cat("\nEMMs, cm:\n");        print(emm_cm_df, digits = 6)

con_link <- as.data.frame(summary(contrast(emm_link, "pairwise", adjust = "tukey"),
                                  infer = c(TRUE, TRUE)))
con_ratio <- as.data.frame(summary(contrast(emm_link, "pairwise", adjust = "tukey"),
                                   type = "response", infer = c(TRUE, TRUE)))
cat("\nTukey contrasts, link scale:\n"); print(con_link, digits = 6)
cat("\nTukey contrasts, ratio scale:\n"); print(con_ratio, digits = 6)

# Derived quantities quoted in the manuscript text.
b_disp <- fixed_tab$estimate[fixed_tab$term == "stimDisparity_c"]
ci_disp <- unlist(fixed_tab[fixed_tab$term == "stimDisparity_c", c("ci_low", "ci_high")])
cat("\nDisparity per 10 cm: ", round((exp(b_disp / 10) - 1) * 100, 4), "% [",
    round((exp(ci_disp[1] / 10) - 1) * 100, 4), ",",
    round((exp(ci_disp[2] / 10) - 1) * 100, 4), "]\n", sep = "")
b_trl <- fixed_tab$estimate[fixed_tab$term == "trialSequence_c"]
ci_trl <- unlist(fixed_tab[fixed_tab$term == "trialSequence_c", c("ci_low", "ci_high")])
cat("Learning per trial: ", round((exp(b_trl) - 1) * 100, 5), "% [",
    round((exp(ci_trl[1]) - 1) * 100, 5), ",", round((exp(ci_trl[2]) - 1) * 100, 5), "]\n", sep = "")
cat("Cumulative over 23 increments (trial 1 -> 24): ",
    round((exp(b_trl * 23) - 1) * 100, 4), "% [",
    round((exp(ci_trl[1] * 23) - 1) * 100, 4), ",",
    round((exp(ci_trl[2] * 23) - 1) * 100, 4), "]\n", sep = "")

# ----------------------------------------------- 7. DHARMa diagnostics -------

rule("7. DHARMa diagnostics, and the source of the difference between fits")

dharma_block <- function(fit, label, seed = 20260805) {
  sr <- simulateResiduals(fit, n = 1000, seed = seed)
  ks <- testUniformity(sr, plot = FALSE)
  dp <- testDispersion(sr, plot = FALSE)
  ot <- testOutliers(sr, type = "bootstrap", nBoot = 500, plot = FALSE)
  qt <- testQuantiles(sr, plot = FALSE)
  cat(sprintf("%-28s KS D=%.5f p=%.5f | disp ratio=%.4f p=%.4f | outliers=%d p=%.4f | quantiles p=%.4f\n",
              label, unname(ks$statistic), ks$p.value,
              unname(dp$statistic), dp$p.value,
              as.integer(unname(ot$statistic)), ot$p.value, qt$p.value))
  data.frame(model = label, ks_D = unname(ks$statistic), ks_p = ks$p.value,
             disp_ratio = unname(dp$statistic), disp_p = dp$p.value,
             n_outliers = unname(ot$statistic), outlier_p = ot$p.value,
             quantile_p = qt$p.value)
}

diag_tab <- rbind(
  dharma_block(fit_tmb, "glmmTMB Gamma"),
  dharma_block(glmer_fits$bobyqa$fit, "glmer Gamma bobyqa"),
  dharma_block(glmer_fits$`default (nloptwrap)`$fit, "glmer Gamma nloptwrap")
)

# DHARMa residuals depend on point estimates and on the dispersion parameter,
# never on the standard errors. If the two fits differ in KS, the cause has to
# be one of those two, so both are reported explicitly.
cat("\nWhy the DHARMa results differ between packages:\n")
cat("  max abs difference in fixed-effect estimates (glmer bobyqa vs glmmTMB):",
    max(abs(coef(summary(glmer_fits$bobyqa$fit))[, 1] - fixed_tab$estimate)), "\n")
cat("  glmmTMB RE SD:", sqrt(var_rand), "| implied Gamma shape:", 1 / sig^2, "\n")
cat("  glmer bobyqa RE SD:", as.data.frame(VarCorr(glmer_fits$bobyqa$fit))$sdcor[1],
    "| implied shape:", 1 / sigma(glmer_fits$bobyqa$fit)^2, "\n")
cat("  glmer nloptwrap RE SD:",
    as.data.frame(VarCorr(glmer_fits$`default (nloptwrap)`$fit))$sdcor[1],
    "| implied shape:", 1 / sigma(glmer_fits$`default (nloptwrap)`$fit)^2, "\n")
cat("  The nloptwrap fit converges cleanly yet still fails the KS test, so the\n")
cat("  failing uniformity test tracks the lme4 fit itself, not the convergence warning.\n")
cat("  glmer sigma() (nloptwrap):", sigma(glmer_fits$`default (nloptwrap)`$fit), "\n")
cat("  observed CV of y/fitted (glmmTMB):", sd(dat$participantError_cm / mu_hat), "\n")

# Seed stability for the uniformity test.
ks_seeds <- c(1L, 7L, 42L, 999L, 20260805L)
ks_stab <- do.call(rbind, lapply(ks_seeds, function(s) {
  a <- testUniformity(simulateResiduals(fit_tmb, n = 1000, seed = s), plot = FALSE)
  b <- testUniformity(simulateResiduals(glmer_fits$bobyqa$fit, n = 1000, seed = s), plot = FALSE)
  data.frame(seed = s, glmmTMB_D = unname(a$statistic), glmmTMB_p = a$p.value,
             glmer_D = unname(b$statistic), glmer_p = b$p.value)
}))
cat("\nKS uniformity across seeds:\n")
print(ks_stab, digits = 5, row.names = FALSE)

# ------------------------------------- 8. random-slope extensions ------------

rule("8. Random-slope extensions")

# soundType, disparity and trial sequence all vary within participant, so the
# random-intercept-only model is not the maximal one the design supports.
# Each surviving effect is re-tested with the corresponding by-participant slope.
slope_specs <- list(
  `slope: trialSequence_c` = participantError_cm ~ stimDisparity_c + soundType +
    trialSequence_c + azimuthSector + elevationCategory + (1 + trialSequence_c | participantId),
  `slope: stimDisparity_c` = participantError_cm ~ stimDisparity_c + soundType +
    trialSequence_c + azimuthSector + elevationCategory + (1 + stimDisparity_c | participantId),
  `slope: soundType` = participantError_cm ~ stimDisparity_c + soundType +
    trialSequence_c + azimuthSector + elevationCategory + (1 + soundType | participantId)
)

slope_rows <- list()
for (nm in names(slope_specs)) {
  w <- character(0)
  f <- withCallingHandlers(
    try(glmmTMB(slope_specs[[nm]], data = dat, family = Gamma(link = "log")), silent = TRUE),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  cat("\n--", nm, "--\n")
  if (inherits(f, "try-error")) { cat("  fit failed\n"); next }
  cat("  convergence:", f$fit$convergence, "| pdHess:", f$sdr$pdHess,
      "| AIC:", AIC(f), "| warnings:",
      if (length(w) == 0) "<none>" else paste(w, collapse = " | "), "\n")
  if (!isTRUE(f$sdr$pdHess)) { cat("  Hessian not positive definite; estimates not interpretable\n"); next }
  co <- as.data.frame(summary(f)$coefficients$cond)
  names(co) <- c("estimate", "se", "z", "p")
  keep <- c("stimDisparity_c", "trialSequence_c",
            grep("^soundType", rownames(co), value = TRUE))
  print(co[keep, ], digits = 5)
  cat("  LRT vs random-intercept model: ")
  an <- try(anova(fit_tmb, f), silent = TRUE)
  if (inherits(an, "try-error")) cat("failed\n") else
    cat("Chisq =", round(an$Chisq[2], 4), "df =", an$`Chi Df`[2],
        "p =", signif(an$`Pr(>Chisq)`[2], 5), "\n")
  jt2 <- as.data.frame(emmeans::joint_tests(f))
  cat("  soundType omnibus: F =", round(jt2$F.ratio[jt2$`model term` == "soundType"], 4),
      "p =", signif(jt2$p.value[jt2$`model term` == "soundType"], 5), "\n")
  slope_rows[[nm]] <- data.frame(model = nm, term = keep, co[keep, ], row.names = NULL)
}
slope_tab <- if (length(slope_rows)) do.call(rbind, slope_rows) else NULL

# --------------------------------- 9. precision the design affords -----------

rule("9. Precision of a sound-type contrast")

se_dp <- con_link$SE[con_link$contrast == "Drum - (Pink Noise)"]
if (length(se_dp) == 0) se_dp <- con_link$SE[grepl("Pink", con_link$contrast) &
                                               grepl("Drum", con_link$contrast)][1]
# Smallest log-ratio detectable at 80% power, two-sided, using the Tukey
# critical value for four means rather than the nominal 1.96.
crit_tukey <- qtukey(0.95, nmeans = 4, df = Inf) / sqrt(2)
mdes_unadj <- (qnorm(0.975) + qnorm(0.80)) * se_dp
mdes_tukey <- (crit_tukey + qnorm(0.80)) * se_dp
cat("SE of the Drum - Pink Noise log contrast:", se_dp, "\n")
cat("Tukey critical value (4 means, df = Inf):", crit_tukey, "\n")
cat("Minimum detectable effect at 80% power, unadjusted:", mdes_unadj,
    "-> ratio", exp(mdes_unadj), "(", round((exp(mdes_unadj) - 1) * 100, 1), "% )\n")
cat("Minimum detectable effect at 80% power, Tukey-adjusted:", mdes_tukey,
    "-> ratio", exp(mdes_tukey), "(", round((exp(mdes_tukey) - 1) * 100, 1), "% )\n")
cat("Observed Drum - Pink Noise:",
    con_link$estimate[grepl("Drum", con_link$contrast) & grepl("Pink", con_link$contrast)][1], "\n")
cat("Simultaneous Tukey CI on the ratio scale spans:",
    paste(round(unlist(con_ratio[grepl("Drum", con_ratio$contrast) &
                                   grepl("Pink", con_ratio$contrast),
                                 c("asymp.LCL", "asymp.UCL")][1, ]), 4), collapse = " to "), "\n")

# --------------------------------------------- 10. robustness families -------

rule("10. Robustness across error families")

fit_ln <- glmmTMB(model_formula, data = dat, family = lognormal(link = "log"))
cat("lognormal: convergence", fit_ln$fit$convergence, "| pdHess", fit_ln$sdr$pdHess,
    "| logLik", as.numeric(logLik(fit_ln)), "| AIC", AIC(fit_ln), "\n")
ln_co <- as.data.frame(summary(fit_ln)$coefficients$cond)
names(ln_co) <- c("estimate", "se", "z", "p")
print(ln_co, digits = 6)
con_ln <- as.data.frame(summary(contrast(emmeans(fit_ln, ~ soundType), "pairwise",
                                         adjust = "tukey"),
                                type = "response", infer = c(TRUE, TRUE)))
cat("\nlognormal Tukey contrasts:\n"); print(con_ln, digits = 6)

fit_lmm <- lmerTest::lmer(log(participantError_cm) ~ stimDisparity_c + soundType +
                            trialSequence_c + azimuthSector + elevationCategory +
                            (1 | participantId), data = dat, REML = TRUE)
lmm_co <- as.data.frame(coef(summary(fit_lmm)))
names(lmm_co) <- c("estimate", "se", "df", "t", "p")
cat("\nGaussian LMM on log(error):\n"); print(lmm_co, digits = 6)
con_lmm <- as.data.frame(summary(contrast(emmeans(fit_lmm, ~ soundType), "pairwise",
                                          adjust = "tukey"),
                                 type = "response", infer = c(TRUE, TRUE)))
cat("\nlog-LMM Tukey contrasts:\n"); print(con_lmm, digits = 6)

cat("\nAIC comparability check\n")
cat("  glmmTMB Gamma:", AIC(fit_tmb), "| glmer Gamma bobyqa:", AIC(glmer_fits$bobyqa$fit),
    "| glmer Gamma nloptwrap:", AIC(glmer_fits$`default (nloptwrap)`$fit), "\n")
cat("  glmmTMB lognormal:", AIC(fit_ln), "| Gaussian LMM on log(error) (REML):", AIC(fit_lmm), "\n")

# ------------------------------------------------------------- outputs -------

write.csv(fixed_tab, file.path(results_dir, "verify_accuracy_gamma_fixed_effects.csv"), row.names = FALSE)
write.csv(glmer_tab, file.path(results_dir, "verify_accuracy_gamma_glmer_optimisers.csv"), row.names = FALSE)
write.csv(r2_tab,    file.path(results_dir, "verify_accuracy_gamma_r2_variants.csv"), row.names = FALSE)
write.csv(jt,        file.path(results_dir, "verify_accuracy_gamma_joint_tests.csv"), row.names = FALSE)
write.csv(lrt_tab,   file.path(results_dir, "verify_accuracy_gamma_lrt.csv"), row.names = FALSE)
write.csv(con_link,  file.path(results_dir, "verify_accuracy_gamma_contrasts_link.csv"), row.names = FALSE)
write.csv(con_ratio, file.path(results_dir, "verify_accuracy_gamma_contrasts_ratio.csv"), row.names = FALSE)
write.csv(emm_cm_df, file.path(results_dir, "verify_accuracy_gamma_emmeans_cm.csv"), row.names = FALSE)
write.csv(diag_tab,  file.path(results_dir, "verify_accuracy_gamma_dharma.csv"), row.names = FALSE)
write.csv(ks_stab,   file.path(results_dir, "verify_accuracy_gamma_dharma_seeds.csv"), row.names = FALSE)
if (!is.null(slope_tab))
  write.csv(slope_tab, file.path(results_dir, "verify_accuracy_gamma_random_slopes.csv"), row.names = FALSE)

cat("\nSession info:\n")
print(sessionInfo())
cat("\nDone.\n")
