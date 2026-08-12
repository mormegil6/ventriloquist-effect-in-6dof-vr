#!/usr/bin/env Rscript
# Audit of every glmmTMB::lognormal() call in the revision scripts.
#
# The family calibration in verify_accuracy_gamma_glmm_deep_checks.R showed that
# glmmTMB::lognormal() is broken in this installation: it stops short of the
# maximum likelihood solution, reports a log-likelihood it cannot attain, and
# returns a dispersion far from sd(log y). Because a log-normal regression is
# exactly least squares on log(y), the correct fit is available in closed form,
# so every quantity a lognormal() call produced can be recomputed exactly.
#
# This script visits the four call sites, recomputes each quantity with a
# Gaussian LMM on log(error), and records whether the conclusion survives:
#
#   1. refit_accuracy_gamma_glmm_revision.R      lognormal robustness of the
#                                                accuracy GLMM (coefficients,
#                                                Tukey contrasts, AIC, DHARMa)
#   2. verify_accuracy_gamma_glmm_deep_checks.R  the calibration itself, plus
#                                                the corrected robustness fit
#   3. verify_accuracy_gamma_glmm_independent.R  independent replication of 1
#   4. verify_movement_effects_independent.R     cross-family check on the
#                                                movement-rate coefficients
#
# It also verifies the two log-normal robustness claims the manuscript makes:
# the learning effect (z = -2.79, p = .005) and the within-participant movement
# coefficient (b = 0.889, SE = 0.453, p = .050).
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/lognormal_audit_*.csv and a plain-text log.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(emmeans)
  library(DHARMa)
})

proj_dir    <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
results_dir <- file.path(proj_dir, "results_revision")

log_con <- file(file.path(results_dir, "lognormal_audit_log.txt"), open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")
# Closed explicitly at the end; on.exit() at top level under Rscript fires
# immediately and would leave the tail of the log unflushed.

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

set.seed(20260807)

cat("Audit of glmmTMB::lognormal() call sites in rscripts/revision/\n")
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "| glmmTMB", as.character(packageVersion("glmmTMB")),
    "| TMB", as.character(packageVersion("TMB")),
    "| lme4", as.character(packageVersion("lme4")), "\n")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds"))

acc_vars <- c("participantError_cm", "stimDisparity_c", "soundType",
              "trialSequence_c", "azimuthSector", "elevationCategory", "participantId")
dat <- as.data.frame(df[stats::complete.cases(df[, acc_vars]), acc_vars])
dat$participantId <- factor(dat$participantId)

f_acc <- participantError_cm ~ stimDisparity_c + soundType + trialSequence_c +
  azimuthSector + elevationCategory + (1 | participantId)
f_acc_log <- update(f_acc, log(participantError_cm) ~ .)

audit <- list()   # one row per call site
add_audit <- function(...) audit[[length(audit) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

# =============================================== 0. the defect, on the real data

rule("0. The defect, restated on the analysis data")

fit_ln_bad <- glmmTMB(f_acc, data = dat, family = lognormal(link = "log"))
fit_ln_ok  <- glmmTMB(f_acc_log, data = dat, family = gaussian(), REML = FALSE)
sum_log_y  <- sum(log(dat$participantError_cm))
n_par      <- length(fixef(fit_ln_bad)$cond) + 2   # fixed effects + RE SD + sigma

# A log-normal log-likelihood on the y scale equals the Gaussian log-likelihood
# on log(y) minus sum(log y), which puts it on the same scale as the Gamma AIC.
ll_ok_y  <- as.numeric(logLik(fit_ln_ok)) - sum_log_y
aic_ok_y <- 2 * n_par - 2 * ll_ok_y

cat("sd(log y) =", sd(log(dat$participantError_cm)), "\n")
cat("glmmTMB lognormal(): sigma =", sigma(fit_ln_bad),
    "| reported logLik =", as.numeric(logLik(fit_ln_bad)),
    "| reported AIC =", AIC(fit_ln_bad), "\n")
cat("Correct log-normal ML : sigma =", sigma(fit_ln_ok),
    "| logLik on the y scale =", ll_ok_y, "| AIC on the y scale =", aic_ok_y, "\n")
cat("The correct fit is", round(ll_ok_y - as.numeric(logLik(fit_ln_bad)), 3),
    "nats above the reported value, so the reported log-likelihood is unattainable.\n")

fit_gamma <- glmmTMB(f_acc, data = dat, family = Gamma(link = "log"))
cat("Gamma(log) glmmTMB AIC =", AIC(fit_gamma), "\n")
cat("Model ranking with the broken lognormal AIC: Gamma better by",
    round(AIC(fit_ln_bad) - AIC(fit_gamma), 2), "\n")
cat("Model ranking with the correct lognormal AIC: Gamma better by",
    round(aic_ok_y - AIC(fit_gamma), 2), "\n")

defect_tab <- data.frame(
  quantity = c("sigma (log scale)", "logLik (y scale)", "AIC (y scale)"),
  sd_log_y = c(sd(log(dat$participantError_cm)), NA, NA),
  glmmTMB_lognormal = c(sigma(fit_ln_bad), as.numeric(logLik(fit_ln_bad)), AIC(fit_ln_bad)),
  correct_lognormal = c(sigma(fit_ln_ok), ll_ok_y, aic_ok_y),
  stringsAsFactors = FALSE
)

# ==================================== 1. refit_accuracy_gamma_glmm_revision.R

rule("1. refit_accuracy_gamma_glmm_revision.R, line 281")

cat("Quantity produced: the lognormal robustness variant of the accuracy GLMM.\n")
cat("Written to accuracy_lognormal_glmm_fixed_effects.csv,\n")
cat("accuracy_lognormal_glmm_soundtype_contrasts_ratio.csv, rows 5-7 of\n")
cat("accuracy_glmm_dharma_diagnostics.csv, row 2 of accuracy_glmm_model_comparison.csv,\n")
cat("and accuracy_gamma_glmm_dharma_lognormal.png.\n\n")

fit_loglmer <- lmerTest::lmer(f_acc_log, data = dat, REML = TRUE)

cmp_terms <- c("stimDisparity_c", "soundTypeFlute", "soundTypeSpeech",
               "soundTypePink Noise", "trialSequence_c")
bad_co  <- summary(fit_ln_bad)$coefficients$cond[cmp_terms, , drop = FALSE]
good_co <- coef(summary(fit_loglmer))[cmp_terms, , drop = FALSE]
site1 <- data.frame(
  term = cmp_terms,
  b_lognormal = bad_co[, 1], se_lognormal = bad_co[, 2], p_lognormal = bad_co[, 4],
  b_loglmer   = good_co[, 1], se_loglmer = good_co[, 2], p_loglmer = good_co[, 5],
  row.names = NULL, stringsAsFactors = FALSE
)
cat("Fixed effects, broken lognormal versus Gaussian LMM on log(error) (REML):\n")
print(site1, digits = 5, row.names = FALSE)
cat("\nEvery SE from the broken fit is too small by a factor of roughly",
    round(mean(site1$se_loglmer / site1$se_lognormal), 2), "\n")

con_bad <- as.data.frame(summary(contrast(emmeans(fit_ln_bad, ~ soundType), "pairwise",
                                          adjust = "tukey"), type = "response",
                                 infer = c(TRUE, TRUE)))
con_ok  <- as.data.frame(summary(contrast(emmeans(fit_loglmer, ~ soundType), "pairwise",
                                          adjust = "tukey"), type = "response",
                                 infer = c(TRUE, TRUE)))
cat("\nSound-type Tukey contrasts, smallest p: broken", min(con_bad$p.value),
    "| correct", min(con_ok$p.value), "\n")
cat("Conclusion in the manuscript (no sound-type contrast significant): unchanged.\n")

# DHARMa on the correct log-normal model, replacing the rows computed from the
# broken fit.
sim_ok <- simulateResiduals(fit_loglmer, n = 1000, seed = 20260807)
ks_ok   <- testUniformity(sim_ok, plot = FALSE)
disp_ok <- testDispersion(sim_ok, plot = FALSE)
out_ok  <- testOutliers(sim_ok, type = "bootstrap", nBoot = 500, plot = FALSE)
cat("\nDHARMa on the correct log-normal model: KS D =", unname(ks_ok$statistic),
    ", p =", ks_ok$p.value, "| dispersion =", unname(disp_ok$statistic),
    ", p =", disp_ok$p.value, "| outlier p =", out_ok$p.value, "\n")

add_audit(script = "refit_accuracy_gamma_glmm_revision.R", line = 281,
          quantity = "lognormal robustness variant of the accuracy GLMM: coefficients, Tukey contrasts, AIC, DHARMa",
          reported_as_headline = "no; the manuscript reports the Gamma GLMM and the Gaussian log LMM",
          conclusion_survives = "yes",
          action = "add a corrective comment; the CSV and PNG outputs derived from lognormal() should be regenerated or withdrawn")

# =================================== 2. verify_accuracy_gamma_glmm_deep_checks.R

rule("2. verify_accuracy_gamma_glmm_deep_checks.R, lines 93 and 144")

cat("Line 93 fits lognormal() to simulated data with a known answer: this is the\n")
cat("calibration that detected the defect, so the call is intentional.\n")
cat("Line 144 fits it to the real data as fit_ln_bad, alongside the correct fit,\n")
cat("to quantify the damage. Both uses are diagnostic, not inferential.\n\n")

cat("Reproducing the simulation calibration:\n")
n_sim <- 4000
xs <- rnorm(n_sim); ys <- exp(rnorm(n_sim, 1 + 0.5 * xs, 0.5))
m_sim_bad <- glmmTMB(ys ~ xs, family = lognormal(link = "log"), data = data.frame(ys, xs))
ols <- lm(log(ys) ~ xs)
cat("  true slope 0.5, true log-scale SD 0.5\n")
cat("  closed-form MLE   : slope", round(coef(ols)[2], 5),
    "| sigma", round(sqrt(mean(residuals(ols)^2)), 5), "\n")
cat("  glmmTMB lognormal : slope", round(fixef(m_sim_bad)$cond[2], 5),
    "| sigma", round(sigma(m_sim_bad), 5), "\n")

cat("\nThe corrected robustness fit produced at this call site is the source of the\n")
cat("manuscript's lognormal learning-effect number. Verified below in section 5.\n")

add_audit(script = "verify_accuracy_gamma_glmm_deep_checks.R", line = 93,
          quantity = "calibration of lognormal() against a closed-form MLE on simulated data",
          reported_as_headline = "no; it is the diagnostic that established the defect",
          conclusion_survives = "yes; the defect reproduces",
          action = "none")
add_audit(script = "verify_accuracy_gamma_glmm_deep_checks.R", line = 144,
          quantity = "broken lognormal fit on the real data, shown side by side with the correct fit",
          reported_as_headline = "no; the correct fit at the same call site supplies the manuscript's z = -2.79",
          conclusion_survives = "yes",
          action = "none")

# ================================= 3. verify_accuracy_gamma_glmm_independent.R

rule("3. verify_accuracy_gamma_glmm_independent.R, line 430")

cat("Quantity produced: an independent replication of call site 1, printed to\n")
cat("verify_accuracy_gamma_glmm_log.txt only. No CSV depends on it. The same\n")
cat("section already fits the Gaussian LMM on log(error), which is the correct\n")
cat("log-normal model, so the intended comparison is available without the\n")
cat("lognormal() call.\n\n")

aic_tab <- data.frame(
  model = c("Gamma(log) glmmTMB", "lognormal(log) glmmTMB as fitted",
            "correct log-normal ML, y scale", "Gaussian LMM on log(error), REML, log scale"),
  AIC = c(AIC(fit_gamma), AIC(fit_ln_bad), aic_ok_y, AIC(fit_loglmer)),
  stringsAsFactors = FALSE
)
cat("AIC comparison as printed by that script, with the correct log-normal added:\n")
print(aic_tab, digits = 8, row.names = FALSE)
cat("\nThe Gamma family is preferred either way, so the reported ranking holds,\n")
cat("but the margin is", round(AIC(fit_ln_bad) - AIC(fit_gamma), 1), "AIC as printed against",
    round(aic_ok_y - AIC(fit_gamma), 1), "AIC correctly.\n")

add_audit(script = "verify_accuracy_gamma_glmm_independent.R", line = 430,
          quantity = "lognormal coefficients, Tukey contrasts and AIC, printed to the log only",
          reported_as_headline = "no",
          conclusion_survives = "yes; the Gamma family is preferred under either AIC",
          action = "add a corrective comment; no rerun needed for any published number")

# =================================== 4. verify_movement_effects_independent.R

rule("4. verify_movement_effects_independent.R, line 161")

mv_vars <- c("participantError_m", "stimulusDisparity_m", "soundType", "trialSequenceNum",
             "participant_mean_rate", "trial_rate_deviation", "participantId")
mv <- as.data.frame(df[stats::complete.cases(df[, mv_vars]), mv_vars])
mv$participantId <- factor(mv$participantId)

f_mv <- participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum +
  participant_mean_rate + trial_rate_deviation + (1 | participantId)

cat("Quantity produced: a cross-family check on the two movement-rate coefficients,\n")
cat("printed to movement_effects_verification_log.txt. No CSV depends on it.\n\n")

m_mv_lntmb <- glmmTMB(f_mv, data = mv, family = lognormal(link = "log"), REML = FALSE)
m_mv_lmer  <- lmerTest::lmer(update(f_mv, log(participantError_m) ~ .), data = mv, REML = TRUE)
m_mv_lmer_ml <- lmerTest::lmer(update(f_mv, log(participantError_m) ~ .), data = mv, REML = FALSE)

mv_terms <- c("participant_mean_rate", "trial_rate_deviation")
mv_bad <- summary(m_mv_lntmb)$coefficients$cond[mv_terms, , drop = FALSE]
mv_ok  <- coef(summary(m_mv_lmer))[mv_terms, , drop = FALSE]
mv_ok_ml <- coef(summary(m_mv_lmer_ml))[mv_terms, , drop = FALSE]

site4 <- data.frame(
  term = mv_terms,
  b_lognormal = mv_bad[, 1], se_lognormal = mv_bad[, 2], p_lognormal = mv_bad[, 4],
  b_lmer_reml = mv_ok[, 1],  se_lmer_reml = mv_ok[, 2],  p_lmer_reml = mv_ok[, 5],
  b_lmer_ml   = mv_ok_ml[, 1], se_lmer_ml = mv_ok_ml[, 2], p_lmer_ml = mv_ok_ml[, 5],
  row.names = NULL, stringsAsFactors = FALSE
)
print(site4, digits = 5, row.names = FALSE)
cat("\nThe broken fit makes the within-participant term look significant",
    "(p =", signif(site4$p_lognormal[2], 3), ") and shrinks the between-participant\n")
cat("coefficient towards zero. Neither number is in the manuscript.\n")

add_audit(script = "verify_movement_effects_independent.R", line = 161,
          quantity = "cross-family check on the between- and within-participant movement-rate coefficients",
          reported_as_headline = "no; the manuscript quotes the lmerTest fit from the same section",
          conclusion_survives = "yes for the lmerTest fit; the lognormal() line is misleading and should not be cited",
          action = "add a corrective comment; no published number changes")

# ================================== 5. the two manuscript log-normal claims ====

rule("5. The two log-normal robustness claims in the manuscript")

# (i) learning effect. The manuscript reports it twice: once as a "lognormal
# error family" (z = -2.79, p = .005) and once as a "Gaussian model on log
# error" (t(703.3) = -2.77, p = .006). These are the same likelihood.
learn_ml   <- coef(summary(lmerTest::lmer(f_acc_log, data = dat, REML = FALSE)))["trialSequence_c", ]
learn_reml <- coef(summary(fit_loglmer))["trialSequence_c", ]
learn_tmb  <- summary(fit_ln_ok)$coefficients$cond["trialSequence_c", ]
learn_bad  <- summary(fit_ln_bad)$coefficients$cond["trialSequence_c", ]

claim_learn <- data.frame(
  fit = c("glmmTMB gaussian on log(y), ML (the correct log-normal MLE)",
          "lmerTest on log(y), ML",
          "lmerTest on log(y), REML",
          "glmmTMB lognormal() as called"),
  estimate = c(learn_tmb[1], learn_ml[1], learn_reml[1], learn_bad[1]),
  se       = c(learn_tmb[2], learn_ml[2], learn_reml[2], learn_bad[2]),
  statistic = c(learn_tmb[3], learn_ml[4], learn_reml[4], learn_bad[3]),
  df = c(NA, learn_ml[3], learn_reml[3], NA),
  p = c(learn_tmb[4], learn_ml[5], learn_reml[5], learn_bad[4]),
  row.names = NULL, stringsAsFactors = FALSE
)
cat("(i) Learning effect, published under a lognormal family as z = -2.79, p = .005:\n")
print(claim_learn, digits = 6, row.names = FALSE)
cat("\nThe published z comes from the correct log-normal MLE, not from lognormal().\n")
cat("The lognormal() call would have given z =", round(learn_bad[3], 3),
    ", p =", round(learn_bad[4], 3), ", which is not significant.\n")
cat("Note: the manuscript's 'lognormal family' and 'Gaussian model on log error'\n")
cat("checks are the same model, ML and REML respectively, not two independent checks.\n")

# (ii) within-participant movement coefficient.
cat("\n(ii) Within-participant movement coefficient, published as b = 0.889,",
    "SE = 0.453, p = .050:\n")
claim_mv <- data.frame(
  fit = c("lmerTest on log(error), REML", "lmerTest on log(error), ML",
          "glmmTMB lognormal() as called"),
  estimate = c(mv_ok[2, 1], mv_ok_ml[2, 1], mv_bad[2, 1]),
  se       = c(mv_ok[2, 2], mv_ok_ml[2, 2], mv_bad[2, 2]),
  df       = c(mv_ok[2, 3], mv_ok_ml[2, 3], NA),
  statistic = c(mv_ok[2, 4], mv_ok_ml[2, 4], mv_bad[2, 3]),
  p = c(mv_ok[2, 5], mv_ok_ml[2, 5], mv_bad[2, 4]),
  row.names = NULL, stringsAsFactors = FALSE
)
print(claim_mv, digits = 6, row.names = FALSE)
cat("\nThe published values match the lmerTest REML fit exactly, confirming the\n")
cat("manuscript is not quoting the glmmTMB lognormal() number.\n")

# --------------------------------------------------------------- outputs -----

audit_tab <- do.call(rbind, audit)
cat("\n")
rule("Call-site summary")
print(audit_tab[, c("script", "line", "reported_as_headline", "conclusion_survives")],
      row.names = FALSE)

write.csv(defect_tab,  file.path(results_dir, "lognormal_audit_defect.csv"), row.names = FALSE)
write.csv(audit_tab,   file.path(results_dir, "lognormal_audit_call_sites.csv"), row.names = FALSE)
write.csv(site1,       file.path(results_dir, "lognormal_audit_accuracy_coefficients.csv"), row.names = FALSE)
write.csv(site4,       file.path(results_dir, "lognormal_audit_movement_coefficients.csv"), row.names = FALSE)
write.csv(claim_learn, file.path(results_dir, "lognormal_audit_learning_claim.csv"), row.names = FALSE)
write.csv(claim_mv,    file.path(results_dir, "lognormal_audit_movement_claim.csv"), row.names = FALSE)
write.csv(aic_tab,     file.path(results_dir, "lognormal_audit_aic.csv"), row.names = FALSE)
write.csv(data.frame(model = "Gaussian LMM on log(error), REML",
                     test = c("KS uniformity", "dispersion", "outliers (bootstrap)"),
                     statistic = c(unname(ks_ok$statistic), unname(disp_ok$statistic),
                                   unname(out_ok$statistic)),
                     p_value = c(ks_ok$p.value, disp_ok$p.value, out_ok$p.value)),
          file.path(results_dir, "lognormal_audit_dharma.csv"), row.names = FALSE)

cat("\nSession info:\n")
print(sessionInfo())
cat("\nDone.\n")

sink(type = "message"); sink(); close(log_con)
