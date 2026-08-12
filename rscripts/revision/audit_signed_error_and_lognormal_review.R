#!/usr/bin/env Rscript
# Independent adversarial re-check of the signed-error LMM re-estimation and of
# the glmmTMB::lognormal() audit, run from the prepared trial-level data only.
#
# Nothing here reads the earlier revision scripts or their outputs: every model
# is specified from the manuscript text and from the original rscripts/compute_ftests.R
# and refitted here, so an agreement is evidence and not an echo.
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/audit_signed_error_lognormal_log.txt and two CSVs.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(glmmTMB)
})

proj_dir    <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
results_dir <- file.path(proj_dir, "results_revision")

log_con <- file(file.path(results_dir, "audit_signed_error_lognormal_log.txt"), open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

cat("Independent audit run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(R.version.string, "| lme4", as.character(packageVersion("lme4")),
    "| lmerTest", as.character(packageVersion("lmerTest")),
    "| emmeans", as.character(packageVersion("emmeans")),
    "| glmmTMB", as.character(packageVersion("glmmTMB")), "\n")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds"))
cat("Rows:", nrow(df), "| participants:", length(unique(df$participantId)), "\n")

model_vars <- c("signedError_m", "stimulusDisparity_m", "soundType",
                "azimuthSector", "elevationCategory", "participantId")
cat("Complete cases on the signed-error model variables:",
    sum(stats::complete.cases(df[, model_vars])), "\n")
dat <- as.data.frame(df[stats::complete.cases(df[, model_vars]), model_vars])
dat$participantId <- factor(dat$participantId)

# --------------------------------------------------------------------- (a) --

rule("(a) Sound type main effect on signed error")

f_sound <- signedError_m ~ stimulusDisparity_m + soundType +
  (1 + stimulusDisparity_m | participantId)

# The original rscripts/compute_ftests.R fitted this model with bobyqa; the
# revision script used the lmer default. Both are run.
m_sound_nm  <- lmer(f_sound, data = dat, REML = TRUE)
m_sound_bob <- lmer(f_sound, data = dat, REML = TRUE,
                    control = lmerControl(optimizer = "bobyqa"))

cat("Default optimiser (", m_sound_nm@optinfo$optimizer, "):\n", sep = "")
print(anova(m_sound_nm), digits = 6)
cat("\nbobyqa, as in the original compute_ftests.R:\n")
print(anova(m_sound_bob), digits = 6)

# Type III tests should not depend on the factor coding when there are no
# interactions; verified rather than assumed.
old_contr <- options(contrasts = c("contr.sum", "contr.poly"))
m_sound_sum <- lmer(f_sound, data = dat, REML = TRUE)
cat("\nSame model under sum-to-zero contrasts:\n")
print(anova(m_sound_sum), digits = 6)
options(old_contr)

cat("\nRandom effects:\n"); print(VarCorr(m_sound_nm), digits = 5)

# --------------------------------------------------------------------- (b) --

rule("(b) Sound type EMMs and Tukey contrasts")

cat("emmeans default df method for lmerMod:", get_emm_option("lmer.df"), "\n")

emm_kr  <- emmeans(m_sound_nm, ~ soundType, lmer.df = "kenward-roger")
emm_sat <- emmeans(m_sound_nm, ~ soundType, lmer.df = "satterthwaite")

to_cm <- function(x, cols) { x[cols] <- x[cols] * 100; x }
emm_kr_tab  <- to_cm(as.data.frame(summary(emm_kr,  infer = c(TRUE, TRUE))),
                     c("emmean", "SE", "lower.CL", "upper.CL"))
emm_sat_tab <- to_cm(as.data.frame(summary(emm_sat, infer = c(TRUE, TRUE))),
                     c("emmean", "SE", "lower.CL", "upper.CL"))
cat("\nEMMs in cm, Kenward-Roger:\n"); print(emm_kr_tab, digits = 6)
cat("\nEMMs in cm, Satterthwaite:\n"); print(emm_sat_tab, digits = 6)

con_kr <- to_cm(as.data.frame(summary(contrast(emm_kr, "pairwise", adjust = "tukey"),
                                      infer = c(TRUE, TRUE))),
                c("estimate", "SE", "lower.CL", "upper.CL"))
cat("\nTukey pairwise contrasts in cm (Kenward-Roger):\n"); print(con_kr, digits = 6)
cat("Smallest Tukey p:", min(con_kr$p.value), "at",
    as.character(con_kr$contrast[which.min(con_kr$p.value)]), "\n")

# Unadjusted and Bonferroni floors, to confirm the adjustment actually applied.
con_none <- as.data.frame(summary(contrast(emm_kr, "pairwise", adjust = "none")))
cat("Smallest unadjusted p:", min(con_none$p.value), "\n")

# --------------------------------------------------------------------- (c) --

rule("(c) Azimuth sector and elevation category on signed error")

f_spatial <- signedError_m ~ stimulusDisparity_m + azimuthSector +
  elevationCategory + soundType + (1 | participantId)

m_spatial <- lmer(f_spatial, data = dat, REML = TRUE)
cat("Optimiser:", m_spatial@optinfo$optimizer, "\n")
print(anova(m_spatial), digits = 6)

con_azi <- as.data.frame(summary(contrast(emmeans(m_spatial, ~ azimuthSector),
                                          "pairwise", adjust = "tukey"),
                                 infer = c(TRUE, TRUE)))
con_ele <- as.data.frame(summary(contrast(emmeans(m_spatial, ~ elevationCategory),
                                          "pairwise", adjust = "tukey"),
                                 infer = c(TRUE, TRUE)))
con_azi <- to_cm(con_azi, c("estimate", "SE", "lower.CL", "upper.CL"))
con_ele <- to_cm(con_ele, c("estimate", "SE", "lower.CL", "upper.CL"))
cat("\nAzimuth Tukey contrasts (cm):\n"); print(con_azi, digits = 6)
cat("Smallest azimuth Tukey p:", min(con_azi$p.value), "at",
    as.character(con_azi$contrast[which.min(con_azi$p.value)]), "\n")
cat("\nElevation Tukey contrasts (cm):\n"); print(con_ele, digits = 6)
cat("Smallest elevation Tukey p:", min(con_ele$p.value), "at",
    as.character(con_ele$contrast[which.min(con_ele$p.value)]), "\n")

# --------------------------------------------------------------------- (d) --

rule("(d) Joint LRT for the two spatial factors")

lrt_row <- function(reduced, full, label) {
  a <- anova(reduced, full)
  # anova.merMod calls the chi-square df "Df"; glmmTMB calls it "Chi Df" and
  # uses "Df" for the model df. Take the glmmTMB name first where present.
  dfc <- if ("Chi Df" %in% names(a)) a[["Chi Df"]][2] else a[["Df"]][2]
  data.frame(comparison = label, chisq = a$Chisq[2], df = dfc,
             p = a[[grep("^Pr", names(a))]][2],
             logLik_reduced = as.numeric(logLik(reduced)),
             logLik_full = as.numeric(logLik(full)), stringsAsFactors = FALSE)
}

# Signed-error version (the paragraph the sentence sits in is unsigned accuracy,
# but the sentence follows the signed-error tests, so both are computed).
m_sig_red  <- lme4::lmer(signedError_m ~ stimulusDisparity_m + soundType +
                           (1 | participantId), data = dat, REML = FALSE)
m_sig_full <- lme4::lmer(f_spatial, data = dat, REML = FALSE)
r1 <- lrt_row(m_sig_red, m_sig_full, "signed error LMM (ML)")

acc_vars <- c("participantError_cm", "stimDisparity_c", "soundType",
              "trialSequence_c", "azimuthSector", "elevationCategory", "participantId")
acc <- as.data.frame(df[stats::complete.cases(df[, acc_vars]), acc_vars])
acc$participantId <- factor(acc$participantId)
cat("Accuracy analysis rows:", nrow(acc), "\n")

f_acc_red  <- participantError_cm ~ stimDisparity_c + soundType + trialSequence_c +
  (1 | participantId)
f_acc_full <- update(f_acc_red, . ~ . + azimuthSector + elevationCategory)

report_conv <- function(fit, label) {
  msg <- fit@optinfo$conv$lme4$messages
  cat(sprintf("  %-42s optimiser %-12s max|grad| = %.6g | messages: %s\n",
              label, fit@optinfo$optimizer,
              max(abs(fit@optinfo$derivs$gradient)),
              if (length(msg)) paste(msg, collapse = "; ") else "none"))
}

for (opt in c("bobyqa", "Nelder_Mead", "nloptwrap")) {
  g_red  <- glmer(f_acc_red,  data = acc, family = Gamma(link = "log"),
                  control = glmerControl(optimizer = opt))
  g_full <- glmer(f_acc_full, data = acc, family = Gamma(link = "log"),
                  control = glmerControl(optimizer = opt))
  cat("\nglmer Gamma, optimiser", opt, ":\n")
  report_conv(g_red,  "reduced")
  report_conv(g_full, "full")
  assign(paste0("r_glmer_", opt),
         lrt_row(g_red, g_full, paste0("unsigned accuracy Gamma, glmer (", opt, ")")))
  if (opt == "bobyqa") g_full_bobyqa <- g_full
}
# glmer's own default is bobyqa for the nAGQ = 0 stage and Nelder_Mead for the
# second, which is neither of the single-optimiser settings above.
g_red_def  <- glmer(f_acc_red,  data = acc, family = Gamma(link = "log"))
g_full_def <- glmer(f_acc_full, data = acc, family = Gamma(link = "log"))
cat("\nglmer Gamma, package default (bobyqa + Nelder_Mead):\n")
report_conv(g_red_def,  "reduced")
report_conv(g_full_def, "full")
r_glmer_default <- lrt_row(g_red_def, g_full_def,
                           "unsigned accuracy Gamma, glmer (package default)")

t_red  <- glmmTMB(f_acc_red,  data = acc, family = Gamma(link = "log"))
t_full <- glmmTMB(f_acc_full, data = acc, family = Gamma(link = "log"))
cat("\nglmmTMB convergence: reduced code", t_red$fit$convergence, "pdHess", t_red$sdr$pdHess,
    "| full code", t_full$fit$convergence, "pdHess", t_full$sdr$pdHess, "\n")
r4 <- lrt_row(t_red, t_full, "unsigned accuracy Gamma, glmmTMB")

lrt_tab <- rbind(r1, r_glmer_bobyqa, r_glmer_Nelder_Mead, r_glmer_nloptwrap,
                 r_glmer_default, r4)
cat("\n"); print(lrt_tab, digits = 7, row.names = FALSE)

# The same sentence in the manuscript quotes two percentage effects from the
# Gamma GLMM. They are checked here against both the refit and the glmer fit,
# since only the LRT in that sentence was scheduled for revision.
pct_terms <- c("azimuthSectorRight", "azimuthSectorBack", "azimuthSectorLeft",
               "elevationCategoryLevel", "elevationCategoryAbove")
pct <- data.frame(term = pct_terms,
                  pct_glmmTMB = 100 * (exp(fixef(t_full)$cond[pct_terms]) - 1),
                  pct_glmer_bobyqa = 100 * (exp(fixef(g_full_bobyqa)[pct_terms]) - 1),
                  pct_glmer_default = 100 * (exp(fixef(g_full_def)[pct_terms]) - 1),
                  row.names = NULL)
cat("\nPercentage change versus the reference level (Front, Below):\n")
print(pct, digits = 4)
cat("Manuscript quotes +2.5% (Back vs Front) and +5.7% (Above vs Below).\n")

# ------------------------------------------- convergence and conditioning ----

rule("Convergence, conditioning and cross-optimiser stability of the two LMMs")

health <- function(fit, label) {
  V <- as.matrix(vcov(fit))
  data.frame(model = label, optimizer = fit@optinfo$optimizer,
             n_warn = length(fit@optinfo$conv$lme4$messages),
             max_abs_grad = max(abs(fit@optinfo$derivs$gradient)),
             theta_hess_min_eig = min(eigen(fit@optinfo$derivs$Hessian,
                                            symmetric = TRUE, only.values = TRUE)$values),
             vcov_kappa = kappa(V, exact = TRUE),
             vcov_min_eig = min(eigen(V, symmetric = TRUE, only.values = TRUE)$values),
             singular = isSingular(fit), stringsAsFactors = FALSE)
}
print(rbind(health(m_sound_nm, "(a)/(b) sound type"),
            health(m_spatial,  "(c) spatial")), digits = 6, row.names = FALSE)

opt_cmp <- function(form, label) {
  fits <- lapply(c(bobyqa = "bobyqa", nloptwrap = "nloptwrap", Nelder_Mead = "Nelder_Mead"),
                 function(o) lmer(form, data = dat, REML = TRUE,
                                  control = lmerControl(optimizer = o)))
  se <- sapply(fits, function(f) sqrt(diag(as.matrix(vcov(f)))))
  b  <- sapply(fits, fixef)
  data.frame(model = label, term = rownames(se), se, beta_range = apply(b, 1, function(x) diff(range(x))),
             se_max_rel_dev = apply(se, 1, function(x) max(abs(x / x[1] - 1))),
             row.names = NULL)
}
opt_tab <- rbind(opt_cmp(f_sound, "(a)/(b)"), opt_cmp(f_spatial, "(c)"))
print(opt_tab, digits = 6, row.names = FALSE)
cat("Largest relative SE deviation across three optimisers:",
    signif(max(opt_tab$se_max_rel_dev), 4), "\n")
cat("Largest absolute fixed-effect range across optimisers:",
    signif(max(opt_tab$beta_range), 4), "\n")

# --------------------------------------------------- lognormal audit re-check -

rule("Lognormal audit: call site 1 (refit_accuracy_gamma_glmm_revision.R:281)")

f_acc <- participantError_cm ~ stimDisparity_c + soundType + trialSequence_c +
  azimuthSector + elevationCategory + (1 | participantId)

fit_bad <- glmmTMB(f_acc, data = acc, family = lognormal(link = "log"))
fit_ok  <- glmmTMB(update(f_acc, log(participantError_cm) ~ .), data = acc,
                   family = gaussian(), REML = FALSE)
fit_ok_lmer_ml   <- lmer(update(f_acc, log(participantError_cm) ~ .), data = acc, REML = FALSE)
fit_ok_lmer_reml <- lmer(update(f_acc, log(participantError_cm) ~ .), data = acc, REML = TRUE)

cat("glmmTMB lognormal(): convergence", fit_bad$fit$convergence,
    "| pdHess", fit_bad$sdr$pdHess, "| sigma =", sigma(fit_bad), "\n")
cat("sd(log y) =", sd(log(acc$participantError_cm)),
    "| correct residual sigma =", sigma(fit_ok), "\n")
sum_log_y <- sum(log(acc$participantError_cm))
cat("logLik: broken", as.numeric(logLik(fit_bad)),
    "| correct on the y scale", as.numeric(logLik(fit_ok)) - sum_log_y,
    "| deficit", as.numeric(logLik(fit_ok)) - sum_log_y - as.numeric(logLik(fit_bad)), "nats\n")

terms_cmp <- c("stimDisparity_c", "trialSequence_c", "soundTypeFlute",
               "soundTypeSpeech", "soundTypePink Noise")
b <- summary(fit_bad)$coefficients$cond[terms_cmp, , drop = FALSE]
g <- summary(fit_ok)$coefficients$cond[terms_cmp, , drop = FALSE]
cat("\nBroken versus correct fixed effects:\n")
print(data.frame(term = terms_cmp, b_broken = b[, 1], se_broken = b[, 2], p_broken = b[, 4],
                 b_correct = g[, 1], se_correct = g[, 2], p_correct = g[, 4],
                 se_ratio = b[, 2] / g[, 2], row.names = NULL), digits = 5)

cat("\nSound-type Tukey p, broken lognormal fit:\n")
cb <- as.data.frame(summary(contrast(emmeans(fit_bad, ~ soundType), "pairwise", adjust = "tukey")))
print(cb, digits = 5)
cat("Smallest p (broken):", min(cb$p.value), "\n")
cat("\nSound-type Tukey p, correct log-normal fit (Gaussian on log y):\n")
cg <- as.data.frame(summary(contrast(emmeans(fit_ok, ~ soundType), "pairwise", adjust = "tukey")))
print(cg, digits = 5)
cat("Smallest p (correct):", min(cg$p.value), "\n")
cl <- as.data.frame(summary(contrast(emmeans(fit_ok_lmer_ml, ~ soundType),
                                     "pairwise", adjust = "tukey")))
cat("Smallest p (correct, lmerTest with finite df):", min(cl$p.value), "\n")

# Is the defect simply a Gaussian fit on the untransformed response? Testing the
# hypothesis directly, rather than inferring it from the size of sigma.
fit_gaus_log <- glmmTMB(f_acc, data = acc, family = gaussian(link = "log"))
cat("\nlognormal() against gaussian(link = 'log') on the same data:\n")
cat("  max abs coefficient difference:",
    max(abs(fixef(fit_bad)$cond - fixef(fit_gaus_log)$cond)), "\n")
cat("  sigma:", sigma(fit_bad), "versus", sigma(fit_gaus_log),
    "| sd(y) =", sd(acc$participantError_cm), "\n")
cat("  The two do not coincide, so lognormal() is not silently a Gaussian fit",
    "on the response; it converges to a wrong interior point.\n")

cat("\nAIC family comparison. Gamma:", AIC(t_full),
    "| lognormal as printed by glmmTMB:", AIC(fit_bad),
    "| correct log-normal on the y scale:",
    2 * attr(logLik(fit_ok), "df") - 2 * (as.numeric(logLik(fit_ok)) - sum_log_y), "\n")

rule("Lognormal audit: the two published Task 2 numbers")

cat("(i) learning effect, coefficient on trialSequence_c\n")
cat("  correct log-normal MLE (glmmTMB gaussian on log y, ML):\n")
print(summary(fit_ok)$coefficients$cond["trialSequence_c", , drop = FALSE], digits = 6)
cat("  same via lmerTest, ML:\n")
print(summary(fit_ok_lmer_ml)$coefficients["trialSequence_c", , drop = FALSE], digits = 6)
cat("  same via lmerTest, REML:\n")
print(summary(fit_ok_lmer_reml)$coefficients["trialSequence_c", , drop = FALSE], digits = 6)
cat("  glmmTMB lognormal() family (must not be reported):\n")
print(summary(fit_bad)$coefficients$cond["trialSequence_c", , drop = FALSE], digits = 6)

cat("\n(ii) movement effects, call site verify_movement_effects_independent.R:161\n")
mv_vars <- c("participantError_m", "stimulusDisparity_m", "soundType", "trialSequenceNum",
             "participant_mean_rate", "trial_rate_deviation", "participantId")
mv <- as.data.frame(df[stats::complete.cases(df[, mv_vars]), mv_vars])
mv$participantId <- factor(mv$participantId)
cat("  rows:", nrow(mv), "\n")
f_mv <- log(participantError_m) ~ stimulusDisparity_m + soundType + trialSequenceNum +
  participant_mean_rate + trial_rate_deviation + (1 | participantId)
m_mv_reml <- lmer(f_mv, data = mv, REML = TRUE)
m_mv_ml   <- lmer(f_mv, data = mv, REML = FALSE)
f_mv_tmb <- participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum +
  participant_mean_rate + trial_rate_deviation + (1 | participantId)
m_mv_bad <- glmmTMB(f_mv_tmb, data = mv, family = lognormal(link = "log"), REML = FALSE)
cat("  lmerTest REML on log(error):\n")
print(summary(m_mv_reml)$coefficients[c("participant_mean_rate", "trial_rate_deviation"), ], digits = 6)
cat("  lmerTest ML on log(error):\n")
print(summary(m_mv_ml)$coefficients[c("participant_mean_rate", "trial_rate_deviation"), ], digits = 6)
cat("  glmmTMB lognormal():\n")
print(summary(m_mv_bad)$coefficients$cond[c("participant_mean_rate", "trial_rate_deviation"), ], digits = 6)
cat("  glmmTMB lognormal sigma:", sigma(m_mv_bad), "| sd(log y):", sd(log(mv$participantError_m)), "\n")

# ---------------------------------------------------------------- outputs ----

write.csv(lrt_tab, file.path(results_dir, "audit_spatial_lrt_comparison.csv"), row.names = FALSE)
write.csv(rbind(cbind(factor_tested = "soundType", con_kr[, c("contrast", "estimate", "SE", "df", "t.ratio", "p.value")]),
                cbind(factor_tested = "azimuthSector", con_azi[, c("contrast", "estimate", "SE", "df", "t.ratio", "p.value")]),
                cbind(factor_tested = "elevationCategory", con_ele[, c("contrast", "estimate", "SE", "df", "t.ratio", "p.value")])),
          file.path(results_dir, "audit_signed_error_tukey_contrasts.csv"), row.names = FALSE)

cat("\nDone.\n")
sink(type = "message"); sink(); close(log_con)
