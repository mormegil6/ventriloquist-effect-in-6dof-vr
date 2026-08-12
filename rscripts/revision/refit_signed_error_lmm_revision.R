#!/usr/bin/env Rscript
# Re-estimation of the three signed-error (visual bias) Gaussian LMM results.
#
# These are the only Analysis 2 / Analysis 4 results that were not refitted
# during the revision, so the repository would not otherwise reproduce them.
# They come from the same lme4 pipeline whose Gamma GLMMs turned out to have
# degenerate variance-covariance matrices, so each fit is checked here for
# convergence and for a well-conditioned Hessian, and every standard error is
# compared across two optimisers.
#
#   (a) sound type main effect on signed error
#   (b) estimated marginal means and Tukey contrasts for sound type
#   (c) azimuth sector and elevation category effects on signed error
#   (d) joint likelihood ratio test for adding the two spatial factors
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/signed_error_lmm_*.csv and a plain-text log.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

proj_dir    <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
results_dir <- file.path(proj_dir, "results_revision")

log_con <- file(file.path(results_dir, "signed_error_lmm_refit_log.txt"), open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")
# The sink is closed explicitly at the end of the script. on.exit() is not used:
# registered at top level under Rscript it fires immediately, which leaves the
# connection open and the tail of the log unflushed.

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

cat("Signed-error LMMs: re-estimation for the revision\n")
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "| lme4", as.character(packageVersion("lme4")),
    "| lmerTest", as.character(packageVersion("lmerTest")),
    "| emmeans", as.character(packageVersion("emmeans")), "\n")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds"))
model_vars <- c("signedError_m", "stimulusDisparity_m", "soundType",
                "azimuthSector", "elevationCategory", "participantId")
dat <- as.data.frame(df[stats::complete.cases(df[, model_vars]), model_vars])
dat$participantId <- factor(dat$participantId)
cat("Trials:", nrow(dat), "| participants:", nlevels(dat$participantId), "\n")

# Convergence and conditioning of a fitted LMM. A degenerate Hessian is what
# broke the Gamma GLMMs, so both the random-effects optimum and the fixed-effect
# vcov are inspected explicitly rather than trusted.
fit_health <- function(fit, label) {
  opt   <- fit@optinfo
  grad  <- max(abs(opt$derivs$gradient))
  hess  <- opt$derivs$Hessian
  V     <- as.matrix(vcov(fit))
  data.frame(
    model            = label,
    optimizer        = opt$optimizer,
    n_warnings       = length(opt$conv$lme4$messages),
    warnings         = if (length(opt$conv$lme4$messages))
                         paste(opt$conv$lme4$messages, collapse = "; ") else "none",
    max_abs_gradient = grad,
    theta_hessian_min_eigen = min(eigen(hess, symmetric = TRUE, only.values = TRUE)$values),
    vcov_condition_number   = kappa(V, exact = TRUE),
    vcov_min_eigen          = min(eigen(V, symmetric = TRUE, only.values = TRUE)$values),
    is_singular      = isSingular(fit),
    stringsAsFactors = FALSE
  )
}

# Refit under a second optimiser and report the largest relative disagreement
# in the fixed-effect standard errors. Under a sound Hessian this is ~1e-6;
# in the degenerate Gamma fits it reached a factor of 240.
optimiser_agreement <- function(formula, data, label) {
  fits <- list(
    bobyqa    = lmer(formula, data = data, REML = TRUE,
                     control = lmerControl(optimizer = "bobyqa")),
    nloptwrap = lmer(formula, data = data, REML = TRUE,
                     control = lmerControl(optimizer = "nloptwrap"))
  )
  se <- sapply(fits, function(f) sqrt(diag(as.matrix(vcov(f)))))
  b  <- sapply(fits, fixef)
  out <- data.frame(
    model     = label,
    term      = rownames(se),
    se_bobyqa = se[, "bobyqa"],
    se_nloptwrap = se[, "nloptwrap"],
    se_ratio  = se[, "nloptwrap"] / se[, "bobyqa"],
    beta_abs_diff = abs(b[, "nloptwrap"] - b[, "bobyqa"]),
    row.names = NULL, stringsAsFactors = FALSE
  )
  attr(out, "fits") <- fits
  out
}

# ================================================ (a) sound type, signed error

rule("(a) Sound type main effect on signed error")

f_sound <- signedError_m ~ stimulusDisparity_m + soundType +
  (1 + stimulusDisparity_m | participantId)

m_sound <- lmer(f_sound, data = dat, REML = TRUE)

anova_sound <- as.data.frame(anova(m_sound))
anova_sound <- data.frame(model = "sound type", term = rownames(anova_sound),
                          anova_sound, row.names = NULL)
cat("Type III F tests (Satterthwaite):\n")
print(anova_sound, digits = 6)

cat("\nRandom effects:\n"); print(VarCorr(m_sound), digits = 4)

# ================================================ (b) EMMs and Tukey contrasts

rule("(b) Estimated marginal means and Tukey contrasts for sound type")

emm_sound <- emmeans(m_sound, ~ soundType)
emm_tab <- as.data.frame(summary(emm_sound, infer = c(TRUE, TRUE)))
# Convert the metre-scale response to centimetres, as reported in the paper.
cm_cols <- c("emmean", "SE", "lower.CL", "upper.CL")
emm_tab[cm_cols] <- emm_tab[cm_cols] * 100
cat("EMMs in cm (Kenward-Roger df):\n"); print(emm_tab, digits = 6)

con_sound <- as.data.frame(summary(contrast(emm_sound, "pairwise", adjust = "tukey"),
                                   infer = c(TRUE, TRUE)))
con_sound[c("estimate", "SE", "lower.CL", "upper.CL")] <-
  con_sound[c("estimate", "SE", "lower.CL", "upper.CL")] * 100
cat("\nTukey pairwise contrasts in cm:\n"); print(con_sound, digits = 6)
cat("\nSmallest Tukey p:", min(con_sound$p.value), "\n")

# ====================================== (c) spatial configuration, signed error

rule("(c) Azimuth sector and elevation category on signed error")

f_spatial <- signedError_m ~ stimulusDisparity_m + azimuthSector +
  elevationCategory + soundType + (1 | participantId)

m_spatial <- lmer(f_spatial, data = dat, REML = TRUE)

anova_spatial <- as.data.frame(anova(m_spatial))
anova_spatial <- data.frame(model = "spatial configuration", term = rownames(anova_spatial),
                            anova_spatial, row.names = NULL)
cat("Type III F tests (Satterthwaite):\n")
print(anova_spatial, digits = 6)

con_azi <- as.data.frame(summary(contrast(emmeans(m_spatial, ~ azimuthSector),
                                          "pairwise", adjust = "tukey"),
                                 infer = c(TRUE, TRUE)))
con_ele <- as.data.frame(summary(contrast(emmeans(m_spatial, ~ elevationCategory),
                                          "pairwise", adjust = "tukey"),
                                 infer = c(TRUE, TRUE)))
for (nm in c("estimate", "SE", "lower.CL", "upper.CL")) {
  con_azi[[nm]] <- con_azi[[nm]] * 100
  con_ele[[nm]] <- con_ele[[nm]] * 100
}
cat("\nAzimuth Tukey contrasts in cm:\n"); print(con_azi, digits = 6)
cat("Smallest azimuth Tukey p:", min(con_azi$p.value), "\n")
cat("\nElevation Tukey contrasts in cm:\n"); print(con_ele, digits = 6)
cat("Smallest elevation Tukey p:", min(con_ele$p.value), "\n")

# ============================================ (d) joint LRT for spatial factors

rule("(d) Joint likelihood ratio test for adding the spatial factors")

# The manuscript places this test in the unsigned-accuracy paragraph, so the
# Gamma GLMM version is the one that carries the published number. Both the
# signed-error and the unsigned-accuracy versions are computed here: the first
# because it belongs with (c), the second to identify which fit produced 2.12.
lrt_row <- function(reduced, full, label) {
  a <- anova(reduced, full)
  data.frame(comparison = label,
             chisq = a$Chisq[2], df = a$Df[2], p = a$`Pr(>Chisq)`[2],
             logLik_reduced = as.numeric(logLik(reduced)),
             logLik_full = as.numeric(logLik(full)),
             stringsAsFactors = FALSE)
}

m_sig_red  <- lme4::lmer(signedError_m ~ stimulusDisparity_m + soundType +
                           (1 | participantId), data = dat, REML = FALSE)
m_sig_full <- lme4::lmer(f_spatial, data = dat, REML = FALSE)
lrt_signed <- lrt_row(m_sig_red, m_sig_full, "signed error LMM, ML")

# Unsigned accuracy. The published value predates the glmmTMB refit, so both
# the original degenerate glmer fit and the trustworthy glmmTMB fit are run.
acc_vars <- c("participantError_cm", "stimDisparity_c", "soundType",
              "trialSequence_c", "azimuthSector", "elevationCategory", "participantId")
acc <- as.data.frame(df[stats::complete.cases(df[, acc_vars]), acc_vars])
acc$participantId <- factor(acc$participantId)
f_acc_red  <- participantError_cm ~ stimDisparity_c + soundType + trialSequence_c +
  (1 | participantId)
f_acc_full <- update(f_acc_red, . ~ . + azimuthSector + elevationCategory)

g_red  <- glmer(f_acc_red,  data = acc, family = Gamma(link = "log"),
                control = glmerControl(optimizer = "bobyqa"))
g_full <- glmer(f_acc_full, data = acc, family = Gamma(link = "log"),
                control = glmerControl(optimizer = "bobyqa"))
lrt_glmer <- lrt_row(g_red, g_full, "unsigned accuracy Gamma GLMM, glmer (degenerate)")

suppressPackageStartupMessages(library(glmmTMB))
t_red  <- glmmTMB(f_acc_red,  data = acc, family = Gamma(link = "log"))
t_full <- glmmTMB(f_acc_full, data = acc, family = Gamma(link = "log"))
a_tmb <- anova(t_red, t_full)
lrt_tmb <- data.frame(comparison = "unsigned accuracy Gamma GLMM, glmmTMB",
                      chisq = a_tmb$Chisq[2], df = a_tmb$`Chi Df`[2],
                      p = a_tmb$`Pr(>Chisq)`[2],
                      logLik_reduced = as.numeric(logLik(t_red)),
                      logLik_full = as.numeric(logLik(t_full)),
                      stringsAsFactors = FALSE)

lrt_tab <- rbind(lrt_signed, lrt_glmer, lrt_tmb)
print(lrt_tab, digits = 6, row.names = FALSE)

# ================================ convergence, conditioning, optimiser stability

rule("Convergence, Hessian conditioning and cross-optimiser stability")

health <- rbind(fit_health(m_sound,   "(a)/(b) sound type"),
                fit_health(m_spatial, "(c) spatial configuration"))
print(health[, c("model", "optimizer", "n_warnings", "max_abs_gradient",
                 "theta_hessian_min_eigen", "vcov_condition_number",
                 "vcov_min_eigen", "is_singular")], digits = 5, row.names = FALSE)
if (any(health$n_warnings > 0)) {
  cat("\nConvergence messages:\n"); print(health[, c("model", "warnings")], row.names = FALSE)
}

opt_tab <- rbind(optimiser_agreement(f_sound,   dat, "(a)/(b) sound type"),
                 optimiser_agreement(f_spatial, dat, "(c) spatial configuration"))
cat("\nFixed-effect standard errors, bobyqa versus nloptwrap:\n")
print(opt_tab, digits = 6, row.names = FALSE)
cat("\nLargest relative SE disagreement across optimisers:",
    signif(max(abs(opt_tab$se_ratio - 1)), 4), "\n")
cat("Largest absolute coefficient disagreement:",
    signif(max(opt_tab$beta_abs_diff), 4), "\n")

# ------------------------------------------------------------------- outputs --

write.csv(rbind(anova_sound, anova_spatial),
          file.path(results_dir, "signed_error_lmm_ftests.csv"), row.names = FALSE)
write.csv(emm_tab,  file.path(results_dir, "signed_error_lmm_soundtype_emmeans_cm.csv"), row.names = FALSE)
write.csv(con_sound, file.path(results_dir, "signed_error_lmm_soundtype_contrasts_cm.csv"), row.names = FALSE)
write.csv(rbind(cbind(factor_tested = "azimuthSector", con_azi),
                cbind(factor_tested = "elevationCategory", con_ele)),
          file.path(results_dir, "signed_error_lmm_spatial_contrasts_cm.csv"), row.names = FALSE)
write.csv(lrt_tab,  file.path(results_dir, "signed_error_lmm_spatial_lrt.csv"), row.names = FALSE)
write.csv(health,   file.path(results_dir, "signed_error_lmm_convergence.csv"), row.names = FALSE)
write.csv(opt_tab,  file.path(results_dir, "signed_error_lmm_optimiser_stability.csv"), row.names = FALSE)

cat("\nSession info:\n")
print(sessionInfo())
cat("\nDone.\n")

sink(type = "message"); sink(); close(log_con)
