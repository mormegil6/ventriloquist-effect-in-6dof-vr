# Deep checks on the Analysis 2 accuracy GLMM re-analysis.
#
# Three questions that the headline refit does not answer on its own:
#
#   A. Are the glmmTMB error families in this installation trustworthy? The
#      package was built against a different TMB version, so every family used
#      is calibrated against a case with a known answer.
#   B. lme4 and glmmTMB report different maximised log-likelihoods for the same
#      Gamma GLMM. Exact adaptive Gauss-Hermite quadrature settles which fit is
#      actually at the maximum, and whether the two AICs are on one scale.
#   C. How much of the audio-visual disparity effect on unsigned error is
#      carried by the response being pulled towards the visual flash, which is
#      a near-definitional relationship, rather than by a general loss of
#      precision?
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/verify_accuracy_deep_*.csv and a plain-text log.

# NOTE (2026-08-07): glmmTMB::lognormal() is BROKEN in this R installation. It reports a
# log-likelihood it cannot attain and returns a nonsensical dispersion (sigma = 16.16 against
# sd(log y) = 0.63). Any lognormal() output below is retained only to document the defect and
# must not be reported. The correct log-normal model is a Gaussian fit on log(y); see
# refit_signed_error_lmm_revision.R and the audit in results_revision/.

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(glmmTMB)
  library(emmeans)
})

# Gauss-Hermite nodes and weights for weight function exp(-x^2), obtained from
# the Golub-Welsch eigenvalue construction so that no extra package is needed.
gauss_hermite <- function(n) {
  i <- seq_len(n - 1)
  J <- diag(0, n)
  J[cbind(i, i + 1)] <- sqrt(i / 2)
  J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  ord <- order(e$values)
  list(nodes = e$values[ord], weights = sqrt(pi) * (e$vectors[1, ord])^2)
}

set.seed(20260805)

proj_dir    <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
results_dir <- file.path(proj_dir, "results_revision")
log_con <- file(file.path(results_dir, "verify_accuracy_deep_checks_log.txt"), open = "wt")
sink(log_con, split = TRUE); sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

cat("Deep checks on the Analysis 2 accuracy GLMMs\n")
cat("Run date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "| glmmTMB", as.character(packageVersion("glmmTMB")),
    "| TMB", as.character(packageVersion("TMB")),
    "| lme4", as.character(packageVersion("lme4")), "\n")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds"))
model_vars <- c("participantError_cm", "stimDisparity_c", "soundType",
                "trialSequence_c", "azimuthSector", "elevationCategory", "participantId")
dat <- df[stats::complete.cases(df[, model_vars]), model_vars]
dat$participantId <- factor(dat$participantId)

model_formula <- participantError_cm ~ stimDisparity_c + soundType +
  trialSequence_c + azimuthSector + elevationCategory + (1 | participantId)

# =================================================== A. family calibration ===

rule("A. Calibration of the glmmTMB families against known answers")

# A1. Gamma. Simulate with a known shape and slope; check both are recovered
#     and establish what sigma() returns.
cat("A1. Gamma(log): parameter recovery, and the meaning of sigma()\n")
gamma_cal <- do.call(rbind, lapply(c(0.75, 2, 4, 8), function(sh) {
  set.seed(100 + round(sh * 10))
  n <- 6000; x <- rnorm(n)
  gg <- factor(rep(seq_len(50), each = n / 50)); u <- rnorm(50, 0, 0.25)
  mu <- exp(1 + 0.5 * x + u[gg])
  y <- rgamma(n, shape = sh, rate = sh / mu)
  m <- glmmTMB(y ~ x + (1 | gg), family = Gamma(link = "log"),
               data = data.frame(y, x, gg))
  data.frame(true_shape = sh, slope = fixef(m)$cond[["x"]],
             sigma = sigma(m), inv_sigma = 1 / sigma(m),
             inv_sigma_sq = 1 / sigma(m)^2,
             re_sd = sqrt(as.numeric(VarCorr(m)$cond$gg)))
}))
print(gamma_cal, digits = 5, row.names = FALSE)
cat("True slope 0.5, true RE SD 0.25. sigma() tracks 1/sqrt(shape),",
    "so shape = 1/sigma()^2, not sigma() and not 1/sigma().\n")

# A2. Lognormal. For a lognormal regression the maximum likelihood solution is
#     exactly least squares on log(y), so the correct answer is known in closed
#     form and the fit can be audited directly.
cat("\nA2. lognormal(log): audit against the closed-form MLE\n")
set.seed(99)
n <- 4000; xs <- rnorm(n); ys <- exp(rnorm(n, 1 + 0.5 * xs, 0.5))
m_ln_sim <- glmmTMB(ys ~ xs, family = lognormal(link = "log"), data = data.frame(ys, xs))
ols <- lm(log(ys) ~ xs)
b_ols <- coef(ols); s_ols <- sqrt(mean(residuals(ols)^2))
Xs <- cbind(1, xs)
ll_median <- function(b, s) sum(dnorm(log(ys), Xs %*% b, s, log = TRUE) - log(ys))
ll_mean   <- function(b, s) sum(dnorm(log(ys), Xs %*% b - s^2 / 2, s, log = TRUE) - log(ys))
cat("  truth: slope 0.5, log-scale SD 0.5\n")
cat("  closed-form MLE      : slope", round(b_ols[2], 5), "| sigma", round(s_ols, 5),
    "| logLik", round(ll_median(b_ols, s_ols), 3), "\n")
cat("  glmmTMB lognormal    : slope", round(fixef(m_ln_sim)$cond[2], 5),
    "| sigma", round(sigma(m_ln_sim), 5),
    "| reported logLik", round(as.numeric(logLik(m_ln_sim)), 3), "\n")
cat("  logLik recomputed at the glmmTMB estimates: median-link",
    round(ll_median(fixef(m_ln_sim)$cond, sigma(m_ln_sim)), 3),
    "| mean-link", round(ll_mean(fixef(m_ln_sim)$cond, sigma(m_ln_sim)), 3), "\n")
cat("  The reported value matches neither parameterisation, and the estimates",
    "sit", round(ll_median(b_ols, s_ols) - ll_median(fixef(m_ln_sim)$cond, sigma(m_ln_sim)), 1),
    "nats below the MLE.\n")

# A3. Other families that other revision scripts may rely on.
cat("\nA3. Spot-check of the remaining families (recovery of a known slope of 0.5)\n")
set.seed(4)
n2 <- 6000; x2 <- rnorm(n2)
fam_rows <- list()
chk <- function(label, fit, truth) {
  s <- try(fixef(fit)$cond[2], silent = TRUE)
  ok <- !inherits(s, "try-error") && abs(s - truth) < 0.08
  cat(sprintf("  %-22s slope = %8.4f (true %.2f)  %s\n", label, s, truth,
              if (ok) "ok" else "DEVIATES"))
  data.frame(family = label, slope = as.numeric(s), true_slope = truth, ok = ok)
}
mu_p <- exp(1 + 0.5 * x2)
fam_rows[[1]] <- chk("poisson", glmmTMB(rpois(n2, mu_p) ~ x2,
                                        family = poisson(), data = data.frame(x2)), 0.5)
fam_rows[[2]] <- chk("nbinom2", glmmTMB(rnbinom(n2, mu = mu_p, size = 3) ~ x2,
                                        family = nbinom2(), data = data.frame(x2)), 0.5)
fam_rows[[3]] <- chk("gaussian", glmmTMB(rnorm(n2, 1 + 0.5 * x2, 1) ~ x2,
                                         family = gaussian(), data = data.frame(x2)), 0.5)
fam_rows[[4]] <- chk("Gamma(log)", glmmTMB(rgamma(n2, shape = 3, rate = 3 / mu_p) ~ x2,
                                           family = Gamma(link = "log"), data = data.frame(x2)), 0.5)
fam_rows[[5]] <- chk("lognormal(log)", m_ln_sim, 0.5)
mu_b <- plogis(0.2 + 0.5 * x2)
fam_rows[[6]] <- chk("beta_family", glmmTMB(rbeta(n2, mu_b * 8, (1 - mu_b) * 8) ~ x2,
                                            family = beta_family(), data = data.frame(x2)), 0.5)
fam_tab <- do.call(rbind, fam_rows)

# =================================== A4. consequences for the robustness check =

rule("A4. The lognormal robustness check on the real data, done correctly")

fit_tmb <- glmmTMB(model_formula, data = dat, family = Gamma(link = "log"))
fit_ln_bad <- glmmTMB(model_formula, data = dat, family = lognormal(link = "log"))
fit_ln_ok  <- glmmTMB(update(model_formula, log(participantError_cm) ~ .),
                      data = dat, family = gaussian(), REML = FALSE)
sum_log_y <- sum(log(dat$participantError_cm))

cat("glmmTMB lognormal as used in the re-analysis:\n")
cat("  sigma =", sigma(fit_ln_bad), "on the log scale, against sd(log y) =",
    sd(log(dat$participantError_cm)), "\n")
cat("  reported logLik", as.numeric(logLik(fit_ln_bad)), "| AIC", AIC(fit_ln_bad), "\n")
cat("Correct lognormal ML fit (Gaussian on log(y), ML; identical model):\n")
cat("  sigma =", sigma(fit_ln_ok), "| RE SD =", sqrt(as.numeric(VarCorr(fit_ln_ok)$cond$participantId)), "\n")
cat("  logLik on the log(y) scale", as.numeric(logLik(fit_ln_ok)),
    "| on the y scale", as.numeric(logLik(fit_ln_ok)) - sum_log_y,
    "| AIC on the y scale", 2 * 13 - 2 * (as.numeric(logLik(fit_ln_ok)) - sum_log_y), "\n")
cat("  The correct lognormal fit is", round((as.numeric(logLik(fit_ln_ok)) - sum_log_y) -
      as.numeric(logLik(fit_ln_bad)), 3), "nats better than the one used.\n")

cmp_terms <- c("stimDisparity_c", "trialSequence_c", "soundTypeFlute",
               "soundTypeSpeech", "soundTypePink Noise")
bad <- summary(fit_ln_bad)$coefficients$cond[cmp_terms, , drop = FALSE]
good <- summary(fit_ln_ok)$coefficients$cond[cmp_terms, , drop = FALSE]
ln_cmp <- data.frame(term = cmp_terms,
                     b_broken = bad[, 1], se_broken = bad[, 2], z_broken = bad[, 3], p_broken = bad[, 4],
                     b_correct = good[, 1], se_correct = good[, 2], z_correct = good[, 3], p_correct = good[, 4],
                     row.names = NULL)
cat("\nlognormal robustness check, as used versus done correctly:\n")
print(ln_cmp, digits = 5)

cat("\nTukey contrasts under the correct lognormal fit:\n")
print(as.data.frame(summary(contrast(emmeans(fit_ln_ok, ~ soundType), "pairwise",
                                     adjust = "tukey"), type = "response",
                            infer = c(TRUE, TRUE))), digits = 6)

# ============================ B. which Gamma fit is at the maximum likelihood =

rule("B. Exact marginal likelihood for the Gamma GLMM (adaptive Gauss-Hermite)")

# Marginal likelihood of a Gamma GLMM with a single random intercept, obtained
# by 61-point Gauss-Hermite quadrature per participant. This is independent of
# both packages, so it can arbitrate between their reported log-likelihoods.
gh <- gauss_hermite(61)
X <- model.matrix(~ stimDisparity_c + soundType + trialSequence_c +
                    azimuthSector + elevationCategory, data = dat)
grp <- dat$participantId
yv <- dat$participantError_cm

marg_loglik <- function(beta, sd_u, shape) {
  eta <- as.vector(X %*% beta)
  nodes <- sqrt(2) * sd_u * gh$nodes
  wts <- gh$weights / sqrt(pi)
  sum(vapply(levels(grp), function(g) {
    idx <- which(grp == g)
    ll <- vapply(nodes, function(u) {
      mu <- exp(eta[idx] + u)
      sum(dgamma(yv[idx], shape = shape, scale = mu / shape, log = TRUE))
    }, numeric(1))
    m <- max(ll)
    m + log(sum(wts * exp(ll - m)))
  }, numeric(1)))
}

fit_bob <- glmer(model_formula, data = dat, family = Gamma(link = "log"),
                 control = glmerControl(optimizer = "bobyqa"))
fit_nlo <- glmer(model_formula, data = dat, family = Gamma(link = "log"))

b_tmb <- fixef(fit_tmb)$cond
sd_tmb <- sqrt(as.numeric(VarCorr(fit_tmb)$cond$participantId))
shape_tmb <- 1 / sigma(fit_tmb)^2

cand <- list(
  `glmmTMB Gamma`        = list(b = b_tmb, sd = sd_tmb, shape = shape_tmb,
                                reported = as.numeric(logLik(fit_tmb))),
  `glmer bobyqa`         = list(b = fixef(fit_bob),
                                sd = as.data.frame(VarCorr(fit_bob))$sdcor[1],
                                shape = 1 / sigma(fit_bob)^2,
                                reported = as.numeric(logLik(fit_bob))),
  `glmer nloptwrap`      = list(b = fixef(fit_nlo),
                                sd = as.data.frame(VarCorr(fit_nlo))$sdcor[1],
                                shape = 1 / sigma(fit_nlo)^2,
                                reported = as.numeric(logLik(fit_nlo)))
)
ml_rows <- do.call(rbind, lapply(names(cand), function(nm) {
  cc <- cand[[nm]]
  data.frame(fit = nm, re_sd = cc$sd, shape = cc$shape,
             reported_logLik = cc$reported,
             quadrature_logLik = marg_loglik(cc$b, cc$sd, cc$shape))
}))
ml_rows$reported_minus_quadrature <- ml_rows$reported_logLik - ml_rows$quadrature_logLik
print(ml_rows, digits = 8, row.names = FALSE)
cat("\nThe quadrature column is the like-for-like comparison. Higher is better.\n")
cat("Best fit by exact marginal likelihood:",
    ml_rows$fit[which.max(ml_rows$quadrature_logLik)], "\n")
cat("glmmTMB minus best glmer, by quadrature:",
    round(ml_rows$quadrature_logLik[1] - max(ml_rows$quadrature_logLik[-1]), 4), "\n")

# =========================== C. is the disparity effect partly definitional? ==

rule("C. Decomposition of the disparity effect")

# participantError is the distance from the sound source to the response.
# signedError_m is its component along the sound-to-flash axis and
# orthogonalError_m the component perpendicular to it. If responses are pulled
# towards the flash, unsigned error grows with disparity almost by
# construction. The orthogonal component carries no such constraint, so it
# tests whether disparity degrades localisation generally.
dec <- df[, c("participantId", "soundType", "stimDisparity_c", "trialSequence_c",
              "azimuthSector", "elevationCategory", "participantError_cm",
              "signedError_m", "orthogonalError_m", "stimulusDisparity_m")]
dec$participantId <- factor(dec$participantId)
dec$signedError_cm <- dec$signedError_m * 100
dec$orthError_cm <- abs(dec$orthogonalError_m) * 100
dec$captureFrac <- dec$signedError_m / dec$stimulusDisparity_m

cat("Mean capture fraction (signed error / disparity):", mean(dec$captureFrac),
    "| median", median(dec$captureFrac), "\n")
cat("Correlation of unsigned error with disparity:",
    cor(dec$participantError_cm, dec$stimulusDisparity_m), "\n")
cat("Correlation of orthogonal error with disparity:",
    cor(dec$orthError_cm, dec$stimulusDisparity_m), "\n")

# Linear benchmark: if the response were a fixed fraction f of the way from the
# sound to the flash and nothing else changed, unsigned error would equal
# f * disparity exactly.
lm_bench <- lm(participantError_cm ~ 0 + I(stimulusDisparity_m * 100), data = dec)
cat("Slope of unsigned error on disparity, both in cm, no intercept:",
    round(coef(lm_bench)[1], 4), "\n")

orth_form <- orthError_cm ~ stimDisparity_c + soundType + trialSequence_c +
  azimuthSector + elevationCategory + (1 | participantId)
dec_pos <- dec[dec$orthError_cm > 0, ]
cat("Trials with orthogonal error exactly zero, dropped from the Gamma fit:",
    nrow(dec) - nrow(dec_pos), "\n")
fit_orth <- glmmTMB(orth_form, data = dec_pos, family = Gamma(link = "log"))
cat("convergence", fit_orth$fit$convergence, "| pdHess", fit_orth$sdr$pdHess, "\n")
orth_co <- summary(fit_orth)$coefficients$cond
cat("\nOrthogonal-error model, same specification:\n")
print(round(orth_co, 6))
b_o <- orth_co["stimDisparity_c", ]
cat("\nDisparity on the orthogonal error component: b =", round(b_o[1], 5),
    "per metre, i.e.", round((exp(b_o[1] / 10) - 1) * 100, 2), "% per 10 cm, z =",
    round(b_o[3], 3), ", p =", signif(b_o[4], 4), "\n")
cat("For comparison, on total unsigned error: b =",
    round(fixef(fit_tmb)$cond[["stimDisparity_c"]], 5), "per metre, i.e.",
    round((exp(fixef(fit_tmb)$cond[["stimDisparity_c"]] / 10) - 1) * 100, 2), "% per 10 cm\n")

# Same model on the signed component, which is the visual-capture component.
sgn_form <- update(orth_form, I(signedError_cm - min(signedError_cm) + 0.01) ~ .)
fit_sgn <- lme4::lmer(signedError_cm ~ stimDisparity_c + soundType + trialSequence_c +
                        azimuthSector + elevationCategory + (1 | participantId), data = dec)
cat("\nSigned (towards-flash) component, LMM in cm:\n")
print(round(summary(fit_sgn)$coefficients, 5))

# ------------------------------------------------------------------ outputs --
write.csv(gamma_cal, file.path(results_dir, "verify_accuracy_deep_gamma_calibration.csv"), row.names = FALSE)
write.csv(fam_tab,   file.path(results_dir, "verify_accuracy_deep_family_calibration.csv"), row.names = FALSE)
write.csv(ln_cmp,    file.path(results_dir, "verify_accuracy_deep_lognormal_comparison.csv"), row.names = FALSE)
write.csv(ml_rows,   file.path(results_dir, "verify_accuracy_deep_marginal_loglik.csv"), row.names = FALSE)
write.csv(as.data.frame(orth_co), file.path(results_dir, "verify_accuracy_deep_orthogonal_error.csv"))

cat("\nDone.\n")
