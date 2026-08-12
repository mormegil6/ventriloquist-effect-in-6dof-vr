#!/usr/bin/env Rscript
# Re-estimation of the movement-accuracy analysis (revision).
#
# The published between- and within-subject movement coefficients came from a
# Gamma glmer that failed to converge with a degenerate Hessian, producing
# implausible z values (-1097 and 80.11). This script:
#   1. reproduces the degenerate glmer fit for the record,
#   2. refits the same within/between decomposition with glmmTMB,
#   3. repeats the fit on the cm-scaled response for comparability with the text,
#   4. triangulates the between-subject effect three ways (GLMM, participant-level
#      correlation, participant-level bootstrap),
#   5. tests whether the within-subject effect survives difficulty covariates,
#   6. runs the likelihood-ratio test for adding the two movement terms.
#
# Deterministic: set.seed() is used for the bootstrap and residual simulations.

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(glmmTMB)
  library(DHARMa)
})

proj <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir <- file.path(proj, "results_revision")
log_path <- file.path(res_dir, "movement_effects_refit_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

cat("Movement effects re-estimation\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "| glmmTMB", as.character(packageVersion("glmmTMB")),
    "| lme4", as.character(packageVersion("lme4")), "\n\n")

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

adf <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
adf <- as.data.frame(adf)
adf$participantId <- factor(adf$participantId)

cat("n trials =", nrow(adf), "| n participants =", nlevels(adf$participantId), "\n")
cat("movement_rate (m/s): mean", round(mean(adf$movement_rate), 4),
    "SD", round(sd(adf$movement_rate), 4),
    "range", paste(round(range(adf$movement_rate), 4), collapse = " to "), "\n")
pm <- adf$participant_mean_rate[!duplicated(adf$participantId)]
cat("participant_mean_rate (m/s, n = 31): mean", round(mean(pm), 4),
    "SD", round(sd(pm), 4),
    "range", paste(round(range(pm), 4), collapse = " to "), "\n")
cat("trial_rate_deviation (m/s): SD", round(sd(adf$trial_rate_deviation), 4), "\n\n")

# Helper: coefficient table with Wald 95% CI and effect per +0.1 m/s.
coef_table <- function(model, terms, model_label, per = 0.1) {
  cf <- if (inherits(model, "glmmTMB")) summary(model)$coefficients$cond else summary(model)$coefficients
  out <- lapply(terms, function(tm) {
    est <- cf[tm, 1]; se <- cf[tm, 2]; z <- cf[tm, 3]; p <- cf[tm, 4]
    lo <- est - 1.96 * se; hi <- est + 1.96 * se
    data.frame(
      model = model_label, term = tm,
      estimate = est, se = se, z = z, p = p, ci_lo = lo, ci_hi = hi,
      pct_per_0.1ms = (exp(per * est) - 1) * 100,
      pct_lo = (exp(per * lo) - 1) * 100,
      pct_hi = (exp(per * hi) - 1) * 100,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

fmt_row <- function(r) {
  sprintf(paste0("  %-22s beta = %8.4f  SE = %7.4f  z = %7.2f  p = %.4f  ",
                 "95%% CI [%.4f, %.4f]\n    -> %+.1f%% error per +0.1 m/s ",
                 "(95%% CI %+.1f%% to %+.1f%%)\n"),
          r$term, r$estimate, r$se, r$z, r$p, r$ci_lo, r$ci_hi,
          r$pct_per_0.1ms, r$pct_lo, r$pct_hi)
}

mv_terms <- c("participant_mean_rate", "trial_rate_deviation")

f_full <- participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum +
  participant_mean_rate + trial_rate_deviation + (1 | participantId)

# ---------------------------------------------------------------------------
# 1. Reproduce the published (degenerate) glmer fit
# ---------------------------------------------------------------------------

cat(strrep("=", 78), "\n1. PUBLISHED FIT REPRODUCED: glmer Gamma(log), bobyqa\n",
    strrep("=", 78), "\n", sep = "")

glmer_warnings <- character(0)
m_glmer <- withCallingHandlers(
  glmer(f_full, data = adf, family = Gamma(link = "log"),
        control = glmerControl(optimizer = "bobyqa")),
  warning = function(w) {
    glmer_warnings <<- c(glmer_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
print(round(summary(m_glmer)$coefficients, 4))
cat("\nglmer warnings captured:\n")
if (length(glmer_warnings)) {
  for (w in glmer_warnings) cat("  - ", w, "\n", sep = "")
} else {
  cat("  (none)\n")
}
cat("Random-intercept SD:", round(as.numeric(VarCorr(m_glmer)$participantId), 6), "\n")
cat("Singular fit:", isSingular(m_glmer), "\n")
cat("logLik:", as.numeric(logLik(m_glmer)), "\n\n")

# ---------------------------------------------------------------------------
# 2. glmmTMB refit, metre-scaled response (primary model)
# ---------------------------------------------------------------------------

cat(strrep("=", 78), "\n2. PRIMARY REFIT: glmmTMB Gamma(log), response in metres\n",
    strrep("=", 78), "\n", sep = "")

m_tmb <- glmmTMB(f_full, data = adf, family = Gamma(link = "log"))
print(summary(m_tmb))

cat("\nConvergence code:", m_tmb$fit$convergence,
    "| message:", m_tmb$fit$message,
    "| positive-definite Hessian:", m_tmb$sdr$pdHess, "\n")
cat("Max absolute gradient at the optimum:",
    max(abs(m_tmb$obj$gr(m_tmb$fit$par))), "\n\n")

tab_m <- coef_table(m_tmb, mv_terms, "glmmTMB_metres")
cat("Movement coefficients (metre scale):\n")
for (i in seq_len(nrow(tab_m))) cat(fmt_row(tab_m[i, ]))

# Between-subject effect on the observed spread of participant mean rates.
b_between <- tab_m$estimate[tab_m$term == "participant_mean_rate"]
sd_pm <- sd(pm)
cat(sprintf("\n  Between-subject effect per +1 SD of participant mean rate (%.4f m/s): %+.1f%%\n",
            sd_pm, (exp(b_between * sd_pm) - 1) * 100))
b_within <- tab_m$estimate[tab_m$term == "trial_rate_deviation"]
sd_dev <- sd(adf$trial_rate_deviation)
cat(sprintf("  Within-subject effect per +1 SD of trial rate deviation (%.4f m/s): %+.1f%%\n\n",
            sd_dev, (exp(b_within * sd_dev) - 1) * 100))

# Residual diagnostics for the primary model.
set.seed(20260805)
sim_res <- simulateResiduals(m_tmb, n = 1000)
cat("DHARMa uniformity KS p =", signif(testUniformity(sim_res, plot = FALSE)$p.value, 4),
    "| dispersion p =", signif(testDispersion(sim_res, plot = FALSE)$p.value, 4), "\n\n")

# ---------------------------------------------------------------------------
# 3. cm-scaled response, for comparability with the published text
# ---------------------------------------------------------------------------

cat(strrep("=", 78), "\n3. cm-SCALED RESPONSE (glmmTMB Gamma(log))\n",
    strrep("=", 78), "\n", sep = "")

f_cm <- update(f_full, participantError_cm ~ .)
m_tmb_cm <- glmmTMB(f_cm, data = adf, family = Gamma(link = "log"))
cat("Convergence code:", m_tmb_cm$fit$convergence,
    "| positive-definite Hessian:", m_tmb_cm$sdr$pdHess, "\n")
tab_cm <- coef_table(m_tmb_cm, mv_terms, "glmmTMB_cm")
for (i in seq_len(nrow(tab_cm))) cat(fmt_row(tab_cm[i, ]))
cat(sprintf("\n  Intercept shift m -> cm: %.4f (log(100) = %.4f)\n",
            fixef(m_tmb_cm)$cond[["(Intercept)"]] - fixef(m_tmb)$cond[["(Intercept)"]],
            log(100)))
cat("  Max abs difference in movement slopes between scales:",
    signif(max(abs(tab_cm$estimate - tab_m$estimate)), 3), "\n\n")

# ---------------------------------------------------------------------------
# 4. Triangulating the between-subject effect (31 clusters)
# ---------------------------------------------------------------------------

cat(strrep("=", 78), "\n4. TRIANGULATION OF THE BETWEEN-SUBJECT EFFECT (n = 31)\n",
    strrep("=", 78), "\n", sep = "")

pdat <- adf %>%
  group_by(participantId) %>%
  summarise(
    mean_rate = mean(movement_rate),
    mean_error_m = mean(participantError_m),
    median_error_m = median(participantError_m),
    mean_error_cm = mean(participantError_cm),
    .groups = "drop"
  )
cat("Participant-level table: n =", nrow(pdat), "\n\n")

ct_p <- cor.test(pdat$mean_rate, pdat$mean_error_m, method = "pearson")
ct_s <- suppressWarnings(cor.test(pdat$mean_rate, pdat$mean_error_m, method = "spearman"))
ct_p_med <- cor.test(pdat$mean_rate, pdat$median_error_m, method = "pearson")

cat(sprintf("(2a) Pearson  r = %+.4f, 95%% CI [%+.4f, %+.4f], t(%d) = %.3f, p = %.4f\n",
            ct_p$estimate, ct_p$conf.int[1], ct_p$conf.int[2],
            ct_p$parameter, ct_p$statistic, ct_p$p.value))
cat(sprintf("(2b) Spearman rho = %+.4f, S = %.1f, p = %.4f\n",
            ct_s$estimate, ct_s$statistic, ct_s$p.value))
cat(sprintf("     (sensitivity: Pearson vs participant MEDIAN error r = %+.4f, p = %.4f)\n\n",
            ct_p_med$estimate, ct_p_med$p.value))

# Bootstrap over participants.
set.seed(20260805)
B <- 5000
idx <- seq_len(nrow(pdat))
boot <- replicate(B, {
  s <- sample(idx, replace = TRUE)
  x <- pdat$mean_rate[s]; y <- pdat$mean_error_m[s]
  if (sd(x) == 0 || sd(y) == 0) return(rep(NA_real_, 4))
  c(cor(x, y),
    suppressWarnings(cor(x, y, method = "spearman")),
    unname(coef(lm(y ~ x))[2]),
    unname(coef(lm(log(y) ~ x))[2]))
})
boot_r <- boot[1, ]; boot_rho <- boot[2, ]
boot_slope <- boot[3, ]; boot_logslope <- boot[4, ]
cat("Bootstrap: B =", B, "resamples of participants;",
    sum(is.na(boot_r)), "degenerate resamples dropped\n")
q_r <- quantile(boot_r, c(.025, .975), na.rm = TRUE)
q_rho <- quantile(boot_rho, c(.025, .975), na.rm = TRUE)
q_sl <- quantile(boot_slope, c(.025, .975), na.rm = TRUE)
cat(sprintf("(3a) Pearson r  = %+.4f, bootstrap percentile 95%% CI [%+.4f, %+.4f]\n",
            ct_p$estimate, q_r[1], q_r[2]))
cat(sprintf("(3b) Spearman rho = %+.4f, bootstrap percentile 95%% CI [%+.4f, %+.4f]\n",
            ct_s$estimate, q_rho[1], q_rho[2]))
obs_slope <- unname(coef(lm(mean_error_m ~ mean_rate, data = pdat))[2])
cat(sprintf("(3c) OLS slope = %+.4f m error per 1 m/s, bootstrap 95%% CI [%+.4f, %+.4f]\n",
            obs_slope, q_sl[1], q_sl[2]))
cat(sprintf("     -> %+.2f cm error per +0.1 m/s (95%% CI %+.2f to %+.2f cm)\n",
            obs_slope * 0.1 * 100, q_sl[1] * 0.1 * 100, q_sl[2] * 0.1 * 100))
cat(sprintf("     bootstrap share of resamples with r < 0: %.3f\n",
            mean(boot_r < 0, na.rm = TRUE)))

# Log-scale participant-level slope: directly comparable to the GLMM
# between-subject coefficient, which is also on the log-error scale.
q_ls <- quantile(boot_logslope, c(.025, .975), na.rm = TRUE)
obs_logslope <- unname(coef(lm(log(mean_error_m) ~ mean_rate, data = pdat))[2])
cat(sprintf("(3d) Participant-level slope on log error = %+.4f, bootstrap 95%% CI [%+.4f, %+.4f]\n",
            obs_logslope, q_ls[1], q_ls[2]))
cat(sprintf("     -> %+.1f%% error per +0.1 m/s (95%% CI %+.1f%% to %+.1f%%); GLMM gave %+.4f\n\n",
            (exp(0.1 * obs_logslope) - 1) * 100,
            (exp(0.1 * q_ls[1]) - 1) * 100, (exp(0.1 * q_ls[2]) - 1) * 100,
            b_between))

# ---------------------------------------------------------------------------
# 5. Does the within-subject effect survive difficulty covariates?
# ---------------------------------------------------------------------------

cat(strrep("=", 78), "\n5. WITHIN-SUBJECT EFFECT WITH DIFFICULTY COVARIATES\n",
    strrep("=", 78), "\n", sep = "")

f_ang <- update(f_full, . ~ . + angDisp_mean_deg)
m_ang <- glmmTMB(f_ang, data = adf, family = Gamma(link = "log"))
cat("Model A: + angDisp_mean_deg | convergence", m_ang$fit$convergence,
    "| pdHess", m_ang$sdr$pdHess, "\n")
tab_ang <- coef_table(m_ang, mv_terms, "glmmTMB_plus_angDisp")
for (i in seq_len(nrow(tab_ang))) cat(fmt_row(tab_ang[i, ]))
cf_ang <- summary(m_ang)$coefficients$cond
cat(sprintf("  angDisp_mean_deg: beta = %.5f, SE = %.5f, z = %.2f, p = %.4g\n\n",
            cf_ang["angDisp_mean_deg", 1], cf_ang["angDisp_mean_deg", 2],
            cf_ang["angDisp_mean_deg", 3], cf_ang["angDisp_mean_deg", 4]))

# Response time is the other obvious difficulty proxy: slow trials are hard
# trials, and movement rate is path length divided by duration.
f_rt <- update(f_ang, . ~ . + responseTime_s)
m_rt <- glmmTMB(f_rt, data = adf, family = Gamma(link = "log"))
cat("Model B: + angDisp_mean_deg + responseTime_s | convergence", m_rt$fit$convergence,
    "| pdHess", m_rt$sdr$pdHess, "\n")
tab_rt <- coef_table(m_rt, mv_terms, "glmmTMB_plus_angDisp_rt")
for (i in seq_len(nrow(tab_rt))) cat(fmt_row(tab_rt[i, ]))
cf_rt <- summary(m_rt)$coefficients$cond
cat(sprintf("  responseTime_s: beta = %.5f, SE = %.5f, z = %.2f, p = %.4g\n",
            cf_rt["responseTime_s", 1], cf_rt["responseTime_s", 2],
            cf_rt["responseTime_s", 3], cf_rt["responseTime_s", 4]))
cat(sprintf("\n  trial_rate_deviation across specifications: base %.4f -> +angDisp %.4f (%+.1f%%) -> +angDisp+RT %.4f (%+.1f%%)\n\n",
            b_within,
            tab_ang$estimate[tab_ang$term == "trial_rate_deviation"],
            100 * (tab_ang$estimate[tab_ang$term == "trial_rate_deviation"] / b_within - 1),
            tab_rt$estimate[tab_rt$term == "trial_rate_deviation"],
            100 * (tab_rt$estimate[tab_rt$term == "trial_rate_deviation"] / b_within - 1)))

# ---------------------------------------------------------------------------
# 6. Likelihood-ratio test for the two movement terms
# ---------------------------------------------------------------------------

cat(strrep("=", 78), "\n6. LIKELIHOOD-RATIO TEST: base vs base + 2 movement terms\n",
    strrep("=", 78), "\n", sep = "")

f_base <- participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum +
  (1 | participantId)
m_base_tmb <- glmmTMB(f_base, data = adf, family = Gamma(link = "log"))
cat("Base model convergence:", m_base_tmb$fit$convergence,
    "| pdHess:", m_base_tmb$sdr$pdHess, "\n")
lrt <- anova(m_base_tmb, m_tmb)
print(lrt)
chi <- lrt$Chisq[2]; df_lrt <- lrt$`Chi Df`[2]; p_lrt <- lrt$`Pr(>Chisq)`[2]
cat(sprintf("\nglmmTMB LRT: chi2(%d) = %.3f, p = %.4g  (logLik %.3f -> %.3f, increase %.3f)\n",
            df_lrt, chi, p_lrt, as.numeric(logLik(m_base_tmb)),
            as.numeric(logLik(m_tmb)),
            as.numeric(logLik(m_tmb)) - as.numeric(logLik(m_base_tmb))))

# The published LRT used glmer; reproduce it to document the degeneracy.
m_base_glmer <- suppressWarnings(
  glmer(f_base, data = adf, family = Gamma(link = "log"),
        control = glmerControl(optimizer = "bobyqa"))
)
lrt_glmer <- suppressWarnings(anova(m_base_glmer, m_glmer))
cat("\nOriginal-style glmer LRT for the record:\n")
print(lrt_glmer)
cat(sprintf("glmer logLik: base %.3f -> movement %.3f (change %.3f)\n\n",
            as.numeric(logLik(m_base_glmer)), as.numeric(logLik(m_glmer)),
            as.numeric(logLik(m_glmer)) - as.numeric(logLik(m_base_glmer))))

# Separate single-term LRTs so each movement effect gets its own test.
m_no_between <- glmmTMB(update(f_full, . ~ . - participant_mean_rate),
                        data = adf, family = Gamma(link = "log"))
m_no_within <- glmmTMB(update(f_full, . ~ . - trial_rate_deviation),
                       data = adf, family = Gamma(link = "log"))
lrt_between <- anova(m_no_between, m_tmb)
lrt_within <- anova(m_no_within, m_tmb)
cat(sprintf("LRT participant_mean_rate: chi2(%d) = %.3f, p = %.4g\n",
            lrt_between$`Chi Df`[2], lrt_between$Chisq[2], lrt_between$`Pr(>Chisq)`[2]))
cat(sprintf("LRT trial_rate_deviation:  chi2(%d) = %.3f, p = %.4g\n\n",
            lrt_within$`Chi Df`[2], lrt_within$Chisq[2], lrt_within$`Pr(>Chisq)`[2]))

# ---------------------------------------------------------------------------
# 7. The undecomposed (pooled) movement_rate model
#    This is the model that produced the degenerate LRT reported in the paper
#    (chi2 = 0, p = 1, with the log-likelihood decreasing when a term was added)
#    and the non-significant beta = -0.443, SE = 0.406, z = -1.09.
# ---------------------------------------------------------------------------

cat(strrep("=", 78), "\n7. UNDECOMPOSED movement_rate MODEL (source of the chi2 = 0 LRT)\n",
    strrep("=", 78), "\n", sep = "")

f_pooled <- update(f_base, . ~ . + movement_rate)
m_pool_glmer <- suppressWarnings(
  glmer(f_pooled, data = adf, family = Gamma(link = "log"),
        control = glmerControl(optimizer = "bobyqa"))
)
cat("glmer pooled model, movement_rate row:\n")
print(round(summary(m_pool_glmer)$coefficients["movement_rate", , drop = FALSE], 4))
lrt_pool_glmer <- suppressWarnings(anova(m_base_glmer, m_pool_glmer))
print(lrt_pool_glmer)
cat(sprintf("glmer logLik: base %.3f -> pooled %.3f (change %.3f; a nested model cannot lose logLik)\n\n",
            as.numeric(logLik(m_base_glmer)), as.numeric(logLik(m_pool_glmer)),
            as.numeric(logLik(m_pool_glmer)) - as.numeric(logLik(m_base_glmer))))

m_pool_tmb <- glmmTMB(f_pooled, data = adf, family = Gamma(link = "log"))
cat("glmmTMB pooled model | convergence", m_pool_tmb$fit$convergence,
    "| pdHess", m_pool_tmb$sdr$pdHess, "\n")
tab_pool <- coef_table(m_pool_tmb, "movement_rate", "glmmTMB_pooled_movement_rate")
cat(fmt_row(tab_pool[1, ]))
lrt_pool <- anova(m_base_tmb, m_pool_tmb)
cat(sprintf("glmmTMB LRT (1 df): chi2 = %.3f, p = %.4g (logLik %.3f -> %.3f)\n\n",
            lrt_pool$Chisq[2], lrt_pool$`Pr(>Chisq)`[2],
            as.numeric(logLik(m_base_tmb)), as.numeric(logLik(m_pool_tmb))))

# ---------------------------------------------------------------------------
# 8. What is the between-subject effect actually about?
#    movement_rate is path length divided by trajectory duration, so a high rate
#    can mean more movement or a shorter trial. The paper's claim is worded as
#    "moved more", which is path length, not rate.
# ---------------------------------------------------------------------------

cat(strrep("=", 78), "\n8. RATE VERSUS AMOUNT AT THE PARTICIPANT LEVEL\n",
    strrep("=", 78), "\n", sep = "")

pdat2 <- adf %>%
  group_by(participantId) %>%
  summarise(
    mean_rate = mean(movement_rate),
    mean_path = mean(total_path_length),
    mean_duration = mean(trajectory_duration),
    mean_rt = mean(responseTime_s),
    mean_error_m = mean(participantError_m),
    .groups = "drop"
  )

corr_row <- function(x, xlab) {
  ct <- cor.test(pdat2[[x]], pdat2$mean_error_m)
  cat(sprintf("  %-28s r = %+.3f, 95%% CI [%+.3f, %+.3f], p = %.4g\n",
              xlab, ct$estimate, ct$conf.int[1], ct$conf.int[2], ct$p.value))
  data.frame(predictor = x, label = xlab, r = unname(ct$estimate),
             ci_lo = ct$conf.int[1], ci_hi = ct$conf.int[2], p_value = ct$p.value,
             stringsAsFactors = FALSE)
}
cat("Participant-level correlations with mean unsigned error (n = 31):\n")
rate_amount <- rbind(
  corr_row("mean_rate", "mean movement rate (m/s)"),
  corr_row("mean_path", "mean path length (m)"),
  corr_row("mean_duration", "mean trajectory duration (s)"),
  corr_row("mean_rt", "mean response time (s)")
)
cat(sprintf("\n  Correlation of mean rate with mean path length: r = %+.3f\n",
            cor(pdat2$mean_rate, pdat2$mean_path)))
cat(sprintf("  Correlation of mean rate with mean duration:    r = %+.3f\n\n",
            cor(pdat2$mean_rate, pdat2$mean_duration)))

# Rate is path over duration, so decompose the log-error relationship into the
# numerator and the denominator to see which side carries the association.
m_decomp <- lm(log(mean_error_m) ~ log(mean_path) + log(mean_duration), data = pdat2)
cat("Participant-level decomposition, lm(log mean error ~ log path + log duration):\n")
print(round(summary(m_decomp)$coefficients, 4))
cat(sprintf("  adjusted R2 = %.3f\n", summary(m_decomp)$adj.r.squared))
# If the two coefficients are equal and opposite, only the ratio (the rate) matters.
cc <- coef(m_decomp)[c("log(mean_path)", "log(mean_duration)")]
vv <- vcov(m_decomp)[c("log(mean_path)", "log(mean_duration)"),
                     c("log(mean_path)", "log(mean_duration)")]
s <- sum(cc); s_se <- sqrt(sum(vv))
cat(sprintf("  Wald test that the two coefficients are equal and opposite: sum = %.4f, SE = %.4f, t(%d) = %.3f, p = %.3f\n\n",
            s, s_se, df.residual(m_decomp), s / s_se,
            2 * pt(-abs(s / s_se), df.residual(m_decomp))))

# ---------------------------------------------------------------------------
# Machine-readable output
# ---------------------------------------------------------------------------

coefs_out <- rbind(tab_m, tab_cm, tab_ang, tab_rt, tab_pool)
glmer_tab <- summary(m_glmer)$coefficients[mv_terms, , drop = FALSE]
coefs_out <- rbind(
  data.frame(
    model = "glmer_published_degenerate", term = mv_terms,
    estimate = glmer_tab[, 1], se = glmer_tab[, 2], z = glmer_tab[, 3], p = glmer_tab[, 4],
    ci_lo = glmer_tab[, 1] - 1.96 * glmer_tab[, 2],
    ci_hi = glmer_tab[, 1] + 1.96 * glmer_tab[, 2],
    pct_per_0.1ms = (exp(0.1 * glmer_tab[, 1]) - 1) * 100,
    pct_lo = (exp(0.1 * (glmer_tab[, 1] - 1.96 * glmer_tab[, 2])) - 1) * 100,
    pct_hi = (exp(0.1 * (glmer_tab[, 1] + 1.96 * glmer_tab[, 2])) - 1) * 100,
    stringsAsFactors = FALSE
  ),
  coefs_out
)
rownames(coefs_out) <- NULL
write.csv(coefs_out, file.path(res_dir, "movement_effects_coefficients.csv"), row.names = FALSE)

triangulation <- data.frame(
  approach = c("GLMM between-subject (glmmTMB)",
               "Participant-level Pearson",
               "Participant-level Spearman",
               "Participant-level Pearson (bootstrap)",
               "Participant-level Spearman (bootstrap)",
               "Participant-level OLS slope (bootstrap)",
               "Participant-level log-error slope (bootstrap)"),
  statistic = c("beta (log scale, per 1 m/s)", "r", "rho", "r", "rho",
                "slope (m error per 1 m/s)", "slope (log error per 1 m/s)"),
  value = c(b_between, ct_p$estimate, ct_s$estimate, ct_p$estimate,
            ct_s$estimate, obs_slope, obs_logslope),
  ci_lo = c(tab_m$ci_lo[tab_m$term == "participant_mean_rate"],
            ct_p$conf.int[1], NA, q_r[1], q_rho[1], q_sl[1], q_ls[1]),
  ci_hi = c(tab_m$ci_hi[tab_m$term == "participant_mean_rate"],
            ct_p$conf.int[2], NA, q_r[2], q_rho[2], q_sl[2], q_ls[2]),
  p_value = c(summary(m_tmb)$coefficients$cond["participant_mean_rate", 4],
              ct_p$p.value, ct_s$p.value, NA, NA, NA, NA),
  n = c(nrow(adf), rep(nrow(pdat), 6)),
  stringsAsFactors = FALSE
)
write.csv(triangulation, file.path(res_dir, "movement_effects_triangulation.csv"), row.names = FALSE)

lrt_out <- data.frame(
  test = c("glmmTMB: +2 movement terms", "glmer (published style): +2 movement terms",
           "glmmTMB: participant_mean_rate", "glmmTMB: trial_rate_deviation",
           "glmmTMB: pooled movement_rate", "glmer (published): pooled movement_rate"),
  chisq = c(chi, lrt_glmer$Chisq[2], lrt_between$Chisq[2], lrt_within$Chisq[2],
            lrt_pool$Chisq[2], lrt_pool_glmer$Chisq[2]),
  df = c(df_lrt, lrt_glmer$Df[2], lrt_between$`Chi Df`[2], lrt_within$`Chi Df`[2],
         lrt_pool$`Chi Df`[2], lrt_pool_glmer$Df[2]),
  p_value = c(p_lrt, lrt_glmer$`Pr(>Chisq)`[2], lrt_between$`Pr(>Chisq)`[2],
              lrt_within$`Pr(>Chisq)`[2], lrt_pool$`Pr(>Chisq)`[2],
              lrt_pool_glmer$`Pr(>Chisq)`[2]),
  logLik_reduced = c(as.numeric(logLik(m_base_tmb)), as.numeric(logLik(m_base_glmer)),
                     as.numeric(logLik(m_no_between)), as.numeric(logLik(m_no_within)),
                     as.numeric(logLik(m_base_tmb)), as.numeric(logLik(m_base_glmer))),
  logLik_full = c(as.numeric(logLik(m_tmb)), as.numeric(logLik(m_glmer)),
                  as.numeric(logLik(m_tmb)), as.numeric(logLik(m_tmb)),
                  as.numeric(logLik(m_pool_tmb)), as.numeric(logLik(m_pool_glmer))),
  stringsAsFactors = FALSE
)
write.csv(lrt_out, file.path(res_dir, "movement_effects_lrt.csv"), row.names = FALSE)

write.csv(pdat2, file.path(res_dir, "movement_effects_participant_level.csv"), row.names = FALSE)
write.csv(rate_amount, file.path(res_dir, "movement_effects_rate_vs_amount.csv"), row.names = FALSE)

cat("Written:\n")
cat("  movement_effects_coefficients.csv\n  movement_effects_triangulation.csv\n")
cat("  movement_effects_lrt.csv\n  movement_effects_participant_level.csv\n")
cat("  movement_effects_rate_vs_amount.csv\n")
cat("  movement_effects_refit_log.txt\n")
