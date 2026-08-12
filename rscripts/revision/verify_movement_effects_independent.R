# Independent verification of the movement-rate refit reported for the TVCG revision.
#
# Purpose: re-derive every headline number of the "movement-refit" analysis from the
# prepared trial-level dataset, without reusing the original analysis script, and add
# checks the original did not perform: cross-optimiser and cross-package agreement of
# the variance-covariance matrix, an lme4 logLik comparability check, small-cluster
# inference for the between-participant term, informativeness of the within-participant
# null, and leverage of individual participants on the participant-level correlation.
#
# Run with R 4.4 (arm64). Outputs go to results_revision/.

# NOTE (2026-08-07): glmmTMB::lognormal() is BROKEN in this R installation. It reports a
# log-likelihood it cannot attain and returns a nonsensical dispersion (sigma = 16.16 against
# sd(log y) = 0.63). Any lognormal() output below is retained only to document the defect and
# must not be reported. The correct log-normal model is a Gaussian fit on log(y); see
# refit_signed_error_lmm_revision.R and the audit in results_revision/.

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(glmmTMB)
  library(DHARMa)
  library(dplyr)
})

set.seed(20260805)

res_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations/results_revision"
d <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
d <- as.data.frame(d)

log_path <- file.path(res_dir, "movement_effects_verification_log.txt")
con <- file(log_path, open = "wt")
sink(con, split = TRUE)
sink(con, type = "message")

hdr <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

hdr("0. DATA INTEGRITY")
cat("trials:", nrow(d), " participants:", length(unique(d$participantId)), "\n")
cat("complete cases on model variables:",
    sum(complete.cases(d[, c("participantError_m", "stimulusDisparity_m", "soundType",
                             "trialSequenceNum", "participant_mean_rate",
                             "trial_rate_deviation", "angDisp_mean_deg", "responseTime_s")])), "\n")
cat("all errors strictly positive:", all(d$participantError_m > 0), "\n")
chk <- d %>% group_by(participantId) %>%
  mutate(pm = mean(movement_rate), dev = movement_rate - pm) %>% ungroup()
cat("max |participant_mean_rate - recomputed group mean|:", max(abs(chk$pm - chk$participant_mean_rate)), "\n")
cat("max |trial_rate_deviation - recomputed deviation| :", max(abs(chk$dev - chk$trial_rate_deviation)), "\n")
cat("max |movement_rate - path/duration|:", max(abs(d$movement_rate - d$total_path_length / d$trajectory_duration)), "\n")
cat("trials per participant: min", min(table(d$participantId)), " max", max(table(d$participantId)), "\n")

pm_sd  <- sd(tapply(d$movement_rate, d$participantId, mean))
dev_sd <- sd(d$trial_rate_deviation)
cat("SD of participant mean rate (between, n=31):", round(pm_sd, 5), "m/s\n")
cat("SD of trial rate deviation  (within, n=744):", round(dev_sd, 5), "m/s\n")
cat("range of participant mean rate:", round(range(tapply(d$movement_rate, d$participantId, mean)), 4), "\n")

f_primary <- participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum +
  participant_mean_rate + trial_rate_deviation + (1 | participantId)
f_base <- participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum + (1 | participantId)
f_pooled <- participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum +
  movement_rate + (1 | participantId)

# ---------------------------------------------------------------- 1. glmer reproduction
hdr("1. REPRODUCTION OF THE PUBLISHED (DEGENERATE) glmer FIT")
w <- NULL
m_glmer <- withCallingHandlers(
  glmer(f_primary, data = d, family = Gamma(link = "log"),
        control = glmerControl(optimizer = "bobyqa")),
  warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
)
cat("warnings emitted:\n"); if (length(w)) cat(paste0("  ", w, collapse = "\n"), "\n") else cat("  none\n")
print(round(summary(m_glmer)$coefficients, 6))
cat("isSingular:", isSingular(m_glmer), "\n")
cat("range of fixed-effect SEs:", range(sqrt(diag(vcov(m_glmer)))), "\n")
cat("logLik:", as.numeric(logLik(m_glmer)), "\n")

# Variance vs SD. as.numeric(VarCorr(fit)$participantId) is the VARIANCE, not the SD;
# reading it as an SD makes an ordinary variance component look like a collapsed one.
vc <- as.data.frame(VarCorr(m_glmer))
cat("\nrandom-intercept VARIANCE (glmer):", vc$vcov[vc$grp == "participantId"], "\n")
cat("random-intercept SD       (glmer):", vc$sdcor[vc$grp == "participantId"], "\n")
cat("residual sigma            (glmer):", sigma(m_glmer), "\n")

# Does a different lme4 optimiser give the same answer? A well-conditioned optimum should.
opt_tab <- do.call(rbind, lapply(
  list(list(lab = "bobyqa",      f = function() m_glmer),
       list(lab = "Nelder_Mead", f = function() glmer(f_primary, data = d, family = Gamma(link = "log"),
              control = glmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 1e5)))),
       list(lab = "nAGQ=0",      f = function() glmer(f_primary, data = d, family = Gamma(link = "log"),
              nAGQ = 0, control = glmerControl(optimizer = "bobyqa")))),
  function(o) {
    fit <- suppressWarnings(o$f())
    v <- as.data.frame(VarCorr(fit))
    data.frame(optimizer = o$lab,
               logLik = as.numeric(logLik(fit)),
               RE_SD = v$sdcor[v$grp == "participantId"],
               beta_between = unname(fixef(fit)["participant_mean_rate"]),
               SE_between = unname(sqrt(diag(vcov(fit)))["participant_mean_rate"]),
               beta_within = unname(fixef(fit)["trial_rate_deviation"]),
               SE_within = unname(sqrt(diag(vcov(fit)))["trial_rate_deviation"]))
  }))
cat("\nlme4 across optimisers (same data, same formula):\n")
print(opt_tab, row.names = FALSE, digits = 6)

# ------------------------------------------------------------- 2. glmmTMB primary refit
hdr("2. glmmTMB PRIMARY REFIT")
m <- glmmTMB(f_primary, data = d, family = Gamma(link = "log"), REML = FALSE)
s <- summary(m)
print(round(s$coefficients$cond, 6))
cat("\nconvergence code:", m$fit$convergence, " message:", m$fit$message, "\n")
cat("positive-definite Hessian:", m$sdr$pdHess, "\n")
cat("max abs gradient (sdr$gradient.fixed):", max(abs(m$sdr$gradient.fixed)), "\n")
cat("random-intercept SD:", attr(VarCorr(m)$cond$participantId, "stddev"),
    " variance:", VarCorr(m)$cond$participantId[1], "\n")
cat("dispersion sigma():", sigma(m), " (glmer:", sigma(m_glmer), ")\n")
cat("logLik:", as.numeric(logLik(m)), " df:", attr(logLik(m), "df"), "\n")
cat("vcov condition number (fixed effects):", kappa(vcov(m)$cond, exact = TRUE), "\n")

ci <- confint(m, parm = "beta_", method = "wald")
print(round(ci, 6))

cf <- s$coefficients$cond
b_btw <- cf["participant_mean_rate", 1]; se_btw <- cf["participant_mean_rate", 2]
b_wth <- cf["trial_rate_deviation", 1];  se_wth <- cf["trial_rate_deviation", 2]

cat("\nInterpretable scale:\n")
pct <- function(b, mult) 100 * (exp(b * mult) - 1)
cat("  between, per +0.1 m/s:", round(pct(b_btw, .1), 3), "% CI",
    round(pct(b_btw - 1.96 * se_btw, .1), 3), "to", round(pct(b_btw + 1.96 * se_btw, .1), 3), "\n")
cat("  between, per +1 SD (", round(pm_sd, 4), " m/s):", round(pct(b_btw, pm_sd), 3), "% CI",
    round(pct(b_btw - 1.96 * se_btw, pm_sd), 3), "to", round(pct(b_btw + 1.96 * se_btw, pm_sd), 3), "\n")
cat("  within,  per +0.1 m/s:", round(pct(b_wth, .1), 3), "% CI",
    round(pct(b_wth - 1.96 * se_wth, .1), 3), "to", round(pct(b_wth + 1.96 * se_wth, .1), 3), "\n")
cat("  within,  per +1 SD (", round(dev_sd, 4), " m/s):", round(pct(b_wth, dev_sd), 3), "% CI",
    round(pct(b_wth - 1.96 * se_wth, dev_sd), 3), "to", round(pct(b_wth + 1.96 * se_wth, dev_sd), 3), "\n")
cat("  0.1 m/s expressed in between-SD units:", round(0.1 / pm_sd, 3), "SD\n")
cat("  0.1 m/s expressed in within-SD units :", round(0.1 / dev_sd, 3), "SD\n")

# --------------------------------------------------- 3. cross-optimiser / cross-family
hdr("3. CROSS-OPTIMISER AND CROSS-FAMILY STABILITY OF THE glmmTMB FIT")
m_bfgs <- glmmTMB(f_primary, data = d, family = Gamma(link = "log"), REML = FALSE,
                  control = glmmTMBControl(optimizer = optim,
                                           optArgs = list(method = "BFGS")))
cat("BFGS: beta_between", fixef(m_bfgs)$cond["participant_mean_rate"],
    " SE", sqrt(diag(vcov(m_bfgs)$cond))["participant_mean_rate"],
    " RE SD", attr(VarCorr(m_bfgs)$cond$participantId, "stddev"),
    " logLik", as.numeric(logLik(m_bfgs)), "\n")

m_reml <- glmmTMB(f_primary, data = d, family = Gamma(link = "log"), REML = TRUE)
cat("REML: beta_between", fixef(m_reml)$cond["participant_mean_rate"],
    " SE", sqrt(diag(vcov(m_reml)$cond))["participant_mean_rate"], "\n")

# Family robustness. The Gamma log link models log E[Y]; a log-normal LMM models E[log Y].
# If the Gamma assumption holds these give the same slopes, so disagreement is diagnostic.
m_ln <- lmer(log(participantError_m) ~ stimulusDisparity_m + soundType + trialSequenceNum +
               participant_mean_rate + trial_rate_deviation + (1 | participantId),
             data = d, REML = TRUE)
cat("\nlog-normal LMM on log(error), lmerTest Satterthwaite (REML):\n")
print(round(summary(m_ln)$coefficients[c("participant_mean_rate", "trial_rate_deviation"), ], 5))
m_ln_ml <- lmer(log(participantError_m) ~ stimulusDisparity_m + soundType + trialSequenceNum +
                  participant_mean_rate + trial_rate_deviation + (1 | participantId),
                data = d, REML = FALSE)
cat("same, ML:\n")
print(round(summary(m_ln_ml)$coefficients[c("participant_mean_rate", "trial_rate_deviation"), ], 5))
m_lntmb <- glmmTMB(f_primary, data = d, family = lognormal(link = "log"), REML = FALSE)
cat("glmmTMB lognormal family:\n")
print(round(summary(m_lntmb)$coefficients$cond[c("participant_mean_rate", "trial_rate_deviation"), ], 5))

# Are the two lme4 / glmmTMB Gamma log-likelihoods on the same scale?
m_base_tmb   <- glmmTMB(f_base, data = d, family = Gamma(link = "log"), REML = FALSE)
m_base_glmer <- glmer(f_base, data = d, family = Gamma(link = "log"),
                      control = glmerControl(optimizer = "bobyqa"))
cat("\nlogLik comparability check (same base model, same data):\n")
cat("  glmmTMB base logLik:", as.numeric(logLik(m_base_tmb)),
    " df:", attr(logLik(m_base_tmb), "df"), "\n")
cat("  glmer   base logLik:", as.numeric(logLik(m_base_glmer)),
    " df:", attr(logLik(m_base_glmer), "df"), "\n")
cat("  difference (glmer - glmmTMB):", as.numeric(logLik(m_base_glmer)) - as.numeric(logLik(m_base_tmb)), "\n")
cat("  glmer Gamma shape (1/sigma^2):", 1 / sigma(m_base_glmer)^2,
    "  glmmTMB shape:", sigma(m_base_tmb), "\n")

# ------------------------------------------------------------------------ 4. DHARMa
hdr("4. RESIDUAL DIAGNOSTICS (DHARMa, 1000 sims)")
set.seed(20260805)
sim <- simulateResiduals(m, n = 1000)
u <- testUniformity(sim); disp <- testDispersion(sim, plot = FALSE); out <- testOutliers(sim, plot = FALSE)
cat("uniformity KS: D =", u$statistic, " p =", u$p.value, "\n")
cat("dispersion   : stat =", disp$statistic, " p =", disp$p.value, "\n")
cat("outliers     : p =", out$p.value, "  (observed", out$statistic,
    "of 744, expected", round(744 * out$null.value, 2), ")\n")
cat("quantile dev : p =", testQuantiles(sim, plot = FALSE)$p.value, "\n")

# ---------------------------------------------------------------------------- 5. LRTs
hdr("5. LIKELIHOOD RATIO TESTS")
m_btw_only <- glmmTMB(update(f_base, . ~ . + participant_mean_rate), data = d,
                      family = Gamma(link = "log"), REML = FALSE)
m_wth_only <- glmmTMB(update(f_base, . ~ . + trial_rate_deviation), data = d,
                      family = Gamma(link = "log"), REML = FALSE)
lrt <- function(a, b, lab) {
  ch <- 2 * (as.numeric(logLik(b)) - as.numeric(logLik(a)))
  df <- attr(logLik(b), "df") - attr(logLik(a), "df")
  cat(sprintf("%-42s chi2(%d) = %.4f, p = %.6g  [logLik %.4f -> %.4f]\n",
              lab, df, ch, pchisq(ch, df, lower.tail = FALSE), logLik(a), logLik(b)))
  c(chisq = ch, df = df, p = pchisq(ch, df, lower.tail = FALSE))
}
l_both <- lrt(m_base_tmb, m, "base -> +both movement terms")
l_btw  <- lrt(m_wth_only, m, "within-only -> +between (between alone)")
l_wth  <- lrt(m_btw_only, m, "between-only -> +within (within alone)")
l_btw_m <- lrt(m_base_tmb, m_btw_only, "base -> +between only (marginal)")
l_wth_m <- lrt(m_base_tmb, m_wth_only, "base -> +within only (marginal)")

hdr("5b. THE PUBLISHED DEGENERATE glmer LRT (pooled movement_rate)")
m_pool_glmer <- glmer(f_pooled, data = d, family = Gamma(link = "log"),
                      control = glmerControl(optimizer = "bobyqa"))
cat("glmer base   logLik:", as.numeric(logLik(m_base_glmer)), "\n")
cat("glmer pooled logLik:", as.numeric(logLik(m_pool_glmer)), "\n")
cat("delta logLik on ADDING a parameter:", as.numeric(logLik(m_pool_glmer)) - as.numeric(logLik(m_base_glmer)), "\n")
print(anova(m_base_glmer, m_pool_glmer))
m_pool_tmb <- glmmTMB(f_pooled, data = d, family = Gamma(link = "log"), REML = FALSE)
print(round(summary(m_pool_tmb)$coefficients$cond["movement_rate", , drop = FALSE], 6))
lrt(m_base_tmb, m_pool_tmb, "glmmTMB base -> +pooled movement_rate")

# ---------------------------------------------------------------- 6. covariate control
hdr("6. DIFFICULTY-CONTROL MODELS")
m_ad  <- glmmTMB(update(f_primary, . ~ . + angDisp_mean_deg), data = d,
                 family = Gamma(link = "log"), REML = FALSE)
m_adr <- glmmTMB(update(f_primary, . ~ . + angDisp_mean_deg + responseTime_s), data = d,
                 family = Gamma(link = "log"), REML = FALSE)
cat("+angDisp:\n"); print(round(summary(m_ad)$coefficients$cond[
  c("participant_mean_rate", "trial_rate_deviation", "angDisp_mean_deg"), ], 6))
cat("\n+angDisp+responseTime:\n"); print(round(summary(m_adr)$coefficients$cond[
  c("participant_mean_rate", "trial_rate_deviation", "angDisp_mean_deg", "responseTime_s"), ], 6))

# ------------------------------------------------------------ 7. participant level
hdr("7. PARTICIPANT-LEVEL TRIANGULATION (n = 31)")
pl <- d %>% group_by(participantId) %>%
  summarise(mean_rate = mean(movement_rate),
            mean_err  = mean(participantError_m),
            med_err   = median(participantError_m),
            mean_path = mean(total_path_length),
            mean_dur  = mean(trajectory_duration),
            mean_rt   = mean(responseTime_s), .groups = "drop")

ct <- cor.test(pl$mean_rate, pl$mean_err)
cat("Pearson  r =", ct$estimate, " CI", ct$conf.int, " t(", ct$parameter, ") =", ct$statistic,
    " p =", ct$p.value, "\n")
st <- suppressWarnings(cor.test(pl$mean_rate, pl$mean_err, method = "spearman"))
cat("Spearman rho =", st$estimate, " S =", st$statistic, " p =", st$p.value, "\n")
cat("Pearson with participant MEDIAN error: r =", cor(pl$mean_rate, pl$med_err), "\n")

fit_log <- lm(log(mean_err) ~ mean_rate, data = pl)
cat("\nparticipant-level OLS, log mean error ~ mean rate:\n")
print(round(summary(fit_log)$coefficients, 5))
fit_raw <- lm(mean_err ~ mean_rate, data = pl)
cat("raw-scale slope:", coef(fit_raw)[2], "m per 1 m/s =",
    100 * coef(fit_raw)[2] / 10, "cm per 0.1 m/s\n")

# leverage / influence on the correlation
cat("\nleave-one-participant-out range of Pearson r:\n")
loo <- sapply(seq_len(nrow(pl)), function(i) cor(pl$mean_rate[-i], pl$mean_err[-i]))
cat("  min", min(loo), " max", max(loo), " most influential:", pl$participantId[which.max(abs(loo - ct$estimate))], "\n")
cat("max Cook's D (participant-level OLS on log error):", max(cooks.distance(fit_log)), "\n")
cat("participant with max mean rate:", pl$participantId[which.max(pl$mean_rate)],
    " rate =", max(pl$mean_rate), " (next:", sort(pl$mean_rate, decreasing = TRUE)[2], ")\n")

# bootstrap over participants
B <- 5000
set.seed(20260805)
bs <- replicate(B, {
  i <- sample(nrow(pl), replace = TRUE)
  x <- pl$mean_rate[i]; y <- pl$mean_err[i]
  if (sd(x) == 0 || sd(y) == 0) return(c(NA, NA, NA))
  c(cor(x, y), suppressWarnings(cor(x, y, method = "spearman")), coef(lm(log(y) ~ x))[2])
})
cat("\nbootstrap (5000 participant resamples, seed 20260805):\n")
cat("  Pearson  CI:", quantile(bs[1, ], c(.025, .975), na.rm = TRUE), "\n")
cat("  Spearman CI:", quantile(bs[2, ], c(.025, .975), na.rm = TRUE), "\n")
cat("  log-slope CI:", quantile(bs[3, ], c(.025, .975), na.rm = TRUE), "\n")
cat("  degenerate resamples dropped:", sum(is.na(bs[1, ])), "\n")
cat("  share of resamples with r < 0:", mean(bs[1, ] < 0, na.rm = TRUE),
    " (count with r >= 0:", sum(bs[1, ] >= 0, na.rm = TRUE), ")\n")
# The share is close enough to 1 that printing it at three decimals rounds to 1.000.
# Repeat across seeds so the reported value is not a rounding artefact.
cat("  share across seeds:\n")
for (s in c(1, 2, 42, 999, 2026)) {
  set.seed(s)
  rr <- replicate(B, { i <- sample(nrow(pl), replace = TRUE)
                       if (sd(pl$mean_rate[i]) == 0) NA else cor(pl$mean_rate[i], pl$mean_err[i]) })
  cat(sprintf("    seed %5d: share r<0 = %.4f (%d resamples with r >= 0)\n",
              s, mean(rr < 0, na.rm = TRUE), sum(rr >= 0, na.rm = TRUE)))
}

# permutation test of the between effect (exact-ish, respects 31 clusters)
set.seed(20260805)
perm <- replicate(2e5, cor(sample(pl$mean_rate), pl$mean_err))
cat("\npermutation p (two-sided, 200000 perms) for participant-level r:",
    (1 + sum(abs(perm) >= abs(ct$estimate))) / 200001, "\n")

# ------------------------------------------------------- 8. rate vs amount vs duration
hdr("8. RATE VS AMOUNT VS DURATION (participant level)")
for (v in c("mean_path", "mean_dur", "mean_rt")) {
  cc <- cor.test(pl[[v]], pl$mean_err)
  cat(sprintf("%-10s vs mean error: r = %+.4f  CI [%+.3f, %+.3f]  p = %.4f\n",
              v, cc$estimate, cc$conf.int[1], cc$conf.int[2], cc$p.value))
}
cat("cor(mean_path, mean_dur) =", cor(pl$mean_path, pl$mean_dur), "\n")
cat("cor(mean_rate, mean_dur) =", cor(pl$mean_rate, pl$mean_dur), "\n")
cat("cor(mean_rate, mean_path) =", cor(pl$mean_rate, pl$mean_path), "\n")

fit_pd <- lm(log(mean_err) ~ log(mean_path) + log(mean_dur), data = pl)
cat("\nlog mean error ~ log path + log duration:\n")
print(round(summary(fit_pd)$coefficients, 5))
cat("adj R2:", summary(fit_pd)$adj.r.squared, "\n")
L <- c(0, 1, 1)
est <- sum(L * coef(fit_pd)); se <- sqrt(drop(t(L) %*% vcov(fit_pd) %*% L))
cat("Wald test b_path + b_dur = 0: sum =", est, " SE =", se, " t(", df.residual(fit_pd), ") =",
    est / se, " p =", 2 * pt(abs(est / se), df.residual(fit_pd), lower.tail = FALSE), "\n")
fit_ratio <- lm(log(mean_err) ~ log(mean_rate), data = pl)
cat("log mean error ~ log mean rate: adj R2 =", summary(fit_ratio)$adj.r.squared,
    " slope =", coef(fit_ratio)[2], "\n")

# ----------------------------------------------------------- 9. within-null diagnostics
hdr("9. IS THE WITHIN-PARTICIPANT NULL INFORMATIVE?")
cat("within-term CI on the log scale:", b_wth - 1.96 * se_wth, "to", b_wth + 1.96 * se_wth, "\n")
cat("smallest between-style effect the within CI still excludes (log scale):\n")
cat("  upper bound of |effect| ruled out below:", b_wth - 1.96 * se_wth, "\n")
cat("  the between estimate is", b_btw, "-> is it inside the within CI?",
    b_btw > b_wth - 1.96 * se_wth && b_btw < b_wth + 1.96 * se_wth, "\n")
# formal test that between and within slopes are equal (Hausman-style contrast)
V <- vcov(m)$cond
L2 <- rep(0, nrow(V)); names(L2) <- rownames(V)
L2["participant_mean_rate"] <- 1; L2["trial_rate_deviation"] <- -1
dif <- sum(L2 * fixef(m)$cond); sed <- sqrt(drop(t(L2) %*% V %*% L2))
cat("between - within contrast:", dif, " SE =", sed, " z =", dif / sed,
    " p =", 2 * pnorm(abs(dif / sed), lower.tail = FALSE), "\n")
# within-subject range actually observed
cat("trial_rate_deviation range:", range(d$trial_rate_deviation),
    " IQR:", IQR(d$trial_rate_deviation), "\n")

# within-participant correlation, computed directly per participant
wc <- d %>% group_by(participantId) %>%
  summarise(r = cor(movement_rate, participantError_m), .groups = "drop")
cat("per-participant within-subject correlations: mean r =", mean(wc$r),
    " median =", median(wc$r), " n positive =", sum(wc$r > 0), "/", nrow(wc), "\n")
print(t.test(wc$r))

# ----------------------------------------------------------------- 10. cm scale check
hdr("10. cm-SCALED RESPONSE")
m_cm <- glmmTMB(update(f_primary, participantError_cm ~ .), data = d,
                family = Gamma(link = "log"), REML = FALSE)
cat("max abs slope difference:", max(abs(fixef(m_cm)$cond[-1] - fixef(m)$cond[-1])), "\n")
cat("intercept shift:", fixef(m_cm)$cond[1] - fixef(m)$cond[1], " log(100) =", log(100), "\n")

# ------------------------------------------------------------------------ save outputs
hdr("11. WRITING CSV OUTPUT")
coef_out <- rbind(
  data.frame(model = "glmmTMB primary", term = rownames(cf), cf, check.names = FALSE),
  data.frame(model = "glmer bobyqa (degenerate)", term = rownames(summary(m_glmer)$coefficients),
             summary(m_glmer)$coefficients, check.names = FALSE)
)
names(coef_out)[3:6] <- c("estimate", "std_error", "statistic", "p_value")
write.csv(coef_out, file.path(res_dir, "movement_effects_verification_coefficients.csv"), row.names = FALSE)

summ <- data.frame(
  quantity = c("beta_between", "se_between", "z_between", "p_between",
               "beta_within", "se_within", "z_within", "p_within",
               "RE_SD_glmmTMB", "RE_SD_glmer", "logLik_glmmTMB", "logLik_glmer",
               "pearson_r", "pearson_p", "spearman_rho", "spearman_p",
               "chi2_both", "chi2_between", "chi2_within",
               "path_r", "path_p", "dur_r"),
  value = c(b_btw, se_btw, cf["participant_mean_rate", 3], cf["participant_mean_rate", 4],
            b_wth, se_wth, cf["trial_rate_deviation", 3], cf["trial_rate_deviation", 4],
            attr(VarCorr(m)$cond$participantId, "stddev"),
            vc$sdcor[vc$grp == "participantId"],
            as.numeric(logLik(m)), as.numeric(logLik(m_glmer)),
            ct$estimate, ct$p.value, st$estimate, st$p.value,
            l_both["chisq"], l_btw["chisq"], l_wth["chisq"],
            cor.test(pl$mean_path, pl$mean_err)$estimate,
            cor.test(pl$mean_path, pl$mean_err)$p.value,
            cor(pl$mean_dur, pl$mean_err))
)
write.csv(summ, file.path(res_dir, "movement_effects_verification_summary.csv"), row.names = FALSE)
write.csv(pl, file.path(res_dir, "movement_effects_verification_participant_level.csv"), row.names = FALSE)

cat("\nsessionInfo:\n"); print(sessionInfo()$otherPkgs[["glmmTMB"]]$Version)
sink(type = "message"); sink(); close(con)
