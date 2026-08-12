#!/usr/bin/env Rscript
# Follow-up adversarial checks on the Reviewer 1 primary-model extensions.
#
#  A. sound type is manipulated WITHIN participants (6 trials per participant per
#     sound type), so the fixed-effects-only interaction test is anti-conservative.
#     Refit with by-participant random slopes for sound type and for the
#     disparity x sound type interaction and check whether the null survives and
#     how much wider the intervals become.
#  B. What slope differences are actually compatible with the data? The correct
#     bound is the largest upper confidence limit across the pairwise contrasts,
#     not the half-width of the widest interval.
#  C. Optimiser robustness of the three-term random-effects fit M2rs, which is the
#     fit most at risk of a degenerate variance-covariance estimate.
#  D. Profile likelihood on M0 to check how well the reported random-effect
#     parameters are identified.

library(lme4)
library(lmerTest)
library(emmeans)

set.seed(20260805)
emm_options(lmer.df = "satterthwaite")

root <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res  <- file.path(root, "results_revision")
d    <- readRDS(file.path(res, "analysis_df_revision.rds"))

con <- file(file.path(res, "verify_reviewer1_followups_log.txt"), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 5e5))

rule("A. IS SOUND TYPE WITHIN PARTICIPANT?")
print(table(d$participantId, d$soundType)[1:5, ])
cat("trials per participant per sound type (unique):",
    paste(sort(unique(as.vector(table(d$participantId, d$soundType)))), collapse = ","), "\n")

m1 <- lmer(signedError_m ~ stimulusDisparity_m * soundType +
             (1 + stimulusDisparity_m | participantId),
           data = d, REML = TRUE, control = ctrl)

fit_try <- function(form, label) {
  msgs <- character(0)
  fit <- withCallingHandlers(
    tryCatch(lmer(form, data = d, REML = TRUE, control = ctrl), error = function(e) e),
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  cat("\n---", label, "---\n")
  if (inherits(fit, "error")) { cat("FIT FAILED:", conditionMessage(fit), "\n"); return(NULL) }
  cat("singular:", isSingular(fit), "  n warnings/messages:", length(msgs), "\n")
  if (length(msgs)) cat("  ", paste(trimws(msgs), collapse = " || "), "\n")
  print(anova(fit, type = 3))
  fit
}

m1_rs_st <- fit_try(signedError_m ~ stimulusDisparity_m * soundType +
                      (1 + stimulusDisparity_m + soundType | participantId),
                    "M1 + by-participant random slope for soundType")
m1_rs_ix <- fit_try(signedError_m ~ stimulusDisparity_m * soundType +
                      (1 + stimulusDisparity_m * soundType | participantId),
                    "M1 + maximal (disparity x soundType) random slopes")

if (!is.null(m1_rs_st)) {
  cat("\nREML LRT M1 vs M1 + random soundType slope:\n")
  print(anova(m1, m1_rs_st, refit = FALSE))
  cat("\nPer-sound-type slopes under the random-soundType-slope model:\n")
  ets <- summary(emtrends(m1_rs_st, ~ soundType, var = "stimulusDisparity_m"),
                 infer = c(TRUE, TRUE))
  print(ets)
  cat("CI widths:", paste(round(ets$upper.CL - ets$lower.CL, 4), collapse = ", "), "\n")
  cat("\nPairwise contrasts under that model (Tukey):\n")
  print(summary(pairs(emtrends(m1_rs_st, ~ soundType, var = "stimulusDisparity_m")),
                infer = c(TRUE, TRUE)))
}

rule("B. WHAT SLOPE DIFFERENCES ARE ACTUALLY COMPATIBLE WITH THE DATA?")
et  <- emtrends(m1, ~ soundType, var = "stimulusDisparity_m")
pw_t <- summary(pairs(et), infer = c(TRUE, TRUE))
pw_u <- summary(pairs(et, adjust = "none"), infer = c(TRUE, TRUE))
pooled <- 0.3465964   # M0 capture slope

bnd <- data.frame(
  contrast   = pw_t$contrast,
  estimate   = pw_t$estimate,
  tukey_lo   = pw_t$lower.CL, tukey_hi = pw_t$upper.CL,
  tukey_width = pw_t$upper.CL - pw_t$lower.CL,
  unadj_lo   = pw_u$lower.CL, unadj_hi = pw_u$upper.CL,
  unadj_width = pw_u$upper.CL - pw_u$lower.CL)
bnd$tukey_hi_pct_of_pooled <- 100 * bnd$tukey_hi / pooled
bnd$unadj_hi_pct_of_pooled <- 100 * bnd$unadj_hi / pooled
print(bnd, row.names = FALSE, digits = 4)
cat("\nwidest Tukey pairwise CI width:", max(bnd$tukey_width),
    " (half-width", max(bnd$tukey_width) / 2, ")\n")
cat("widest UNADJUSTED pairwise CI width:", max(bnd$unadj_width),
    " (half-width", max(bnd$unadj_width) / 2, ")\n")
cat("largest slope difference still inside a Tukey CI (max upper CL):",
    max(bnd$tukey_hi), "=", 100 * max(bnd$tukey_hi) / pooled, "% of the pooled slope\n")
cat("largest slope difference still inside an UNADJUSTED CI (max upper CL):",
    max(bnd$unadj_hi), "=", 100 * max(bnd$unadj_hi) / pooled, "% of the pooled slope\n")
write.csv(bnd, file.path(res, "verify_reviewer1_slope_difference_bounds.csv"), row.names = FALSE)

rule("C. OPTIMISER ROBUSTNESS OF M2rs (three random-effect terms)")
f2rs <- signedError_m ~ stimulusDisparity_m * trialSequence_c + soundType +
  (1 + stimulusDisparity_m + trialSequence_c | participantId)
opts <- c("bobyqa", "Nelder_Mead", "nloptwrap", "nlminbwrap")
tab <- do.call(rbind, lapply(opts, function(o) {
  fit <- lmer(f2rs, data = d, REML = TRUE,
              control = lmerControl(optimizer = o, optCtrl = if (o == "bobyqa")
                list(maxfun = 5e5) else list()))
  s <- summary(fit)$coefficients; v <- as.data.frame(VarCorr(fit))
  data.frame(optimizer = o, logLik = as.numeric(logLik(fit)),
             b_disp = s["stimulusDisparity_m", "Estimate"],
             se_disp = s["stimulusDisparity_m", "Std. Error"],
             b_ix = s["stimulusDisparity_m:trialSequence_c", "Estimate"],
             p_ix = s["stimulusDisparity_m:trialSequence_c", "Pr(>|t|)"],
             sd_trial = v$sdcor[v$var1 == "trialSequence_c" & is.na(v$var2)],
             singular = isSingular(fit))
}))
print(tab, row.names = FALSE)
cat("\nspread in logLik:", diff(range(tab$logLik)),
    " spread in SD_trial:", diff(range(tab$sd_trial)), "\n")

rule("D. PROFILE LIKELIHOOD ON M0 RANDOM EFFECTS")
m0 <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
             (1 + stimulusDisparity_m | participantId),
           data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
pr <- confint(m0, method = "profile", oldNames = FALSE)
w  <- confint(m0, method = "Wald")
cat("Profile 95% CI:\n"); print(pr)
cat("\nWald 95% CI (fixed effects only):\n"); print(w)

cat("\nParametric bootstrap CI for the capture slope and the RE parameters (500 sims):\n")
bootci <- tryCatch(confint(m0, method = "boot", nsim = 500, seed = 20260805,
                           parm = c("sd_(Intercept)|participantId",
                                    "cor_stimulusDisparity_m.(Intercept)|participantId",
                                    "sd_stimulusDisparity_m|participantId",
                                    "stimulusDisparity_m"),
                           oldNames = FALSE),
                   error = function(e) paste("boot failed:", conditionMessage(e)))
print(bootci)

sink(type = "message"); sink(); close(con)
