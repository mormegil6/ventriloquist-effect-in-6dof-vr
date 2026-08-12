# verify_tolerance_covariance.R
#
# Focused check on the random-effect covariance of the primary signed-error LMM,
# because the 10th/90th-percentile listener profiles in the tolerance lookup table
# depend entirely on the random-intercept SD and the intercept-slope correlation.
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/verify_tolerance_covariance_log.txt

suppressPackageStartupMessages({ library(lme4); library(dplyr) })

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir  <- file.path(proj_dir, "results_revision")
con <- file(file.path(res_dir, "verify_tolerance_covariance_log.txt"), open = "wt")
sink(con, split = TRUE); sink(con, type = "message")
on.exit({ sink(type = "message"); sink(); close(con) }, add = TRUE)
say <- function(...) cat(..., "\n", sep = "")

d <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
ctrl <- lmerControl(optimizer = "bobyqa")

# Full model (as published) and two nested restrictions of the covariance matrix.
m_full <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                 (1 + stimulusDisparity_m | participantId),
               data = d, REML = FALSE, control = ctrl)
m_diag <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                 (1 | participantId) + (0 + stimulusDisparity_m | participantId),
               data = d, REML = FALSE, control = ctrl)
m_slop <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                 (0 + stimulusDisparity_m | participantId),
               data = d, REML = FALSE, control = ctrl)

say("ML log-likelihoods: full ", round(logLik(m_full), 4),
    " | diagonal (rho = 0) ", round(logLik(m_diag), 4),
    " | slope only (tau0 = 0, rho = 0) ", round(logLik(m_slop), 4))
say("\nLRT full vs diagonal (rho = 0, 1 df):")
print(anova(m_diag, m_full))
say("\nLRT diagonal vs slope-only (tau0 = 0, 1 df, boundary):")
print(anova(m_slop, m_diag))
say("\nLRT full vs slope-only (2 df):")
print(anova(m_slop, m_full))

# The 15 cm prediction under each covariance assumption, for the 10th and 90th
# percentile listener. This is the number the manuscript paragraph relies on.
z <- qnorm(0.90)
prof15 <- function(fit, label) {
  fx <- fixef(fit)
  b0 <- unname(fx[1]) + mean(c(0, fx[grep("^soundType", names(fx))]))
  b1 <- unname(fx[2])
  vc <- VarCorr(fit)
  if (length(vc) == 1 && ncol(vc[[1]]) == 2) {
    sd <- attr(vc$participantId, "stddev"); rho <- attr(vc$participantId, "correlation")[1, 2]
    tau0 <- unname(sd[1]); tau1 <- unname(sd[2])
  } else if (length(vc) == 2) {
    tau0 <- unname(attr(vc[[which(sapply(vc, function(v) rownames(v)[1] == "(Intercept)"))]], "stddev")[1])
    tau1 <- unname(attr(vc[[which(sapply(vc, function(v) rownames(v)[1] != "(Intercept)"))]], "stddev")[1])
    rho <- 0
  } else {
    tau0 <- 0; tau1 <- unname(attr(vc[[1]], "stddev")[1]); rho <- 0
  }
  b0d <- if (tau1 > 0) rho * (tau0 / tau1) * (z * tau1) else 0
  for (x in c(0.15, 0.70)) {
    lo <- (b0 - b0d + (b1 - z * tau1) * x) * 100
    hi <- (b0 + b0d + (b1 + z * tau1) * x) * 100
    say(sprintf("%-28s %3.0f cm: p10 = %6.2f cm, p90 = %6.2f cm, ratio = %.2f",
                label, x * 100, lo, hi, hi / lo))
  }
  say(sprintf("   tau0 = %.5f, tau1 = %.5f, rho = %+.4f", tau0, tau1, rho))
}
say("")
prof15(m_full, "full (published)")
prof15(m_diag, "diagonal (rho = 0)")
prof15(m_slop, "slope only")

# Observed per-participant displacement at the small end of the range, as a
# model-free check on the claim that listeners are indistinguishable at 15 cm.
small <- d %>% filter(stimulusDisparity_m <= 0.22) %>%
  group_by(participantId) %>%
  summarise(n = n(), mean_signed_cm = mean(signedError_m) * 100,
            mean_disp_cm = mean(stimulusDisparity_m) * 100, .groups = "drop")
say("\nobserved per-participant mean signed error on trials with disparity <= 22 cm")
say("participants ", nrow(small), ", trials ", sum(small$n))
say("quantiles of per-participant mean signed error (cm): ",
    paste(round(quantile(small$mean_signed_cm, c(0, .1, .5, .9, 1)), 2), collapse = " "))
say("SD across participants (cm): ", round(sd(small$mean_signed_cm), 2))
say("one-way ANOVA of signed error on participant, disparity <= 22 cm:")
sub <- d %>% filter(stimulusDisparity_m <= 0.22)
print(anova(lm(signedError_m ~ participantId, data = sub)))

print(sessionInfo())
