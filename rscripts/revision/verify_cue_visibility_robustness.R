# Robustness of the cue-visibility (head-orientation proxy) result:
# random-effects structure, likelihood-ratio test for a by-participant slope,
# and sensitivity to the 55 deg field-of-view threshold.
#
# Outputs: results_revision/verify_cue_robustness_log.txt
#          results_revision/verify_cue_robustness.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

set.seed(20260805)
options(digits = 10)

root <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res  <- file.path(root, "results_revision")
con  <- file(file.path(res, "verify_cue_robustness_log.txt"), open = "wt")
say  <- function(...) { t <- paste0(...); cat(t, "\n", sep = ""); cat(t, "\n", sep = "", file = con) }
show_df <- function(x, digits = 6) for (l in capture.output(
  print(as.data.frame(x), digits = digits, row.names = FALSE))) say(l)

d <- readRDS(file.path(res, "analysis_df_revision.rds")) %>%
  mutate(signedError_cm = 100 * signedError_m, bias_pct = 100 * ventriloquistBias)
s <- readRDS(file.path(res, "trajectory_samples_angular.rds"))
ctrl <- lmerControl(optimizer = "bobyqa")

centre <- function(dat, v) {
  dat$.x <- dat[[v]]
  dat %>% group_by(participantId) %>% mutate(.pm = mean(.x)) %>% ungroup() %>%
    mutate(frac_within = .x - .pm, frac_between = .pm - mean(.pm))
}

# ---- (a) random-effects structure -----------------------------------------
say("=== (a) SENSITIVITY OF THE WITHIN-PARTICIPANT CUE EFFECT TO RE STRUCTURE ===")
dd <- centre(d, "frac_time_cue_within_55deg")
specs <- list(
  "(1|pid)  [as reported]"          = signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType + (1 | participantId),
  "(1+disparity|pid)  [primary RE]" = signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType + (1 + stimulusDisparity_m | participantId),
  "(1+frac_within|pid)  [maximal]"  = signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType + (1 + frac_within | participantId),
  "(1+frac_within+disparity|pid)"   = signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType + (1 + frac_within + stimulusDisparity_m | participantId),
  "(1|pid) + duration"              = signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType + scale(responseTime_s) + (1 | participantId)
)
rows <- list()
for (nm in names(specs)) {
  w <- character(0)
  m <- withCallingHandlers(lmer(specs[[nm]], dd, REML = TRUE, control = ctrl),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") },
    message = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleMessage") })
  cf <- coef(summary(m))
  rows[[nm]] <- data.frame(
    spec = nm,
    within_est_per10pp = cf["frac_within", "Estimate"] / 10,
    within_se_per10pp  = cf["frac_within", "Std. Error"] / 10,
    within_df = cf["frac_within", "df"], within_p = cf["frac_within", "Pr(>|t|)"],
    between_est_per10pp = cf["frac_between", "Estimate"] / 10,
    between_p = cf["frac_between", "Pr(>|t|)"],
    singular = isSingular(m),
    warnings = if (length(w)) paste(w, collapse = " | ") else "none",
    stringsAsFactors = FALSE)
}
tbl <- bind_rows(rows)
show_df(tbl, 5)
write.csv(tbl, file.path(res, "verify_cue_robustness.csv"), row.names = FALSE)

say("\n-- LRT: is a by-participant slope for frac_within supported? (ML) --")
m0 <- lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType +
             (1 | participantId), dd, REML = FALSE, control = ctrl)
m1 <- lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType +
             (1 + frac_within | participantId), dd, REML = FALSE, control = ctrl)
show_df(cbind(model = c("(1|pid)", "(1+frac_within|pid)"), as.data.frame(anova(m0, m1))), 5)

say("\n-- sign test on the 31 per-participant correlations --")
pp <- dd %>% group_by(participantId) %>%
  summarise(r = cor(frac_time_cue_within_55deg, signedError_cm), .groups = "drop")
bt <- binom.test(sum(pp$r > 0), nrow(pp))
say("positive in ", sum(pp$r > 0), "/", nrow(pp), "   binom.test p = ", format(bt$p.value),
    "   median r = ", format(median(pp$r)))
say("one-sample t on Fisher-z of per-participant r: p = ",
    format(t.test(atanh(pp$r))$p.value))

# ---- (b) threshold sensitivity --------------------------------------------
say("\n=== (b) SENSITIVITY TO THE 55 DEG THRESHOLD ===")
thr_tbl <- list()
for (th in c(20, 30, 40, 55, 70, 90)) {
  ft <- s %>% group_by(participantId, trialSequenceNum) %>%
    summarise(f = sum(w_time * (cueEcc_deg <= th)) / sum(w_time), .groups = "drop")
  dj <- d %>% left_join(ft, by = c("participantId", "trialSequenceNum"))
  ctp <- cor.test(dj$f, dj$signedError_cm)
  dj2 <- centre(dj, "f")
  mi <- lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType +
               (1 | participantId), dj2, REML = TRUE, control = ctrl)
  ms <- suppressWarnings(suppressMessages(
    lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m + soundType +
           (1 + frac_within | participantId), dj2, REML = TRUE, control = ctrl)))
  thr_tbl[[as.character(th)]] <- data.frame(
    threshold_deg = th,
    pooled_frac_samples = mean(s$cueEcc_deg <= th),
    mean_trial_frac = mean(dj$f),
    pooled_r = unname(ctp$estimate), pooled_p = ctp$p.value,
    lmm_ri_within_per10pp = coef(summary(mi))["frac_within", "Estimate"] / 10,
    lmm_ri_p = coef(summary(mi))["frac_within", "Pr(>|t|)"],
    lmm_rs_within_per10pp = coef(summary(ms))["frac_within", "Estimate"] / 10,
    lmm_rs_p = coef(summary(ms))["frac_within", "Pr(>|t|)"])
}
show_df(bind_rows(thr_tbl), 5)

# ---- (c) cue vs source eccentricity, like-for-like -------------------------
say("\n=== (c) CUE vs SOURCE ECCENTRICITY, MATCHED METRICS ===")
say("pooled over samples   : cue ", format(100 * mean(s$cueEcc_deg <= 55)), "%   source ",
    format(100 * mean(s$srcEcc_deg <= 55)), "%")
say("time-weighted pooled  : cue ",
    format(100 * sum(s$w_time * (s$cueEcc_deg <= 55)) / sum(s$w_time)), "%   source ",
    format(100 * sum(s$w_time * (s$srcEcc_deg <= 55)) / sum(s$w_time)), "%")
say("per-trial mean of the trial fractions : cue ",
    format(100 * mean(d$frac_time_cue_within_55deg)), "%   source ",
    format(100 * mean(d$frac_time_src_within_55deg)), "%")
pd <- d$frac_time_cue_within_55deg - d$frac_time_src_within_55deg
tt <- t.test(pd)
say("paired difference (cue - source), per trial: M = ", format(mean(pd)),
    "  t(", tt$parameter, ") = ", format(unname(tt$statistic)), "  p = ", format(tt$p.value))
ppd <- d %>% group_by(participantId) %>%
  summarise(dd = mean(frac_time_cue_within_55deg - frac_time_src_within_55deg), .groups = "drop")
tt2 <- t.test(ppd$dd)
say("participant-level paired difference: M = ", format(mean(ppd$dd)),
    "  t(", tt2$parameter, ") = ", format(unname(tt2$statistic)), "  p = ", format(tt2$p.value),
    "  positive in ", sum(ppd$dd > 0), "/", nrow(ppd))
say("NOTE: cue and source are on average ", format(mean(d$stimulusDisparity_m)),
    " m apart, so their eccentricities are not independent measures.")

close(con)
