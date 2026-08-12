# Follow-up verification: (a) optimizer stability of the outlier-sensitivity
# refits, (b) dimensional comparison of the orthogonal and on-axis error
# components, (c) session-level accounting of trial time.
#
# Outputs: results_revision/verify_convergence_log.txt
#          results_revision/verify_allfit_signed3sd.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

set.seed(20260805)
options(digits = 10)

root <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res  <- file.path(root, "results_revision")
con  <- file(file.path(res, "verify_convergence_log.txt"), open = "wt")
say  <- function(...) { t <- paste0(...); cat(t, "\n", sep = ""); cat(t, "\n", sep = "", file = con) }
show_df <- function(x, digits = 6) for (l in capture.output(
  print(as.data.frame(x), digits = digits, row.names = FALSE))) say(l)

d <- readRDS(file.path(res, "analysis_df_revision.rds")) %>%
  mutate(signedError_cm = 100 * signedError_m)
s <- readRDS(file.path(res, "trajectory_samples_angular.rds"))

cut_unsigned <- mean(d$participantError_cm) + 3 * sd(d$participantError_cm)
cut_lo <- mean(d$signedError_cm) - 3 * sd(d$signedError_cm)
cut_hi <- mean(d$signedError_cm) + 3 * sd(d$signedError_cm)

subsets <- list(
  full       = d,
  rt_ge_3    = filter(d, responseTime_s >= 3),
  err_3sd    = filter(d, participantError_cm <= cut_unsigned),
  signed_3sd = filter(d, signedError_cm >= cut_lo, signedError_cm <= cut_hi),
  both       = filter(d, responseTime_s >= 3, participantError_cm <= cut_unsigned)
)
form <- signedError_m ~ stimulusDisparity_m + soundType +
  (1 + stimulusDisparity_m | participantId)

# ---- (a) allFit across every available optimizer, per subset ---------------
say("=== (a) OPTIMIZER STABILITY OF THE PRIMARY REFITS ===")
rows <- list()
for (nm in names(subsets)) {
  f0 <- lmer(form, data = subsets[[nm]], REML = TRUE,
             control = lmerControl(optimizer = "bobyqa"))
  af <- suppressWarnings(suppressMessages(allFit(f0, verbose = FALSE)))
  ll <- summary(af)$llik
  ok <- vapply(af, function(x) inherits(x, "merMod"), logical(1))
  for (o in names(af)[ok]) {
    fx <- af[[o]]
    cf <- coef(summary(as(fx, "lmerModLmerTest")))["stimulusDisparity_m", ]
    vc <- as.data.frame(VarCorr(fx))
    rows[[paste(nm, o)]] <- data.frame(
      subset = nm, optimizer = o, logLik_REML = as.numeric(ll[o]),
      est = cf[["Estimate"]], se = cf[["Std. Error"]],
      df = cf[["df"]], p = cf[["Pr(>|t|)"]],
      sd_int = vc$sdcor[vc$var1 == "(Intercept)" & is.na(vc$var2)],
      sd_slope = vc$sdcor[vc$var1 == "stimulusDisparity_m" & is.na(vc$var2)],
      cor_int_slope = vc$sdcor[!is.na(vc$var2)],
      singular = isSingular(fx),
      stringsAsFactors = FALSE)
  }
}
allfit_tbl <- bind_rows(rows) %>%
  group_by(subset) %>% mutate(dLogLik = logLik_REML - max(logLik_REML)) %>% ungroup()
show_df(allfit_tbl, 7)
write.csv(allfit_tbl, file.path(res, "verify_allfit_signed3sd.csv"), row.names = FALSE)

say("\n-- spread of the disparity slope and its SE within each subset --")
show_df(allfit_tbl %>% group_by(subset) %>%
          summarise(n_opt = n(), n_singular = sum(singular),
                    est_min = min(est), est_max = max(est),
                    se_min = min(se), se_max = max(se),
                    df_min = min(df), df_max = max(df),
                    logLik_range = max(logLik_REML) - min(logLik_REML),
                    .groups = "drop"), 7)

# Which solution has the higher REML criterion in the signed_3sd subset?
say("\n-- signed_3sd: bobyqa vs Nelder_Mead, deviance and correlation parameter --")
d3 <- subsets$signed_3sd
fb <- lmer(form, d3, REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
fn <- suppressWarnings(suppressMessages(
  lmer(form, d3, REML = TRUE, control = lmerControl(optimizer = "Nelder_Mead"))))
say("bobyqa      REML logLik = ", format(as.numeric(logLik(fb))),
    "   theta = ", paste(format(getME(fb, "theta")), collapse = ", "))
say("Nelder_Mead REML logLik = ", format(as.numeric(logLik(fn))),
    "   theta = ", paste(format(getME(fn, "theta")), collapse = ", "))
say("bobyqa corr(int,slope)      = ", format(as.data.frame(VarCorr(fb))$sdcor[3]))
say("Nelder_Mead corr(int,slope) = ", format(as.data.frame(VarCorr(fn))$sdcor[3]))
say("bobyqa RE SDs      : ", paste(format(as.data.frame(VarCorr(fb))$sdcor[1:2]), collapse = ", "))
say("Nelder_Mead RE SDs : ", paste(format(as.data.frame(VarCorr(fn))$sdcor[1:2]), collapse = ", "))

# also check the full-data fit for the same pathology
say("\n-- correlation parameter in every subset (bobyqa) --")
for (nm in names(subsets)) {
  f <- lmer(form, subsets[[nm]], REML = TRUE, control = lmerControl(optimizer = "bobyqa"))
  vc <- as.data.frame(VarCorr(f))
  say(sprintf("%-11s sd_int=%.5f sd_slope=%.5f corr=%+.4f singular=%s",
              nm, vc$sdcor[1], vc$sdcor[2], vc$sdcor[3], isSingular(f)))
}

# ---- (b) dimensional comparison of error components ------------------------
say("\n=== (b) ORTHOGONAL vs ON-AXIS: DIMENSIONAL CHECK ===")
oe <- 100 * d$orthogonalError_m
se <- d$signedError_cm
say("orthogonal is the norm of a 2-D residual; signed is a 1-D scalar.")
say("mean orth / mean |signed|            = ", format(mean(oe) / mean(abs(se))))
say("expected ratio under isotropic noise = ", format(sqrt(pi / 2) / sqrt(2 / pi)),
    "  (= pi/2)")
say("observed / isotropic expectation     = ", format((mean(oe) / mean(abs(se))) / (pi / 2)))
say("per-axis RMS, off-axis  = ", format(sqrt(mean(oe^2) / 2)), " cm")
say("per-axis RMS, on-axis   = ", format(sqrt(mean(se^2))), " cm")
say("on-axis SD about mean   = ", format(sd(se)), " cm  (mean = ", format(mean(se)), " cm)")
say("on-axis residual SD after removing disparity effect:")
lm_ax <- lm(se ~ stimulusDisparity_m, data = d)
say("  ", format(summary(lm_ax)$sigma), " cm")

# ---- (c) session-level accounting -----------------------------------------
say("\n=== (c) SESSION TIME ACCOUNTING ===")
sess <- s %>% group_by(participantId) %>%
  summarise(span = max(timestamp) - min(timestamp),
            t_first = min(timestamp), .groups = "drop") %>%
  left_join(d %>% group_by(participantId) %>%
              summarise(sum_rt = sum(responseTime_s),
                        sum_traj = sum(trajectory_duration), .groups = "drop"),
            by = "participantId") %>%
  mutate(pct_rt = 100 * sum_rt / span, pct_traj = 100 * sum_traj / span) %>%
  arrange(pct_rt)
show_df(sess, 7)
say("sessions with < 95% of span covered by trial time: ", sum(sess$pct_rt < 95))
say("min pct_rt = ", format(min(sess$pct_rt)), " (", sess$participantId[1], ")")
say("unaccounted time in that session = ", format(sess$span[1] - sess$sum_rt[1]), " s")

# is the shortfall a within-session gap or a leading/trailing offset?
p_bad <- sess$participantId[1]
runs_bad <- s %>% filter(participantId == p_bad) %>%
  group_by(run_id) %>%
  summarise(t0 = min(timestamp), t1 = max(timestamp), n = n(), .groups = "drop") %>%
  arrange(t0) %>% mutate(gap_before = t0 - lag(t1), dur = t1 - t0)
say("\n-- run structure of ", p_bad, " --")
show_df(runs_bad, 7)
say("sum of run durations = ", format(sum(runs_bad$dur)),
    " s vs span ", format(sess$span[1]), " s")

close(con)
