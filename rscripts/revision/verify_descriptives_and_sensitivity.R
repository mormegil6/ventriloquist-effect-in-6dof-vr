# Independent verification of the revision descriptives and outlier-sensitivity
# analyses. Recomputes every headline number from the prepared trial-level data
# frame and the per-sample trajectory file, without reusing the original script.
#
# Inputs : results_revision/analysis_df_revision.rds
#          results_revision/trajectory_samples_angular.rds
# Outputs: results_revision/verify_descriptives_log.txt
#          results_revision/verify_*.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lme4)
  library(lmerTest)
})

set.seed(20260805)
options(digits = 10)

root <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res  <- file.path(root, "results_revision")
logf <- file.path(res, "verify_descriptives_log.txt")

con <- file(logf, open = "wt")
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = con)
}
show_df <- function(x, digits = 6) {
  out <- capture.output(print(as.data.frame(x), digits = digits, row.names = FALSE))
  for (l in out) say(l)
}

d <- readRDS(file.path(res, "analysis_df_revision.rds"))
s <- readRDS(file.path(res, "trajectory_samples_angular.rds"))

say("=== VERIFICATION RUN ===")
say("R ", getRversion(), "  lme4 ", as.character(packageVersion("lme4")),
    "  lmerTest ", as.character(packageVersion("lmerTest")))
say("trials: ", nrow(d), "   participants: ", n_distinct(d$participantId))
say("samples: ", nrow(s), "   runs: ", n_distinct(paste(s$participantId, s$run_id)))

# ---------------------------------------------------------------- (1) distance
say("\n=== (1) LISTENING DISTANCE ===")
wq <- function(x, w, p) {
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  approx(cw, x, xout = p, method = "linear", rule = 2, ties = "ordered")$y
}
dist_all <- s$dist_source_m
say("per-sample mean   = ", format(mean(dist_all)))
say("per-sample SD     = ", format(sd(dist_all)))
say("per-sample quart  = ", paste(format(quantile(dist_all, c(0, .25, .5, .75, 1))), collapse = " "))
say("time-wtd mean     = ", format(sum(s$w_time * dist_all) / sum(s$w_time)))
say("time-wtd median   = ", format(wq(dist_all, s$w_time, 0.5)))
n05  <- sum(dist_all < 0.5); n0574 <- sum(dist_all < 0.574)
say("n < 0.5 m   = ", n05,  "  (", format(100 * n05 / nrow(s)), "%)")
say("n < 0.574 m = ", n0574, "  (", format(100 * n0574 / nrow(s)), "%)")
say("  [<= instead of < : ", sum(dist_all <= 0.5), " / ", sum(dist_all <= 0.574), "]")
say("time-wtd frac <0.5   = ", format(100 * sum(s$w_time * (dist_all < 0.5)) / sum(s$w_time)))
say("time-wtd frac <0.574 = ", format(100 * sum(s$w_time * (dist_all < 0.574)) / sum(s$w_time)))
say("n samples > 3 m = ", sum(dist_all > 3), "   max = ", format(max(dist_all)))

trial_dist <- s %>%
  group_by(participantId, trialSequenceNum) %>%
  summarise(tw_mean = sum(w_time * dist_source_m) / sum(w_time),
            dmin    = min(dist_source_m),
            donset  = dist_source_m[which.min(sample_idx + 1e6 * (run_id != min(run_id)))],
            any05   = any(dist_source_m < 0.5),
            any0574 = any(dist_source_m < 0.574),
            any010  = any(dist_source_m < 0.10),
            .groups = "drop")
say("per-trial time-wtd mean dist: M = ", format(mean(trial_dist$tw_mean)),
    "  SD = ", format(sd(trial_dist$tw_mean)), "  Mdn = ", format(median(trial_dist$tw_mean)))
say("per-trial min dist: M = ", format(mean(trial_dist$dmin)),
    "  Mdn = ", format(median(trial_dist$dmin)))
say("dist at trial onset: M = ", format(mean(trial_dist$donset)),
    "  Mdn = ", format(median(trial_dist$donset)))
say("trials with any sample <0.5 m   : ", sum(trial_dist$any05), "/", nrow(trial_dist),
    " (", format(100 * mean(trial_dist$any05)), "%)")
say("trials with any sample <0.574 m : ", sum(trial_dist$any0574), " (",
    format(100 * mean(trial_dist$any0574)), "%)")
say("trials with any sample <0.10 m  : ", sum(trial_dist$any010), " (",
    format(100 * mean(trial_dist$any010)), "%)")

pp_min <- s %>% group_by(participantId) %>%
  summarise(dmin_cm = 100 * min(dist_source_m), .groups = "drop")
say("per-participant min dist (cm): M = ", format(mean(pp_min$dmin_cm)),
    "  SD = ", format(sd(pp_min$dmin_cm)), "  Mdn = ", format(median(pp_min$dmin_cm)),
    "  range ", format(min(pp_min$dmin_cm)), "-", format(max(pp_min$dmin_cm)))
say("participants reaching <10 cm: ", sum(pp_min$dmin_cm < 10), "/", nrow(pp_min))

# ---------------------------------------------------------- (2) cue visibility
say("\n=== (2) CUE VISIBILITY ===")
say("cueEcc per-sample: M = ", format(mean(s$cueEcc_deg)), "  SD = ", format(sd(s$cueEcc_deg)),
    "  quart = ", paste(format(quantile(s$cueEcc_deg, c(0, .25, .5, .75, 1))), collapse = " "))
for (th in c(30, 55, 90)) {
  say("pooled samples with cue within ", th, " deg: ",
      format(100 * mean(s$cueEcc_deg <= th)), "%   (time-wtd ",
      format(100 * sum(s$w_time * (s$cueEcc_deg <= th)) / sum(s$w_time)), "%)")
}
say("pooled samples with SOURCE within 55 deg: ", format(100 * mean(s$srcEcc_deg <= 55)), "%")
say("time-wtd pooled SOURCE within 55 deg   : ",
    format(100 * sum(s$w_time * (s$srcEcc_deg <= 55)) / sum(s$w_time)), "%")

fc <- d$frac_time_cue_within_55deg; fs <- d$frac_time_src_within_55deg
say("per-trial frac cue<55: M = ", format(mean(fc)), " SD = ", format(sd(fc)),
    " Mdn = ", format(median(fc)), " Q1 = ", format(quantile(fc, .25)),
    " Q3 = ", format(quantile(fc, .75)), " range ", format(min(fc)), "-", format(max(fc)))
say("per-trial frac src<55: M = ", format(mean(fs)), " Mdn = ", format(median(fs)))
say("trials with cue outside 55 deg throughout: ", sum(fc == 0))
say("trials with cue inside 55 deg >90% of time: ", sum(fc > 0.90))

# pooled correlations (trial-level, ignoring participant clustering)
d <- d %>% mutate(signedError_cm = 100 * signedError_m,
                  bias_pct = 100 * ventriloquistBias)
cor_tbl <- bind_rows(
  broom_ct <- {
    ct <- cor.test(fc, d$signedError_cm)
    data.frame(pair = "frac_cue55 vs signedError_cm", method = "pearson",
               est = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2],
               stat = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value)
  },
  { ct <- suppressWarnings(cor.test(fc, d$signedError_cm, method = "spearman"))
    data.frame(pair = "frac_cue55 vs signedError_cm", method = "spearman",
               est = unname(ct$estimate), lo = NA, hi = NA,
               stat = unname(ct$statistic), df = NA, p = ct$p.value) },
  { ct <- cor.test(fc, d$bias_pct)
    data.frame(pair = "frac_cue55 vs bias_pct", method = "pearson",
               est = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2],
               stat = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value) },
  { ct <- suppressWarnings(cor.test(fc, d$bias_pct, method = "spearman"))
    data.frame(pair = "frac_cue55 vs bias_pct", method = "spearman",
               est = unname(ct$estimate), lo = NA, hi = NA,
               stat = unname(ct$statistic), df = NA, p = ct$p.value) },
  { ct <- cor.test(d$cueEcc_mean_deg, d$signedError_cm)
    data.frame(pair = "cueEcc_mean vs signedError_cm", method = "pearson",
               est = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2],
               stat = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value) },
  { ct <- cor.test(fc, d$stimulusDisparity_m)
    data.frame(pair = "frac_cue55 vs disparity", method = "pearson",
               est = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2],
               stat = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value) },
  { ct <- cor.test(fc, d$responseTime_s)
    data.frame(pair = "frac_cue55 vs responseTime_s", method = "pearson",
               est = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2],
               stat = unname(ct$statistic), df = unname(ct$parameter), p = ct$p.value) }
)
say("\n-- pooled trial-level correlations --")
show_df(cor_tbl)
write.csv(cor_tbl, file.path(res, "verify_cue_correlations.csv"), row.names = FALSE)

# per-participant correlations
pp_cor <- d %>% group_by(participantId) %>%
  summarise(r = cor(frac_time_cue_within_55deg, signedError_cm), .groups = "drop")
say("per-participant r: positive in ", sum(pp_cor$r > 0), "/", nrow(pp_cor),
    ", median r = ", format(median(pp_cor$r)))

# within/between decomposition (Mundlak)
d <- d %>% group_by(participantId) %>%
  mutate(fc_pm = mean(frac_time_cue_within_55deg)) %>% ungroup() %>%
  mutate(frac_within  = frac_time_cue_within_55deg - fc_pm,
         frac_between = fc_pm - mean(fc_pm))

ctrl <- lmerControl(optimizer = "bobyqa")
m_cue <- lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m +
                soundType + (1 | participantId), data = d, REML = TRUE, control = ctrl)
say("\n-- cue-visibility LMM on signed error (cm) --")
show_df(cbind(term = rownames(coef(summary(m_cue))), as.data.frame(coef(summary(m_cue)))), 6)
say("isSingular = ", isSingular(m_cue))

m_bias <- lmer(bias_pct ~ frac_within + frac_between + soundType +
                 (1 | participantId), data = d, REML = TRUE, control = ctrl)
say("\n-- cue-visibility LMM on bias proportion (%) --")
show_df(cbind(term = rownames(coef(summary(m_bias))), as.data.frame(coef(summary(m_bias)))), 6)

# sensitivity: add trial duration; add random slope
m_cue_dur <- lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m +
                    soundType + scale(responseTime_s) + (1 | participantId),
                  data = d, REML = TRUE, control = ctrl)
say("\n-- same model + scaled trial duration --")
show_df(cbind(term = rownames(coef(summary(m_cue_dur))), as.data.frame(coef(summary(m_cue_dur)))), 6)

m_cue_rs <- lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m +
                   soundType + (1 + frac_within | participantId),
                 data = d, REML = TRUE, control = ctrl)
say("\n-- same model + random slope for frac_within --")
show_df(cbind(term = rownames(coef(summary(m_cue_rs))), as.data.frame(coef(summary(m_cue_rs)))), 6)
say("isSingular(random slope model) = ", isSingular(m_cue_rs))

# same, but with the primary model's random-slope-for-disparity structure
m_cue_ps <- lmer(signedError_cm ~ frac_within + frac_between + stimulusDisparity_m +
                   soundType + (1 + stimulusDisparity_m | participantId),
                 data = d, REML = TRUE, control = ctrl)
say("\n-- same model with the PRIMARY random-effects structure --")
show_df(cbind(term = rownames(coef(summary(m_cue_ps))), as.data.frame(coef(summary(m_cue_ps)))), 6)

# ------------------------------------------------------- (3) inter-trial gaps
say("\n=== (3) INTER-TRIAL TIMING ===")
runs <- s %>% group_by(participantId, run_id) %>%
  summarise(t_first = min(timestamp), t_last = max(timestamp),
            trialSeq = first(trialSequenceNum), n = n(), .groups = "drop") %>%
  arrange(participantId, t_first)
say("n runs = ", nrow(runs), "  n sessions = ", n_distinct(runs$participantId))
gaps <- runs %>% group_by(participantId) %>%
  mutate(gap = t_first - lag(t_last)) %>% ungroup() %>% filter(!is.na(gap))
say("n gaps = ", nrow(gaps))
say("gap: M = ", format(mean(gaps$gap)), "  Mdn = ", format(median(gaps$gap)),
    "  IQR ", format(quantile(gaps$gap, .25)), "-", format(quantile(gaps$gap, .75)),
    "  range ", format(min(gaps$gap)), "-", format(max(gaps$gap)))
say("n gaps > 0.2 s = ", sum(gaps$gap > 0.2), "   n negative = ", sum(gaps$gap < 0))

say("responseTime_s: Mdn = ", format(median(d$responseTime_s)),
    "  IQR ", format(quantile(d$responseTime_s, .25)), "-", format(quantile(d$responseTime_s, .75)),
    "  range ", format(min(d$responseTime_s)), "-", format(max(d$responseTime_s)),
    "  M = ", format(mean(d$responseTime_s)), "  SD = ", format(sd(d$responseTime_s)))
say("  p05/p95 = ", paste(format(quantile(d$responseTime_s, c(.05, .95))), collapse = " / "))
say("  n < 3 s = ", sum(d$responseTime_s < 3), "   n > 120 s = ", sum(d$responseTime_s > 120))

sess <- s %>% group_by(participantId) %>%
  summarise(span = max(timestamp) - min(timestamp), .groups = "drop") %>%
  left_join(d %>% group_by(participantId) %>%
              summarise(sum_rt = sum(responseTime_s),
                        sum_traj = sum(trajectory_duration), .groups = "drop"),
            by = "participantId") %>%
  mutate(pct_rt = 100 * sum_rt / span, pct_traj = 100 * sum_traj / span)
say("session span: Mdn = ", format(median(sess$span)), " s  range ",
    format(min(sess$span)), "-", format(max(sess$span)), " s")
say("  = ", format(median(sess$span) / 60), " min; range ",
    format(min(sess$span) / 60), "-", format(max(sess$span) / 60), " min")
say("sum(responseTime)/span: Mdn = ", format(median(sess$pct_rt)), "%  range ",
    format(min(sess$pct_rt)), "-", format(max(sess$pct_rt)), "%")
say("sum(trajectory_duration)/span: Mdn = ", format(median(sess$pct_traj)), "%")
write.csv(sess, file.path(res, "verify_session_spans.csv"), row.names = FALSE)

# ------------------------------------------------------- (4) orthogonal error
say("\n=== (4) ORTHOGONAL ERROR ===")
oe <- 100 * d$orthogonalError_m
say("orth (cm): M = ", format(mean(oe)), " SD = ", format(sd(oe)), " Mdn = ", format(median(oe)),
    " Q1 = ", format(quantile(oe, .25)), " Q3 = ", format(quantile(oe, .75)),
    " range ", format(min(oe)), "-", format(max(oe)))
se_cm <- d$signedError_cm; ue_cm <- d$participantError_cm
say("ratio of means vs mean signed  = ", format(mean(oe) / mean(se_cm)))
say("ratio of means vs mean |signed| = ", format(mean(oe) / mean(abs(se_cm))),
    "   (mean |signed| = ", format(mean(abs(se_cm))), ")")
say("ratio of means vs mean unsigned = ", format(mean(oe) / mean(ue_cm)))
say("per-trial orth/unsigned: M = ", format(mean(oe / ue_cm)), " Mdn = ", format(median(oe / ue_cm)))
say("orth > |signed| in ", sum(oe > abs(se_cm)), "/", nrow(d),
    " (", format(100 * mean(oe > abs(se_cm))), "%)")
pp_oe <- d %>% group_by(participantId) %>% summarise(m = mean(100 * orthogonalError_m), .groups = "drop")
say("per-participant mean orth: range ", format(min(pp_oe$m)), "-", format(max(pp_oe$m)),
    "  Mdn = ", format(median(pp_oe$m)))

# dimensionality-corrected comparison: orthogonal is a 2-D norm, signed is 1-D
say("RMS signed (1-D)          = ", format(sqrt(mean(se_cm^2))), " cm")
say("RMS orthogonal (2-D norm) = ", format(sqrt(mean(oe^2))), " cm")
say("RMS orthogonal per axis   = ", format(sqrt(mean(oe^2) / 2)), " cm")
say("SD of signed (bias removed) = ", format(sd(se_cm)), " cm")

# NaN check on the sqrt form
radicand <- d$participantError_m^2 - d$signedError_m^2
say("n trials with negative radicand (sqrt form -> NaN): ", sum(radicand < 0))
say("min radicand = ", format(min(radicand)))
sqrt_form <- sqrt(pmax(radicand, 0))
say("max |vector form - sqrt form| = ", format(max(abs(d$orthogonalError_m - sqrt_form))), " m")

# ------------------------------------------------------ (5) outlier sensitivity
say("\n=== (5) OUTLIER SENSITIVITY ===")
cut_unsigned <- mean(d$participantError_cm) + 3 * sd(d$participantError_cm)
cut_lo <- mean(se_cm) - 3 * sd(se_cm); cut_hi <- mean(se_cm) + 3 * sd(se_cm)
say("unsigned cut = ", format(cut_unsigned), " cm -> ", sum(d$participantError_cm > cut_unsigned), " trials")
say("signed cuts  = ", format(cut_lo), " to ", format(cut_hi), " cm -> ",
    sum(se_cm < cut_lo | se_cm > cut_hi), " trials")
say("RT < 3 s -> ", sum(d$responseTime_s < 3), " trials")

extreme <- d %>% arrange(responseTime_s) %>%
  select(participantId, trialSequenceNum, soundType, stimulusDisparity_m,
         responseTime_s, participantError_cm, signedError_m, ventriloquistBias,
         n_samples, total_path_length) %>% head(3)
say("\n-- 3 shortest trials --"); show_df(extreme)
biggest <- d %>% arrange(desc(participantError_cm)) %>%
  select(participantId, trialSequenceNum, participantError_cm, responseTime_s) %>% head(5)
say("-- 5 largest errors --"); show_df(biggest)
say("-- most negative signed error --")
show_df(d %>% arrange(signedError_cm) %>%
          select(participantId, trialSequenceNum, signedError_cm) %>% head(3))

subsets <- list(
  full        = d,
  rt_ge_3     = filter(d, responseTime_s >= 3),
  err_3sd     = filter(d, participantError_cm <= cut_unsigned),
  signed_3sd  = filter(d, se_cm >= cut_lo & se_cm <= cut_hi),
  both        = filter(d, responseTime_s >= 3, participantError_cm <= cut_unsigned)
)

fit_one <- function(dat, opt) {
  w <- character(0)
  f <- withCallingHandlers(
    lmer(signedError_m ~ stimulusDisparity_m + soundType +
           (1 + stimulusDisparity_m | participantId),
         data = dat, REML = TRUE, control = lmerControl(optimizer = opt)),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") },
    message = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleMessage") })
  cf <- coef(summary(f))["stimulusDisparity_m", ]
  vc <- as.data.frame(VarCorr(f))
  list(fit = f,
       row = data.frame(
         optimizer = opt, n = nobs(f),
         est = cf[["Estimate"]], se = cf[["Std. Error"]], df = cf[["df"]],
         t = cf[["t value"]], p = cf[["Pr(>|t|)"]],
         resid_sd = sigma(f),
         sd_int = vc$sdcor[vc$grp == "participantId" & vc$var1 == "(Intercept)" & is.na(vc$var2)],
         sd_slope = vc$sdcor[vc$grp == "participantId" & vc$var1 == "stimulusDisparity_m" & is.na(vc$var2)],
         corr_is = vc$sdcor[!is.na(vc$var2)],
         singular = isSingular(f),
         rcond_vcov = rcond(as.matrix(vcov(f))),
         rcond_hess = tryCatch(rcond(f@optinfo$derivs$Hessian), error = function(e) NA_real_),
         min_eig_hess = tryCatch(min(eigen(f@optinfo$derivs$Hessian, only.values = TRUE)$values),
                                 error = function(e) NA_real_),
         warnings = if (length(w)) paste(w, collapse = " | ") else "none",
         stringsAsFactors = FALSE))
}

opts <- c("bobyqa", "Nelder_Mead", "nloptwrap")
sens <- list()
for (nm in names(subsets)) {
  for (o in opts) {
    r <- fit_one(subsets[[nm]], o)
    sens[[paste(nm, o)]] <- cbind(subset = nm, r$row)
  }
}
sens <- bind_rows(sens)
base_se <- sens$se[sens$subset == "full" & sens$optimizer == "bobyqa"]
base_est <- sens$est[sens$subset == "full" & sens$optimizer == "bobyqa"]
sens <- sens %>% mutate(pct_change = 100 * (est - base_est) / base_est,
                        shift_in_base_se = (est - base_est) / base_se)
say("\n-- primary-model refits across subsets and optimizers --")
show_df(sens %>% select(subset, optimizer, n, est, se, df, t, p, pct_change,
                        shift_in_base_se, singular, rcond_vcov, min_eig_hess, warnings), 6)
write.csv(sens, file.path(res, "verify_outlier_sensitivity.csv"), row.names = FALSE)

# Wald CIs and profile CIs on the bobyqa full fit
f_full <- fit_one(d, "bobyqa")$fit
ci_w <- confint(f_full, parm = "stimulusDisparity_m", method = "Wald")
say("full-fit Wald CI: ", format(ci_w[1]), " to ", format(ci_w[2]))
ci_p <- suppressMessages(suppressWarnings(
  confint(f_full, parm = "stimulusDisparity_m", method = "profile")))
say("full-fit profile CI: ", format(ci_p[1]), " to ", format(ci_p[2]))
say("anova F test: ")
show_df(cbind(term = rownames(anova(f_full)), as.data.frame(anova(f_full))), 6)

# Do the subset CIs actually overlap the full-data estimate?
sb <- sens %>% filter(optimizer == "bobyqa") %>%
  mutate(lo = est - 1.96 * se, hi = est + 1.96 * se,
         covers_full_est = base_est >= lo & base_est <= hi)
say("\n-- Wald CIs per subset (bobyqa) --")
show_df(sb %>% select(subset, n, est, lo, hi, covers_full_est), 6)

say("\n=== END ===")
close(con)
