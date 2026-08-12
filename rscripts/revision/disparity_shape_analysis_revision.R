# Shape of the audiovisual bias as a function of stimulus disparity
#
# Question: does the ventriloquist-type bias scale proportionally with the
# sound-to-flash disparity across the tested range (0.15-0.70 m), or does it
# saturate / break down beyond some point?
#
# Three converging approaches, all on the signed error (component of the
# response displacement along the sound-to-flash axis):
#   (1) GAMM smooth vs linear   - is the disparity effect curved?
#   (2) segmented (broken-stick) LMM with a grid search over the breakpoint
#   (3) direct proportionality test on the bias ratio (signed error / disparity)
# Plus a description of the tested range in egocentric angular terms.
#
# Input : results_revision/analysis_df_revision.rds  (built by build_analysis_df_revision.R)
# Output: results_revision/disparity_shape_*.csv and disparity_shape_log.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
  library(lme4)
  library(lmerTest)
  library(quantreg)
})

set.seed(20260805)

proj_dir    <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
results_dir <- file.path(proj_dir, "results_revision")
log_path    <- file.path(results_dir, "disparity_shape_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

cat("Disparity-response shape analysis\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "| mgcv", as.character(packageVersion("mgcv")),
    "| lme4", as.character(packageVersion("lme4")),
    "| quantreg", as.character(packageVersion("quantreg")), "\n\n")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds"))
df <- df %>%
  mutate(participantId = factor(participantId),
         soundType     = factor(soundType))
df <- as.data.frame(df)

cat("N trials:", nrow(df), " N participants:", nlevels(df$participantId), "\n")
cat("Disparity: ", length(unique(df$stimulusDisparity_m)), " distinct values, range ",
    sprintf("%.3f-%.3f m", min(df$stimulusDisparity_m), max(df$stimulusDisparity_m)), "\n\n", sep = "")

# ---------------------------------------------------------------------------
# 0. Descriptives: binned means and medians of signed error and bias ratio
# ---------------------------------------------------------------------------

cat("== 0. Binned descriptives ==\n")

bin_summary <- function(data, bin_var, label) {
  data %>%
    group_by(bin = .data[[bin_var]]) %>%
    summarise(n              = n(),
              disp_mean_m    = mean(stimulusDisparity_m),
              signed_mean_cm = 100 * mean(signedError_m),
              signed_sd_cm   = 100 * sd(signedError_m),
              signed_med_cm  = 100 * median(signedError_m),
              bias_mean_pct  = 100 * mean(ventriloquistBias),
              bias_sd_pct    = 100 * sd(ventriloquistBias),
              bias_med_pct   = 100 * median(ventriloquistBias),
              bias_iqr_pct   = 100 * IQR(ventriloquistBias),
              .groups = "drop") %>%
    mutate(binning = label, .before = 1)
}

df <- df %>%
  mutate(disp_quintile = cut(stimulusDisparity_m,
                             breaks = quantile(stimulusDisparity_m, seq(0, 1, 0.2)),
                             include.lowest = TRUE))

bins_tbl <- bind_rows(
  bin_summary(df, "disparityRange", "published 3-bin"),
  bin_summary(df, "disp_quintile",  "quintile")
)
print(as.data.frame(bins_tbl), digits = 4)
write.csv(bins_tbl, file.path(results_dir, "disparity_shape_bins.csv"), row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# 1. GAMM: smooth vs linear disparity effect on signed error
# ---------------------------------------------------------------------------

cat("== 1. GAMM smooth vs linear (signed error) ==\n")

rand_spec <- list(participantId = ~ 1 + stimulusDisparity_m)

# REML fits (default) for the smooth's edf and approximate p-value
gamm_smooth <- gamm(signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType,
                    random = rand_spec, data = df)
gamm_linear <- gamm(signedError_m ~ stimulusDisparity_m + soundType,
                    random = rand_spec, data = df)

cat("-- smooth model (REML) --\n")
print(summary(gamm_smooth$gam))
cat("-- linear model (REML) --\n")
print(summary(gamm_linear$gam))

edf_smooth <- summary(gamm_smooth$gam)$s.table[1, "edf"]
p_smooth   <- summary(gamm_smooth$gam)$s.table[1, "p-value"]
f_smooth   <- summary(gamm_smooth$gam)$s.table[1, "F"]

# ML fits for a like-for-like AIC comparison (fixed-effect structures differ)
gamm_smooth_ml <- gamm(signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType,
                       random = rand_spec, data = df, method = "ML")
gamm_linear_ml <- gamm(signedError_m ~ stimulusDisparity_m + soundType,
                       random = rand_spec, data = df, method = "ML")

aic_gamm <- data.frame(
  model  = c("GAMM s(disparity, k=6)", "LMM linear disparity"),
  df     = c(attr(logLik(gamm_smooth_ml$lme), "df"), attr(logLik(gamm_linear_ml$lme), "df")),
  logLik = c(as.numeric(logLik(gamm_smooth_ml$lme)), as.numeric(logLik(gamm_linear_ml$lme))),
  AIC    = c(AIC(gamm_smooth_ml$lme), AIC(gamm_linear_ml$lme)),
  BIC    = c(BIC(gamm_smooth_ml$lme), BIC(gamm_linear_ml$lme))
)
aic_gamm$dAIC <- aic_gamm$AIC - min(aic_gamm$AIC)
cat("-- ML AIC comparison (gamm lme components) --\n")
print(aic_gamm, digits = 6)

# Cross-check with a pure-GAM formulation (random effects as smooths, REML)
gam_smooth <- gam(signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType +
                    s(participantId, bs = "re") +
                    s(participantId, stimulusDisparity_m, bs = "re"),
                  data = df, method = "REML")
gam_linear <- gam(signedError_m ~ stimulusDisparity_m + soundType +
                    s(participantId, bs = "re") +
                    s(participantId, stimulusDisparity_m, bs = "re"),
                  data = df, method = "REML")
cat("\n-- cross-check: gam() with random-effect smooths --\n")
cat("edf of s(disparity):", round(summary(gam_smooth)$s.table[1, "edf"], 3),
    " p =", signif(summary(gam_smooth)$s.table[1, "p-value"], 3), "\n")
cat("AIC smooth:", round(AIC(gam_smooth), 3), " AIC linear:", round(AIC(gam_linear), 3),
    " dAIC(smooth - linear):", round(AIC(gam_smooth) - AIC(gam_linear), 3), "\n\n")

# Parametric alternatives to a straight line (ML LMMs), as a further curvature check
df <- df %>% mutate(disp2 = stimulusDisparity_m^2,
                    disp_log = log(stimulusDisparity_m),
                    disp_sqrt = sqrt(stimulusDisparity_m))
lmm_lin  <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                   (1 + stimulusDisparity_m | participantId), data = df, REML = FALSE)
lmm_quad <- lmer(signedError_m ~ stimulusDisparity_m + disp2 + soundType +
                   (1 + stimulusDisparity_m | participantId), data = df, REML = FALSE)
lmm_log  <- lmer(signedError_m ~ disp_log + soundType +
                   (1 + stimulusDisparity_m | participantId), data = df, REML = FALSE)
lmm_sqrt <- lmer(signedError_m ~ disp_sqrt + soundType +
                   (1 + stimulusDisparity_m | participantId), data = df, REML = FALSE)

shape_tbl <- data.frame(
  model  = c("linear", "quadratic", "log", "sqrt"),
  df     = sapply(list(lmm_lin, lmm_quad, lmm_log, lmm_sqrt), function(m) attr(logLik(m), "df")),
  logLik = sapply(list(lmm_lin, lmm_quad, lmm_log, lmm_sqrt), function(m) as.numeric(logLik(m))),
  AIC    = sapply(list(lmm_lin, lmm_quad, lmm_log, lmm_sqrt), AIC),
  BIC    = sapply(list(lmm_lin, lmm_quad, lmm_log, lmm_sqrt), BIC)
)
shape_tbl$dAIC <- shape_tbl$AIC - min(shape_tbl$AIC)
cat("-- parametric shape alternatives (ML LMMs, same random structure) --\n")
print(shape_tbl, digits = 6)
cat("quadratic term:\n")
print(summary(lmm_quad)$coefficients["disp2", , drop = FALSE], digits = 4)
cat("LRT linear vs quadratic: chisq =",
    signif(2 * (logLik(lmm_quad) - logLik(lmm_lin)), 4), " df = 1  p =",
    signif(pchisq(2 * as.numeric(logLik(lmm_quad) - logLik(lmm_lin)), 1, lower.tail = FALSE), 4), "\n")
cat("linear model singular fit:", isSingular(lmm_lin), "\n\n")

cat("-- linear LMM fixed effects (REML, lmerTest) --\n")
lmm_lin_reml <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                       (1 + stimulusDisparity_m | participantId), data = df)
print(summary(lmm_lin_reml)$coefficients, digits = 4)
ci_lin <- confint(lmm_lin_reml, parm = "stimulusDisparity_m", method = "profile")
cat("slope 95% profile CI:", sprintf("%.4f to %.4f", ci_lin[1, 1], ci_lin[1, 2]), "\n\n")

write.csv(rbind(
  data.frame(comparison = "GAMM ML", aic_gamm),
  data.frame(comparison = "parametric shapes (ML LMM)",
             model = shape_tbl$model, df = shape_tbl$df, logLik = shape_tbl$logLik,
             AIC = shape_tbl$AIC, BIC = shape_tbl$BIC, dAIC = shape_tbl$dAIC)
), file.path(results_dir, "disparity_shape_model_comparison.csv"), row.names = FALSE)

# Fitted smooth on a prediction grid, for plotting / inspection in the paper
pred_grid <- data.frame(stimulusDisparity_m = seq(0.15, 0.70, by = 0.005),
                        soundType = factor("Drum", levels = levels(df$soundType)))
pr_s <- predict(gamm_smooth$gam, newdata = pred_grid, se.fit = TRUE)
pr_l <- predict(gamm_linear$gam, newdata = pred_grid, se.fit = TRUE)
smooth_fit <- data.frame(stimulusDisparity_m = pred_grid$stimulusDisparity_m,
                         gamm_fit = as.numeric(pr_s$fit), gamm_se = as.numeric(pr_s$se.fit),
                         lin_fit  = as.numeric(pr_l$fit), lin_se  = as.numeric(pr_l$se.fit))
write.csv(smooth_fit, file.path(results_dir, "disparity_shape_gamm_fit.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Segmented (broken-stick) LMM: grid search over the breakpoint
# ---------------------------------------------------------------------------

cat("== 2. Segmented LMM, breakpoint grid search 0.20-0.60 m (1 cm steps) ==\n")

breaks <- seq(0.20, 0.60, by = 0.01)

fit_segmented <- function(c_val, data) {
  data$hinge <- pmax(data$stimulusDisparity_m - c_val, 0)
  lmer(signedError_m ~ stimulusDisparity_m + hinge + soundType +
         (1 + stimulusDisparity_m | participantId), data = data, REML = FALSE)
}

seg_fits <- lapply(breaks, fit_segmented, data = df)
ll_lin   <- as.numeric(logLik(lmm_lin))

seg_profile <- data.frame(
  breakpoint_m = breaks,
  logLik       = sapply(seg_fits, function(m) as.numeric(logLik(m))),
  AIC          = sapply(seg_fits, AIC),
  BIC          = sapply(seg_fits, BIC),
  slope1       = sapply(seg_fits, function(m) fixef(m)["stimulusDisparity_m"]),
  hinge        = sapply(seg_fits, function(m) fixef(m)["hinge"]),
  singular     = sapply(seg_fits, isSingular)
)
seg_profile$slope2   <- seg_profile$slope1 + seg_profile$hinge
seg_profile$LRT_stat <- 2 * (seg_profile$logLik - ll_lin)
print(seg_profile, digits = 5, row.names = FALSE)
write.csv(seg_profile, file.path(results_dir, "disparity_shape_breakpoint_profile.csv"),
          row.names = FALSE)

best_i   <- which.max(seg_profile$logLik)
best_c   <- seg_profile$breakpoint_m[best_i]
best_fit <- seg_fits[[best_i]]
lrt_obs  <- seg_profile$LRT_stat[best_i]

cat("\nBest breakpoint:", best_c, "m\n")
cat("logLik linear:", round(ll_lin, 4), " logLik segmented:",
    round(seg_profile$logLik[best_i], 4), "\n")
cat("AIC linear:", round(AIC(lmm_lin), 3), " AIC segmented:", round(AIC(best_fit), 3),
    " dAIC:", round(AIC(best_fit) - AIC(lmm_lin), 3), "\n")
cat("BIC linear:", round(BIC(lmm_lin), 3), " BIC segmented:", round(BIC(best_fit), 3),
    " dBIC:", round(BIC(best_fit) - BIC(lmm_lin), 3), "\n")
cat("Naive LRT (1 df, breakpoint treated as known): chisq =", round(lrt_obs, 4),
    " p =", signif(pchisq(lrt_obs, 1, lower.tail = FALSE), 4), "\n")
cat("NOTE: the breakpoint is not identified under the null (no hinge), so this\n")
cat("      naive p-value is anticonservative. A parametric bootstrap of the\n")
cat("      maximised LRT follows.\n")
cat("-- segmented model at the best breakpoint --\n")
print(summary(best_fit)$coefficients, digits = 4)
cat("slope before:", round(seg_profile$slope1[best_i], 4),
    " slope after:", round(seg_profile$slope2[best_i], 4), "\n\n")

# Parametric bootstrap: simulate under the linear (null) model, redo the grid
# search on each replicate, and compare the observed maximised LRT to the null
# distribution of maximised LRTs.
n_boot <- 1000
cat("Parametric bootstrap of the maximised LRT, B =", n_boot, "\n")
sim_y <- simulate(lmm_lin, nsim = n_boot, seed = 20260805)

# Singular random-slope fits and lme4 gradient warnings are expected on
# replicates simulated with a small slope variance; they are counted, not printed.
# No RNG is used inside the loop (all responses are pre-simulated), so the result
# does not depend on how the replicates are distributed over cores.
refit_quiet <- function(model, y) {
  n_warn <- 0L
  fit <- withCallingHandlers(
    suppressMessages(refit(model, newresp = y)),
    warning = function(cond) { n_warn <<- n_warn + 1L; invokeRestart("muffleWarning") })
  list(logLik = as.numeric(logLik(fit)), singular = isSingular(fit), n_warn = n_warn)
}

boot_one <- function(b) {
  y   <- sim_y[[b]]
  r0  <- refit_quiet(lmm_lin, y)
  rs  <- lapply(seg_fits, refit_quiet, y = y)
  lls <- vapply(rs, `[[`, numeric(1), "logLik")
  data.frame(replicate = b,
             max_LRT = 2 * (max(lls) - r0$logLik),
             argmax_breakpoint_m = breaks[which.max(lls)],
             n_singular = sum(vapply(rs, `[[`, logical(1), "singular")),
             n_warn = r0$n_warn + sum(vapply(rs, `[[`, integer(1), "n_warn")))
}

n_cores <- max(1L, min(4L, parallel::detectCores() - 1L))
t_start <- Sys.time()
boot_tbl <- do.call(rbind, parallel::mclapply(seq_len(n_boot), boot_one, mc.cores = n_cores))
boot_max <- boot_tbl$max_LRT
boot_c   <- boot_tbl$argmax_breakpoint_m

cat("bootstrap time:", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 2),
    "min on", n_cores, "cores\n")
cat("singular segmented fits across replicates:", sum(boot_tbl$n_singular),
    "of", n_boot * length(breaks), "\n")
cat("lme4 convergence warnings across replicates:", sum(boot_tbl$n_warn),
    "of", n_boot * (length(breaks) + 1), "fits;",
    sum(boot_tbl$n_warn > 0), "replicates affected\n")

p_boot <- (1 + sum(boot_max >= lrt_obs)) / (n_boot + 1)
cat("observed maximised LRT:", round(lrt_obs, 4), "\n")
cat("null distribution quantiles (50/90/95/99%):",
    paste(round(quantile(boot_max, c(0.5, 0.9, 0.95, 0.99)), 3), collapse = ", "), "\n")
cat("bootstrap p-value:", signif(p_boot, 4), "\n")
cat("under a true straight line the grid search still returns a breakpoint;\n")
cat("distribution of the selected breakpoint (fraction at the grid edges 0.20/0.60):",
    round(mean(boot_c %in% range(breaks)), 3), "\n\n")

write.csv(boot_tbl, file.path(results_dir, "disparity_shape_breakpoint_bootstrap.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Direct proportionality test on the bias ratio
# ---------------------------------------------------------------------------

cat("== 3. Proportionality test: ventriloquistBias ~ disparity ==\n")
cat("Exact proportional scaling -> zero slope; saturation -> negative slope.\n\n")

bias_lmm <- lmer(ventriloquistBias ~ stimulusDisparity_m + (1 | participantId), data = df)
print(summary(bias_lmm)$coefficients, digits = 4)
ci_bias <- confint(bias_lmm, parm = c("(Intercept)", "stimulusDisparity_m"),
                   method = "profile")
print(ci_bias, digits = 4)
cat("singular fit:", isSingular(bias_lmm), "\n")
cat("slope (bias units per m):", round(fixef(bias_lmm)["stimulusDisparity_m"], 4),
    " -> ", round(100 * fixef(bias_lmm)["stimulusDisparity_m"], 2),
    "percentage points of bias per 1 m of disparity;\n across the tested 0.55 m range that is ",
    round(100 * fixef(bias_lmm)["stimulusDisparity_m"] * 0.55, 2), " percentage points.\n\n", sep = "")

# Robust version A: symmetric 5% trimming of the bias ratio
lo <- quantile(df$ventriloquistBias, 0.05)
hi <- quantile(df$ventriloquistBias, 0.95)
df_trim <- df %>% filter(ventriloquistBias >= lo, ventriloquistBias <= hi)
cat("Robust A: 5/95 percentile trimming of the bias ratio (cut at",
    sprintf("%.3f and %.3f", lo, hi), "); n =", nrow(df_trim), "of", nrow(df), "\n")
bias_lmm_trim <- lmer(ventriloquistBias ~ stimulusDisparity_m + (1 | participantId), data = df_trim)
print(summary(bias_lmm_trim)$coefficients, digits = 4)
ci_trim <- confint(bias_lmm_trim, parm = "stimulusDisparity_m", method = "profile")
cat("trimmed slope 95% CI:", sprintf("%.4f to %.4f", ci_trim[1, 1], ci_trim[1, 2]), "\n\n")

# Robust version B: quantile (L1) regression with participant fixed effects.
# tau = 0.5 is the median slope; tau = 0.25 / 0.75 show whether the spread of the
# bias ratio changes with disparity.
rq_slope <- function(tau, data) {
  fit <- rq(ventriloquistBias ~ stimulusDisparity_m + participantId, tau = tau, data = data)
  s   <- summary(fit, se = "boot", R = 2000, bsmethod = "xy")
  cf  <- s$coefficients["stimulusDisparity_m", ]
  data.frame(tau = tau, slope = cf[["Value"]], se = cf[["Std. Error"]],
             t = cf[["t value"]], p = cf[["Pr(>|t|)"]],
             ci_low = cf[["Value"]] - 1.96 * cf[["Std. Error"]],
             ci_high = cf[["Value"]] + 1.96 * cf[["Std. Error"]])
}
cat("Robust B: quantile regression (quantreg::rq, participant fixed effects),\n")
cat("          bootstrap SE, R = 2000. Note rq may report a nonunique solution\n")
cat("          for this design; the slope is stable to the reported precision.\n")
rq_tbl <- bind_rows(lapply(c(0.25, 0.5, 0.75), rq_slope, data = df))
print(rq_tbl, digits = 4, row.names = FALSE)
cat("\n")

# Same test on the signed error but with disparity as the offset: an equivalent way
# of asking whether the slope on signed error differs from proportionality is to test
# whether the intercept of the signed-error model is zero (proportional through origin).
cat("Intercept of the signed-error LMM (zero implies a line through the origin,\n")
cat("i.e. exact proportionality with slope = mean bias):\n")
print(summary(lmm_lin_reml)$coefficients["(Intercept)", , drop = FALSE], digits = 4)
ci_int <- confint(lmm_lin_reml, parm = "(Intercept)", method = "profile")
cat("intercept 95% profile CI:", sprintf("%.4f to %.4f m", ci_int[1, 1], ci_int[1, 2]), "\n\n")

bias_slopes <- data.frame(
  method  = c("LMM (1|participant)", "LMM, 5/95 trimmed",
              "rq tau=0.25 (FE)", "rq tau=0.50 (FE)", "rq tau=0.75 (FE)"),
  n       = c(nrow(df), nrow(df_trim), nrow(df), nrow(df), nrow(df)),
  slope   = c(fixef(bias_lmm)["stimulusDisparity_m"],
              fixef(bias_lmm_trim)["stimulusDisparity_m"], rq_tbl$slope),
  ci_low  = c(ci_bias["stimulusDisparity_m", 1], ci_trim[1, 1], rq_tbl$ci_low),
  ci_high = c(ci_bias["stimulusDisparity_m", 2], ci_trim[1, 2], rq_tbl$ci_high)
)
print(bias_slopes, digits = 4, row.names = FALSE)
write.csv(bias_slopes, file.path(results_dir, "disparity_shape_bias_slopes.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Egocentric angular size of the tested disparities
# ---------------------------------------------------------------------------

cat("== 4. Tested range in egocentric angular terms ==\n")

ang_vars <- c("angDisp_onset_deg", "angDisp_mean_deg", "angDisp_median_deg",
              "angDisp_min_deg", "angDisp_max_deg")
ang_tbl <- lapply(ang_vars, function(v) {
  x <- df[[v]]
  data.frame(variable = v, n = sum(!is.na(x)), mean = mean(x), sd = sd(x),
             min = min(x), q05 = quantile(x, 0.05), q25 = quantile(x, 0.25),
             median = median(x), q75 = quantile(x, 0.75), q95 = quantile(x, 0.95),
             max = max(x),
             pct_over_30deg = 100 * mean(x > 30),
             pct_over_15deg = 100 * mean(x > 15))
}) %>% bind_rows()
rownames(ang_tbl) <- NULL
print(ang_tbl, digits = 4)
write.csv(ang_tbl, file.path(results_dir, "disparity_shape_angular_summary.csv"), row.names = FALSE)

cat("\nAngular disparity at onset by published disparity bin:\n")
ang_by_bin <- df %>%
  group_by(disparityRange) %>%
  summarise(n = n(),
            onset_median = median(angDisp_onset_deg),
            onset_q25 = quantile(angDisp_onset_deg, 0.25),
            onset_q75 = quantile(angDisp_onset_deg, 0.75),
            onset_pct_over_30 = 100 * mean(angDisp_onset_deg > 30),
            mean_median = median(angDisp_mean_deg),
            mean_q25 = quantile(angDisp_mean_deg, 0.25),
            mean_q75 = quantile(angDisp_mean_deg, 0.75),
            mean_pct_over_30 = 100 * mean(angDisp_mean_deg > 30),
            .groups = "drop")
print(as.data.frame(ang_by_bin), digits = 4)
write.csv(ang_by_bin, file.path(results_dir, "disparity_shape_angular_by_bin.csv"), row.names = FALSE)

cat("\nCorrelation of metric disparity with angular disparity:\n")
cat("  onset: r =", round(cor(df$stimulusDisparity_m, df$angDisp_onset_deg), 3),
    " Spearman =", round(cor(df$stimulusDisparity_m, df$angDisp_onset_deg, method = "spearman"), 3), "\n")
cat("  mean : r =", round(cor(df$stimulusDisparity_m, df$angDisp_mean_deg), 3),
    " Spearman =", round(cor(df$stimulusDisparity_m, df$angDisp_mean_deg, method = "spearman"), 3), "\n")

# Is the bias curved in angular rather than metric terms? Same smooth/linear check.
cat("\nGAMM on angular disparity at onset (secondary check):\n")
gamm_ang <- gamm(signedError_m ~ s(angDisp_onset_deg, k = 6) + soundType,
                 random = list(participantId = ~ 1), data = df)
print(summary(gamm_ang$gam)$s.table, digits = 4)

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
