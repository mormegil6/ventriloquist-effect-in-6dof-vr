# Independent verification of the disparity-shape analysis
#
# Re-fits every headline model of disparity_shape_analysis_revision.R from the
# prepared trial-level data, and adds the checks that the original script does
# not run:
#   A. data integrity (no silent filtering, ratio definition, bin counts)
#   B. GAMM smooth vs linear, incl. basis / k / select sensitivity
#   C. parametric shape family, and how much curvature the data actually exclude
#   D. segmented grid search, random-effects conditioning, optimizer agreement
#   E. proportionality test, with and without soundType, Wald vs profile CIs
#   F. residual heteroscedasticity (validity of the Gaussian bootstrap null)
#   G. angular disparity descriptives and the shape of the angular smooth
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/disparity_shape_verification_log.txt and *_verify_*.csv

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
log_path    <- file.path(results_dir, "disparity_shape_verification_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

cat("Independent verification of the disparity-shape analysis\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R:", R.version.string, "| mgcv", as.character(packageVersion("mgcv")),
    "| lme4", as.character(packageVersion("lme4")),
    "| quantreg", as.character(packageVersion("quantreg")), "\n\n")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds"))
df <- df %>%
  mutate(participantId = factor(participantId),
         soundType     = factor(soundType))
df <- as.data.frame(df)

# ---------------------------------------------------------------------------
# A. Data integrity
# ---------------------------------------------------------------------------

cat("== A. Data integrity ==\n")
cat("rows:", nrow(df), " participants:", nlevels(df$participantId),
    " soundType levels:", nlevels(df$soundType), "\n")
key_vars <- c("signedError_m", "stimulusDisparity_m", "ventriloquistBias",
              "soundType", "participantId", "angDisp_onset_deg", "angDisp_mean_deg")
cat("NAs in model variables:\n")
print(sapply(df[key_vars], function(x) sum(is.na(x))))
cat("distinct disparity values:", length(unique(df$stimulusDisparity_m)),
    sprintf(" range %.3f-%.3f m\n", min(df$stimulusDisparity_m), max(df$stimulusDisparity_m)))
cat("bias ratio == signedError / disparity?  max abs deviation:",
    max(abs(df$ventriloquistBias - df$signedError_m / df$stimulusDisparity_m)), "\n")
cat("trials per published bin:\n"); print(table(df$disparityRange))
cat("trials per participant: min", min(table(df$participantId)),
    " max", max(table(df$participantId)), "\n")
cat("mean bias %:", round(100 * mean(df$ventriloquistBias), 3),
    " SD:", round(100 * sd(df$ventriloquistBias), 3), "\n")
cat("bias range:", round(range(df$ventriloquistBias), 3), "\n\n")

bin_tbl <- df %>%
  group_by(disparityRange) %>%
  summarise(n = n(), disp_mean = mean(stimulusDisparity_m),
            signed_mean_cm = 100 * mean(signedError_m),
            bias_mean_pct = 100 * mean(ventriloquistBias),
            bias_med_pct = 100 * median(ventriloquistBias),
            bias_sd_pct = 100 * sd(ventriloquistBias), .groups = "drop")
print(as.data.frame(bin_tbl), digits = 5)
cat("\n")

# ---------------------------------------------------------------------------
# B. GAMM smooth vs linear, and sensitivity to the smoothing set-up
# ---------------------------------------------------------------------------

cat("== B. GAMM smooth vs linear ==\n")
rand_spec <- list(participantId = ~ 1 + stimulusDisparity_m)

gamm_smooth <- gamm(signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType,
                    random = rand_spec, data = df)
st <- summary(gamm_smooth$gam)$s.table
cat("REML smooth: edf =", round(st[1, "edf"], 4), " Ref.df =", st[1, "Ref.df"],
    " F =", round(st[1, "F"], 3), " p =", signif(st[1, "p-value"], 4), "\n")

gamm_smooth_ml <- gamm(signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType,
                       random = rand_spec, data = df, method = "ML")
gamm_linear_ml <- gamm(signedError_m ~ stimulusDisparity_m + soundType,
                       random = rand_spec, data = df, method = "ML")
gamm_aic <- data.frame(
  model  = c("GAMM smooth", "GAMM linear"),
  df     = c(attr(logLik(gamm_smooth_ml$lme), "df"), attr(logLik(gamm_linear_ml$lme), "df")),
  logLik = c(as.numeric(logLik(gamm_smooth_ml$lme)), as.numeric(logLik(gamm_linear_ml$lme))),
  AIC    = c(AIC(gamm_smooth_ml$lme), AIC(gamm_linear_ml$lme)),
  BIC    = c(BIC(gamm_smooth_ml$lme), BIC(gamm_linear_ml$lme)))
print(gamm_aic, digits = 7, row.names = FALSE)

gam_smooth <- gam(signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType +
                    s(participantId, bs = "re") +
                    s(participantId, stimulusDisparity_m, bs = "re"),
                  data = df, method = "REML")
gam_linear <- gam(signedError_m ~ stimulusDisparity_m + soundType +
                    s(participantId, bs = "re") +
                    s(participantId, stimulusDisparity_m, bs = "re"),
                  data = df, method = "REML")
cat("gam() cross-check: edf =", round(summary(gam_smooth)$s.table[1, "edf"], 4),
    " AIC smooth", round(AIC(gam_smooth), 3), " AIC linear", round(AIC(gam_linear), 3),
    " dAIC", round(AIC(gam_smooth) - AIC(gam_linear), 3), "\n")

# sensitivity: basis dimension, basis type, and shrinkage smoothers
sens <- lapply(list(
  list(lab = "tp k=6",        f = signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType),
  list(lab = "tp k=10",       f = signedError_m ~ s(stimulusDisparity_m, k = 10) + soundType),
  list(lab = "tp k=20",       f = signedError_m ~ s(stimulusDisparity_m, k = 20) + soundType),
  list(lab = "cr k=10",       f = signedError_m ~ s(stimulusDisparity_m, bs = "cr", k = 10) + soundType),
  list(lab = "ts k=10 (shrink)", f = signedError_m ~ s(stimulusDisparity_m, bs = "ts", k = 10) + soundType)
), function(z) {
  m <- gam(update(z$f, . ~ . + s(participantId, bs = "re") +
                    s(participantId, stimulusDisparity_m, bs = "re")),
           data = df, method = "REML")
  s <- summary(m)$s.table
  data.frame(spec = z$lab, edf = s[1, "edf"], F = s[1, "F"],
             p = s[1, "p-value"], AIC = AIC(m))
}) %>% bind_rows()
cat("\nsmooth sensitivity (gam with RE smooths, REML):\n")
print(sens, digits = 4, row.names = FALSE)
cat("AIC of the matching linear gam:", round(AIC(gam_linear), 3), "\n\n")

# ---------------------------------------------------------------------------
# C. Parametric shape family, and the size of curvature still compatible
# ---------------------------------------------------------------------------

cat("== C. Parametric shapes and curvature that is NOT excluded ==\n")
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
  BIC    = sapply(list(lmm_lin, lmm_quad, lmm_log, lmm_sqrt), BIC),
  singular = sapply(list(lmm_lin, lmm_quad, lmm_log, lmm_sqrt), isSingular))
shape_tbl$dAIC_vs_best <- shape_tbl$AIC - min(shape_tbl$AIC)
print(shape_tbl, digits = 7, row.names = FALSE)
cat("quadratic term (ML fit):\n")
print(summary(lmm_quad)$coefficients["disp2", , drop = FALSE], digits = 5)
lrt_q <- 2 * as.numeric(logLik(lmm_quad) - logLik(lmm_lin))
cat("LRT linear vs quadratic: chisq =", signif(lrt_q, 5),
    " p =", signif(pchisq(lrt_q, 1, lower.tail = FALSE), 5), "\n")

# How large a departure from a straight line is still inside the quadratic CI?
q_ci <- confint(lmm_quad, parm = "disp2", method = "Wald")
x1 <- min(df$stimulusDisparity_m); x2 <- max(df$stimulusDisparity_m); xm <- (x1 + x2) / 2
sag <- function(b2) b2 * (xm^2 - (x1^2 + x2^2) / 2)   # deviation of b2*x^2 from its chord
cat("\nquadratic coefficient 95% Wald CI:", sprintf("%.4f to %.4f", q_ci[1, 1], q_ci[1, 2]), "\n")
cat("implied max deviation from the straight line over 0.15-0.70 m (mid-range sag):\n")
cat("  at CI lower bound:", sprintf("%.2f cm", 100 * sag(q_ci[1, 1])),
    "  at point estimate:", sprintf("%.2f cm", 100 * sag(fixef(lmm_quad)["disp2"])),
    "  at CI upper bound:", sprintf("%.2f cm", 100 * sag(q_ci[1, 2])), "\n")
cat("for reference the fitted linear rise across the range is",
    sprintf("%.2f cm", 100 * fixef(lmm_lin)["stimulusDisparity_m"] * (x2 - x1)), "\n")

# how different are the log/sqrt fits from the line in cm across the range?
pg <- data.frame(stimulusDisparity_m = seq(x1, x2, length.out = 201))
pg$soundType <- factor(levels(df$soundType)[1], levels = levels(df$soundType))
pg$disp2 <- pg$stimulusDisparity_m^2
pg$disp_log <- log(pg$stimulusDisparity_m)
pg$disp_sqrt <- sqrt(pg$stimulusDisparity_m)
pl <- predict(lmm_lin,  newdata = pg, re.form = NA)
pv <- lapply(list(quadratic = lmm_quad, log = lmm_log, sqrt = lmm_sqrt),
             function(m) predict(m, newdata = pg, re.form = NA))
cat("\nmax |fitted - linear fitted| across the range (population level):\n")
for (nm in names(pv)) cat(" ", nm, sprintf("%.2f cm", 100 * max(abs(pv[[nm]] - pl))), "\n")

cat("\nlinear LMM (REML, lmerTest):\n")
lmm_lin_reml <- lmer(signedError_m ~ stimulusDisparity_m + soundType +
                       (1 + stimulusDisparity_m | participantId), data = df)
print(summary(lmm_lin_reml)$coefficients, digits = 5)
ci_lin <- confint(lmm_lin_reml, parm = c("(Intercept)", "stimulusDisparity_m"), method = "profile")
print(ci_lin, digits = 5)
cat("\n")

# ---------------------------------------------------------------------------
# D. Segmented grid search + conditioning and optimizer checks
# ---------------------------------------------------------------------------

cat("== D. Segmented LMM ==\n")
breaks <- seq(0.20, 0.60, by = 0.01)
fit_segmented <- function(c_val, data) {
  data$hinge <- pmax(data$stimulusDisparity_m - c_val, 0)
  lmer(signedError_m ~ stimulusDisparity_m + hinge + soundType +
         (1 + stimulusDisparity_m | participantId), data = data, REML = FALSE)
}
seg_fits <- lapply(breaks, fit_segmented, data = df)
ll_lin <- as.numeric(logLik(lmm_lin))
seg_profile <- data.frame(
  breakpoint_m = breaks,
  logLik = sapply(seg_fits, function(m) as.numeric(logLik(m))),
  AIC = sapply(seg_fits, AIC), BIC = sapply(seg_fits, BIC),
  slope1 = sapply(seg_fits, function(m) fixef(m)["stimulusDisparity_m"]),
  hinge = sapply(seg_fits, function(m) fixef(m)["hinge"]),
  singular = sapply(seg_fits, isSingular))
seg_profile$slope2 <- seg_profile$slope1 + seg_profile$hinge
seg_profile$LRT <- 2 * (seg_profile$logLik - ll_lin)
best_i <- which.max(seg_profile$logLik)
cat("logLik linear:", round(ll_lin, 6), "\n")
cat("profile range:", sprintf("%.4f to %.4f (span %.4f)",
    min(seg_profile$logLik), max(seg_profile$logLik),
    diff(range(seg_profile$logLik))), "\n")
cat("best breakpoint:", seg_profile$breakpoint_m[best_i],
    " logLik", round(seg_profile$logLik[best_i], 6),
    " LRT", round(seg_profile$LRT[best_i], 5), "\n")
cat("second best:", seg_profile$breakpoint_m[order(-seg_profile$logLik)[2]],
    " logLik", round(sort(seg_profile$logLik, decreasing = TRUE)[2], 6), "\n")
cat("dAIC (segmented - linear):", round(AIC(seg_fits[[best_i]]) - AIC(lmm_lin), 4),
    " dBIC:", round(BIC(seg_fits[[best_i]]) - BIC(lmm_lin), 4), "\n")
cat("any singular in grid:", any(seg_profile$singular), "\n")
cat("segmented coefficients at the best breakpoint:\n")
print(summary(seg_fits[[best_i]])$coefficients, digits = 5)
cat("slope before", round(seg_profile$slope1[best_i], 5),
    " slope after", round(seg_profile$slope2[best_i], 5), "\n")
write.csv(seg_profile, file.path(results_dir, "disparity_shape_verify_breakpoint_profile.csv"),
          row.names = FALSE)

# conditioning of the random-effects covariance and of the fixed-effect vcov
cat("\n-- conditioning checks on the linear LMM --\n")
vc <- VarCorr(lmm_lin_reml)$participantId
cat("random-effect SDs:", round(sqrt(diag(vc)), 5),
    " correlation:", round(attr(vc, "correlation")[1, 2], 5), "\n")
cat("RE covariance eigenvalues:", signif(eigen(vc, only.values = TRUE)$values, 5), "\n")
cat("RE covariance condition number:", signif(kappa(vc, exact = TRUE), 5), "\n")
cat("theta (Cholesky):", signif(getME(lmm_lin_reml, "theta"), 5), "\n")
cat("isSingular:", isSingular(lmm_lin_reml), "\n")
V <- as.matrix(vcov(lmm_lin_reml))
cat("fixed-effect vcov eigenvalues:", signif(eigen(V, only.values = TRUE)$values, 4), "\n")
cat("fixed-effect vcov condition number:", signif(kappa(V, exact = TRUE), 5), "\n")

cat("\n-- optimizer agreement (linear LMM, REML) --\n")
opt_cmp <- lapply(c("bobyqa", "Nelder_Mead", "nloptwrap"), function(o) {
  m <- update(lmm_lin_reml, control = lmerControl(optimizer = o))
  data.frame(optimizer = o, logLik = as.numeric(logLik(m)),
             slope = unname(fixef(m)["stimulusDisparity_m"]),
             se_slope = unname(sqrt(diag(vcov(m)))[2]),
             intercept = unname(fixef(m)["(Intercept)"]),
             se_intercept = unname(sqrt(diag(vcov(m)))[1]),
             sd_int = sqrt(VarCorr(m)$participantId[1, 1]),
             sd_slope = sqrt(VarCorr(m)$participantId[2, 2]),
             singular = isSingular(m))
}) %>% bind_rows()
print(opt_cmp, digits = 6, row.names = FALSE)

cat("\n-- optimizer agreement (segmented LMM at the best breakpoint, ML) --\n")
dfb <- df; dfb$hinge <- pmax(dfb$stimulusDisparity_m - seg_profile$breakpoint_m[best_i], 0)
opt_seg <- lapply(c("bobyqa", "Nelder_Mead", "nloptwrap"), function(o) {
  m <- lmer(signedError_m ~ stimulusDisparity_m + hinge + soundType +
              (1 + stimulusDisparity_m | participantId), data = dfb, REML = FALSE,
            control = lmerControl(optimizer = o))
  data.frame(optimizer = o, logLik = as.numeric(logLik(m)),
             slope1 = unname(fixef(m)["stimulusDisparity_m"]),
             hinge = unname(fixef(m)["hinge"]),
             se_hinge = unname(sqrt(diag(vcov(m)))["hinge"]),
             singular = isSingular(m))
}) %>% bind_rows()
print(opt_seg, digits = 6, row.names = FALSE)
cat("\n")

# ---------------------------------------------------------------------------
# E. Proportionality test on the bias ratio
# ---------------------------------------------------------------------------

cat("== E. Proportionality ==\n")
bias_lmm <- lmer(ventriloquistBias ~ stimulusDisparity_m + (1 | participantId), data = df)
print(summary(bias_lmm)$coefficients, digits = 5)
ci_bias <- confint(bias_lmm, parm = c("(Intercept)", "stimulusDisparity_m"), method = "profile")
print(ci_bias, digits = 5)
cat("singular:", isSingular(bias_lmm),
    " participant SD:", round(sqrt(VarCorr(bias_lmm)$participantId[1, 1]), 5), "\n")
cat("slope across 0.55 m:",
    round(100 * fixef(bias_lmm)["stimulusDisparity_m"] * 0.55, 3), "pp;",
    "CI across 0.55 m:", round(100 * ci_bias["stimulusDisparity_m", ] * 0.55, 3), "pp\n")

cat("\nsame model + soundType (the signed-error models include it; this one does not):\n")
bias_lmm_st <- lmer(ventriloquistBias ~ stimulusDisparity_m + soundType + (1 | participantId), data = df)
print(summary(bias_lmm_st)$coefficients["stimulusDisparity_m", , drop = FALSE], digits = 5)
cat("\nsame model + by-participant random slope:\n")
bias_lmm_rs <- lmer(ventriloquistBias ~ stimulusDisparity_m +
                      (1 + stimulusDisparity_m | participantId), data = df)
print(summary(bias_lmm_rs)$coefficients["stimulusDisparity_m", , drop = FALSE], digits = 5)
cat("singular:", isSingular(bias_lmm_rs), "\n")

lo <- quantile(df$ventriloquistBias, 0.05); hi <- quantile(df$ventriloquistBias, 0.95)
df_trim <- df %>% filter(ventriloquistBias >= lo, ventriloquistBias <= hi)
cat("\ntrimmed at", sprintf("%.4f / %.4f", lo, hi), " n =", nrow(df_trim), "\n")
bias_trim <- lmer(ventriloquistBias ~ stimulusDisparity_m + (1 | participantId), data = df_trim)
print(summary(bias_trim)$coefficients, digits = 5)
ci_trim <- confint(bias_trim, parm = "stimulusDisparity_m", method = "profile")
cat("trimmed slope 95% profile CI:", sprintf("%.5f to %.5f", ci_trim[1, 1], ci_trim[1, 2]), "\n")

set.seed(20260805)
rq_slope <- function(tau, data) {
  fit <- rq(ventriloquistBias ~ stimulusDisparity_m + participantId, tau = tau, data = data)
  s <- summary(fit, se = "boot", R = 2000, bsmethod = "xy")
  cf <- s$coefficients["stimulusDisparity_m", ]
  data.frame(tau = tau, slope = cf[["Value"]], se = cf[["Std. Error"]],
             t = cf[["t value"]], p = cf[["Pr(>|t|)"]],
             ci_low = cf[["Value"]] - 1.96 * cf[["Std. Error"]],
             ci_high = cf[["Value"]] + 1.96 * cf[["Std. Error"]])
}
rq_tbl <- bind_rows(lapply(c(0.25, 0.5, 0.75), rq_slope, data = df))
cat("\nquantile regression with participant fixed effects (xy bootstrap SE, R = 2000):\n")
print(rq_tbl, digits = 5, row.names = FALSE)

# bin medians (published 3 bins and quintiles)
df$disp_quintile <- cut(df$stimulusDisparity_m,
                        breaks = quantile(df$stimulusDisparity_m, seq(0, 1, 0.2)),
                        include.lowest = TRUE)
cat("\nbias medians by quintile:\n")
print(as.data.frame(df %>% group_by(disp_quintile) %>%
  summarise(n = n(), med = 100 * median(ventriloquistBias),
            mean = 100 * mean(ventriloquistBias), sd = 100 * sd(ventriloquistBias),
            .groups = "drop")), digits = 4)
cat("\n")

# ---------------------------------------------------------------------------
# F. Is the Gaussian homoscedastic null adequate? (validity of the bootstrap)
# ---------------------------------------------------------------------------

cat("== F. Residual structure of the linear LMM (bootstrap null adequacy) ==\n")
res <- residuals(lmm_lin_reml)
cat("residual SD by disparity tertile:\n")
print(tapply(res, cut(df$stimulusDisparity_m, quantile(df$stimulusDisparity_m, c(0, 1/3, 2/3, 1)),
                      include.lowest = TRUE), sd), digits = 4)
cat("Levene-type test (abs residual ~ disparity, OLS):\n")
print(summary(lm(abs(res) ~ df$stimulusDisparity_m))$coefficients, digits = 5)
cat("residual skewness:", round(mean(scale(res)^3), 4),
    " excess kurtosis:", round(mean(scale(res)^4) - 3, 4), "\n")
cat("Shapiro-Wilk on residuals: W =", round(shapiro.test(res)$statistic, 5),
    " p =", signif(shapiro.test(res)$p.value, 4), "\n\n")

# ---------------------------------------------------------------------------
# G. Angular disparity
# ---------------------------------------------------------------------------

cat("== G. Angular disparity ==\n")
ang_vars <- c("angDisp_onset_deg", "angDisp_mean_deg", "angDisp_min_deg", "angDisp_max_deg")
ang_tbl <- lapply(ang_vars, function(v) {
  x <- df[[v]]
  data.frame(variable = v, mean = mean(x), sd = sd(x), min = min(x),
             q05 = quantile(x, 0.05), q25 = quantile(x, 0.25), median = median(x),
             q75 = quantile(x, 0.75), q95 = quantile(x, 0.95), max = max(x),
             pct_over_30 = 100 * mean(x > 30), pct_over_15 = 100 * mean(x > 15))
}) %>% bind_rows()
rownames(ang_tbl) <- NULL
print(ang_tbl, digits = 5)
cat("\nby published bin (onset):\n")
print(as.data.frame(df %>% group_by(disparityRange) %>%
  summarise(n = n(), onset_med = median(angDisp_onset_deg),
            onset_pct30 = 100 * mean(angDisp_onset_deg > 30),
            mean_med = median(angDisp_mean_deg),
            mean_pct30 = 100 * mean(angDisp_mean_deg > 30), .groups = "drop")), digits = 4)
cat("\ncorrelations: onset r =", round(cor(df$stimulusDisparity_m, df$angDisp_onset_deg), 4),
    " Spearman", round(cor(df$stimulusDisparity_m, df$angDisp_onset_deg, method = "spearman"), 4),
    "| trial mean r =", round(cor(df$stimulusDisparity_m, df$angDisp_mean_deg), 4),
    " Spearman", round(cor(df$stimulusDisparity_m, df$angDisp_mean_deg, method = "spearman"), 4), "\n")

cat("\nsecondary GAMM on onset angular disparity:\n")
gamm_ang <- gamm(signedError_m ~ s(angDisp_onset_deg, k = 6) + soundType,
                 random = list(participantId = ~ 1), data = df)
print(summary(gamm_ang$gam)$s.table, digits = 5)

# what does that curvature look like? fitted smooth on a grid, and the bias ratio
ang_grid <- data.frame(angDisp_onset_deg = seq(5, 100, by = 5),
                       soundType = factor(levels(df$soundType)[1], levels = levels(df$soundType)))
pa <- predict(gamm_ang$gam, newdata = ang_grid, se.fit = TRUE)
ang_fit <- data.frame(angDisp_onset_deg = ang_grid$angDisp_onset_deg,
                      fit_cm = 100 * as.numeric(pa$fit), se_cm = 100 * as.numeric(pa$se.fit))
print(ang_fit, digits = 4, row.names = FALSE)
write.csv(ang_fit, file.path(results_dir, "disparity_shape_verify_angular_smooth.csv"),
          row.names = FALSE)

cat("\nbias ratio vs angular disparity (does capture fall at large angles?):\n")
bias_ang <- lmer(ventriloquistBias ~ angDisp_onset_deg + (1 | participantId), data = df)
print(summary(bias_ang)$coefficients, digits = 5)
bias_ang_m <- lmer(ventriloquistBias ~ angDisp_mean_deg + (1 | participantId), data = df)
print(summary(bias_ang_m)$coefficients, digits = 5)
cat("\nbias ratio median by onset-angle bin:\n")
print(as.data.frame(df %>%
  mutate(abin = cut(angDisp_onset_deg, c(0, 15, 30, 45, 60, 200))) %>%
  group_by(abin) %>%
  summarise(n = n(), disp_mean = mean(stimulusDisparity_m),
            bias_med = 100 * median(ventriloquistBias),
            bias_mean = 100 * mean(ventriloquistBias),
            signed_mean_cm = 100 * mean(signedError_m), .groups = "drop")), digits = 4)

# joint model: is there a residual angular effect once metric disparity is in?
cat("\nsigned error ~ disparity + s(onset angle), both in the model:\n")
gam_joint <- gam(signedError_m ~ stimulusDisparity_m + s(angDisp_onset_deg, k = 6) + soundType +
                   s(participantId, bs = "re") + s(participantId, stimulusDisparity_m, bs = "re"),
                 data = df, method = "REML")
print(summary(gam_joint)$s.table, digits = 5)
print(summary(gam_joint)$p.table, digits = 5)

write.csv(rbind(
  data.frame(section = "gamm ML AIC", gamm_aic),
  data.frame(section = "parametric shapes", model = shape_tbl$model, df = shape_tbl$df,
             logLik = shape_tbl$logLik, AIC = shape_tbl$AIC, BIC = shape_tbl$BIC)),
  file.path(results_dir, "disparity_shape_verify_model_comparison.csv"), row.names = FALSE)
write.csv(rq_tbl, file.path(results_dir, "disparity_shape_verify_rq.csv"), row.names = FALSE)
write.csv(opt_cmp, file.path(results_dir, "disparity_shape_verify_optimizers.csv"), row.names = FALSE)

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
