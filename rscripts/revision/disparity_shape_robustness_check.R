# Robustness of the disparity-shape null results to assumption violations
#
# The published shape analysis assumes a homoscedastic Gaussian linear mixed
# model for the signed error, and reports an exactly linear penalised smooth
# (edf = 1.00). This script tests how much of that depends on assumptions the
# data do not satisfy:
#   1. residual normality and constant variance of the signed-error LMM
#   2. effect of modelling the variance (varPower in disparity)
#   3. influence of the heaviest residuals on the smooth, the quadratic term
#      and the breakpoint search
#   4. leave-one-participant-out stability of the smooth
#   5. heteroscedasticity of the bias ratio and its effect on the slope
#
# Input : results_revision/analysis_df_revision.rds
# Output: results_revision/disparity_shape_robustness_*.csv, *_robustness_log.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
  library(nlme)
  library(lme4)
  library(lmerTest)
})

set.seed(20260805)

proj_dir    <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
results_dir <- file.path(proj_dir, "results_revision")
log_path    <- file.path(results_dir, "disparity_shape_robustness_log.txt")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(); close(log_con) }, add = TRUE)

cat("Robustness of the disparity-shape null results\n")
cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

df <- readRDS(file.path(results_dir, "analysis_df_revision.rds")) %>%
  mutate(participantId = factor(participantId), soundType = factor(soundType)) %>%
  as.data.frame()
df$disp2 <- df$stimulusDisparity_m^2

re_smooths <- ~ . + s(participantId, bs = "re") + s(participantId, stimulusDisparity_m, bs = "re")
gam_shape <- function(data, k = 6) {
  gam(update(signedError_m ~ s(stimulusDisparity_m, k = k) + soundType, re_smooths),
      data = data, method = "REML")
}
lmm_lin_f  <- signedError_m ~ stimulusDisparity_m + soundType + (1 + stimulusDisparity_m | participantId)
lmm_quad_f <- signedError_m ~ stimulusDisparity_m + disp2 + soundType + (1 + stimulusDisparity_m | participantId)

breaks <- seq(0.20, 0.60, by = 0.01)
grid_search <- function(data) {
  m0 <- lmer(lmm_lin_f, data = data, REML = FALSE)
  lls <- sapply(breaks, function(cv) {
    data$hinge <- pmax(data$stimulusDisparity_m - cv, 0)
    as.numeric(logLik(lmer(signedError_m ~ stimulusDisparity_m + hinge + soundType +
                             (1 + stimulusDisparity_m | participantId),
                           data = data, REML = FALSE)))
  })
  c(best_c = breaks[which.max(lls)], max_LRT = 2 * (max(lls) - as.numeric(logLik(m0))))
}

# ---------------------------------------------------------------------------
# 1-2. Residual diagnostics and a variance function
# ---------------------------------------------------------------------------

cat("== 1. Residual diagnostics of the linear LMM ==\n")
m_lin <- lmer(lmm_lin_f, data = df)
r <- residuals(m_lin)
cat("Shapiro-Wilk W =", round(shapiro.test(r)$statistic, 5),
    " p =", signif(shapiro.test(r)$p.value, 4), "\n")
cat("skewness", round(mean(scale(r)^3), 4), " excess kurtosis", round(mean(scale(r)^4) - 3, 4), "\n")
cat("residual SD by disparity tertile:\n")
print(round(tapply(r, cut(df$stimulusDisparity_m, quantile(df$stimulusDisparity_m, c(0, 1/3, 2/3, 1)),
                          include.lowest = TRUE), sd), 5))
cat("abs(residual) ~ disparity slope:\n")
print(summary(lm(abs(r) ~ df$stimulusDisparity_m))$coefficients, digits = 5)

cat("\n== 2. Modelling the variance (varPower in disparity) ==\n")
rs <- list(participantId = ~ 1 + stimulusDisparity_m)
g_homo <- gamm(signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType, random = rs, data = df)
g_hetero <- gamm(signedError_m ~ s(stimulusDisparity_m, k = 6) + soundType, random = rs, data = df,
                 weights = varPower(form = ~ stimulusDisparity_m))
cat("edf homoscedastic:", round(summary(g_homo$gam)$s.table[1, "edf"], 4),
    " edf varPower:", round(summary(g_hetero$gam)$s.table[1, "edf"], 4), "\n")
cat("varPower delta:", round(as.numeric(g_hetero$lme$modelStruct$varStruct[1]), 4),
    " LRT for the variance function (1 df):",
    round(2 * (as.numeric(logLik(g_hetero$lme)) - as.numeric(logLik(g_homo$lme))), 4), "\n")

# ---------------------------------------------------------------------------
# 3. Influence of the heaviest residuals
# ---------------------------------------------------------------------------

cat("\n== 3. Influence of extreme residuals on the shape conclusions ==\n")
z <- residuals(lmer(lmm_lin_f, data = df, REML = FALSE)) / sigma(lmer(lmm_lin_f, data = df, REML = FALSE))
cat("trials with |standardised residual| > 3:", sum(abs(z) > 3),
    " > 4:", sum(abs(z) > 4), "\n")
cat("participants contributing the |z| > 3 trials:",
    paste(sort(unique(as.character(df$participantId[abs(z) > 3]))), collapse = " "), "\n\n")

trim_tbl <- lapply(c(Inf, 4, 3), function(cut_z) {
  d <- df[abs(z) <= cut_z, ]
  gs <- gam_shape(d)
  mq <- lmer(lmm_quad_f, data = d, REML = FALSE)
  ml <- lmer(lmm_lin_f, data = d, REML = FALSE)
  gr <- grid_search(d)
  cq <- summary(mq)$coefficients["disp2", ]
  data.frame(cut_abs_z = cut_z, n = nrow(d),
             edf = summary(gs)$s.table[1, "edf"],
             slope = unname(fixef(ml)["stimulusDisparity_m"]),
             quad_b = cq[["Estimate"]], quad_se = cq[["Std. Error"]], quad_p = cq[["Pr(>|t|)"]],
             best_c = unname(gr["best_c"]), max_LRT = unname(gr["max_LRT"]))
}) %>% bind_rows()
print(trim_tbl, digits = 5, row.names = FALSE)
write.csv(trim_tbl, file.path(results_dir, "disparity_shape_robustness_trimming.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Leave-one-participant-out stability of the smooth
# ---------------------------------------------------------------------------

cat("\n== 4. Leave-one-participant-out edf of the smooth ==\n")
loo <- lapply(levels(df$participantId), function(p) {
  d <- droplevels(df[df$participantId != p, ])
  gs <- gam_shape(d)
  mq <- lmer(lmm_quad_f, data = d, REML = FALSE)
  data.frame(dropped = p, n = nrow(d), edf = summary(gs)$s.table[1, "edf"],
             quad_b = unname(fixef(mq)["disp2"]),
             quad_p = summary(mq)$coefficients["disp2", "Pr(>|t|)"])
}) %>% bind_rows()
cat("edf: min", round(min(loo$edf), 4), " median", round(median(loo$edf), 4),
    " max", round(max(loo$edf), 4), "\n")
cat("participants whose removal pushes edf above 1.5:",
    paste(loo$dropped[loo$edf > 1.5], collapse = " "), "\n")
cat("quadratic p: min", signif(min(loo$quad_p), 4), " max", signif(max(loo$quad_p), 4),
    "; number of leave-one-out fits with p < .05:", sum(loo$quad_p < 0.05), "\n")
print(loo[order(-loo$edf), ][1:5, ], digits = 5, row.names = FALSE)
write.csv(loo, file.path(results_dir, "disparity_shape_robustness_loo.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Bias ratio: heteroscedasticity and the slope estimate
# ---------------------------------------------------------------------------

cat("\n== 5. Bias ratio, homoscedastic vs variance-modelled ==\n")
b0 <- lme(ventriloquistBias ~ stimulusDisparity_m, random = ~ 1 | participantId,
          data = df, method = "REML")
b1 <- lme(ventriloquistBias ~ stimulusDisparity_m, random = ~ 1 | participantId,
          data = df, method = "REML", weights = varPower(form = ~ stimulusDisparity_m),
          control = lmeControl(opt = "optim"))
cat("homoscedastic (as published):\n"); print(summary(b0)$tTable, digits = 5)
cat("varPower in disparity:\n");        print(summary(b1)$tTable, digits = 5)
cat("varPower delta:", round(as.numeric(b1$modelStruct$varStruct[1]), 4), "\n")
print(anova(b0, b1))
ci1 <- intervals(b1, which = "fixed")$fixed
cat("variance-modelled slope 95% CI:",
    sprintf("%.4f to %.4f", ci1["stimulusDisparity_m", 1], ci1["stimulusDisparity_m", 3]), "\n")

bias_tbl <- data.frame(
  model = c("LMM, homoscedastic (published)", "LME, varPower(disparity)"),
  slope = c(fixef(b0)[["stimulusDisparity_m"]], fixef(b1)[["stimulusDisparity_m"]]),
  se = c(summary(b0)$tTable["stimulusDisparity_m", "Std.Error"],
         summary(b1)$tTable["stimulusDisparity_m", "Std.Error"]),
  p = c(summary(b0)$tTable["stimulusDisparity_m", "p-value"],
        summary(b1)$tTable["stimulusDisparity_m", "p-value"]),
  AIC = c(AIC(b0), AIC(b1)))
print(bias_tbl, digits = 5, row.names = FALSE)
write.csv(bias_tbl, file.path(results_dir, "disparity_shape_robustness_bias.csv"),
          row.names = FALSE)

cat("\nDone:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
