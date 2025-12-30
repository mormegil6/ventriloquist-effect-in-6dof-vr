#!/usr/bin/env Rscript
# Statistical Analysis - 6DoF VR Ventriloquist Effect
# Reads raw data and prints statistical analysis outputs

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(MuMIn)
})

options(width = 120)
options(warn = -1)

# ========== Data Import ==========
cat("# 6DoF VR Ventriloquist Effect - Statistical Analysis Results\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

resolve_data_folder <- function() {
  candidates <- c("../data", "data", "../../data")
  for (p in candidates) {
    if (dir.exists(p)) return(p)
  }
  ups <- c(".", "..", "../..", "../../..")
  for (u in ups) {
    p <- file.path(u, "six-dof-correlations-main", "data")
    if (dir.exists(p)) return(p)
  }
  stop("Data folder not found. Run from 'six-dof-correlations-main/rdocs' or set working dir to the project root.")
}

data_folder <- resolve_data_folder()
json_files <- sort(list.files(data_folder, pattern = "\\.json$", full.names = TRUE))
csv_files <- sort(list.files(data_folder, pattern = "\\.csv$", full.names = TRUE))

if (length(json_files) == 0) {
  stop(sprintf("No JSON files found in data folder: %s", normalizePath(data_folder)))
}
if (length(csv_files) == 0) {
  stop(sprintf("No CSV files found in data folder: %s", normalizePath(data_folder)))
}

participant_ids <- paste0("P", sprintf("%03d", seq_along(json_files)))

# Load experiment data
read_experiment <- function(path, pid) {
  df <- as.data.frame(fromJSON(path, flatten = TRUE)$dataEntries)
  df$participantId <- pid
  df
}

experiment_df <- purrr::map2_dfr(json_files, participant_ids, read_experiment)

# Load trajectory data
read_trajectory <- function(path, pid) {
  df <- read_csv(path, show_col_types = FALSE)
  df$participantId <- pid
  df
}

trajectory_df <- purrr::map2_dfr(csv_files, participant_ids, read_trajectory)

cat("## 1. Dataset Overview\n\n")
cat("| Measure | Value |\n")
cat("|---------|-------|\n")
cat(sprintf("| Participants | %d |\n", length(json_files)))
cat(sprintf("| Total trials | %d |\n", nrow(experiment_df)))
cat(sprintf("| Trials per participant | %d |\n", nrow(experiment_df) / length(json_files)))
cat(sprintf("| Trajectory frames | %d |\n\n", nrow(trajectory_df)))

# ========== Data Processing ==========
analysis_df <- experiment_df %>%
  mutate(
    participantError_m = distanceError / 100,
    stimulusDisparity_m = flashDistance,
    soundType = factor(sampleId, levels = 1:4, 
                       labels = c("Drum", "Flute", "Speech", "Pink Noise")),
    disparityRange = cut(flashDistance * 100, breaks = c(0, 30, 40, 70),
                         labels = c("Short (15-30cm)", "Medium (30-40cm)", "Long (40-70cm)"), 
                         include.lowest = TRUE),
    azimuth_adj = ifelse(flashAzimuth > 315, flashAzimuth - 360, flashAzimuth),
    azimuthSector = cut(azimuth_adj, breaks = c(-45, 45, 135, 225, 315),
                        labels = c("Front", "Right", "Back", "Left"), include.lowest = TRUE),
    elevationCategory = cut(flashElevation, breaks = c(-90, -22.5, 22.5, 90),
                            labels = c("Below", "Level", "Above"), include.lowest = TRUE)
  ) %>%
  transmute(
    participantId, trialSequenceNum = level, soundType, disparityRange,
    azimuthSector, elevationCategory,
    stimulusDisparity_m, participantError_m,
    response_x = spherePosition.x, response_y = spherePosition.y, response_z = spherePosition.z,
    sound_x = referencePosition.x, sound_y = referencePosition.y, sound_z = referencePosition.z,
    flash_x = flashSpherePosition.x, flash_y = flashSpherePosition.y, flash_z = flashSpherePosition.z
  )

# Compute signed error (ventriloquist bias)
analysis_df <- analysis_df %>%
  mutate(
    sf_x = flash_x - sound_x, sf_y = flash_y - sound_y, sf_z = flash_z - sound_z,
    sf_len = sqrt(sf_x^2 + sf_y^2 + sf_z^2),
    sf_ux = sf_x/sf_len, sf_uy = sf_y/sf_len, sf_uz = sf_z/sf_len,
    sr_x = response_x - sound_x, sr_y = response_y - sound_y, sr_z = response_z - sound_z,
    signedError_m = sr_x*sf_ux + sr_y*sf_uy + sr_z*sf_uz,
    ventriloquistBias = signedError_m / stimulusDisparity_m
  ) %>%
  select(-sf_x, -sf_y, -sf_z, -sf_len, -sf_ux, -sf_uy, -sf_uz, -sr_x, -sr_y, -sr_z,
         -response_x, -response_y, -response_z, -sound_x, -sound_y, -sound_z,
         -flash_x, -flash_y, -flash_z)

# Compute trajectory metrics
trajectory_summary <- trajectory_df %>%
  rename(trialSequenceNum = level) %>%
  arrange(participantId, trialSequenceNum, timestamp) %>%
  group_by(participantId, trialSequenceNum) %>%
  mutate(
    dx = px - lag(px, default = first(px)),
    dy = py - lag(py, default = first(py)),
    dz = pz - lag(pz, default = first(pz)),
    segment = sqrt(dx^2 + dy^2 + dz^2)
  ) %>%
  summarise(
    total_path_length = sum(segment, na.rm = TRUE),
    trajectory_duration = max(timestamp) - min(timestamp),
    movement_rate = total_path_length / pmax(trajectory_duration, 0.1),
    .groups = "drop"
  )

analysis_df <- left_join(analysis_df, trajectory_summary, by = c("participantId", "trialSequenceNum"))

# Centered predictors for stable GLMMs
analysis_df <- analysis_df %>%
  mutate(
    participantError_cm = participantError_m * 100,
    stimDisparity_c = stimulusDisparity_m - mean(stimulusDisparity_m),
    trialSequence_c = trialSequenceNum - mean(trialSequenceNum)
  )

# ========== Section 2: Raw Data Sample ==========
cat("## 2. Raw Data Sample\n\n")
cat("First 10 observations:\n\n")
cat("```\n")
print(head(analysis_df %>% 
             select(participantId, trialSequenceNum, soundType, disparityRange,
                    stimulusDisparity_m, participantError_m, signedError_m, 
                    ventriloquistBias, total_path_length), 10), 
      row.names = FALSE)
cat("```\n\n")

# ========== Section 3: Summary Statistics ==========
cat("## 3. Descriptive Statistics\n\n")

cat("### 3.1 Overall Performance\n\n")
overall <- analysis_df %>%
  summarise(
    n_trials = n(),
    n_participants = n_distinct(participantId),
    mean_error_cm = mean(participantError_m * 100),
    sd_error_cm = sd(participantError_m * 100),
    median_error_cm = median(participantError_m * 100),
    min_error_cm = min(participantError_m * 100),
    max_error_cm = max(participantError_m * 100),
    mean_signed_cm = mean(signedError_m * 100),
    sd_signed_cm = sd(signedError_m * 100),
    mean_bias = mean(ventriloquistBias),
    sd_bias = sd(ventriloquistBias),
    mean_path_m = mean(total_path_length, na.rm = TRUE),
    sd_path_m = sd(total_path_length, na.rm = TRUE)
  )

cat("| Measure | Value |\n")
cat("|---------|-------|\n")
cat(sprintf("| N (trials) | %d |\n", overall$n_trials))
cat(sprintf("| N (participants) | %d |\n", overall$n_participants))
cat(sprintf("| Mean unsigned error | %.1f cm (SD = %.1f) |\n", overall$mean_error_cm, overall$sd_error_cm))
cat(sprintf("| Median unsigned error | %.1f cm |\n", overall$median_error_cm))
cat(sprintf("| Range unsigned error | %.1f - %.1f cm |\n", overall$min_error_cm, overall$max_error_cm))
cat(sprintf("| Mean signed error | %.1f cm (SD = %.1f) |\n", overall$mean_signed_cm, overall$sd_signed_cm))
cat(sprintf("| Mean ventriloquist bias | %.1f%% (SD = %.1f) |\n", overall$mean_bias*100, overall$sd_bias*100))
cat(sprintf("| Mean path length | %.2f m (SD = %.2f) |\n\n", overall$mean_path_m, overall$sd_path_m))

cat("### 3.2 By Sound Type\n\n")
by_sound <- analysis_df %>%
  group_by(soundType) %>%
  summarise(
    n = n(),
    mean_error_cm = mean(participantError_m * 100),
    sd_error_cm = sd(participantError_m * 100),
    mean_signed_cm = mean(signedError_m * 100),
    sd_signed_cm = sd(signedError_m * 100),
    mean_bias = mean(ventriloquistBias),
    .groups = "drop"
  )
cat("| Sound Type | n | Mean Error (SD) cm | Mean Signed (SD) cm | Bias % |\n")
cat("|------------|---|-------------------|---------------------|--------|\n")
for (i in 1:nrow(by_sound)) {
  cat(sprintf("| %s | %d | %.1f (%.1f) | %.1f (%.1f) | %.1f |\n",
              by_sound$soundType[i], by_sound$n[i], 
              by_sound$mean_error_cm[i], by_sound$sd_error_cm[i],
              by_sound$mean_signed_cm[i], by_sound$sd_signed_cm[i],
              by_sound$mean_bias[i]*100))
}
cat("\n")

cat("### 3.3 By Disparity Range\n\n")
by_disparity <- analysis_df %>%
  group_by(disparityRange) %>%
  summarise(
    n = n(),
    mean_disparity_cm = mean(stimulusDisparity_m * 100),
    mean_error_cm = mean(participantError_m * 100),
    sd_error_cm = sd(participantError_m * 100),
    mean_signed_cm = mean(signedError_m * 100),
    mean_bias = mean(ventriloquistBias),
    .groups = "drop"
  )
cat("| Disparity Range | n | Mean Disparity | Mean Error (SD) cm | Mean Signed cm | Bias % |\n")
cat("|-----------------|---|----------------|-------------------|----------------|--------|\n")
for (i in 1:nrow(by_disparity)) {
  cat(sprintf("| %s | %d | %.1f cm | %.1f (%.1f) | %.1f | %.1f |\n",
              by_disparity$disparityRange[i], by_disparity$n[i],
              by_disparity$mean_disparity_cm[i],
              by_disparity$mean_error_cm[i], by_disparity$sd_error_cm[i],
              by_disparity$mean_signed_cm[i], by_disparity$mean_bias[i]*100))
}
cat("\n")

cat("### 3.4 By Participant\n\n")
cat("```\n")
by_participant <- analysis_df %>%
  group_by(participantId) %>%
  summarise(
    n_trials = n(),
    mean_error_cm = mean(participantError_m * 100),
    sd_error_cm = sd(participantError_m * 100),
    mean_bias = mean(ventriloquistBias),
    mean_path_m = mean(total_path_length, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(by_participant), row.names = FALSE, digits = 3)
cat("```\n\n")

# ========== Section 4: Correlations ==========
cat("## 4. Correlation Matrix\n\n")
cat("Pearson correlations among key variables:\n\n")
cat("```\n")
cor_vars <- analysis_df %>%
  select(stimulusDisparity_m, participantError_m, signedError_m, 
         ventriloquistBias, total_path_length, trialSequenceNum)
cor_matrix <- cor(cor_vars, use = "complete.obs")
print(round(cor_matrix, 3))
cat("```\n\n")

# ========== Section 5: Main Ventriloquist Effect ==========
cat("## 5. Main Ventriloquist Effect (Signed Error Model)\n\n")

cat("### 5.1 Model Specification\n\n")
cat("| Component | Specification |\n")
cat("|-----------|---------------|\n")
cat("| Response | Signed error (meters) - positive = toward flash |\n")
cat("| Fixed effects | Disparity (m), Sound Type |\n")
cat("| Random effects | Participant intercept + disparity slope |\n")
cat("| Family | Gaussian (LMM) |\n")
cat("| Estimation | REML |\n\n")

# Fit models for comparison
model_signed_intercept <- lmer(signedError_m ~ stimulusDisparity_m + soundType + 
                                 (1 | participantId),
                               data = analysis_df, REML = TRUE,
                               control = lmerControl(optimizer = "bobyqa"))

model_signed <- lmer(signedError_m ~ stimulusDisparity_m + soundType + 
                       (1 + stimulusDisparity_m | participantId),
                     data = analysis_df, REML = TRUE,
                     control = lmerControl(optimizer = "bobyqa"))

cat("### 5.2 Model Comparison: Random Intercepts vs Random Slopes\n\n")
lrt_signed <- anova(model_signed_intercept, model_signed, refit = FALSE)
cat("```\n")
print(lrt_signed)
cat("```\n\n")
cat(sprintf("**Result:** Random slopes model preferred (Chi-sq = %.2f, df = %d, p %s)\n\n",
            lrt_signed$Chisq[2], lrt_signed$Df[2],
            ifelse(lrt_signed$`Pr(>Chisq)`[2] < 0.001, "< .001",
                   sprintf("= %.3f", lrt_signed$`Pr(>Chisq)`[2]))))

cat("### 5.3 Final Model Summary\n\n")
cat("```\n")
print(summary(model_signed))
cat("```\n\n")

# R-squared
r2_signed <- r.squaredGLMM(model_signed)
cat("### 5.4 Variance Explained (R-squared)\n\n")
cat("| Type | R² |\n")
cat("|------|----|\n")
cat(sprintf("| Marginal R² (fixed effects) | %.3f |\n", r2_signed[1, "R2m"]))
cat(sprintf("| Conditional R² (fixed + random) | %.3f |\n\n", r2_signed[1, "R2c"]))

vc <- VarCorr(model_signed)
var_intercept <- as.numeric(attr(vc$participantId, "stddev")["(Intercept)"])^2
var_residual <- sigma(model_signed)^2
icc_signed <- var_intercept / (var_intercept + var_residual)
cat("### 5.5 Intraclass Correlation (ICC)\n\n")
cat(sprintf("ICC (participant) = %.3f\n\n", icc_signed))
cat(sprintf("Interpretation: %.1f%% of variance in signed error is attributable to between-participant differences.\n\n", icc_signed * 100))

cat("### 5.6 Fixed Effects with 95% CI\n\n")
fe <- fixef(model_signed)
ci_signed <- confint(model_signed, parm = "beta_", method = "Wald")

cat("| Parameter | Estimate | 95% CI | Interpretation |\n")
cat("|-----------|----------|--------|----------------|\n")
cat(sprintf("| Intercept | %.4f m | [%.4f, %.4f] | Baseline at 0 disparity, Drum |\n",
            fe["(Intercept)"], ci_signed["(Intercept)", 1], ci_signed["(Intercept)", 2]))
cat(sprintf("| Disparity | %.4f m/m | [%.4f, %.4f] | %.1f%% visual capture |\n",
            fe["stimulusDisparity_m"], ci_signed["stimulusDisparity_m", 1], 
            ci_signed["stimulusDisparity_m", 2], fe["stimulusDisparity_m"]*100))
cat(sprintf("| Flute vs Drum | %.4f m | [%.4f, %.4f] | %+.2f cm |\n",
            fe["soundTypeFlute"], ci_signed["soundTypeFlute", 1], 
            ci_signed["soundTypeFlute", 2], fe["soundTypeFlute"]*100))
cat(sprintf("| Speech vs Drum | %.4f m | [%.4f, %.4f] | %+.2f cm |\n",
            fe["soundTypeSpeech"], ci_signed["soundTypeSpeech", 1], 
            ci_signed["soundTypeSpeech", 2], fe["soundTypeSpeech"]*100))
cat(sprintf("| Pink Noise vs Drum | %.4f m | [%.4f, %.4f] | %+.2f cm |\n\n",
            fe["soundTypePink Noise"], ci_signed["soundTypePink Noise", 1], 
            ci_signed["soundTypePink Noise", 2], fe["soundTypePink Noise"]*100))

cat("### 5.7 Random Effects Variance\n\n")
cat("```\n")
print(VarCorr(model_signed))
cat("```\n\n")

# ========== Section 6: Localization Accuracy ==========
cat("## 6. Localization Accuracy (Gamma GLMM)\n\n")

cat("### 6.1 Model Specification\n\n")
cat("| Component | Specification |\n")
cat("|-----------|---------------|\n")
cat("| Response | Unsigned error (cm) - always positive |\n")
cat("| Fixed effects | Disparity (centered), Sound Type, Trial (centered), Azimuth, Elevation |\n")
cat("| Random effects | Participant intercept |\n")
cat("| Family | Gamma (log link) |\n")
cat("| Rationale | Right-skewed positive data; multiplicative effects |\n\n")

# Fit models for comparison
model_accuracy_reduced <- glmer(participantError_cm ~ stimDisparity_c + soundType + trialSequence_c +
                                  (1 | participantId),
                                data = analysis_df, family = Gamma(link = "log"),
                                control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

model_accuracy <- glmer(participantError_cm ~ stimDisparity_c + soundType + trialSequence_c +
                          azimuthSector + elevationCategory +
                          (1 | participantId),
                        data = analysis_df, family = Gamma(link = "log"),
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

cat("### 6.2 Model Comparison: With vs Without Spatial Factors\n\n")
lrt_spatial <- anova(model_accuracy_reduced, model_accuracy)
cat("```\n")
print(lrt_spatial)
cat("```\n\n")
cat(sprintf("**Result:** Spatial factors %s model fit (Chi-sq = %.2f, df = %d, p %s)\n\n",
            ifelse(lrt_spatial$`Pr(>Chisq)`[2] < 0.05, "significantly improve", "do not improve"),
            lrt_spatial$Chisq[2], lrt_spatial$Df[2],
            ifelse(lrt_spatial$`Pr(>Chisq)`[2] < 0.001, "< .001",
                   sprintf("= %.3f", lrt_spatial$`Pr(>Chisq)`[2]))))

cat("### 6.3 Final Model Summary\n\n")
cat("```\n")
print(summary(model_accuracy))
cat("```\n\n")

# R-squared for GLMM
r2_acc <- tryCatch(r.squaredGLMM(model_accuracy), error = function(e) NULL)
if (!is.null(r2_acc)) {
  cat("### 6.4 Variance Explained (R-squared)\n\n")
  cat("| Type | R² |\n")
  cat("|------|----|\n")
  cat(sprintf("| Marginal R² (fixed effects) | %.3f |\n", r2_acc[1, "R2m"]))
  cat(sprintf("| Conditional R² (fixed + random) | %.3f |\n\n", r2_acc[1, "R2c"]))
}

# Fitted vs observed correlation
analysis_df$fitted_error_cm <- fitted(model_accuracy)
fit_cor <- cor(analysis_df$fitted_error_cm, analysis_df$participantError_cm)
cat("### 6.5 Model Fit: Fitted vs Observed Correlation\n\n")
cat(sprintf("r(fitted, observed) = %.3f\n\n", fit_cor))

cat("### 6.6 Effect Sizes (Multiplicative, Log Link)\n\n")
fe_acc <- fixef(model_accuracy)
cat("| Effect | Per-unit change | Interpretation |\n")
cat("|--------|-----------------|----------------|\n")
cat(sprintf("| Disparity | +10cm → %.1f%% | Larger disparity increases error |\n", 
            (exp(fe_acc["stimDisparity_c"]*0.1)-1)*100))
cat(sprintf("| Trial sequence | per trial → %.2f%% | Learning effect |\n", 
            (exp(fe_acc["trialSequence_c"])-1)*100))
cat(sprintf("| Flute vs Drum | %.1f%% | Sound type effect |\n", (exp(fe_acc["soundTypeFlute"])-1)*100))
cat(sprintf("| Speech vs Drum | %.1f%% | Sound type effect |\n", (exp(fe_acc["soundTypeSpeech"])-1)*100))
cat(sprintf("| Pink Noise vs Drum | %.1f%% | Sound type effect |\n\n", (exp(fe_acc["soundTypePink Noise"])-1)*100))

# ========== Section 7: Sound Type Effects ==========
cat("## 7. Sound Type Effects (Post-Hoc Comparisons)\n\n")

cat("### 7.1 Estimated Marginal Means (Signed Error)\n\n")
cat("```\n")
emm <- emmeans(model_signed, ~ soundType)
print(summary(emm))
cat("```\n\n")

cat("### 7.2 Pairwise Contrasts (Signed Error, Tukey-adjusted)\n\n")
cat("```\n")
print(pairs(emm))
cat("```\n\n")

cat("### 7.3 Estimated Marginal Means (Accuracy, response scale)\n\n")
cat("```\n")
emm_acc <- emmeans(model_accuracy, ~ soundType, type = "response")
print(summary(emm_acc))
cat("```\n\n")

cat("### 7.4 Pairwise Contrasts (Accuracy, Tukey-adjusted)\n\n")
cat("```\n")
print(pairs(emm_acc))
cat("```\n\n")

# ========== Section 8: Movement Effects ==========
cat("## 8. Movement Effects (6DoF Exploration)\n\n")

cat("### 8.1 Movement Decomposition\n\n")
cat("To disentangle individual strategies from trial-level effects:\n")
cat("- **Between-subject (mean_rate)**: participant's average movement rate across all trials\n")
cat("- **Within-subject (rate_dev)**: trial-level deviation from participant's mean\n\n")

# Between/within decomposition
analysis_df <- analysis_df %>%
  group_by(participantId) %>%
  mutate(
    mean_rate = mean(movement_rate, na.rm = TRUE),
    rate_dev = movement_rate - mean_rate
  ) %>%
  ungroup()

model_movement <- glmer(participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum +
                          mean_rate + rate_dev + (1 | participantId),
                        data = analysis_df, family = Gamma(link = "log"),
                        control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

cat("### 8.2 Model Summary\n\n")
cat("```\n")
print(summary(model_movement))
cat("```\n\n")

# R-squared
r2_mov <- tryCatch(r.squaredGLMM(model_movement), error = function(e) NULL)
if (!is.null(r2_mov)) {
  cat("### 8.3 Variance Explained\n\n")
  cat(sprintf("Marginal R² = %.3f | Conditional R² = %.3f\n\n", r2_mov[1, "R2m"], r2_mov[1, "R2c"]))
}

cat("### 8.4 Movement Effect Interpretation\n\n")
fe_mov <- fixef(model_movement)
cat("| Effect | Beta | Per 0.1 m/s | Interpretation |\n")
cat("|--------|------|-------------|----------------|\n")
cat(sprintf("| Between-subject (mean_rate) | %.4f | %.1f%% | %s |\n",
            fe_mov["mean_rate"], (exp(fe_mov["mean_rate"]*0.1)-1)*100,
            ifelse(fe_mov["mean_rate"] < 0, "More movement = lower error", "More movement = higher error")))
cat(sprintf("| Within-subject (rate_dev) | %.4f | %.1f%% | %s |\n\n",
            fe_mov["rate_dev"], (exp(fe_mov["rate_dev"]*0.1)-1)*100,
            ifelse(fe_mov["rate_dev"] > 0, "Extra movement = difficulty marker", "Extra movement = helpful")))

# ========== Section 9: Spatial Configuration ==========
cat("## 9. Spatial Configuration Effects\n\n")

model_spatial <- lmer(signedError_m ~ stimulusDisparity_m + azimuthSector + elevationCategory +
                        soundType + (1 | participantId),
                      data = analysis_df, REML = TRUE)

cat("### 9.1 Model Summary\n\n")
cat("```\n")
print(summary(model_spatial))
cat("```\n\n")

cat("### 9.2 Post-Hoc: Azimuth Sectors\n\n")
cat("```\n")
emm_az <- emmeans(model_spatial, ~ azimuthSector)
print(pairs(emm_az))
cat("```\n\n")

cat("### 9.3 Post-Hoc: Elevation Categories\n\n")
cat("```\n")
emm_el <- emmeans(model_spatial, ~ elevationCategory)
print(pairs(emm_el))
cat("```\n\n")

# ========== Section 10: Individual Differences ==========
cat("## 10. Individual Differences in Ventriloquist Susceptibility\n\n")

re <- ranef(model_signed)$participantId
fe_disp <- fixef(model_signed)["stimulusDisparity_m"]

individual_slopes <- data.frame(
  participantId = rownames(re),
  intercept = fe["(Intercept)"] + re[["(Intercept)"]],
  slope = fe_disp + re$stimulusDisparity_m,
  bias_pct = (fe_disp + re$stimulusDisparity_m) * 100
)

cat("### 10.1 Random Slope Variance Summary\n\n")
cat("| Statistic | Value |\n")
cat("|-----------|-------|\n")
cat(sprintf("| Population mean slope | %.1f%% |\n", fe_disp * 100))
cat(sprintf("| Min individual slope | %.1f%% |\n", min(individual_slopes$bias_pct)))
cat(sprintf("| Max individual slope | %.1f%% |\n", max(individual_slopes$bias_pct)))
cat(sprintf("| SD of individual slopes | %.1f%% |\n", sd(individual_slopes$bias_pct)))
cat(sprintf("| Range | %.1f percentage points |\n\n", max(individual_slopes$bias_pct) - min(individual_slopes$bias_pct)))

cat("### 10.2 Individual Slopes (BLUPs)\n\n")
cat("```\n")
cat("Individual ventriloquist susceptibility (% of disparity captured by vision):\n\n")
print(individual_slopes %>% arrange(desc(bias_pct)), row.names = FALSE)
cat("```\n\n")

# ========== Section 11: Summary of Key Findings ==========
cat("## 11. Summary of Key Findings\n\n")

cat("### Primary Outcome: Ventriloquist Effect\n\n")
cat("| Finding | Statistic | 95% CI |\n")
cat("|---------|-----------|--------|\n")
cat(sprintf("| Visual capture | %.1f%% of disparity | [%.1f%%, %.1f%%] |\n",
            fe_disp * 100, ci_signed["stimulusDisparity_m", 1]*100, ci_signed["stimulusDisparity_m", 2]*100))
cat(sprintf("| Effect per 10cm disparity | %.2f cm toward flash | [%.2f, %.2f] |\n",
            fe_disp * 10, ci_signed["stimulusDisparity_m", 1]*10, ci_signed["stimulusDisparity_m", 2]*10))
cat(sprintf("| Individual range | %.1f%% to %.1f%% | — |\n\n",
            min(individual_slopes$bias_pct), max(individual_slopes$bias_pct)))

cat("### Secondary Outcomes\n\n")
cat("| Outcome | Effect | p-value |\n")
cat("|---------|--------|--------|\n")
sum_signed <- summary(model_signed)$coefficients
cat(sprintf("| Sound type (overall) | See post-hoc | — |\n"))
cat(sprintf("| Learning effect | %.2f%% per trial | %s |\n",
            (exp(fe_acc["trialSequence_c"])-1)*100,
            ifelse(summary(model_accuracy)$coefficients["trialSequence_c", "Pr(>|z|)"] < 0.001, "< .001", 
                   sprintf("%.3f", summary(model_accuracy)$coefficients["trialSequence_c", "Pr(>|z|)"]))))
cat(sprintf("| Movement benefit (between) | %.1f%% per 0.1 m/s | %s |\n\n",
            abs((exp(fe_mov["mean_rate"]*0.1)-1)*100),
            ifelse(abs(fe_mov["mean_rate"]/summary(model_movement)$coefficients["mean_rate", "Std. Error"]) > 1.96, 
                   "< .05", "ns")))

cat("### Model Fit Statistics\n\n")
cat("| Model | Marginal R² | Conditional R² | AIC |\n")
cat("|-------|-------------|----------------|-----|\n")
cat(sprintf("| Signed error (LMM) | %.3f | %.3f | %.1f |\n", 
            r2_signed[1, "R2m"], r2_signed[1, "R2c"], AIC(model_signed)))
if (!is.null(r2_acc)) {
  cat(sprintf("| Accuracy (Gamma GLMM) | %.3f | %.3f | %.1f |\n", 
              r2_acc[1, "R2m"], r2_acc[1, "R2c"], AIC(model_accuracy)))
}
if (!is.null(r2_mov)) {
  cat(sprintf("| Movement (Gamma GLMM) | %.3f | %.3f | %.1f |\n", 
              r2_mov[1, "R2m"], r2_mov[1, "R2c"], AIC(model_movement)))
}
cat("\n")

cat("---\n\n")
cat("*Report generated by statistical_analysis.R*\n")
