#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Independent verification of the experienced-angular-disparity re-analysis.
#
# This script does not source the analysis script. Every quantity is
# recomputed from the prepared trial-level dataset and the trajectory sample
# file, using independently written geometry code, so that the reported
# numbers can be checked rather than repeated.
#
# Sections
#   V0  data, and an independent recomputation of the per-sample angular
#       disparity from head positions and trial geometry
#   V1  independent angular response measure (signed-angle-about-normal form)
#   V2  reproduction of the published metric model
#   V3  all eight reported model fits, coefficients, AIC/BIC, R2
#   V4  likelihood-ratio tests, including symmetric versions that give the
#       added predictor its own random slope
#   V5  conditioning of the primary fit: multi-optimizer refits, vcov
#       condition number, profile and bootstrap intervals
#   V6  precision of the reported null (nominal added to angular)
#   V7  geometric-artefact placebo: synthetic responses generated from the
#       fitted METRIC model, pushed through the same angular pipeline
#   V8  head-distance confound in the angular formulation
#   V9  descriptives, within-trial variation, correlations, diagnostics
#
# Inputs : results_revision/analysis_df_revision.rds
#          results_revision/trajectory_samples_angular.rds
# Outputs: results_revision/verify_angular_*.csv, verify_angular_log.txt
#
# R 4.4 arm64 framework build. Deterministic (set.seed used for the bootstrap,
# the placebo noise and DHARMa).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(lmerTest)
  library(DHARMa)
})

proj_dir <- "/Users/bart/Library/CloudStorage/OneDrive-Personal/Documents/R/PawelPerkowskiProject_6dofCorrelations"
res_dir  <- file.path(proj_dir, "results_revision")

log_con <- file(file.path(res_dir, "verify_angular_log.txt"), open = "wt")
say <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n", sep = "")
  cat(txt, "\n", sep = "", file = log_con)
}
say_df <- function(df, digits = 5) {
  for (l in capture.output(print(as.data.frame(df), digits = digits, row.names = FALSE))) say(l)
}
rule <- function(title) {
  say(""); say(strrep("=", 78)); say(title); say(strrep("=", 78))
}
# Comparison helper: reported value against recomputed value.
checks <- list()
chk <- function(label, reported, recomputed, tol = 5e-4) {
  ok <- is.na(reported) || abs(reported - recomputed) <= tol
  checks[[length(checks) + 1]] <<- data.frame(
    quantity = label, reported = reported, recomputed = recomputed,
    abs_diff = abs(reported - recomputed), tol = tol, matches = ok)
  say(sprintf("  %-62s reported %-14s recomputed %-14s %s", label,
              formatC(reported, format = "g", digits = 6),
              formatC(recomputed, format = "g", digits = 6),
              ifelse(ok, "MATCH", "*** MISMATCH ***")))
  invisible(ok)
}

rule("INDEPENDENT VERIFICATION OF THE ANGULAR DISPARITY RE-ANALYSIS")
say("Run: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", R.version.string)
say("lme4 ", as.character(packageVersion("lme4")),
    " | lmerTest ", as.character(packageVersion("lmerTest")),
    " | glmmTMB ", as.character(packageVersion("glmmTMB")),
    " | DHARMa ", as.character(packageVersion("DHARMa")))

# ---------------------------------------------------------------------------
# V0  data and independent angular disparity
# ---------------------------------------------------------------------------
rule("V0: DATA AND INDEPENDENT RECOMPUTATION OF ANGULAR DISPARITY")

analysis_df <- readRDS(file.path(res_dir, "analysis_df_revision.rds"))
samples     <- readRDS(file.path(res_dir, "trajectory_samples_angular.rds"))

say("Trials ", nrow(analysis_df), " | participants ", n_distinct(analysis_df$participantId),
    " | trajectory samples ", nrow(samples))

geom <- analysis_df %>%
  select(participantId, trialSequenceNum,
         sound_x, sound_y, sound_z, flash_x, flash_y, flash_z,
         response_x, response_y, response_z)

samp <- samples %>%
  select(participantId, trialSequenceNum, sample_idx, w_time, px, py, pz,
         angDisp_stored = angDisp_deg) %>%
  left_join(geom, by = c("participantId", "trialSequenceNum"))

# Independent geometry: unit vectors, angle between them, and the signed
# in-plane angle of the response direction expressed as a rotation about the
# source-cue normal. This is algebraically equivalent to the projection form
# used in the analysis script but is coded from a different starting point.
unitise <- function(x, y, z) {
  n <- sqrt(x^2 + y^2 + z^2)
  list(x = x / n, y = y / n, z = z / n, n = n)
}
us <- unitise(samp$sound_x - samp$px, samp$sound_y - samp$py, samp$sound_z - samp$pz)
uc <- unitise(samp$flash_x - samp$px, samp$flash_y - samp$py, samp$flash_z - samp$pz)
ur <- unitise(samp$response_x - samp$px, samp$response_y - samp$py, samp$response_z - samp$pz)

cross_x <- function(a, b) a$y * b$z - a$z * b$y
cross_y <- function(a, b) a$z * b$x - a$x * b$z
cross_z <- function(a, b) a$x * b$y - a$y * b$x
dot     <- function(a, b) a$x * b$x + a$y * b$y + a$z * b$z

nx <- cross_x(us, uc); ny <- cross_y(us, uc); nz <- cross_z(us, uc)
nn <- sqrt(nx^2 + ny^2 + nz^2)
samp$angDisp_mine <- atan2(nn, dot(us, uc)) * 180 / pi
nx <- nx / nn; ny <- ny / nn; nz <- nz / nn
srx <- cross_x(us, ur); sry <- cross_y(us, ur); srz <- cross_z(us, ur)
samp$angResp_mine <- atan2(srx * nx + sry * ny + srz * nz, dot(us, ur)) * 180 / pi
samp$dist_source_m <- us$n

say("")
say("Per-sample angular disparity, my code vs the stored trajectory column:")
say(sprintf("  max abs difference = %.3g deg over %d samples",
            max(abs(samp$angDisp_mine - samp$angDisp_stored)), nrow(samp)))

wmean <- function(x, w) sum(x * w) / sum(w)

trial <- samp %>%
  group_by(participantId, trialSequenceNum) %>%
  summarise(
    angDisp_onset_mine = angDisp_mine[which.min(sample_idx)],
    angDisp_mean_mine  = wmean(angDisp_mine, w_time),
    angDisp_final_mine = angDisp_mine[which.max(sample_idx)],
    angDisp_min_mine   = min(angDisp_mine),
    angDisp_max_mine   = max(angDisp_mine),
    angResp_onset_mine = angResp_mine[which.min(sample_idx)],
    angResp_mean_mine  = wmean(angResp_mine, w_time),
    angResp_final_mine = angResp_mine[which.max(sample_idx)],
    dist_onset_mine    = dist_source_m[which.min(sample_idx)],
    dist_mean_mine     = wmean(dist_source_m, w_time),
    n_samp             = n(),
    .groups = "drop")

dat <- analysis_df %>%
  left_join(trial, by = c("participantId", "trialSequenceNum")) %>%
  mutate(
    angDisp_range_deg = angDisp_max_deg - angDisp_min_deg,
    angResp_onset_deg = angResp_onset_mine,
    angResp_mean_deg  = angResp_mean_mine,
    angResp_final_deg = angResp_final_mine,
    angDisp_final_deg = angDisp_final_mine,
    stimDisp_cm       = stimulusDisparity_m - mean(stimulusDisparity_m),
    angDisp_onset_cm  = angDisp_onset_deg - mean(angDisp_onset_deg),
    angDisp_mean_cm   = angDisp_mean_deg - mean(angDisp_mean_deg),
    angDisp_final_cm  = angDisp_final_deg - mean(angDisp_final_deg))

say("")
say("Trial-level angular disparity summaries, my code vs the prepared columns:")
say(sprintf("  onset : max abs diff %.3g deg", max(abs(dat$angDisp_onset_mine - dat$angDisp_onset_deg))))
say(sprintf("  mean  : max abs diff %.3g deg", max(abs(dat$angDisp_mean_mine  - dat$angDisp_mean_deg))))
say(sprintf("  min   : max abs diff %.3g deg", max(abs(dat$angDisp_min_mine   - dat$angDisp_min_deg))))
say(sprintf("  max   : max abs diff %.3g deg", max(abs(dat$angDisp_max_mine   - dat$angDisp_max_deg))))
say(sprintf("  onset head-source distance : max abs diff %.3g m",
            max(abs(dat$dist_onset_mine - dat$dist_at_onset_m))))
say("")
say("Sanity of the angular response definition (should hold by construction):")
say(sprintf("  correlation of angResp_onset with signedError_m / dist_at_onset: r = %.4f",
            cor(dat$angResp_onset_deg, dat$signedError_m / dat$dist_at_onset_m)))
say(sprintf("  correlation of angDisp_onset with stimulusDisparity / dist_at_onset: r = %.4f",
            cor(dat$angDisp_onset_deg, dat$stimulusDisparity_m / dat$dist_at_onset_m)))

# ---------------------------------------------------------------------------
# V1  descriptives of the derived measures
# ---------------------------------------------------------------------------
rule("V1: DESCRIPTIVES OF THE DERIVED MEASURES")
desc <- data.frame(
  measure = c("angDisp_onset_deg", "angDisp_mean_deg", "angDisp_final_deg",
              "angResp_onset_deg", "angResp_mean_deg", "angResp_final_deg",
              "dist_at_onset_m", "signedError_m"),
  mean = sapply(list(dat$angDisp_onset_deg, dat$angDisp_mean_deg, dat$angDisp_final_deg,
                     dat$angResp_onset_deg, dat$angResp_mean_deg, dat$angResp_final_deg,
                     dat$dist_at_onset_m, dat$signedError_m), mean),
  sd = sapply(list(dat$angDisp_onset_deg, dat$angDisp_mean_deg, dat$angDisp_final_deg,
                   dat$angResp_onset_deg, dat$angResp_mean_deg, dat$angResp_final_deg,
                   dat$dist_at_onset_m, dat$signedError_m), sd),
  median = sapply(list(dat$angDisp_onset_deg, dat$angDisp_mean_deg, dat$angDisp_final_deg,
                       dat$angResp_onset_deg, dat$angResp_mean_deg, dat$angResp_final_deg,
                       dat$dist_at_onset_m, dat$signedError_m), median),
  min = sapply(list(dat$angDisp_onset_deg, dat$angDisp_mean_deg, dat$angDisp_final_deg,
                    dat$angResp_onset_deg, dat$angResp_mean_deg, dat$angResp_final_deg,
                    dat$dist_at_onset_m, dat$signedError_m), min),
  max = sapply(list(dat$angDisp_onset_deg, dat$angDisp_mean_deg, dat$angDisp_final_deg,
                    dat$angResp_onset_deg, dat$angResp_mean_deg, dat$angResp_final_deg,
                    dat$dist_at_onset_m, dat$signedError_m), max))
say_df(desc)
write.csv(desc, file.path(res_dir, "verify_angular_descriptives.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# V2  published metric model
# ---------------------------------------------------------------------------
rule("V2: PUBLISHED METRIC MODEL, REPRODUCED")

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
fitq <- function(fml, reml, data_sym = quote(dat))
  suppressWarnings(do.call("lmer", list(formula = fml, data = data_sym,
                                        REML = reml, control = quote(ctrl))))

A0pub <- fitq(signedError_m ~ stimulusDisparity_m + soundType +
                (1 + stimulusDisparity_m | participantId), TRUE)
s <- summary(A0pub)$coefficients
say_df(cbind(term = rownames(s), as.data.frame(s)))
ciw <- confint(A0pub, method = "Wald", parm = "stimulusDisparity_m")
say(sprintf("Wald 95%% CI for the disparity slope: [%.4f, %.4f]", ciw[1], ciw[2]))
chk("A0pub slope",        0.34660, unname(s["stimulusDisparity_m", "Estimate"]), 1e-4)
chk("A0pub SE",           0.05351, unname(s["stimulusDisparity_m", "Std. Error"]), 1e-4)
chk("A0pub t",            6.477,   unname(s["stimulusDisparity_m", "t value"]), 5e-3)
chk("A0pub CI lo",        0.2417,  ciw[1], 5e-4)
chk("A0pub CI hi",        0.4515,  ciw[2], 5e-4)
chk("A0pub intercept (m)", 0.0086739, unname(s["(Intercept)", "Estimate"]), 1e-5)
chk("A0pub intercept p",  0.6876,  unname(s["(Intercept)", "Pr(>|t|)"]), 5e-4)

# The centred fit, to test the claim that the published intercept of 0.009 m is
# only reproducible from the uncentred model.
A0c <- fitq(signedError_m ~ stimDisp_cm + soundType +
              (1 + stimDisp_cm | participantId), TRUE)
sc <- summary(A0c)$coefficients
say("")
say(sprintf("Centred fit intercept = %.6f m (p = %.4f); uncentred intercept = %.6f m",
            sc["(Intercept)", "Estimate"], sc["(Intercept)", "Pr(>|t|)"],
            s["(Intercept)", "Estimate"]))
say(sprintf("Slope identical across centring: %.6f vs %.6f; logLik %.5f vs %.5f",
            sc["stimDisp_cm", "Estimate"], s["stimulusDisparity_m", "Estimate"],
            as.numeric(logLik(A0c)), as.numeric(logLik(A0pub))))

# ---------------------------------------------------------------------------
# V3  the eight reported fits
# ---------------------------------------------------------------------------
rule("V3: THE REPORTED MODEL SET")

f_corr <- function(resp, pred) as.formula(sprintf("%s ~ %s + soundType + (1 + %s | participantId)", resp, pred, pred))
f_zero <- function(resp, pred) as.formula(sprintf("%s ~ %s + soundType + (1 + %s || participantId)", resp, pred, pred))

fit_ok <- function(f) is.null(f@optinfo$conv$lme4$messages) && !isSingular(f, tol = 1e-4)
fit_msg <- function(f) paste(f@optinfo$conv$lme4$messages, collapse = " | ")

say("Random-effects ladder, independently rerun (correlated form first):")
for (resp in c("signedError_m", "angResp_onset_deg")) {
  for (pred in c("stimDisp_cm", "angDisp_onset_cm", "angDisp_mean_cm")) {
    mr <- fitq(f_corr(resp, pred), TRUE); mm <- fitq(f_corr(resp, pred), FALSE)
    say(sprintf("  %-18s ~ %-17s correlated : REML sing=%-5s ML sing=%-5s  reml_corr=%+.3f  %s",
                resp, pred, isSingular(mr, 1e-4), isSingular(mm, 1e-4),
                attr(VarCorr(mr)$participantId, "correlation")[1, 2],
                paste(unique(c(fit_msg(mr), fit_msg(mm))), collapse = " ")))
    zr <- fitq(f_zero(resp, pred), TRUE); zm <- fitq(f_zero(resp, pred), FALSE)
    say(sprintf("  %-18s ~ %-17s zero-corr  : REML sing=%-5s ML sing=%-5s  %s",
                resp, pred, isSingular(zr, 1e-4), isSingular(zm, 1e-4),
                paste(unique(c(fit_msg(zr), fit_msg(zm))), collapse = " ")))
  }
}

A0 <- fitq(f_corr("signedError_m", "stimDisp_cm"), TRUE)
A1 <- fitq(f_corr("signedError_m", "angDisp_onset_cm"), TRUE)
A2 <- fitq(f_corr("signedError_m", "angDisp_mean_cm"), TRUE)
B0 <- fitq(f_zero("angResp_onset_deg", "stimDisp_cm"), TRUE)
B1 <- fitq(f_zero("angResp_onset_deg", "angDisp_onset_cm"), TRUE)
B2 <- fitq(f_zero("angResp_onset_deg", "angDisp_mean_cm"), TRUE)
B2m <- fitq(f_zero("angResp_mean_deg", "angDisp_mean_cm"), TRUE)

A0ml <- fitq(f_corr("signedError_m", "stimDisp_cm"), FALSE)
A1ml <- fitq(f_corr("signedError_m", "angDisp_onset_cm"), FALSE)
A2ml <- fitq(f_corr("signedError_m", "angDisp_mean_cm"), FALSE)
B0ml <- fitq(f_zero("angResp_onset_deg", "stimDisp_cm"), FALSE)
B1ml <- fitq(f_zero("angResp_onset_deg", "angDisp_onset_cm"), FALSE)
B2ml <- fitq(f_zero("angResp_onset_deg", "angDisp_mean_cm"), FALSE)
B2mml <- fitq(f_zero("angResp_mean_deg", "angDisp_mean_cm"), FALSE)

# Nakagawa/Johnson R2. The random component is the mean over observations of the
# variance contributed by the random effects, computed here from the full sparse
# Z matrix and the Gp block pointers rather than from mmList, so that it is an
# independent implementation of the same quantity.
r2_nj <- function(fit) {
  vf <- as.numeric(var(as.vector(model.matrix(fit) %*% fixef(fit))))
  vc <- VarCorr(fit)
  Z  <- as.matrix(getME(fit, "Z"))     # n x sum(q_i * nlev_i), level-major
  Gp <- getME(fit, "Gp")
  vr <- 0
  for (i in seq_along(vc)) {
    S  <- as.matrix(vc[[i]])
    q  <- nrow(S)
    cols <- (Gp[i] + 1):Gp[i + 1]
    Zb <- Z[, cols, drop = FALSE]
    nl <- length(cols) / q
    # within a block lme4 stores the q components consecutively for each level
    Sbig <- kronecker(diag(nl), S)
    vr <- vr + mean(rowSums((Zb %*% Sbig) * Zb))
  }
  ve <- sigma(fit)^2
  c(R2m = vf / (vf + vr + ve), R2c = (vf + vr) / (vf + vr + ve))
}
# cross-check against the mmList form used in the analysis script
r2_nj_mm <- function(fit) {
  vf <- as.numeric(var(as.vector(model.matrix(fit) %*% fixef(fit))))
  vc <- VarCorr(fit); mm <- getME(fit, "mmList"); vr <- 0
  for (i in seq_along(mm)) {
    M <- as.matrix(mm[[i]]); S <- as.matrix(vc[[i]])
    vr <- vr + mean(rowSums((M %*% S) * M))
  }
  ve <- sigma(fit)^2
  c(R2m = vf / (vf + vr + ve), R2c = (vf + vr) / (vf + vr + ve))
}

mods <- list(A0 = list(A0, A0ml), A1 = list(A1, A1ml), A2 = list(A2, A2ml),
             B0 = list(B0, B0ml), B1 = list(B1, B1ml), B2 = list(B2, B2ml),
             B2m = list(B2m, B2mml), A0pub = list(A0pub, fitq(
               signedError_m ~ stimulusDisparity_m + soundType +
                 (1 + stimulusDisparity_m | participantId), FALSE)))

comp <- bind_rows(lapply(names(mods), function(nm) {
  fr <- mods[[nm]][[1]]; fm <- mods[[nm]][[2]]
  tn <- attr(terms(fr), "term.labels")[1]
  sm <- summary(fr)$coefficients
  an <- anova(fr)
  r2 <- r2_nj(fr)
  data.frame(model = nm, response = as.character(formula(fr))[2], predictor = tn,
             n = nobs(fr), estimate = sm[tn, "Estimate"], se = sm[tn, "Std. Error"],
             df = sm[tn, "df"], t = sm[tn, "t value"], p = sm[tn, "Pr(>|t|)"],
             F_value = an[tn, "F value"], F_df2 = an[tn, "DenDF"],
             ci_lo = sm[tn, "Estimate"] - 1.959964 * sm[tn, "Std. Error"],
             ci_hi = sm[tn, "Estimate"] + 1.959964 * sm[tn, "Std. Error"],
             AIC_ml = AIC(fm), BIC_ml = BIC(fm), logLik_ml = as.numeric(logLik(fm)),
             R2m = unname(r2["R2m"]), R2c = unname(r2["R2c"]),
             R2m_mmList = unname(r2_nj_mm(fr)["R2m"]),
             R2c_mmList = unname(r2_nj_mm(fr)["R2c"]),
             singular = isSingular(fr, 1e-4))
}))
say(""); say_df(comp)
write.csv(comp, file.path(res_dir, "verify_angular_model_comparison.csv"), row.names = FALSE)

say("")
say("Headline coefficient checks:")
g <- function(m, col) comp[comp$model == m, col]
chk("B1 slope (deg/deg)",  0.33590, g("B1", "estimate"), 2e-3)
chk("B1 SE",               0.06600, g("B1", "se"), 2e-3)
chk("B1 Wald CI lo",       0.2066,  g("B1", "ci_lo"), 3e-3)
chk("B1 Wald CI hi",       0.4651,  g("B1", "ci_hi"), 3e-3)
chk("B1 F",                25.93,   g("B1", "F_value"), 0.2)
chk("B2m slope",           0.42936, g("B2m", "estimate"), 2e-3)
chk("B2m CI lo",           0.2977,  g("B2m", "ci_lo"), 3e-3)
chk("B2m CI hi",           0.5611,  g("B2m", "ci_hi"), 3e-3)
chk("B2 slope",            0.32298, g("B2", "estimate"), 2e-3)
chk("B0 slope (deg/m)",    24.386,  g("B0", "estimate"), 0.05)
chk("A1 slope (m/deg)",    0.0012080, g("A1", "estimate"), 2e-5)
chk("A2 slope (m/deg)",    0.0027400, g("A2", "estimate"), 2e-5)
say("")
say("Fit statistics:")
chk("A0 AIC",   -510.473, g("A0", "AIC_ml"), 0.05)
chk("A1 AIC",   -447.261, g("A1", "AIC_ml"), 0.05)
chk("A2 AIC",   -508.038, g("A2", "AIC_ml"), 0.05)
chk("B0 AIC",   6889.209, g("B0", "AIC_ml"), 0.2)
chk("B1 AIC",   6810.081, g("B1", "AIC_ml"), 0.2)
chk("B2 AIC",   6847.3,   g("B2", "AIC_ml"), 0.3)
chk("B1 vs B0 dAIC", 79.13, g("B0", "AIC_ml") - g("B1", "AIC_ml"), 0.3)
chk("A0 R2m",   0.08173, g("A0", "R2m"), 1e-3)
chk("A0 R2c",   0.21080, g("A0", "R2c"), 1e-3)
chk("B1 R2m",   0.09898, g("B1", "R2m"), 1e-3)
chk("B1 R2c",   0.18335, g("B1", "R2c"), 1e-3)
chk("B0 R2m",   0.02492, g("B0", "R2m"), 1e-3)
chk("A1 R2m",   0.02671, g("A1", "R2m"), 1e-3)
chk("A2 R2m",   0.08895, g("A2", "R2m"), 1e-3)
chk("B2m R2m",  0.1403,  g("B2m", "R2m"), 1.5e-3)

# ---------------------------------------------------------------------------
# V4  likelihood-ratio tests
# ---------------------------------------------------------------------------
rule("V4: LIKELIHOOD-RATIO TESTS")

lrt <- function(label, resp, base, add, corr) {
  ff <- if (corr) f_corr else f_zero
  red <- ff(resp, base)
  full <- as.formula(sub(paste0("~ ", base, " \\+"), paste0("~ ", base, " + ", add, " +"),
                         paste(deparse(red), collapse = "")))
  a <- anova(fitq(red, FALSE), fitq(full, FALSE))
  say(sprintf("  %-64s chi2(%d) = %8.4f  p = %.4g", label, a$Df[2], a$Chisq[2], a$`Pr(>Chisq)`[2]))
  data.frame(comparison = label, chisq = a$Chisq[2], df = a$Df[2], p = a$`Pr(>Chisq)`[2])
}
say("As specified in the analysis (random slope for the BASE predictor only):")
lr <- bind_rows(
  lrt("A: onset angular added to nominal",        "signedError_m", "stimDisp_cm", "angDisp_onset_cm", TRUE),
  lrt("A: time-mean angular added to nominal",    "signedError_m", "stimDisp_cm", "angDisp_mean_cm", TRUE),
  lrt("A: nominal added to onset angular",        "signedError_m", "angDisp_onset_cm", "stimDisp_cm", TRUE),
  lrt("B: onset angular added to nominal",        "angResp_onset_deg", "stimDisp_cm", "angDisp_onset_cm", FALSE),
  lrt("B: nominal added to onset angular",        "angResp_onset_deg", "angDisp_onset_cm", "stimDisp_cm", FALSE),
  lrt("B: time-mean angular added to onset angular", "angResp_onset_deg", "angDisp_onset_cm", "angDisp_mean_cm", FALSE))

chk("LRT B: angular added to nominal chi2", 48.4194, lr$chisq[lr$comparison == "B: onset angular added to nominal"], 0.05)
chk("LRT B: nominal added to angular chi2", 0.10872, lr$chisq[lr$comparison == "B: nominal added to onset angular"], 0.01)
chk("LRT A: time-mean added to nominal chi2", 12.34507, lr$chisq[lr$comparison == "A: time-mean angular added to nominal"], 0.05)
chk("LRT A: nominal added to onset angular chi2", 55.93277, lr$chisq[lr$comparison == "A: nominal added to onset angular"], 0.05)
chk("LRT A: onset angular added to nominal chi2", 0.095, lr$chisq[lr$comparison == "A: onset angular added to nominal"], 0.01)

say("")
say("Symmetric versions: both predictors in the fixed part AND both given a random")
say("slope, so that the test of the added predictor is not anti-conservative from an")
say("omitted random slope. Reduced model drops only the fixed effect.")
sym_lrt <- function(label, resp, p1, p2) {
  full <- as.formula(sprintf("%s ~ %s + %s + soundType + (1 + %s + %s || participantId)",
                             resp, p1, p2, p1, p2))
  red  <- as.formula(sprintf("%s ~ %s + soundType + (1 + %s + %s || participantId)",
                             resp, p1, p1, p2))
  mf <- fitq(full, FALSE); mr <- fitq(red, FALSE)
  a <- anova(mr, mf)
  say(sprintf("  %-64s chi2(%d) = %8.4f  p = %.4g  [sing full=%s red=%s]",
              label, a$Df[2], a$Chisq[2], a$`Pr(>Chisq)`[2],
              isSingular(mf, 1e-4), isSingular(mr, 1e-4)))
  data.frame(comparison = label, chisq = a$Chisq[2], df = a$Df[2], p = a$`Pr(>Chisq)`[2])
}
lr_sym <- bind_rows(
  sym_lrt("B sym: onset angular added to nominal", "angResp_onset_deg", "stimDisp_cm", "angDisp_onset_cm"),
  sym_lrt("B sym: nominal added to onset angular", "angResp_onset_deg", "angDisp_onset_cm", "stimDisp_cm"),
  sym_lrt("A sym: time-mean angular added to nominal", "signedError_m", "stimDisp_cm", "angDisp_mean_cm"),
  sym_lrt("A sym: onset angular added to nominal", "signedError_m", "stimDisp_cm", "angDisp_onset_cm"))
write.csv(bind_rows(lr, lr_sym), file.path(res_dir, "verify_angular_lrt.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# V5  conditioning of the primary fit
# ---------------------------------------------------------------------------
rule("V5: CONDITIONING AND INTERVAL STABILITY OF THE PRIMARY FIT B1")

say("Fixed-effect vcov of B1:")
V <- as.matrix(vcov(B1))
ev <- eigen(V, symmetric = TRUE)$values
say(sprintf("  dimension %d, eigenvalues %s", nrow(V),
            paste(formatC(ev, format = "e", digits = 3), collapse = ", ")))
say(sprintf("  condition number = %.4g, min eigenvalue = %.4g, all positive = %s",
            max(ev) / min(ev), min(ev), all(ev > 0)))
say(sprintf("  random-effect theta = %s, sigma = %.4f",
            paste(sprintf("%.5f", getME(B1, "theta")), collapse = ", "), sigma(B1)))

f_b1 <- f_zero("angResp_onset_deg", "angDisp_onset_cm")

say("")
say("Multi-optimizer refit of B1 (each optimizer run explicitly):")
opts <- c("bobyqa", "Nelder_Mead", "nlminbwrap", "nloptwrap")
af_tab <- bind_rows(lapply(opts, function(nm) {
  cc <- lmerControl(optimizer = nm)
  z <- try(suppressWarnings(lmer(f_b1, data = dat, REML = TRUE, control = cc)), silent = TRUE)
  if (inherits(z, "try-error")) return(data.frame(optimizer = nm, estimate = NA_real_,
                                                  se = NA_real_, logLik = NA_real_, singular = NA))
  sm <- summary(z)$coefficients
  data.frame(optimizer = nm, estimate = sm["angDisp_onset_cm", "Estimate"],
             se = sm["angDisp_onset_cm", "Std. Error"],
             logLik = as.numeric(logLik(z)), singular = isSingular(z, 1e-4))
}))
say_df(af_tab)
af_tab <- af_tab[!is.na(af_tab$estimate), ]
say(sprintf("  spread across optimizers: slope range %.6f, SE range %.6f",
            diff(range(af_tab$estimate)), diff(range(af_tab$se))))

say("")
say("Profile and bootstrap intervals for the B1 slope:")
prof <- suppressWarnings(confint(B1, parm = "angDisp_onset_cm", method = "profile"))
say(sprintf("  profile 95%% CI [%.4f, %.4f]", prof[1, 1], prof[1, 2]))
chk("B1 profile CI lo", 0.2060, prof[1, 1], 3e-3)
chk("B1 profile CI hi", 0.4681, prof[1, 2], 3e-3)

set.seed(4242)
pids <- unique(dat$participantId)
boot <- vapply(seq_len(1000), function(i) {
  take <- sample(pids, length(pids), replace = TRUE)
  bd <- bind_rows(lapply(seq_along(take), function(j) {
    z <- dat[dat$participantId == take[j], , drop = FALSE]
    z$participantId <- paste0("B", j); z }))
  m <- try(suppressWarnings(lmer(f_b1, data = bd, REML = TRUE, control = ctrl)), silent = TRUE)
  if (inherits(m, "try-error")) NA_real_ else unname(fixef(m)["angDisp_onset_cm"])
}, numeric(1))
bci <- quantile(boot, c(.025, .975), na.rm = TRUE)
say(sprintf("  cluster bootstrap: %d/1000 succeeded, mean %.4f, SD %.4f, 95%% CI [%.4f, %.4f], %d resamples <= 0",
            sum(!is.na(boot)), mean(boot, na.rm = TRUE), sd(boot, na.rm = TRUE),
            bci[1], bci[2], sum(boot <= 0, na.rm = TRUE)))

say("")
say("Student-t heavy-tail refits via glmmTMB (checking the reported non-convergence):")
t_try <- function(label, cl) {
  m <- try(suppressWarnings(eval(cl)), silent = TRUE)
  if (inherits(m, "try-error")) { say(sprintf("  %-34s ERROR", label)); return(NULL) }
  pd <- isTRUE(m$sdr$pdHess)
  cf <- tryCatch(summary(m)$coefficients$cond["angDisp_onset_cm", ], error = function(e) rep(NA_real_, 4))
  say(sprintf("  %-34s pdHess = %-5s beta = %.4f  SE = %s", label, pd, cf[1],
              formatC(cf[2], format = "g", digits = 4)))
  data.frame(check = label, pdHess = pd, estimate = cf[1], se = cf[2])
}
tt <- bind_rows(
  t_try("student_t random slope", quote(glmmTMB::glmmTMB(
    angResp_onset_deg ~ angDisp_onset_cm + soundType + (1 | participantId) +
      (0 + angDisp_onset_cm | participantId), data = dat,
    family = glmmTMB::t_family(link = "identity")))),
  t_try("student_t random intercept only", quote(glmmTMB::glmmTMB(
    angResp_onset_deg ~ angDisp_onset_cm + soundType + (1 | participantId), data = dat,
    family = glmmTMB::t_family(link = "identity")))),
  t_try("gaussian glmmTMB, same RE as B1", quote(glmmTMB::glmmTMB(
    angResp_onset_deg ~ angDisp_onset_cm + soundType + (1 | participantId) +
      (0 + angDisp_onset_cm | participantId), data = dat))))
write.csv(tt, file.path(res_dir, "verify_angular_heavytail.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# V6  precision of the reported null
# ---------------------------------------------------------------------------
rule("V6: HOW INFORMATIVE IS THE 'NOMINAL ADDS NOTHING' NULL?")

both_B <- fitq(as.formula(paste("angResp_onset_deg ~ angDisp_onset_cm + stimDisp_cm + soundType +",
                                "(1 + angDisp_onset_cm || participantId)")), TRUE)
sb <- summary(both_B)$coefficients
say_df(cbind(term = rownames(sb), as.data.frame(sb)))
est <- sb["stimDisp_cm", "Estimate"]; se <- sb["stimDisp_cm", "Std. Error"]
say(sprintf("Nominal slope in the joint angular model: %.3f deg/m, 95%% CI [%.3f, %.3f]",
            est, est - 1.96 * se, est + 1.96 * se))
say(sprintf("For comparison, the nominal slope ALONE on this response (B0) is %.3f deg/m.",
            g("B0", "estimate")))
say("The interval is wide relative to that marginal effect, so the null is an absence")
say("of evidence at this precision, not evidence that nominal separation is irrelevant.")

# ---------------------------------------------------------------------------
# V7  geometric-artefact placebo
# ---------------------------------------------------------------------------
rule("V7: PLACEBO -- ARE THE ANGULAR RESULTS A GEOMETRIC RESTATEMENT OF THE METRIC FIT?")

say("Both sides of the angular formulation are evaluated from the same onset head")
say("pose, so theta_R is approximately signedError/d and theta_D approximately")
say("stimulusDisparity/d with the SAME d. The test below asks how much of the")
say("angular result survives that shared geometry.")
say("")
say("Synthetic responses are constructed from the observed geometry with an")
say("in-line displacement taken from the fitted METRIC model (population level,")
say("no random effects), so that by construction the only truth in the data is a")
say("constant metric capture fraction. The angular pipeline is then rerun on them.")

u_sf <- with(dat, {
  dx <- flash_x - sound_x; dy <- flash_y - sound_y; dz <- flash_z - sound_z
  n <- sqrt(dx^2 + dy^2 + dz^2); cbind(dx / n, dy / n, dz / n) })
fitted_metric <- predict(A0pub, re.form = NA)          # population-level metric prediction

onset_head <- dat %>% select(participantId, trialSequenceNum) %>%
  left_join(samp %>% group_by(participantId, trialSequenceNum) %>%
              slice_min(sample_idx, n = 1) %>%
              select(participantId, trialSequenceNum, px, py, pz) %>% ungroup(),
            by = c("participantId", "trialSequenceNum"))

ang_resp_onset <- function(rx, ry, rz) {
  s <- unitise(dat$sound_x - onset_head$px, dat$sound_y - onset_head$py, dat$sound_z - onset_head$pz)
  cc <- unitise(dat$flash_x - onset_head$px, dat$flash_y - onset_head$py, dat$flash_z - onset_head$pz)
  r <- unitise(rx - onset_head$px, ry - onset_head$py, rz - onset_head$pz)
  ax <- cross_x(s, cc); ay <- cross_y(s, cc); az <- cross_z(s, cc)
  an <- sqrt(ax^2 + ay^2 + az^2); ax <- ax / an; ay <- ay / an; az <- az / an
  atan2(cross_x(s, r) * ax + cross_y(s, r) * ay + cross_z(s, r) * az, dot(s, r)) * 180 / pi
}
# check the helper reproduces the observed angular response
say(sprintf("  helper check on the real responses: max abs diff %.3g deg",
            max(abs(ang_resp_onset(dat$response_x, dat$response_y, dat$response_z) -
                      dat$angResp_onset_deg))))

set.seed(20260806)
# residual orthogonal displacement of the real response, kept so the synthetic
# responses have realistic off-axis scatter
proj <- with(dat, cbind(response_x - sound_x, response_y - sound_y, response_z - sound_z))
along <- rowSums(proj * u_sf)
orth  <- proj - along * u_sf

make_synth <- function(inline) {
  p <- inline * u_sf + orth
  data.frame(x = dat$sound_x + p[, 1], y = dat$sound_y + p[, 2], z = dat$sound_z + p[, 3])
}
placebos <- list(
  "P1 metric model truth (b = 0.3466 x disparity)" = make_synth(fitted_metric),
  "P2 constant in-line shift (mean signed error)"  = make_synth(rep(mean(dat$signedError_m), nrow(dat))),
  "P3 zero in-line shift (no capture at all)"      = make_synth(rep(0, nrow(dat))))

pl_rows <- bind_rows(lapply(names(placebos), function(nm) {
  sy <- placebos[[nm]]
  d2 <- dat
  d2$angResp_onset_deg <- ang_resp_onset(sy$x, sy$y, sy$z)
  d2$signedError_m <- rowSums((cbind(sy$x, sy$y, sy$z) -
                                 cbind(dat$sound_x, dat$sound_y, dat$sound_z)) * u_sf)
  m_ang <- suppressWarnings(lmer(f_b1, data = d2, REML = TRUE, control = ctrl))
  m_nom <- suppressWarnings(lmer(f_zero("angResp_onset_deg", "stimDisp_cm"),
                                 data = d2, REML = TRUE, control = ctrl))
  m_ang_ml <- suppressWarnings(lmer(f_b1, data = d2, REML = FALSE, control = ctrl))
  m_nom_ml <- suppressWarnings(lmer(f_zero("angResp_onset_deg", "stimDisp_cm"),
                                    data = d2, REML = FALSE, control = ctrl))
  sm <- summary(m_ang)$coefficients
  say(sprintf("  %-48s angular slope = %+.4f (SE %.4f, t = %5.2f), dAIC(nominal - angular) = %+.1f, R2m = %.4f",
              nm, sm["angDisp_onset_cm", "Estimate"], sm["angDisp_onset_cm", "Std. Error"],
              sm["angDisp_onset_cm", "t value"], AIC(m_nom_ml) - AIC(m_ang_ml),
              unname(r2_nj(m_ang)["R2m"])))
  data.frame(placebo = nm, slope = sm["angDisp_onset_cm", "Estimate"],
             se = sm["angDisp_onset_cm", "Std. Error"], t = sm["angDisp_onset_cm", "t value"],
             p = sm["angDisp_onset_cm", "Pr(>|t|)"],
             dAIC_nominal_minus_angular = AIC(m_nom_ml) - AIC(m_ang_ml),
             R2m_angular = unname(r2_nj(m_ang)["R2m"]),
             R2m_nominal = unname(r2_nj(m_nom)["R2m"]))
}))
say("")
say("Observed data for comparison:")
say(sprintf("  %-48s angular slope = %+.4f (SE %.4f, t = %5.2f), dAIC(nominal - angular) = %+.1f, R2m = %.4f",
            "OBSERVED", g("B1", "estimate"), g("B1", "se"), g("B1", "t"),
            g("B0", "AIC_ml") - g("B1", "AIC_ml"), g("B1", "R2m")))
write.csv(pl_rows, file.path(res_dir, "verify_angular_placebo.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# V8  head-distance confound
# ---------------------------------------------------------------------------
rule("V8: HEAD-TO-SOURCE DISTANCE AS A SHARED DRIVER")

say(sprintf("  r(dist_at_onset, angDisp_onset) = %+.4f", cor(dat$dist_at_onset_m, dat$angDisp_onset_deg)))
say(sprintf("  r(dist_at_onset, |angResp_onset|) = %+.4f", cor(dat$dist_at_onset_m, abs(dat$angResp_onset_deg))))
say(sprintf("  r(1/dist_at_onset, angDisp_onset) = %+.4f", cor(1 / dat$dist_at_onset_m, dat$angDisp_onset_deg)))

dat$invdist_c <- 1 / dat$dist_at_onset_m - mean(1 / dat$dist_at_onset_m)
m_d1 <- suppressWarnings(lmer(angResp_onset_deg ~ angDisp_onset_cm + invdist_c + soundType +
                                (1 + angDisp_onset_cm || participantId), data = dat,
                              REML = TRUE, control = ctrl))
say("")
say("B1 with inverse onset distance added as a covariate:")
say_df(cbind(term = rownames(summary(m_d1)$coefficients), as.data.frame(summary(m_d1)$coefficients)))

# Restrict to trials where the head is not in the near field at onset.
for (thr in c(0.5, 0.75, 1.0)) {
  d3 <- dat %>% filter(dist_at_onset_m >= thr)
  m3 <- suppressWarnings(lmer(f_b1, data = d3, REML = TRUE, control = ctrl))
  s3 <- summary(m3)$coefficients
  say(sprintf("  onset distance >= %.2f m: n = %3d, slope = %+.4f (SE %.4f, t = %5.2f)",
              thr, nrow(d3), s3["angDisp_onset_cm", "Estimate"],
              s3["angDisp_onset_cm", "Std. Error"], s3["angDisp_onset_cm", "t value"]))
}

# ---------------------------------------------------------------------------
# V9  descriptives, robustness, diagnostics
# ---------------------------------------------------------------------------
rule("V9: WITHIN-TRIAL VARIATION, CORRELATIONS, ROBUSTNESS, DIAGNOSTICS")

rq <- quantile(dat$angDisp_range_deg, c(0, .05, .25, .5, .75, .95, 1))
say(sprintf("Within-trial range: mean %.1f SD %.1f min %.1f Q05 %.1f Q25 %.1f med %.1f Q75 %.1f Q95 %.1f max %.1f",
            mean(dat$angDisp_range_deg), sd(dat$angDisp_range_deg), rq[1], rq[2], rq[3], rq[4], rq[5], rq[6], rq[7]))
chk("within-trial range mean", 91.1, mean(dat$angDisp_range_deg), 0.05)
chk("within-trial range SD",   45.5, sd(dat$angDisp_range_deg), 0.05)
chk("within-trial range median", 87.7, median(dat$angDisp_range_deg), 0.05)
for (thr in c(30, 60, 90, 120))
  say(sprintf("  > %3d deg: %3d trials (%.1f%%)", thr, sum(dat$angDisp_range_deg > thr),
              100 * mean(dat$angDisp_range_deg > thr)))
chk("pct range > 30 deg", 89.9, 100 * mean(dat$angDisp_range_deg > 30), 0.1)
chk("pct range > 90 deg", 48.3, 100 * mean(dat$angDisp_range_deg > 90), 0.1)
chk("median range / onset ratio", 3.36, median(dat$angDisp_range_deg / dat$angDisp_onset_deg), 0.01)
chk("SD of onset angular disparity", 23.1, sd(dat$angDisp_onset_deg), 0.05)

say("")
ang_cols <- c("angDisp_onset_deg", "angDisp_mean_deg", "angDisp_median_deg",
              "angDisp_min_deg", "angDisp_max_deg", "angDisp_range_deg", "angDisp_final_deg")
wc <- function(x, g) x - ave(x, g, FUN = mean)
cor_tab <- bind_rows(lapply(ang_cols, function(cc) {
  ct <- cor.test(dat$stimulusDisparity_m, dat[[cc]])
  data.frame(summary = cc, r = unname(ct$estimate), r2 = unname(ct$estimate)^2,
             r_within = cor(wc(dat$stimulusDisparity_m, dat$participantId),
                            wc(dat[[cc]], dat$participantId)))
}))
say_df(cor_tab)
chk("r nominal vs onset angular", 0.4501, cor_tab$r[1], 1e-3)
chk("r nominal vs mean angular",  0.6572, cor_tab$r[2], 1e-3)
chk("r between the two angular summaries", 0.535,
    cor(dat$angDisp_onset_deg, dat$angDisp_mean_deg), 2e-3)
write.csv(cor_tab, file.path(res_dir, "verify_angular_correlations.csv"), row.names = FALSE)

say("")
say("Per-participant slopes of B1:")
ps <- coef(B1)$participantId[["angDisp_onset_cm"]]
say(sprintf("  min %.3f Q25 %.3f median %.3f Q75 %.3f max %.3f SD %.3f",
            min(ps), quantile(ps, .25), median(ps), quantile(ps, .75), max(ps), sd(ps)))
chk("B1 per-participant slope min", -0.192, min(ps), 2e-3)
chk("B1 per-participant slope max",  0.887, max(ps), 2e-3)
chk("B1 per-participant slope SD",   0.225, sd(ps), 2e-3)
B1_ri <- fitq(as.formula("angResp_onset_deg ~ angDisp_onset_cm + soundType + (1 | participantId)"), TRUE)
a_rs <- anova(B1_ri, B1, refit = FALSE)
say(sprintf("  random-slope LRT: chi2(%d) = %.4f, p = %.4g", a_rs$Df[2], a_rs$Chisq[2], a_rs$`Pr(>Chisq)`[2]))
chk("random-slope LRT chi2", 32.19, a_rs$Chisq[2], 0.05)

say("")
say("Robustness refits of the primary slope:")
rob <- function(label, fml, d) {
  m <- suppressWarnings(lmer(fml, data = d, REML = TRUE, control = ctrl))
  tn <- attr(terms(m), "term.labels")[1]
  s <- summary(m)$coefficients[tn, ]
  say(sprintf("  %-42s n = %3d  beta = %+.4f  SE = %.4f  t = %5.2f  singular = %s",
              label, nobs(m), s["Estimate"], s["Std. Error"], s["t value"], isSingular(m, 1e-4)))
  data.frame(check = label, n = nobs(m), estimate = s["Estimate"], se = s["Std. Error"], row.names = NULL)
}
rob_tab <- bind_rows(
  rob("full sample (B1)", f_b1, dat),
  rob("onset angular disparity <= 90 deg", f_b1, dat %>% filter(angDisp_onset_deg <= 90)),
  rob("onset angular disparity <= 60 deg", f_b1, dat %>% filter(angDisp_onset_deg <= 60)),
  rob("excluding revisited trials", f_b1, dat %>% filter(!trial_has_revisit)),
  rob("final-pose frame", f_zero("angResp_final_deg", "angDisp_final_cm"), dat),
  rob("matched time-weighted frame (B2m)", f_zero("angResp_mean_deg", "angDisp_mean_cm"), dat))
write.csv(rob_tab, file.path(res_dir, "verify_angular_robustness.csv"), row.names = FALSE)
chk("robustness <=90 deg slope", 0.3844, rob_tab$estimate[2], 3e-3)
chk("robustness <=60 deg slope", 0.3790, rob_tab$estimate[3], 3e-3)
chk("robustness no-revisit slope", 0.3371, rob_tab$estimate[4], 3e-3)
chk("robustness final-pose slope", 0.4681, rob_tab$estimate[5], 5e-3)
chk("n trials with onset angle <= 90", 721, rob_tab$n[2], 0)
chk("n trials with onset angle <= 60", 672, rob_tab$n[3], 0)

say("")
say("DHARMa residual diagnostics:")
set.seed(20260805)
simB <- simulateResiduals(B1, n = 1000, seed = 20260805)
ksB <- testUniformity(simB, plot = FALSE); dsB <- testDispersion(simB, plot = FALSE)
simA <- simulateResiduals(A0pub, n = 1000, seed = 20260805)
ksA <- testUniformity(simA, plot = FALSE); dsA <- testDispersion(simA, plot = FALSE)
say(sprintf("  B1     KS D = %.4f p = %.4g | dispersion %.4f p = %.3f", ksB$statistic, ksB$p.value, dsB$statistic, dsB$p.value))
say(sprintf("  A0pub  KS D = %.4f p = %.4g | dispersion %.4f p = %.3f", ksA$statistic, ksA$p.value, dsA$statistic, dsA$p.value))
chk("B1 DHARMa KS D", 0.1342, unname(ksB$statistic), 5e-3)
chk("A0pub DHARMa KS D", 0.0753, unname(ksA$statistic), 5e-3)

# ---------------------------------------------------------------------------
rule("SUMMARY OF CHECKS")
chk_tab <- bind_rows(checks)
say_df(chk_tab)
say("")
say(sprintf("Checks matching: %d of %d", sum(chk_tab$matches), nrow(chk_tab)))
if (any(!chk_tab$matches)) {
  say("Mismatches:")
  say_df(chk_tab[!chk_tab$matches, ])
}
write.csv(chk_tab, file.path(res_dir, "verify_angular_checks.csv"), row.names = FALSE)
say(""); say("Done.")
close(log_con)
