# 6DoF VR Ventriloquist Effect - Statistical Analysis Results

Generated: 2025-12-28 15:21:44 

## 1. Dataset Overview

| Measure | Value |
|---------|-------|
| Participants | 31 |
| Total trials | 744 |
| Trials per participant | 24 |
| Trajectory frames | 251666 |

## 2. Raw Data Sample

First 10 observations:

```
 participantId trialSequenceNum  soundType   disparityRange stimulusDisparity_m participantError_m signedError_m
          P001                1      Flute   Long (40-70cm)                0.40          0.2012292   -0.14826287
          P001                2     Speech   Long (40-70cm)                0.61          0.3156049    0.22471552
          P001                3 Pink Noise   Long (40-70cm)                0.40          0.4371269    0.04484877
          P001                4       Drum Medium (30-40cm)                0.36          0.4964415    0.34577399
          P001                5     Speech Medium (30-40cm)                0.34          0.4911234    0.23153717
          P001                6      Flute   Long (40-70cm)                0.40          0.3913430   -0.13486513
          P001                7 Pink Noise  Short (15-30cm)                0.24          0.5694570    0.10274400
          P001                8 Pink Noise   Long (40-70cm)                0.40          0.6843969    0.48369829
          P001                9 Pink Noise  Short (15-30cm)                0.24          0.9024074   -0.02115637
          P001               10      Flute Medium (30-40cm)                0.38          0.5493741    0.49699437
 ventriloquistBias total_path_length
       -0.37065717          6.181477
        0.36838608          7.406651
        0.11212194         10.535612
        0.96048326          4.186433
        0.68099167          9.682172
       -0.33716282          5.625836
        0.42810000          1.957344
        1.20924571          3.585280
       -0.08815155          3.160918
        1.30787993          2.092318
```

## 3. Descriptive Statistics

### 3.1 Overall Performance

| Measure | Value |
|---------|-------|
| N (trials) | 744 |
| N (participants) | 31 |
| Mean unsigned error | 26.2 cm (SD = 17.0) |
| Median unsigned error | 22.3 cm |
| Range unsigned error | 0.6 - 125.6 cm |
| Mean signed error | 13.7 cm (SD = 18.5) |
| Mean ventriloquist bias | 35.9% (SD = 54.6) |
| Mean path length | 8.36 m (SD = 6.21) |

### 3.2 By Sound Type

| Sound Type | n | Mean Error (SD) cm | Mean Signed (SD) cm | Bias % |
|------------|---|-------------------|---------------------|--------|
| Drum | 186 | 28.4 (18.5) | 14.3 (20.1) | 36.4 |
| Flute | 186 | 25.7 (17.3) | 14.9 (18.4) | 40.1 |
| Speech | 186 | 25.6 (16.4) | 13.8 (19.2) | 37.3 |
| Pink Noise | 186 | 25.1 (15.8) | 11.6 (16.0) | 29.7 |

### 3.3 By Disparity Range

| Disparity Range | n | Mean Disparity | Mean Error (SD) cm | Mean Signed cm | Bias % |
|-----------------|---|----------------|-------------------|----------------|--------|
| Short (15-30cm) | 238 | 22.0 cm | 21.9 (16.0) | 7.9 | 35.0 |
| Medium (30-40cm) | 246 | 35.2 cm | 25.2 (16.3) | 13.6 | 38.3 |
| Long (40-70cm) | 260 | 55.0 cm | 31.0 (17.5) | 19.0 | 34.3 |

### 3.4 By Participant

```
 participantId n_trials mean_error_cm sd_error_cm mean_bias mean_path_m
          P001       24          45.7       20.86    0.5298        5.36
          P002       24          25.2       22.07    0.2927       13.42
          P003       24          46.3       28.50    0.2594        6.03
          P004       24          24.1       13.92    0.3613        5.45
          P005       24          31.7       17.26    0.0685        6.04
          P006       24          42.1       22.24    0.7276        6.18
          P007       24          29.0       19.42    0.6231        7.83
          P008       24          20.4        8.04    0.2820        8.93
          P009       24          20.0       11.04    0.1617        8.68
          P010       24          32.1       13.23    0.3309        7.52
          P011       24          24.2       12.64    0.0347        8.09
          P012       24          33.5       19.12    0.4582        8.26
          P013       24          18.9        8.16    0.3864        5.87
          P014       24          24.5       10.47    0.5399        7.64
          P015       24          23.8       10.99    0.5077        4.86
          P016       24          23.0       11.89    0.4932        4.58
          P017       24          24.9        9.44    0.5581        4.79
          P018       24          21.2        8.17    0.5333        5.52
          P019       24          24.9       13.80    0.4596        8.31
          P020       24          23.5       24.65    0.3797        7.41
          P021       24          20.1        9.87    0.3989        5.40
          P022       24          19.8       10.70    0.1302        8.64
          P023       24          22.3       15.93    0.1303       11.72
          P024       24          16.5        8.00    0.2924        7.35
          P025       24          20.3        8.08    0.1486       15.14
          P026       24          29.4       21.79   -0.0140       21.98
          P027       24          25.4       10.78    0.3987        6.61
          P028       24          37.0       16.11    0.7006        8.40
          P029       24          19.5        8.67    0.3944        6.01
          P030       24          23.8       13.75    0.4744       14.99
          P031       24          18.7       19.66    0.0770       12.23
```

## 4. Correlation Matrix

Pearson correlations among key variables:

```
                    stimulusDisparity_m participantError_m signedError_m ventriloquistBias total_path_length
stimulusDisparity_m               1.000              0.249         0.289             0.010             0.049
participantError_m                0.249              1.000         0.544             0.420             0.097
signedError_m                     0.289              0.544         1.000             0.881             0.004
ventriloquistBias                 0.010              0.420         0.881             1.000            -0.017
total_path_length                 0.049              0.097         0.004            -0.017             1.000
trialSequenceNum                  0.009             -0.086        -0.002            -0.012            -0.067
                    trialSequenceNum
stimulusDisparity_m            0.009
participantError_m            -0.086
signedError_m                 -0.002
ventriloquistBias             -0.012
total_path_length             -0.067
trialSequenceNum               1.000
```

## 5. Main Ventriloquist Effect (Signed Error Model)

### 5.1 Model Specification

| Component | Specification |
|-----------|---------------|
| Response | Signed error (meters) - positive = toward flash |
| Fixed effects | Disparity (m), Sound Type |
| Random effects | Participant intercept + disparity slope |
| Family | Gaussian (LMM) |
| Estimation | REML |

### 5.2 Model Comparison: Random Intercepts vs Random Slopes

```
Data: analysis_df
Models:
model_signed_intercept: signedError_m ~ stimulusDisparity_m + soundType + (1 | participantId)
model_signed: signedError_m ~ stimulusDisparity_m + soundType + (1 + stimulusDisparity_m | participantId)
                       npar     AIC     BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)   
model_signed_intercept    7 -471.67 -439.38 242.83   -485.67                        
model_signed              9 -479.64 -438.13 248.82   -497.64 11.971  2   0.002515 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

**Result:** Random slopes model preferred (Chi-sq = 11.97, df = 2, p = 0.003)

### 5.3 Final Model Summary

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
Formula: signedError_m ~ stimulusDisparity_m + soundType + (1 + stimulusDisparity_m |      participantId)
   Data: analysis_df
Control: lmerControl(optimizer = "bobyqa")

REML criterion at convergence: -497.6

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-4.6991 -0.4712 -0.0219  0.4487  5.4138 

Random effects:
 Groups        Name                Variance Std.Dev. Corr 
 participantId (Intercept)         0.001979 0.04449       
               stimulusDisparity_m 0.036155 0.19015  -0.56
 Residual                          0.026936 0.16412       
Number of obs: 744, groups:  participantId, 31

Fixed effects:
                      Estimate Std. Error         df t value Pr(>|t|)    
(Intercept)          8.674e-03  2.145e-02  5.196e+01   0.404    0.688    
stimulusDisparity_m  3.466e-01  5.351e-02  3.050e+01   6.477 3.43e-07 ***
soundTypeFlute       1.204e-02  1.706e-02  6.888e+02   0.706    0.481    
soundTypeSpeech     -2.502e-04  1.704e-02  6.892e+02  -0.015    0.988    
soundTypePink Noise -2.507e-02  1.704e-02  6.895e+02  -1.471    0.142    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) stmlD_ sndTyF sndTyS
stmlsDsprt_ -0.701                     
soundTypFlt -0.431  0.036              
sondTypSpch -0.411  0.015  0.499       
sndTypPnkNs -0.398  0.001  0.500  0.499
```

### 5.4 Variance Explained (R-squared)

| Type | R² |
|------|----|
| Marginal R² (fixed effects) | 0.082 |
| Conditional R² (fixed + random) | 0.211 |

### 5.5 Intraclass Correlation (ICC)

ICC (participant) = 0.068

Interpretation: 6.8% of variance in signed error is attributable to between-participant differences.

### 5.6 Fixed Effects with 95% CI

| Parameter | Estimate | 95% CI | Interpretation |
|-----------|----------|--------|----------------|
| Intercept | 0.0087 m | [-0.0334, 0.0507] | Baseline at 0 disparity, Drum |
| Disparity | 0.3466 m/m | [0.2417, 0.4515] | 34.7% visual capture |
| Flute vs Drum | 0.0120 m | [-0.0214, 0.0455] | +1.20 cm |
| Speech vs Drum | -0.0003 m | [-0.0337, 0.0332] | -0.03 cm |
| Pink Noise vs Drum | -0.0251 m | [-0.0585, 0.0083] | -2.51 cm |

### 5.7 Random Effects Variance

```
 Groups        Name                Std.Dev. Corr  
 participantId (Intercept)         0.044491       
               stimulusDisparity_m 0.190145 -0.556
 Residual                          0.164122       
```

## 6. Localization Accuracy (Gamma GLMM)

### 6.1 Model Specification

| Component | Specification |
|-----------|---------------|
| Response | Unsigned error (cm) - always positive |
| Fixed effects | Disparity (centered), Sound Type, Trial (centered), Azimuth, Elevation |
| Random effects | Participant intercept |
| Family | Gamma (log link) |
| Rationale | Right-skewed positive data; multiplicative effects |

### 6.2 Model Comparison: With vs Without Spatial Factors

```
Data: analysis_df
Models:
model_accuracy_reduced: participantError_cm ~ stimDisparity_c + soundType + trialSequence_c + (1 | participantId)
model_accuracy: participantError_cm ~ stimDisparity_c + soundType + trialSequence_c + azimuthSector + elevationCategory + (1 | participantId)
                       npar    AIC    BIC logLik -2*log(L)  Chisq Df Pr(>Chisq)
model_accuracy_reduced    8 5822.1 5859.0  -2903    5806.1                     
model_accuracy           13 5830.0 5889.9  -2902    5804.0 2.1239  5     0.8318
```

**Result:** Spatial factors do not improve model fit (Chi-sq = 2.12, df = 5, p = 0.832)

### 6.3 Final Model Summary

```
Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: Gamma  ( log )
Formula: participantError_cm ~ stimDisparity_c + soundType + trialSequence_c +  
    azimuthSector + elevationCategory + (1 | participantId)
   Data: analysis_df
Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))

      AIC       BIC    logLik -2*log(L)  df.resid 
   5830.0    5889.9   -2902.0    5804.0       731 

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.7642 -0.6447 -0.1273  0.4073  7.8058 

Random effects:
 Groups        Name        Variance Std.Dev.
 participantId (Intercept) 0.03805  0.1951  
 Residual                  0.29678  0.5448  
Number of obs: 744, groups:  participantId, 31

Fixed effects:
                         Estimate Std. Error t value Pr(>|z|)    
(Intercept)             3.2305500  0.0035139 919.361  < 2e-16 ***
stimDisparity_c         1.1615783  0.0035435 327.805  < 2e-16 ***
soundTypeFlute         -0.0620169  0.0036190 -17.137  < 2e-16 ***
soundTypeSpeech        -0.0648140  0.0036109 -17.949  < 2e-16 ***
soundTypePink Noise    -0.1007196  0.0035806 -28.129  < 2e-16 ***
trialSequence_c        -0.0092858  0.0020804  -4.463 8.07e-06 ***
azimuthSectorRight     -0.0134232  0.0035173  -3.816 0.000135 ***
azimuthSectorBack       0.0253262  0.0035354   7.164 7.85e-13 ***
azimuthSectorLeft      -0.0004893  0.0036236  -0.135 0.892586    
elevationCategoryLevel  0.0470088  0.0035216  13.349  < 2e-16 ***
elevationCategoryAbove  0.0568219  0.0034972  16.248  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) stmDs_ sndTyF sndTyS sndTPN trlSq_ azmtSR azmtSB azmtSL elvtCL
stmDsprty_c  0.002                                                               
soundTypFlt  0.005 -0.003                                                        
sondTypSpch  0.005 -0.003  0.243                                                 
sndTypPnkNs -0.007  0.003  0.010  0.012                                          
trialSqnc_c -0.003  0.000  0.006  0.005 -0.010                                   
azmthSctrRg  0.004 -0.002 -0.006 -0.007  0.002 -0.004                            
azmthSctrBc  0.004 -0.002 -0.006 -0.008 -0.004 -0.002  0.005                     
azmthSctrLf -0.007  0.003  0.008  0.009  0.237  0.001  0.004 -0.002              
elvtnCtgryL  0.003 -0.002 -0.006 -0.009 -0.001  0.003 -0.002  0.002  0.000       
elvtnCtgryA -0.006  0.002  0.005  0.007 -0.014 -0.006  0.001 -0.003 -0.005  0.003
optimizer (bobyqa) convergence code: 0 (OK)
Model failed to converge with max|grad| = 0.0249302 (tol = 0.002, component 1)

```

### 6.4 Variance Explained (R-squared)

| Type | R² |
|------|----|
| Marginal R² (fixed effects) | 0.096 |
| Conditional R² (fixed + random) | 0.199 |

### 6.5 Model Fit: Fitted vs Observed Correlation

r(fitted, observed) = 0.509

### 6.6 Effect Sizes (Multiplicative, Log Link)

| Effect | Per-unit change | Interpretation |
|--------|-----------------|----------------|
| Disparity | +10cm → 12.3% | Larger disparity increases error |
| Trial sequence | per trial → -0.92% | Learning effect |
| Flute vs Drum | -6.0% | Sound type effect |
| Speech vs Drum | -6.3% | Sound type effect |
| Pink Noise vs Drum | -9.6% | Sound type effect |

## 7. Sound Type Effects (Post-Hoc Comparisons)

### 7.1 Estimated Marginal Means (Signed Error)

```
Registered S3 methods overwritten by 'broom':
  method        from 
  nobs.fitdistr MuMIn
  nobs.multinom MuMIn
 soundType  emmean     SE   df lower.CL upper.CL
 Drum        0.140 0.0162 86.1   0.1078    0.172
 Flute       0.152 0.0162 86.4   0.1198    0.184
 Speech      0.140 0.0162 86.2   0.1076    0.172
 Pink Noise  0.115 0.0162 86.2   0.0827    0.147

Degrees-of-freedom method: kenward-roger 
Confidence level used: 0.95 
```

### 7.2 Pairwise Contrasts (Signed Error, Tukey-adjusted)

```
 contrast            estimate     SE  df t.ratio p.value
 Drum - Flute        -0.01204 0.0171 688  -0.705  0.8950
 Drum - Speech        0.00025 0.0171 689   0.015  1.0000
 Drum - Pink Noise    0.02507 0.0171 689   1.470  0.4564
 Flute - Speech       0.01229 0.0171 695   0.719  0.8894
 Flute - Pink Noise   0.03711 0.0171 688   2.174  0.1315
 Speech - Pink Noise  0.02482 0.0171 693   1.453  0.4665

Degrees-of-freedom method: kenward-roger 
P value adjustment: tukey method for comparing a family of 4 estimates 
```

### 7.3 Estimated Marginal Means (Accuracy, response scale)

```
 soundType  response    SE  df asymp.LCL asymp.UCL
 Drum           26.3 0.110 Inf      26.0      26.5
 Flute          24.7 0.137 Inf      24.4      24.9
 Speech         24.6 0.136 Inf      24.3      24.9
 Pink Noise     23.7 0.133 Inf      23.5      24.0

Results are averaged over the levels of: azimuthSector, elevationCategory 
Confidence level used: 0.95 
Intervals are back-transformed from the log scale 
```

### 7.4 Pairwise Contrasts (Accuracy, Tukey-adjusted)

```
 contrast            ratio      SE  df null z.ratio p.value
 Drum / Flute         1.06 0.00385 Inf    1  17.137  <.0001
 Drum / Speech        1.07 0.00385 Inf    1  17.949  <.0001
 Drum / Pink Noise    1.11 0.00396 Inf    1  28.129  <.0001
 Flute / Speech       1.00 0.00446 Inf    1   0.629  0.9229
 Flute / Pink Noise   1.04 0.00527 Inf    1   7.639  <.0001
 Speech / Pink Noise  1.04 0.00524 Inf    1   7.102  <.0001

Results are averaged over the levels of: azimuthSector, elevationCategory 
P value adjustment: tukey method for comparing a family of 4 estimates 
Tests are performed on the log scale 
```

## 8. Movement Effects (6DoF Exploration)

### 8.1 Movement Decomposition

To disentangle individual strategies from trial-level effects:
- **Between-subject (mean_rate)**: participant's average movement rate across all trials
- **Within-subject (rate_dev)**: trial-level deviation from participant's mean

### 8.2 Model Summary

```
Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 Family: Gamma  ( log )
Formula: participantError_m ~ stimulusDisparity_m + soundType + trialSequenceNum +  
    mean_rate + rate_dev + (1 | participantId)
   Data: analysis_df
Control: glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))

      AIC       BIC    logLik -2*log(L)  df.resid 
  -1041.4    -995.3     530.7   -1061.4       734 

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-1.7561 -0.6386 -0.1241  0.4289  8.1912 

Random effects:
 Groups        Name        Variance Std.Dev.
 participantId (Intercept) 0.02111  0.1453  
 Residual                  0.30155  0.5491  
Number of obs: 744, groups:  participantId, 31

Fixed effects:
                     Estimate Std. Error   t value Pr(>|z|)    
(Intercept)         -0.585319   0.004052  -144.442  < 2e-16 ***
stimulusDisparity_m  1.135869   0.004077   278.603  < 2e-16 ***
soundTypeFlute      -0.080180   0.003998   -20.054  < 2e-16 ***
soundTypeSpeech     -0.070305   0.004036   -17.420  < 2e-16 ***
soundTypePink Noise -0.105581   0.004052   -26.054  < 2e-16 ***
trialSequenceNum    -0.010155   0.001690    -6.008 1.88e-09 ***
mean_rate           -4.615982   0.004207 -1097.099  < 2e-16 ***
rate_dev             0.336911   0.004205    80.113  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) stmlD_ sndTyF sndTyS sndTPN trlSqN men_rt
stmlsDsprt_ -0.002                                          
soundTypFlt -0.004  0.000                                   
sondTypSpch  0.001 -0.001  0.008                            
sndTypPnkNs  0.000 -0.001  0.004  0.001                     
trialSqncNm -0.049 -0.020 -0.010 -0.012 -0.015              
mean_rate   -0.002  0.000 -0.001  0.000  0.001 -0.011       
rate_dev    -0.001  0.001 -0.001  0.001  0.001  0.001 -0.250
optimizer (bobyqa) convergence code: 0 (OK)
Model failed to converge with max|grad| = 0.0194302 (tol = 0.002, component 1)

```

### 8.3 Variance Explained

Marginal R² = 0.180 | Conditional R² = 0.233

### 8.4 Movement Effect Interpretation

| Effect | Beta | Per 0.1 m/s | Interpretation |
|--------|------|-------------|----------------|
| Between-subject (mean_rate) | -4.6160 | -37.0% | More movement = lower error |
| Within-subject (rate_dev) | 0.3369 | 3.4% | Extra movement = difficulty marker |

## 9. Spatial Configuration Effects

### 9.1 Model Summary

```
Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
Formula: signedError_m ~ stimulusDisparity_m + azimuthSector + elevationCategory +      soundType + (1 | participantId)
   Data: analysis_df

REML criterion at convergence: -455.4

Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-4.5383 -0.4895 -0.0360  0.4674  5.2482 

Random effects:
 Groups        Name        Variance Std.Dev.
 participantId (Intercept) 0.003715 0.06095 
 Residual                  0.027825 0.16681 
Number of obs: 744, groups:  participantId, 31

Fixed effects:
                         Estimate Std. Error         df t value Pr(>|t|)    
(Intercept)            -1.972e-03  2.706e-02  4.054e+02  -0.073    0.942    
stimulusDisparity_m     3.611e-01  4.168e-02  7.094e+02   8.664   <2e-16 ***
azimuthSectorRight     -3.818e-03  1.773e-02  7.180e+02  -0.215    0.830    
azimuthSectorBack       2.105e-02  1.817e-02  7.265e+02   1.159    0.247    
azimuthSectorLeft       3.922e-03  1.831e-02  7.213e+02   0.214    0.830    
elevationCategoryLevel -3.158e-04  1.566e-02  7.329e+02  -0.020    0.984    
elevationCategoryAbove  1.795e-03  1.659e-02  7.216e+02   0.108    0.914    
soundTypeFlute          1.341e-02  1.737e-02  7.041e+02   0.772    0.441    
soundTypeSpeech        -2.523e-03  1.743e-02  7.043e+02  -0.145    0.885    
soundTypePink Noise    -2.651e-02  1.732e-02  7.040e+02  -1.531    0.126    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Correlation of Fixed Effects:
            (Intr) stmlD_ azmtSR azmtSB azmtSL elvtCL elvtCA sndTyF sndTyS
stmlsDsprt_ -0.587                                                        
azmthSctrRg -0.304 -0.044                                                 
azmthSctrBc -0.299 -0.019  0.526                                          
azmthSctrLf -0.269 -0.076  0.519  0.523                                   
elvtnCtgryL -0.287  0.010 -0.062 -0.153 -0.087                            
elvtnCtgryA -0.314  0.073 -0.100 -0.097 -0.116  0.543                     
soundTypFlt -0.372  0.041  0.027  0.050  0.074  0.003  0.015              
sondTypSpch -0.348  0.029  0.052  0.025 -0.022 -0.032  0.055  0.494       
sndTypPnkNs -0.334  0.003  0.008  0.017  0.014 -0.004  0.033  0.499  0.499
```

### 9.2 Post-Hoc: Azimuth Sectors

```
 contrast      estimate     SE  df t.ratio p.value
 Front - Right  0.00382 0.0177 718   0.215  0.9965
 Front - Back  -0.02105 0.0182 726  -1.157  0.6541
 Front - Left  -0.00392 0.0183 721  -0.214  0.9965
 Right - Back  -0.02487 0.0175 721  -1.421  0.4868
 Right - Left  -0.00774 0.0177 718  -0.437  0.9720
 Back - Left    0.01713 0.0178 720   0.961  0.7717

Results are averaged over the levels of: elevationCategory, soundType 
Degrees-of-freedom method: kenward-roger 
P value adjustment: tukey method for comparing a family of 4 estimates 
```

### 9.3 Post-Hoc: Elevation Categories

```
 contrast       estimate     SE  df t.ratio p.value
 Below - Level  0.000316 0.0157 733   0.020  0.9998
 Below - Above -0.001795 0.0166 722  -0.108  0.9936
 Level - Above -0.002111 0.0155 725  -0.137  0.9898

Results are averaged over the levels of: azimuthSector, soundType 
Degrees-of-freedom method: kenward-roger 
P value adjustment: tukey method for comparing a family of 3 estimates 
```

## 10. Individual Differences in Ventriloquist Susceptibility

### 10.1 Random Slope Variance Summary

| Statistic | Value |
|-----------|-------|
| Population mean slope | 34.7% |
| Min individual slope | 10.7% |
| Max individual slope | 72.6% |
| SD of individual slopes | 14.9% |
| Range | 61.9 percentage points |

### 10.2 Individual Slopes (BLUPs)

```
Individual ventriloquist susceptibility (% of disparity captured by vision):

 participantId     intercept     slope bias_pct
          P006 -0.0107163668 0.7255882 72.55882
          P007  0.0048607576 0.5850575 58.50575
          P019 -0.0142024806 0.5487607 54.87607
          P015 -0.0062424543 0.5429139 54.29139
          P016 -0.0070811494 0.5341855 53.41855
          P018  0.0088608960 0.4763659 47.63659
          P017  0.0153524194 0.4681239 46.81239
          P028  0.0353053366 0.4519551 45.19551
          P014  0.0129274203 0.4497634 44.97634
          P003 -0.0189221036 0.4291997 42.91997
          P021  0.0017937138 0.4158480 41.58480
          P023 -0.0320942103 0.4045458 40.45458
          P004 -0.0005075882 0.3976172 39.76172
          P030  0.0178152766 0.3681019 36.81019
          P010  0.0092705008 0.3310016 33.10016
          P012  0.0278268912 0.3156815 31.56815
          P031 -0.0245984581 0.2798508 27.98508
          P029  0.0295002309 0.2688603 26.88603
          P001  0.0465277056 0.2687051 26.87051
          P008  0.0141982781 0.2615326 26.15326
          P013  0.0293300450 0.2554847 25.54847
          P027  0.0338055761 0.2488882 24.88882
          P024  0.0166824061 0.2437233 24.37233
          P009  0.0019425081 0.2358394 23.58394
          P020  0.0367790931 0.2269842 22.69842
          P002  0.0254872294 0.2118148 21.18148
          P022  0.0063862704 0.1827425 18.27425
          P026 -0.0131700374 0.1826946 18.26946
          P011 -0.0042559632 0.1698695 16.98695
          P025  0.0133777198 0.1559240 15.59240
          P005  0.0126503567 0.1068654 10.68654
```

## 11. Summary of Key Findings

### Primary Outcome: Ventriloquist Effect

| Finding | Statistic | 95% CI |
|---------|-----------|--------|
| Visual capture | 34.7% of disparity | [24.2%, 45.1%] |
| Effect per 10cm disparity | 3.47 cm toward flash | [2.42, 4.51] |
| Individual range | 10.7% to 72.6% | — |

### Secondary Outcomes

| Outcome | Effect | p-value |
|---------|--------|--------|
| Sound type (overall) | See post-hoc | — |
| Learning effect | -0.92% per trial | < .001 |
| Movement benefit (between) | 37.0% per 0.1 m/s | < .05 |

### Model Fit Statistics

| Model | Marginal R² | Conditional R² | AIC |
|-------|-------------|----------------|-----|
| Signed error (LMM) | 0.082 | 0.211 | -479.6 |
| Accuracy (Gamma GLMM) | 0.096 | 0.199 | 5830.0 |
| Movement (Gamma GLMM) | 0.180 | 0.233 | -1041.4 |

---

*Report generated by statistical_analysis.R*
