[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21907408-blue.svg)](https://doi.org/10.5281/zenodo.21907408) [![R](https://img.shields.io/badge/R-4.4-blue.svg)]() [![tidyverse](https://img.shields.io/badge/tidyverse-2.0-blue.svg)]() [![lme4](https://img.shields.io/badge/lme4-1.1-blue.svg)]() [![lmerTest](https://img.shields.io/badge/lmerTest-3.1-blue.svg)]() [![glmmTMB](https://img.shields.io/badge/glmmTMB-1.9-blue.svg)]() [![emmeans](https://img.shields.io/badge/emmeans-2.0-blue.svg)]() [![DHARMa](https://img.shields.io/badge/DHARMa-0.4-blue.svg)]() [![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

# 6DoF VR Ventriloquist Effect - Supplementary Materials

Supplementary materials for the paper:

***"Ventriloquist Effect in 6-Degrees-of-Freedom Virtual Reality with Object-Based Spatial Audio"***  
Bartłomiej Mróz · Paweł Perkowski · *IEEE Transactions on Visualization and Computer Graphics*, 2026 [major revision under review]

> The repository state corresponding to the originally submitted manuscript is preserved under the tag `tvcg-submission-v1`. Renamed from `ventriloquist-effect-in-6dof-vr-with-hoa` during the revision.

This repository contains:
- Raw experimental data (localization responses + head trajectories)
- The reproducible statistical analysis pipeline
- Machine-readable analysis outputs, including the LaTeX value macros used by the manuscript

For methodology, interpretation, and results discussion, please refer to the main paper.

## Repository Structure

```
.
├── data/
│   ├── LocalizationTestData_*.json   # Trial-level localization responses (31 participants)
│   └── OculusTrajectoryData_*.csv    # Head tracking trajectories (~9 Hz)
├── rdocs/
│   └── statistical_analysis.R        # Original analysis script (as submitted; see tag)
├── rscripts/revision/                # Revision analysis pipeline (current)
│   ├── build_analysis_df_revision.R  # Dataset construction incl. trajectory-derived measures
│   ├── refit_*.R                     # glmmTMB re-estimation of all Gamma GLMMs
│   ├── reviewer1_*.R, disparity_*.R  # Interaction, trial-number, nonlinearity, breakpoint analyses
│   ├── angular_* / verify_*.R        # Experienced angular disparity; independent verifications
│   ├── regenerate_tolerance_lookup.R # Developer tolerance table (paper Table I)
│   └── generate_revision_values.R    # Writes reported statistics into the manuscript as LaTeX macros
├── results/
│   └── analysis_report.md            # Original analysis output (as submitted)
├── results_revision/                 # Machine-readable outputs of the revision pipeline
└── README.md
```

## Data

**Participants**: 31  
**Trials per participant**: 24  
**Total observations**: 744

**Experimental variables**:
- Sound type: 4 levels (Drum, Flute, Speech, Pink Noise)
- Audio-visual disparity: 15–70 cm (quasi-continuous, 56 distinct values)
- Spatial configuration: Azimuth and elevation of visual cue
- 6DoF head movement: Full exploratory movement allowed

**Key measures**:
- `participantError_m`: Unsigned localization error (meters)
- `signedError_m`: Bias toward (+) or away from (−) visual cue
- `ventriloquistBias`: Signed error as proportion of disparity
- Trajectory-derived (revision): experienced angular disparity, head-to-source distance, cue eccentricity, movement rate

## Reproducing the Analysis

### Prerequisites

```r
install.packages(c("jsonlite", "tidyverse", "lme4", "lmerTest", "glmmTMB", "emmeans", "DHARMa", "mgcv"))
```

### Run the Pipeline

```bash
Rscript rscripts/revision/build_analysis_df_revision.R   # builds results_revision/analysis_df_revision.rds
Rscript rscripts/revision/generate_revision_values.R     # regenerates the manuscript's reported values
```

Each script in `rscripts/revision/` is self-contained and documents its inputs and outputs in a header comment; the `verify_*.R` scripts independently re-estimate the headline models. The per-sample trajectory intermediate (`trajectory_samples_angular.rds`, ~17 MB) is not tracked; `build_analysis_df_revision.R` regenerates it from the raw CSVs.

**Note on estimation**: the Gamma GLMMs are fitted with `glmmTMB`; the originally reported `lme4::glmer` Gamma fits did not converge reliably (degenerate variance-covariance matrices) and are retained only under the `tvcg-submission-v1` tag for the record.

## Citation

If you use this data or analysis code, please cite the main paper:

```
[BibTeX citation will be added upon publication]
```

## License

This work is licensed under a [Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: https://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Contact

Bartłomiej Mróz · bartlomiej.mroz@pg.edu.pl · Department of Multimedia Systems, Gdańsk University of Technology · [bmroz.eu](https://bmroz.eu)
