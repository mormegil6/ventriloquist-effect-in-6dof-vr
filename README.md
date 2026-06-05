[![R](https://img.shields.io/badge/R-4.4-blue.svg)]() [![tidyverse](https://img.shields.io/badge/tidyverse-2.0-blue.svg)]() [![lme4](https://img.shields.io/badge/lme4-1.1-blue.svg)]() [![lmerTest](https://img.shields.io/badge/lmerTest-3.1-blue.svg)]() [![emmeans](https://img.shields.io/badge/emmeans-2.0-blue.svg)]() [![MuMIn](https://img.shields.io/badge/MuMIn-1.48-blue.svg)]() [![DHARMa](https://img.shields.io/badge/DHARMa-0.4-blue.svg)]() [![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

# 6DoF VR Ventriloquist Effect - Supplementary Materials

Supplementary materials for the paper:

***"Ventriloquist Effect in 6-Degrees-of-Freedom Virtual Reality with Higher-Order Ambisonics"***  
Bartłomiej Mróz · *IEEE Transactions on Visualization and Computer Graphics*, 2026 [under review]

This repository contains:
- Raw experimental data (localization responses + head trajectories)
- Reproducible statistical analysis script
- Analysis outputs

For methodology, interpretation, and results discussion, please refer to the main paper.

## Repository Structure

```
.
├── data/
│   ├── LocalizationTestData_*.json  # Trial-level localization responses (31 participants)
│   └── OculusTrajectoryData_*.csv   # Head tracking trajectories
├── rdocs/
│   └── statistical_analysis.R       # Reproducible analysis script
├── results/
│   └── analysis_report.md           # Complete analysis output
└── README.md
```

## Data

**Participants**: 31  
**Trials per participant**: 24  
**Total observations**: 744

**Experimental variables**:
- Sound type: 4 levels (Drum, Flute, Speech, Pink Noise)
- Audio-visual disparity: 15–70 cm (continuous)
- Spatial configuration: Azimuth and elevation of visual cue
- 6DoF head movement: Full exploratory movement allowed

**Key measures**:
- `participantError_m`: Unsigned localization error (meters)
- `signedError_m`: Bias toward (+) or away from (−) visual cue
- `ventriloquistBias`: Signed error as proportion of disparity
- `total_path_length`: Head movement during trial (meters)

## Reproducing the Analysis

### Prerequisites

```r
install.packages(c("jsonlite", "tidyverse", "lme4", "lmerTest", "emmeans", "MuMIn", "DHARMa"))
```

### Run the Analysis

```bash
Rscript rdocs/statistical_analysis.R
```

The script performs the complete statistical analysis pipeline:
1. Data import and preprocessing
2. Ventriloquist effect analysis (signed error ~ disparity)
3. Localization accuracy analysis (Gamma GLMM)
4. Sound type comparisons
5. Movement effect analysis
6. Spatial configuration effects

**Output**: Console text with all analysis results (also available in `results/analysis_report.md`)

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
