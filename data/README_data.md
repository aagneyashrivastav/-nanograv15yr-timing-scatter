# Data

## NANOGrav 15-Year Dataset (not included in this repository)

The NANOGrav narrowband residual and timing-model parameter files are
publicly available on Zenodo and must be downloaded separately before
running the analysis scripts.

**Download link:**
https://doi.org/10.5281/zenodo.16051178

### Files needed

After downloading, place the files in the following locations:

```
data/
├── residuals/
│   └── <pulsar>_nb.full.res    (narrowband post-fit residual files)
└── par_files/
    └── <pulsar>.par            (timing-model parameter files)
```

Only the narrowband residual files (`*_nb.full.res`) and their
corresponding `.par` files are used in this analysis. The wideband
files are not required.

## Bayesian Red-Noise Amplitudes (included)

`bayesian_rednoise.csv` contains the published Bayesian posterior
medians and 68% credible intervals for log10(A_RN) from Table 2 of:

> Agazie et al. (The NANOGrav Collaboration), "The NANOGrav 15 yr
> Data Set: Detector Characterization and Noise Budget,"
> The Astrophysical Journal Letters (2023).
> DOI: 10.3847/2041-8213/acda88

Columns:
- `Pulsar`       — pulsar name
- `log10_A_RN`   — posterior median of log10(A_RN)
- `sigma_plus`   — upper 68% credible interval bound (positive)
- `sigma_minus`  — lower 68% credible interval bound (negative)
