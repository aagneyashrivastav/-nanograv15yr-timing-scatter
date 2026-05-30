# Empirical Timing-Scatter Metrics and Red-Noise Correlations
# in the NANOGrav 15-Year Dataset

**Author:** Aagneya Shrivastav, Chaffey College  
**Contact:** aagneyashrivastav@gmail.com

---

## Overview

This repository contains the complete analysis code for the paper
*"Empirical Timing-Scatter Metrics and Red-Noise Correlations in the
NANOGrav 15-Year Dataset"*.

Three empirical timing-scatter statistics — residual RMS (RRMS),
weighted RMS (WRMS), and a normalized excess-scatter statistic (zeta) —
are computed from publicly available NANOGrav narrowband timing residuals
for 23 millisecond pulsars with significant Bayesian red-noise detections.
These are compared against published Bayesian red-noise amplitudes from
Agazie et al. (2023) using Pearson and rank-based correlation analyses
with bootstrap confidence intervals.

---

## Repository Structure

```
nanograv15yr-timing-scatter/
│
├── README.md
│
├── data/
│   ├── README_data.md              <- download instructions for NANOGrav data
│   ├── bayesian_rednoise.csv       <- Table 2 from Agazie et al. (2023)
│   ├── residuals/                  <- place *.full.res files here (not included)
│   └── par_files/                  <- place *.par files here (not included)
│
├── scripts/
│   ├── 01_extract_metrics.py       <- extract RRMS, WRMS, zeta from residuals
│   └── 02_analysis_and_figures.py  <- correlations, bootstrap CIs, figures
│
├── figures/
│   ├── Timing_Scatter_Correlations.png    <- Figure 1 from paper
│   └── Pulsar_Ranking_Comparisons.png     <- Figure 2 from paper
│
└── outputs/
    ├── pulsar_metrics.csv               <- Table II from paper
    └── correlation_statistics.csv       <- all reported r, rho, tau, CI values
```

---

## Data

The NANOGrav 15-year dataset is publicly available at:

> https://doi.org/10.5281/zenodo.16051178

See `data/README_data.md` for full download and setup instructions.

The Bayesian red-noise amplitudes used as the benchmark (`data/bayesian_rednoise.csv`)
are taken directly from Table 2 of the NANOGrav noise budget paper:

> Agazie et al. (2023), DOI: 10.3847/2041-8213/acda88

---

## Requirements

Python 3.9 or later. Install dependencies with:

```bash
pip install numpy pandas matplotlib scipy
```

---

## How to Reproduce the Results

**Step 1:** Download the NANOGrav data (see `data/README_data.md`).

**Step 2:** Extract metrics from residuals:

```bash
python scripts/01_extract_metrics.py
```

This reads the `.res` and `.par` files and writes `outputs/pulsar_metrics.csv`.

**Step 3:** Run the correlation analysis and generate figures:

```bash
python scripts/02_analysis_and_figures.py
```

This reads `outputs/pulsar_metrics.csv` and `data/bayesian_rednoise.csv`,
prints all reported statistics, and saves both figures to `figures/`.

---

## Pre-computed Results

If you do not wish to re-run the full pipeline, pre-computed outputs are
available in `outputs/` and figures in `figures/`. These match exactly
the values reported in the paper.

---

## Citation

If you use this code or data, please cite the associated paper:

> Shrivastav, A. (2026). "Empirical Timing-Scatter Metrics and
> Red-Noise Correlations in the NANOGrav 15-Year Dataset."
> [Journal, when published]

And the NANOGrav datasets used:

> Agazie et al. (2023), "The NANOGrav 15 yr Data Set: Observations
> and Timing of 68 Millisecond Pulsars," ApJL.
> DOI: 10.3847/2041-8213/acda9a

> Agazie et al. (2023), "The NANOGrav 15 yr Data Set: Detector
> Characterization and Noise Budget," ApJL.
> DOI: 10.3847/2041-8213/acda88

> NANOGrav Collaboration (2025), "The NANOGrav 15-Year Data Set,"
> Zenodo. DOI: 10.5281/zenodo.16051178
