# =============================================================================
# 01_extract_metrics.py
#
# Purpose:
#   Extract narrowband timing residuals from the NANOGrav 15-year dataset
#   for the 23 target millisecond pulsars and compute three empirical
#   timing-scatter metrics: RRMS, WRMS, and zeta.
#
#   EFAC and EQUAD uncertainty parameters are also extracted from .par files
#   and saved alongside the metrics.
#
# Inputs:
#   data/residuals/<pulsar>_nb.full.res   (NANOGrav narrowband residual files)
#   data/par_files/<pulsar>.par           (NANOGrav timing-model parameter files)
#
# Output:
#   outputs/pulsar_metrics.csv
#
# Note:
#   The 23 target pulsars are those with significant Bayesian red-noise
#   detections in Agazie et al. (2023), Table 2. The published Bayesian
#   log10(A_RN) values are stored separately in data/bayesian_rednoise.csv
#   and are loaded in script 02.
#
# Author: Aagneya Shrivastav, Chaffey College
# =============================================================================

import numpy as np
import pandas as pd
import glob
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR   = os.path.dirname(SCRIPT_DIR)

# =============================================================================
# TARGET PULSARS
# These 23 pulsars have significant Bayesian red-noise detections reported
# in Agazie et al. (2023), Table 2 of the NANOGrav noise budget paper.
# =============================================================================

TARGET_PULSARS = [
    "B1855+09", "B1937+21", "B1953+29",
    "J0030+0451", "J0437-4715", "J0610-2100", "J0613-0200",
    "J1012+5307", "J1600-3053", "J1614-2230", "J1643-1224",
    "J1705-1903", "J1713+0747", "J1744-1134", "J1745+1017",
    "J1747-4036", "J1802-2124", "J1853+1303", "J1903+0327",
    "J1909-3744", "J1946+3417", "J2145-0750", "J2234+0611"
]

# =============================================================================
# STEP 1 — EXTRACT RESIDUALS AND COMPUTE METRICS
# =============================================================================

print("\n" + "="*60)
print("STEP 1 — Extracting residuals and computing metrics")
print("="*60 + "\n")

results = []

residual_files = glob.glob(
    os.path.join(ROOT_DIR, "data", "residuals", "*_nb.full.res")
)

print(f"Found {len(residual_files)} residual files.\n")

for f in residual_files:

    # Extract pulsar name from filename
    pulsar_name = os.path.basename(f).split('_')[0]

    # Skip pulsars not in the target list
    if pulsar_name not in TARGET_PULSARS:
        continue

    try:
        # Load residual file
        # Column 0 = MJD
        # Column 3 = Post-fit timing residual (microseconds)
        # Column 4 = TOA measurement uncertainty (microseconds)
        raw = np.genfromtxt(f, usecols=(0, 3, 4))

        mjd          = raw[:, 0]
        residuals    = raw[:, 1]
        uncertainties = raw[:, 2]

        # Remove invalid rows
        valid = (
            np.isfinite(residuals) &
            np.isfinite(uncertainties) &
            (uncertainties > 0)
        )
        residuals     = residuals[valid]
        uncertainties = uncertainties[valid]

        if len(residuals) < 5:
            print(f"  SKIPPED {pulsar_name} — too few valid observations")
            continue

        # ------------------------------------------------------------------
        # RRMS: residual root-mean-square
        # Equation (1) from paper
        # ------------------------------------------------------------------
        rrms = np.sqrt(np.mean(residuals**2))

        # ------------------------------------------------------------------
        # WRMS: inverse-variance weighted RMS
        # Equation (2) from paper
        # ------------------------------------------------------------------
        weights = 1.0 / uncertainties**2
        wrms    = np.sqrt(np.sum(weights * residuals**2) / np.sum(weights))

        # ------------------------------------------------------------------
        # Zeta: normalized excess-scatter statistic
        # Equation (3) from paper
        # ------------------------------------------------------------------
        median_sigma = np.median(uncertainties)
        zeta         = wrms / median_sigma

        results.append({
            "Pulsar":               pulsar_name,
            "RRMS_us":              rrms,
            "WRMS_us":              wrms,
            "Median_Sigma_us":      median_sigma,
            "Zeta":                 zeta,
            "N_Observations":       int(np.sum(valid))
        })

        print(f"  Computed: {pulsar_name}")

    except Exception as e:
        print(f"  ERROR — {pulsar_name}: {e}")

metrics_df = pd.DataFrame(results)
print(f"\n{len(metrics_df)} pulsars processed.\n")

# =============================================================================
# STEP 2 — EXTRACT EFAC AND EQUAD FROM PAR FILES
# =============================================================================

print("="*60)
print("STEP 2 — Extracting EFAC and EQUAD from .par_files")
print("="*60 + "\n")

noise_params = []

par_files = glob.glob(
    os.path.join(ROOT_DIR, "data", "par_files", "*.par")
)

# Exclude NoRedNoise variants
par_files = [f for f in par_files if "NoRedNoise" not in f]

print(f"Found {len(par_files)} usable .par_files.\n")

for f in par_files:

    pulsar_name = os.path.basename(f).split('.')[0].split('_')[0]

    if pulsar_name not in TARGET_PULSARS:
        continue

    efacs  = []
    equads = []

    with open(f, 'r') as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == "EFAC":
                efacs.append(float(parts[-1]))
            if parts[0] == "EQUAD":
                equads.append(float(parts[-1]))

    noise_params.append({
        "Pulsar":    pulsar_name,
        # Mean EFAC across all backends
        "EFAC_mean": np.mean(efacs)  if efacs  else np.nan,
        # Max EQUAD: worst-case instrumental noise floor across backends
        "EQUAD_max": np.max(equads)  if equads else np.nan,
        "N_Backends": len(efacs)
    })

    print(f"  Extracted: {pulsar_name}  "
          f"({len(efacs)} backends)")

noise_df = pd.DataFrame(noise_params).drop_duplicates(
    subset="Pulsar", keep="first"
)

# =============================================================================
# STEP 3 — MERGE AND SAVE
# =============================================================================

print("\n" + "="*60)
print("STEP 3 — Merging and saving")
print("="*60 + "\n")

output_df = pd.merge(metrics_df, noise_df, on="Pulsar", how="left")
output_df = output_df.sort_values("Pulsar").reset_index(drop=True)

os.makedirs(os.path.join(ROOT_DIR, "outputs"), exist_ok=True)
output_df.to_csv(
    os.path.join(ROOT_DIR, "outputs", "pulsar_metrics.csv"), index=False
)

print("Saved: outputs/pulsar_metrics.csv")
print(f"\nColumns: {output_df.columns.tolist()}")
print(f"\nFinal table ({len(output_df)} pulsars):")
print(output_df[["Pulsar", "RRMS_us", "WRMS_us", "Zeta"]].to_string(index=False))
print("\nDone.")
