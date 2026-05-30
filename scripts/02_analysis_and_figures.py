# =============================================================================
# 02_analysis_and_figures.py
#
# Purpose:
#   Load the empirical timing-scatter metrics from script 01, merge with
#   published Bayesian red-noise amplitudes from Agazie et al. (2023),
#   compute correlation statistics with bootstrap confidence intervals,
#   and generate Figures 1 and 2 from the paper.
#
# Inputs:
#   outputs/pulsar_metrics.csv          (produced by 01_extract_metrics.py)
#   data/bayesian_rednoise.csv          (Table 2, Agazie et al. 2023)
#
# Outputs:
#   figures/Timing_Scatter_Correlations.png   (Figure 1)
#   figures/Pulsar_Ranking_Comparisons.png    (Figure 2)
#   outputs/correlation_statistics.csv        (all reported statistics)
#
# Author: Aagneya Shrivastav, Chaffey College
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress, spearmanr, pearsonr, rankdata, kendalltau
import os

# =============================================================================
# LOAD DATA
# =============================================================================

metrics_df  = pd.read_csv("outputs/pulsar_metrics.csv")
bayesian_df = pd.read_csv("data/bayesian_rednoise.csv")

# Merge on Pulsar name
df = pd.merge(metrics_df, bayesian_df, on="Pulsar", how="inner")
df = df.sort_values("Pulsar").reset_index(drop=True)

print(f"\nLoaded {len(df)} pulsars.\n")

# Extract arrays
pulsars   = df["Pulsar"].values
rrms      = df["RRMS_us"].values
wrms      = df["WRMS_us"].values
zeta      = df["Zeta"].values
rnamp     = df["log10_A_RN"].values
err_upper = df["sigma_minus"].values   # negative: downward uncertainty
err_lower = df["sigma_plus"].values    # positive: upward uncertainty

# Asymmetric error bars for matplotlib: [downward, upward]
yerr = [np.abs(err_upper), err_lower]

# =============================================================================
# BOOTSTRAP CONFIDENCE INTERVALS
# =============================================================================

def bootstrap_ci(x, y, stat_func, n_boot=10000, ci=95, seed=42):
    """
    Estimate bootstrap confidence interval for a correlation statistic.

    Parameters
    ----------
    x, y      : arrays of equal length
    stat_func : function(x, y) -> scalar statistic
    n_boot    : number of bootstrap iterations
    ci        : confidence level (percent)
    seed      : random seed for reproducibility

    Returns
    -------
    (lower, upper) confidence interval bounds
    """
    rng  = np.random.default_rng(seed)
    n    = len(x)
    boot = [stat_func(x[idx := rng.integers(0, n, n)],
                      y[idx]) for _ in range(n_boot)]
    lo   = np.percentile(boot, (100 - ci) / 2)
    hi   = np.percentile(boot, 100 - (100 - ci) / 2)
    return lo, hi

# =============================================================================
# COMPUTE AND PRINT ALL CORRELATION STATISTICS
# =============================================================================

metrics_dict = {"zeta": zeta, "RRMS": rrms, "WRMS": wrms}

print(f"{'Metric':<6} {'r':>6} {'r²':>6} {'p':>10} "
      f"{'Pearson 95% CI':>22} {'ρ':>6} "
      f"{'Spearman 95% CI':>22} {'τ':>6} "
      f"{'Slope':>8} {'Intercept':>10}")
print("─" * 115)

stats_store  = {}
stats_rows   = []

for name, x in metrics_dict.items():

    slope, intercept, r, p, _ = linregress(x, rnamp)
    r2      = r ** 2
    rho, _  = spearmanr(x, rnamp)
    tau, _  = kendalltau(x, rnamp)

    p_ci = bootstrap_ci(x, rnamp, lambda a, b: pearsonr(a, b)[0])
    s_ci = bootstrap_ci(x, rnamp, lambda a, b: spearmanr(a, b)[0])

    print(f"{name:<6} {r:>6.3f} {r2:>6.3f} {p:>10.4g} "
          f"  [{p_ci[0]:>6.3f}, {p_ci[1]:>6.3f}]  "
          f"{rho:>6.3f}  [{s_ci[0]:>6.3f}, {s_ci[1]:>6.3f}] "
          f"{tau:>6.3f}  {slope:>8.4f}  {intercept:>10.4f}")

    stats_store[name] = dict(
        r=r, r2=r2, p=p, slope=slope, intercept=intercept,
        p_ci=p_ci, s_ci=s_ci, rho=rho, tau=tau
    )
    stats_rows.append({
        "Metric": name, "r": r, "r2": r2, "p": p,
        "Pearson_CI_lo": p_ci[0], "Pearson_CI_hi": p_ci[1],
        "rho": rho, "Spearman_CI_lo": s_ci[0], "Spearman_CI_hi": s_ci[1],
        "tau": tau, "Slope": slope, "Intercept": intercept
    })

# Save statistics table
stats_out = pd.DataFrame(stats_rows)
os.makedirs("outputs", exist_ok=True)
stats_out.to_csv("outputs/correlation_statistics.csv", index=False)
print("\nSaved: outputs/correlation_statistics.csv")

# =============================================================================
# RANKING TABLE
# =============================================================================

rnamp_rank = rankdata(rnamp).astype(int)
zeta_rank  = rankdata(zeta).astype(int)
rrms_rank  = rankdata(rrms).astype(int)
wrms_rank  = rankdata(wrms).astype(int)

print(f"\n\nRANKING TABLE (1 = most stable / smallest scatter)")
print(f"{'Pulsar':<14} {'A_RN Rank':>10} {'ζ Rank':>8} "
      f"{'RRMS Rank':>10} {'WRMS Rank':>10}")
print("─" * 56)

for i in np.argsort(rnamp_rank):
    print(f"{pulsars[i]:<14} {rnamp_rank[i]:>10} {zeta_rank[i]:>8} "
          f"{rrms_rank[i]:>10} {wrms_rank[i]:>10}")

# =============================================================================
# FIGURE 1 — CORRELATION PLOTS
# =============================================================================

plt.rcParams.update({
    'font.size':      14,
    'axes.labelsize': 20,
    'axes.titlesize': 22
})

fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)

plot_configs = [
    {
        'x':      zeta,
        'name':   'zeta',
        'title':  r'$\zeta$ vs $\log_{10}(A_\mathrm{RN})$',
        'xlabel': r'$\zeta$',
        'ylabel': r'$\log_{10}(A_\mathrm{RN})$'
    },
    {
        'x':      rrms,
        'name':   'RRMS',
        'title':  r'$\mathrm{RRMS}$ vs $\log_{10}(A_\mathrm{RN})$',
        'xlabel': r'$\mathrm{RRMS}\ (\mu s)$',
        'ylabel': r'$\log_{10}(A_\mathrm{RN})$'
    },
    {
        'x':      wrms,
        'name':   'WRMS',
        'title':  r'$\mathrm{WRMS}$ vs $\log_{10}(A_\mathrm{RN})$',
        'xlabel': r'$\mathrm{Weighted\ RMS}\ (\mu s)$',
        'ylabel': r'$\log_{10}(A_\mathrm{RN})$'
    },
]

for config, ax in zip(plot_configs, axes):

    x = config['x']
    s = stats_store[config['name']]

    # Scatter with asymmetric error bars
    ax.errorbar(
        x, rnamp, yerr=yerr,
        fmt='o', color='steelblue', ecolor='gray',
        elinewidth=1.2, capsize=3,
        markeredgecolor='black', markersize=7,
        alpha=0.85, zorder=2
    )

    # OLS regression line
    x_line = np.linspace(x.min(), x.max(), 200)
    ax.plot(
        x_line, s['slope'] * x_line + s['intercept'],
        color='crimson', linestyle='--', linewidth=2, zorder=1
    )

    # Statistics annotation box
    txt = (
        f"Pearson $r = {s['r']:.3f}$\n"
        f"95% CI $[{s['p_ci'][0]:.3f},\\ {s['p_ci'][1]:.3f}]$\n"
        f"$R^2 = {s['r2']:.3f}$\n"
        f"$p = {s['p']:.4g}$\n"
        f"Slope $= {s['slope']:.3f}$"
    )
    props = dict(boxstyle='round', facecolor='white',
                 alpha=0.95, edgecolor='black')
    ax.text(0.95, 0.05, txt, transform=ax.transAxes,
            fontsize=13, va='bottom', ha='right', bbox=props)

    ax.set_title(config['title'], pad=15)
    ax.set_xlabel(config['xlabel'])
    ax.set_ylabel(config['ylabel'])
    ax.grid(True, linestyle=':', alpha=0.8)

os.makedirs("figures", exist_ok=True)
plt.savefig('figures/Timing_Scatter_Correlations.png',
            format='png', bbox_inches='tight', dpi=300)
plt.show()
print("\nSaved: figures/Timing_Scatter_Correlations.png")

# =============================================================================
# FIGURE 2 — RANKING COMPARISON PLOTS
# =============================================================================

fig2, axes2 = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)

rank_configs = [
    {
        'rank':   zeta_rank,
        'title':  r'$\zeta$ Rank vs $\log_{10}(A_\mathrm{RN})$ Rank',
        'ylabel': r'Predicted $\zeta$ Rank'
    },
    {
        'rank':   rrms_rank,
        'title':  r'RRMS Rank vs $\log_{10}(A_\mathrm{RN})$ Rank',
        'ylabel': r'Predicted RRMS Rank'
    },
    {
        'rank':   wrms_rank,
        'title':  r'WRMS Rank vs $\log_{10}(A_\mathrm{RN})$ Rank',
        'ylabel': r'Predicted Weighted RMS Rank'
    },
]

for config, ax in zip(rank_configs, axes2):

    mr      = config['rank']
    rho, _  = spearmanr(rnamp_rank, mr)
    tau, _  = kendalltau(rnamp_rank, mr)

    ax.scatter(rnamp_rank, mr,
               color='darkorange', edgecolor='black',
               s=70, alpha=0.85, zorder=2)

    lims = [0, 24]
    ax.plot(lims, lims, 'k--', linewidth=1.5, alpha=0.5, zorder=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    txt   = f"Spearman $\\rho = {rho:.3f}$\nKendall $\\tau = {tau:.3f}$"
    props = dict(boxstyle='round', facecolor='white',
                 alpha=0.95, edgecolor='black')
    ax.text(0.05, 0.95, txt, transform=ax.transAxes,
            fontsize=14, va='top', ha='left', bbox=props)

    ax.set_title(config['title'], pad=15)
    ax.set_xlabel(r'$\log_{10}(A_\mathrm{RN})$ Rank (Benchmark)')
    ax.set_ylabel(config['ylabel'])
    ax.grid(True, linestyle=':', alpha=0.8)

plt.savefig('figures/Pulsar_Ranking_Comparisons.png',
            format='png', bbox_inches='tight', dpi=300)
plt.show()
print("Saved: figures/Pulsar_Ranking_Comparisons.png")
print("\nAll done.")
