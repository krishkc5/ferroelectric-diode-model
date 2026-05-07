"""Per-mechanism breakdown of the fitted model from the wave3_refit run.

Reads 10nm_figures/wave3_refit_20260507.csv (pf_region=fe, mean 0.702) and plots
the contributions of each mechanism vs V for one representative sweep, overlaid
on the experimental measurement.

The CSV has these per-mechanism columns: model_pf_a, model_thermionic_a,
model_background_a, plus model_total_a. Tunneling/TAT/SCLC aren't broken out
in the fit-script CSV; we'll mark this gap on the figure.
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(
    os.path.dirname(os.path.dirname(OUT_DIR)),
    "10nm_figures",
    "wave3_refit_20260507.csv",
)


def load_csv(path):
    sweeps = defaultdict(lambda: defaultdict(list))
    with open(path) as h:
        for row in csv.DictReader(h):
            name = row["sweep_name"]
            for k, v in row.items():
                if k == "sweep_name":
                    continue
                sweeps[name][k].append(float(v))
    return {name: {k: np.array(v) for k, v in cols.items()} for name, cols in sweeps.items()}


def main():
    sweeps = load_csv(CSV_PATH)
    # Pick one sweep — the cleanest data range
    target = "1st2KB_B64-T51"
    s = sweeps[target]
    v = s["voltage_v"]
    j_exp = s["experimental_current_a"]
    j_total = s["model_total_a"]
    j_pf = s["model_pf_a"]
    j_te = s["model_thermionic_a"]
    j_bg = s["model_background_a"]

    # Magnitude floor for log scale
    floor = 1e-15

    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.5))

    # Left: model breakdown vs experimental
    ax = axes[0]
    ax.semilogy(v, np.maximum(np.abs(j_exp), floor), "k.", ms=4, label="Measured", zorder=10)
    ax.semilogy(v, np.maximum(np.abs(j_total), floor), "C0-", lw=2.0, label="Model total")
    ax.semilogy(v, np.maximum(np.abs(j_te), floor), "C2-", lw=1.0, alpha=0.85, label="TE")
    ax.semilogy(v, np.maximum(np.abs(j_pf), floor), "C1-", lw=1.0, alpha=0.85, label="PF")
    ax.semilogy(v, np.maximum(np.abs(j_bg), floor), "C3-", lw=1.0, alpha=0.85, label="Background")
    ax.set_ylim(1e-15, 1e-7)
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("Applied voltage  [V]")
    ax.set_ylabel("|I|  [A]")
    ax.set_title(f"{target}: per-mechanism breakdown\n(fitted pf_region=fe; mean |Δlog10 I|=0.706)")
    ax.legend(fontsize=9, loc="lower right")
    ax.grid(True, alpha=0.3, which="both")

    # Right: residual vs V — model − measured in log10 space
    ax = axes[1]
    log_resid = np.log10(np.maximum(np.abs(j_total), floor)) - np.log10(np.maximum(np.abs(j_exp), floor))
    ax.plot(v, log_resid, "C0-", lw=1.5, label=r"$\log_{10}|J_{model}| - \log_{10}|J_{meas}|$")
    ax.axhline(0, color="k", lw=0.6, alpha=0.5)
    ax.axhspan(-0.3, 0.3, color="C2", alpha=0.1, label="±0.3 decade band")
    ax.axvline(0, color="gray", lw=0.5)
    ax.set_xlabel("Applied voltage  [V]")
    ax.set_ylabel(r"$\Delta \log_{10} I$")
    ax.set_title(f"Residual vs V — where does the fit miss?")
    ax.legend(fontsize=9, loc="upper right")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-2.5, 2.5)

    fig.tight_layout()
    out_path = os.path.join(OUT_DIR, "wave3_mechanism_breakdown.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[wrote] {out_path}")

    # Console summary: dominant mechanism by V region
    print("\nDominant mechanism by V (sweep:", target, "):")
    print(f"{'V':>6}  {'|Jexp|':>10}  {'|Jtot|':>10}  {'TE':>10}  {'PF':>10}  {'BG':>10}  dominant")
    for i in range(0, len(v), max(len(v) // 12, 1)):
        comps = {"TE": abs(j_te[i]), "PF": abs(j_pf[i]), "BG": abs(j_bg[i])}
        dom = max(comps, key=comps.get)
        print(
            f"{v[i]:>+6.2f}  {abs(j_exp[i]):>10.2e}  {abs(j_total[i]):>10.2e}  "
            f"{abs(j_te[i]):>10.2e}  {abs(j_pf[i]):>10.2e}  {abs(j_bg[i]):>10.2e}  {dom}"
        )


if __name__ == "__main__":
    main()
