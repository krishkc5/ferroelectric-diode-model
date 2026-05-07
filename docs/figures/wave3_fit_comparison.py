"""Wave-3 fit comparison: three configurations side-by-side for one representative sweep.

Compares:
  (1) Wave-3 default (no Tung; Tung sigmas hard-zero) — cost 909, mean 0.702
  (2) Wave-3 + Tung only                                 — cost 907, mean 0.703
  (3) Wave-3 + Tung + Tunneling + TAT                    — cost 3113, mean 1.235

The point of this figure is to show that adding more DOF doesn't rescue the fit;
the residual gap is structural and needs a different physical addition (constrained
TAT trap distribution, or a Tung-only model with proper bounds).

Sweep B64-T51 is used for all three panels.
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict

import matplotlib.pyplot as plt
import numpy as np

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = OUT_DIR
TENNM_DIR = os.path.join(os.path.dirname(os.path.dirname(OUT_DIR)), "10nm_figures")

CONFIGS = [
    ("Wave-3 default\n(no Tung; TE+PF+BG)", "wave3_refit_20260507.csv", "C0"),
    ("Wave-3 + Tung only\n(TE inhomogeneity, no tunneling/TAT)", "wave3_refit_tung_v2_20260507.csv", "C1"),
    ("Wave-3 + Tung + Tunneling + TAT\n(fit failed to converge)", "wave3_refit_tung_tat_20260507.csv", "C3"),
]
TARGET_SWEEP = "1st2KB_B64-T51"


def load_sweep(csv_path, target):
    sweeps = defaultdict(lambda: defaultdict(list))
    with open(csv_path) as h:
        for row in csv.DictReader(h):
            name = row["sweep_name"]
            for k, v in row.items():
                if k == "sweep_name":
                    continue
                sweeps[name][k].append(float(v))
    s = sweeps[target]
    return {k: np.array(v) for k, v in s.items()}


def per_sweep_mean_log_error(j_exp, j_model):
    mask = (np.abs(j_exp) > 0) & (np.abs(j_model) > 0)
    if mask.sum() == 0:
        return float("nan")
    return float(np.mean(np.abs(np.log10(np.abs(j_model[mask])) - np.log10(np.abs(j_exp[mask])))))


def main():
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0), sharey=True, sharex=True)

    floor = 1e-15
    for ax, (label, csv_name, color) in zip(axes, CONFIGS):
        path = os.path.join(TENNM_DIR, csv_name)
        if not os.path.exists(path):
            ax.text(0.5, 0.5, f"missing:\n{csv_name}", ha="center", va="center", transform=ax.transAxes)
            continue
        s = load_sweep(path, TARGET_SWEEP)
        v = s["voltage_v"]
        j_exp = s["experimental_current_a"]
        j_total = s["model_total_a"]
        err = per_sweep_mean_log_error(j_exp, j_total)
        ax.semilogy(v, np.maximum(np.abs(j_exp), floor), "k.", ms=3.5, label="Measured", zorder=10, alpha=0.85)
        ax.semilogy(v, np.maximum(np.abs(j_total), floor), color=color, lw=2.0, label="Model")
        ax.axvline(0, color="gray", lw=0.5)
        ax.set_xlabel("Applied voltage  [V]")
        ax.set_title(f"{label}\n|Δlog10 I| = {err:.3f}", fontsize=10)
        ax.set_ylim(1e-15, 1e-7)
        ax.grid(True, alpha=0.3, which="both")
        ax.legend(fontsize=9, loc="lower right")

    axes[0].set_ylabel("|I|  [A]")
    fig.suptitle(
        f"Wave-3 fit comparison on sweep {TARGET_SWEEP} — Tung alone gives no improvement; "
        f"enabling tunneling/TAT makes it worse (fit fails to converge)",
        fontsize=11,
    )
    fig.tight_layout()

    out_path = os.path.join(FIG_DIR, "wave3_fit_comparison.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[wrote] {out_path}")


if __name__ == "__main__":
    main()
