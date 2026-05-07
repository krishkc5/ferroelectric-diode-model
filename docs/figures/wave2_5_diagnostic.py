"""Diagnostic plots for the Wave-2.5 screening-saturation work.

Generates two figures:
  wave2_5_screening_before.png  — linear Tsymbal-Kohlstedt formula, blown up at AlScN scales
  wave2_5_barriers_before.png   — top/bottom Schottky barrier vs polarization, with clamp annotated

Usage (from the worktree root):
    uv run --with matplotlib --with scipy python docs/figures/wave2_5_diagnostic.py
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants

# Stack-aware imports — match tests/test_stack_aware_transport.py convention.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
_SRC = os.path.join(_REPO, "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from materials import Materials  # noqa: E402

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def linear_dv_screen(sigma_uc_cm2: np.ndarray, lambda_tf_nm: float, kappa_inf: float) -> np.ndarray:
    """ΔV_screen [V] from the linear Tsymbal-Kohlstedt formula (Eq. 4)."""
    sigma_si = sigma_uc_cm2 * 1e-2  # uC/cm^2 → C/m^2
    lambda_si = lambda_tf_nm * 1e-9
    return sigma_si * lambda_si / (constants.epsilon_0 * kappa_inf)


def figure_screening_before():
    """Plot ΔV_screen vs σ_pol on the linear formula. Mark BTO/BFO/AlScN scales."""
    sigma_grid = np.linspace(0, 150, 400)  # uC/cm^2
    # Use Ti/Al parameters from materials.py — both have screening_len = 0.045 nm
    # and (after Wave 2 edits) electrode-side optical κ defaulting to 1 (Drude tail).
    lambda_tf_nm = 0.045
    kappa_inf_metal = 1.0

    dv_lin = linear_dv_screen(sigma_grid, lambda_tf_nm, kappa_inf_metal)

    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    ax.plot(sigma_grid, dv_lin, "C0-", lw=2.0, label=r"Linear T–K: $\Delta V = \sigma\lambda_{TF}/(\varepsilon_0 \kappa_{\infty,m})$")

    # Reference scales (typical literature P_r values)
    ref_scales = [
        ("BaTiO$_3$", 26.0, "C2"),
        ("BiFeO$_3$", 60.0, "C1"),
        ("AlScN (Sc$≈0.3$)", 113.0, "C3"),
    ]
    for label, sigma, color in ref_scales:
        dv = linear_dv_screen(np.array([sigma]), lambda_tf_nm, kappa_inf_metal)[0]
        ax.axvline(sigma, color=color, ls=":", alpha=0.6)
        ax.scatter([sigma], [dv], s=70, color=color, zorder=5)
        ax.annotate(
            f"{label}\n$\\sigma$={sigma:.0f} $\\mu$C/cm$^2$\n$\\Delta V$={dv:.2f} V",
            xy=(sigma, dv),
            xytext=(8, 8),
            textcoords="offset points",
            fontsize=9,
            color=color,
        )

    # Physical pinning floor reference
    ax.axhspan(0, 1.0, color="gray", alpha=0.12, zorder=0,
               label="Fermi-level pinning floor (~0.5–1 V)")
    ax.axhline(1.0, color="gray", ls="--", lw=1.0, alpha=0.7)

    ax.set_xlabel(r"Polarization-bound surface charge $|\sigma_{pol}|$  [$\mu$C/cm$^2$]")
    ax.set_ylabel(r"Predicted screening dipole $\Delta V_{screen}$  [V]")
    ax.set_title(
        "BEFORE Path A: linear Thomas–Fermi screening blows past physical pinning floor at AlScN scales\n"
        f"($\\lambda_{{TF}}$={lambda_tf_nm*10:.1f} Å, $\\kappa_{{\\infty,metal}}$={kappa_inf_metal:.1f})"
    )
    ax.set_xlim(0, 150)
    ax.set_ylim(0, max(8.0, dv_lin.max() * 1.05))
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out_path = os.path.join(OUT_DIR, "wave2_5_screening_before.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[wrote] {out_path}")
    return dv_lin


def figure_barriers_before():
    """Plot top/bottom barrier vs P, showing where the clamp at 0 eV kicks in."""
    # Values from src/materials.py (atomic units in code, eV here for plotting)
    wf_ti_ev = 4.33
    wf_al_ev = 4.28
    chi_alscn_ev = 2.0  # post-Wave-2 corrected value

    p_grid = np.linspace(-130, 130, 400)  # uC/cm^2
    sigma_top = -p_grid  # Eq. 1
    sigma_bot = +p_grid  # Eq. 2

    base_top = wf_ti_ev - chi_alscn_ev
    base_bot = wf_al_ev - chi_alscn_ev

    lambda_tf_nm = 0.045
    kappa_inf_metal = 1.0

    # Eq. 5/6 with leading minus
    dv_top = -linear_dv_screen(sigma_top, lambda_tf_nm, kappa_inf_metal)  # eV (numerically)
    dv_bot = -linear_dv_screen(sigma_bot, lambda_tf_nm, kappa_inf_metal)

    barrier_top_unclamped = base_top + dv_top
    barrier_bot_unclamped = base_bot + dv_bot
    barrier_top_clamped = np.maximum(barrier_top_unclamped, 0.0)
    barrier_bot_clamped = np.maximum(barrier_bot_unclamped, 0.0)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0), sharey=True)
    for ax, (label, base, dphi_unclamped, dphi_clamped) in zip(
        axes,
        [
            ("Top (Ti / AlScN)", base_top, barrier_top_unclamped, barrier_top_clamped),
            ("Bottom (Al / AlScN)", base_bot, barrier_bot_unclamped, barrier_bot_clamped),
        ],
    ):
        ax.plot(p_grid, dphi_unclamped, "C3-", lw=1.5, alpha=0.7, label="unclamped (Eq. 5/6)")
        ax.plot(p_grid, dphi_clamped, "C0-", lw=2.5, label=r"$\max(\varphi, 0)$ clamp (current code)")
        ax.axhline(base, color="k", ls=":", lw=1.0, alpha=0.6, label=f"base $W-\\chi$={base:.2f} eV")
        ax.axhline(0, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.axvline(113, color="gray", ls=":", lw=0.8, alpha=0.4)
        ax.axvline(-113, color="gray", ls=":", lw=0.8, alpha=0.4)
        ax.text(113, ax.get_ylim()[0] if ax.get_ylim()[0] > 0 else -8, "+P$_r$", ha="center", fontsize=8, color="gray")
        ax.text(-113, ax.get_ylim()[0] if ax.get_ylim()[0] > 0 else -8, "−P$_r$", ha="center", fontsize=8, color="gray")
        ax.set_xlabel("Polarization $P$  [$\\mu$C/cm$^2$]")
        ax.set_title(label)
        ax.legend(fontsize=9, loc="upper center")
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Effective Schottky barrier  [eV]")
    fig.suptitle(
        "BEFORE Path A: barrier vs polarization — clamp at 0 eV is load-bearing on one branch",
        fontsize=12,
    )
    fig.tight_layout()

    out_path = os.path.join(OUT_DIR, "wave2_5_barriers_before.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[wrote] {out_path}")


if __name__ == "__main__":
    print("Generating BEFORE diagnostic figures for Path A (saturating screening)...")
    figure_screening_before()
    figure_barriers_before()
    print("Done.")
