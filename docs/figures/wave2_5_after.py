"""AFTER Path A: saturating screening + barrier-vs-P plots.

Generates:
  wave2_5_screening_after.png  — overlay linear vs tanh-saturated, with BTO/BFO/AlScN markers
  wave2_5_barriers_after.png   — top/bottom barrier vs polarization, no longer hitting clamp
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as constants

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
_SRC = os.path.join(_REPO, "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def linear_dv(sigma_uc_cm2, lambda_tf_nm, kappa_inf):
    sigma_si = sigma_uc_cm2 * 1e-2
    lambda_si = lambda_tf_nm * 1e-9
    return sigma_si * lambda_si / (constants.epsilon_0 * kappa_inf)


def saturated_dv(sigma_uc_cm2, lambda_tf_nm, kappa_inf, dv_pin_v):
    """Eq. 4': ΔV_screen = ΔV_pin · tanh(ΔV_lin / ΔV_pin)."""
    dv_lin = linear_dv(sigma_uc_cm2, lambda_tf_nm, kappa_inf)
    return dv_pin_v * np.tanh(dv_lin / dv_pin_v)


def figure_screening_after():
    sigma_grid = np.linspace(0, 150, 400)
    lambda_tf_nm = 0.045
    kappa_inf_metal = 1.0
    dv_pin = 1.0  # V

    dv_lin = linear_dv(sigma_grid, lambda_tf_nm, kappa_inf_metal)
    dv_sat = saturated_dv(sigma_grid, lambda_tf_nm, kappa_inf_metal, dv_pin)

    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    ax.plot(sigma_grid, dv_lin, "C0--", lw=1.6, alpha=0.6, label=r"Linear T–K (Eq. 4)  [unphysical at high $\sigma$]")
    ax.plot(sigma_grid, dv_sat, "C0-", lw=2.5,
            label=r"Eq. 4′ saturating: $\Delta V_{pin} \cdot \tanh(\Delta V_{lin}/\Delta V_{pin})$, $\Delta V_{pin}=1.0$ V")

    ref_scales = [
        ("BaTiO$_3$", 26.0, "C2"),
        ("BiFeO$_3$", 60.0, "C1"),
        ("AlScN (Sc$≈0.3$)", 113.0, "C3"),
    ]
    for label, sigma, color in ref_scales:
        dv_l = linear_dv(np.array([sigma]), lambda_tf_nm, kappa_inf_metal)[0]
        dv_s = saturated_dv(np.array([sigma]), lambda_tf_nm, kappa_inf_metal, dv_pin)[0]
        ax.axvline(sigma, color=color, ls=":", alpha=0.5)
        ax.scatter([sigma], [dv_s], s=70, color=color, zorder=5, edgecolor="black", lw=0.5)
        ax.annotate(
            f"{label}\n$\\Delta V_{{lin}}$={dv_l:.2f} V → $\\Delta V_{{sat}}$={dv_s:.2f} V",
            xy=(sigma, dv_s),
            xytext=(8, 14),
            textcoords="offset points",
            fontsize=8.5,
            color=color,
        )

    ax.axhspan(0, dv_pin, color="gray", alpha=0.10, zorder=0)
    ax.axhline(dv_pin, color="gray", ls="--", lw=1.0, alpha=0.7,
               label=f"Fermi-level pinning floor $\\Delta V_{{pin}}={dv_pin:.1f}$ V")

    ax.set_xlabel(r"Polarization-bound surface charge $|\sigma_{pol}|$  [$\mu$C/cm$^2$]")
    ax.set_ylabel(r"Screening dipole $\Delta V_{screen}$  [V]")
    ax.set_title(
        "AFTER Path A: tanh-saturated screening recovers BTO/BFO linear regime AND caps AlScN at the pinning floor"
    )
    ax.set_xlim(0, 150)
    ax.set_ylim(0, max(8.0, dv_lin.max() * 1.05))
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    out_path = os.path.join(OUT_DIR, "wave2_5_screening_after.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[wrote] {out_path}")


def figure_barriers_after():
    wf_ti_ev = 4.33
    wf_al_ev = 4.28
    chi_alscn_ev = 2.0

    p_grid = np.linspace(-130, 130, 400)
    sigma_top = -p_grid  # Eq. 1
    sigma_bot = +p_grid  # Eq. 2

    base_top = wf_ti_ev - chi_alscn_ev
    base_bot = wf_al_ev - chi_alscn_ev

    lambda_tf_nm = 0.045
    kappa_inf_metal = 1.0
    dv_pin = 1.0

    # Linear (BEFORE) — Eq. 5 with leading minus
    barrier_top_lin = base_top + (-linear_dv(sigma_top, lambda_tf_nm, kappa_inf_metal))
    barrier_bot_lin = base_bot + (-linear_dv(sigma_bot, lambda_tf_nm, kappa_inf_metal))
    # Saturated (AFTER)
    barrier_top_sat = base_top + (-saturated_dv(sigma_top, lambda_tf_nm, kappa_inf_metal, dv_pin))
    barrier_bot_sat = base_bot + (-saturated_dv(sigma_bot, lambda_tf_nm, kappa_inf_metal, dv_pin))

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.0), sharey=True)
    for ax, (label, base, lin, sat) in zip(
        axes,
        [
            ("Top (Ti / AlScN)", base_top, barrier_top_lin, barrier_top_sat),
            ("Bottom (Al / AlScN)", base_bot, barrier_bot_lin, barrier_bot_sat),
        ],
    ):
        ax.plot(p_grid, np.maximum(lin, 0.0), "C3-", lw=1.5, alpha=0.5,
                label=r"BEFORE: linear + $\max(\varphi,0)$ clamp")
        ax.plot(p_grid, sat, "C0-", lw=2.5, label="AFTER: Eq. 4′ saturated (no clamp needed)")
        ax.axhline(base, color="k", ls=":", lw=1.0, alpha=0.6, label=f"base $W-\\chi$={base:.2f} eV")
        ax.axhline(0, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.axvline(113, color="gray", ls=":", lw=0.8, alpha=0.4)
        ax.axvline(-113, color="gray", ls=":", lw=0.8, alpha=0.4)
        ax.text(113, -0.6, "+P$_r$", ha="center", fontsize=8, color="gray")
        ax.text(-113, -0.6, "−P$_r$", ha="center", fontsize=8, color="gray")
        ax.set_xlabel("Polarization $P$  [$\\mu$C/cm$^2$]")
        ax.set_title(label)
        ax.legend(fontsize=9, loc="upper center")
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel("Effective Schottky barrier  [eV]")
    axes[0].set_ylim(-1.0, 9.5)
    fig.suptitle(
        "AFTER Path A: barrier vs polarization — both branches stay physical, clamp is no longer load-bearing",
        fontsize=12,
    )
    fig.tight_layout()

    out_path = os.path.join(OUT_DIR, "wave2_5_barriers_after.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[wrote] {out_path}")


if __name__ == "__main__":
    print("Generating AFTER diagnostic figures for Path A (saturating screening)...")
    figure_screening_after()
    figure_barriers_after()
    print("Done.")
