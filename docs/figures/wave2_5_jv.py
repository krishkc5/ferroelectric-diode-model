"""End-to-end J-V comparison: BEFORE Path A (linear screening, clamp load-bearing)
vs AFTER Path A (Eq. 4' tanh saturation).

Sweeps V from -5 to +5 V at fixed P = +Pr (P_up) and -Pr (P_down) for the
Ti / HfOx(0nm) / AlScN(10nm) / Al stack and plots |J| on log scale.

The "BEFORE" curves use a very large fermi_level_pinning_v (effectively linear);
the "AFTER" curves use the new default 1.0 V.
"""
from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))
_SRC = os.path.join(_REPO, "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
sys.path.insert(0, os.path.join(_REPO, "tests"))

from atomicunits import AtomicUnits  # noqa: E402
from materials import Materials  # noqa: E402
from preisach import Ferroelectric as PreisachFerroelectric  # noqa: E402
from transport_fed import FerroelectricDiode  # noqa: E402
from transport_potential import Potential  # noqa: E402
from transport_solver import SelfConsistentSolver  # noqa: E402
from transport import TransportEvaluator, TransportModelParameters  # noqa: E402

# Reuse the test file's helpers — they implement the exact pinned-state pattern.
from test_stack_aware_transport import build_stack, force_polarization_all, solve_at  # noqa: E402

OUT_DIR = os.path.dirname(os.path.abspath(__file__))


def jv_sweep(voltages, params, sign):
    diode, potential, solver, fe_model = build_stack(il_nm=0.0, fe_nm=10.0)
    force_polarization_all(fe_model, sign)
    pin_states = sign * np.ones(fe_model.num_domains, dtype=np.int8)
    j_arr = np.empty_like(voltages, dtype=float)
    for i, v_app in enumerate(voltages):
        state = solve_at(solver, fe_model, pin_states, float(v_app))
        evaluator = TransportEvaluator(diode, potential, parameters=params)
        out = evaluator.evaluate(state)
        j_arr[i] = out.total_a_per_m2
    return j_arr


def main():
    voltages = np.linspace(-5.0, 5.0, 51)

    params_after = TransportModelParameters()  # default fermi_level_pinning_*_v = 1.0
    params_before = TransportModelParameters(
        fermi_level_pinning_top_v=1.0e6,
        fermi_level_pinning_bot_v=1.0e6,
    )

    j_up_before = jv_sweep(voltages, params_before, +1)
    j_down_before = jv_sweep(voltages, params_before, -1)
    j_up_after = jv_sweep(voltages, params_after, +1)
    j_down_after = jv_sweep(voltages, params_after, -1)

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), sharey=True)
    floor = 1e-2

    for ax, (label, j_up, j_down) in zip(
        axes,
        [
            ("BEFORE Path A: linear screening (clamp load-bearing)", j_up_before, j_down_before),
            ("AFTER Path A: Eq. 4′ tanh saturation", j_up_after, j_down_after),
        ],
    ):
        ax.semilogy(voltages, np.maximum(np.abs(j_up), floor), "C0-", lw=2.0, label="$P=+P_r$ (up)")
        ax.semilogy(voltages, np.maximum(np.abs(j_down), floor), "C3-", lw=2.0, label="$P=-P_r$ (down)")
        ax.axvline(0, color="gray", lw=0.6, alpha=0.5)
        ax.set_xlabel("Applied voltage $V_{app}$  [V]")
        ax.set_title(label, fontsize=11)
        ax.grid(True, alpha=0.3, which="both")
        ax.legend(fontsize=10, loc="lower right")

    axes[0].set_ylabel(r"$|J|$  [A/m$^2$]")
    fig.suptitle(
        "End-to-end J–V: Ti / HfO$_x$(0 nm) / AlScN(10 nm) / Al — Path A brings electroresistance into the realistic range",
        fontsize=11,
    )
    fig.tight_layout()

    out_path = os.path.join(OUT_DIR, "wave2_5_jv_compare.png")
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"[wrote] {out_path}")

    # Console table
    def at(v):
        return int(np.argmin(np.abs(voltages - v)))
    print()
    print(f"{'V_app':>6}  {'BEFORE J_up':>14}  {'BEFORE J_down':>14}  {'log10 |ratio|':>14}  | "
          f"{'AFTER J_up':>14}  {'AFTER J_down':>14}  {'log10 |ratio|':>14}")
    for v in [-5, -2, 0, 2, 5]:
        i = at(v)
        b_ratio = np.log10(max(abs(j_up_before[i]), 1e-30) / max(abs(j_down_before[i]), 1e-30))
        a_ratio = np.log10(max(abs(j_up_after[i]), 1e-30) / max(abs(j_down_after[i]), 1e-30))
        print(f"{v:>+6}  {j_up_before[i]:>+14.3e}  {j_down_before[i]:>+14.3e}  {abs(b_ratio):>14.3f}  | "
              f"{j_up_after[i]:>+14.3e}  {j_down_after[i]:>+14.3e}  {abs(a_ratio):>14.3f}")


if __name__ == "__main__":
    print("Generating end-to-end J-V comparison (Path A vs linear screening)...")
    main()
    print("Done.")
