"""Stack-aware transport tests for the post-C-COUPLE physics contract.

These tests target the polarization-barrier coupling, image-force kappa_inf,
PF Frenkel pi-prefactor, and the removal of abs() / voltage_sign sign-stripping
in src/transport.py. They are written against the contract in
docs/polarization_barrier_coupling.md (Eq. 1-15) and the punch list in
docs/transport_coupling_audit.md.

Pre-C-COUPLE: T2, T3 and T4 are expected to FAIL.  T1 and T5 may pass or fail
depending on the current state of the code.

Run with:

    uv run python tests/test_stack_aware_transport.py

Each test prints the numbers it asserts on so failures are diagnosable.
"""

from __future__ import annotations

import math
import os
import sys
import traceback

# Make src/ importable: project layout uses bare imports inside src/.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_SRC = os.path.join(_REPO, "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np  # noqa: E402
import scipy.constants as constants  # noqa: E402

from atomicunits import AtomicUnits  # noqa: E402
from materials import Materials  # noqa: E402
from preisach import Ferroelectric as PreisachFerroelectric  # noqa: E402
from transport_fed import FerroelectricDiode  # noqa: E402
from transport_potential import Potential  # noqa: E402
from transport_solver import SelfConsistentSolver  # noqa: E402
from transport import TransportEvaluator, TransportModelParameters  # noqa: E402


# ---------------------------------------------------------------------------
# Stack construction helper.  The Wave-2 spec stack: Ti(top)/HfOx(IL)/AlScN/Al.
# ---------------------------------------------------------------------------
def build_stack(*, il_nm: float, fe_nm: float, num_domains: int = 400, seed: int = 0):
    """Return (diode, potential, solver, fe_model) for Ti/HfOx/AlScN/Al."""
    fe_model = PreisachFerroelectric(
        num_domains=num_domains,
        c_a_mean=1.54,
        c_a_std=0.04,
        seed=seed,
    )
    diode = FerroelectricDiode(
        insulator_thickness=AtomicUnits.nm_to_bohr(il_nm),
        fe_thickness=AtomicUnits.nm_to_bohr(fe_nm),
        dead_layer_thickness=0.0,
        top_electrode=Materials.titanium_electrode,
        bottom_electrode=Materials.aluminum_electrode,
        insulator=Materials.hfo2,
        ferroelectric=Materials.alscn,
        fe_model=fe_model,
    )
    potential = Potential(diode, 0.0)
    solver = SelfConsistentSolver(
        ferroelectric=fe_model,
        potential=potential,
        max_iter=2000,
        threshold=0.5,
    )
    return diode, potential, solver, fe_model


def force_polarization_zero(fe_model: PreisachFerroelectric, seed: int = 12345) -> None:
    """Force net polarization ~ 0 by deterministic half-up / half-down assignment.

    Half +1, half -1 with a deterministic shuffle so that the *signed* mean of
    s_n * P_s,n cancels to floating-point noise.  Pure 50/50 random would leave
    a non-trivial residual; we explicitly pair each domain with its sign-mate.
    """
    n = fe_model.num_domains
    rng = np.random.default_rng(seed)
    states = np.empty(n, dtype=np.int8)
    states[: n // 2] = 1
    states[n // 2 :] = -1
    rng.shuffle(states)
    fe_model.set_states(states)


def force_polarization_all(fe_model: PreisachFerroelectric, sign: int) -> None:
    if sign not in (-1, 1):
        raise ValueError("sign must be +/-1")
    fe_model.set_states(sign * np.ones(fe_model.num_domains, dtype=np.int8))


def solve_at(solver: SelfConsistentSolver, fe_model: PreisachFerroelectric,
             pin_states, voltage_v: float):
    """Solve the self-consistent state at voltage_v with the domain states pinned.

    The solver mutates fe_model.state_values during bisection; we pin them again
    after the solve so the final BiasPointState reflects exactly the requested
    polarization configuration.
    """
    fe_model.set_states(pin_states)
    voltage_au = AtomicUnits.convert_volts(voltage_v)
    # Build state directly without letting bisection flip domains: use
    # _build_state_snapshot with a mocked final_eval whose 'states' is pin_states.
    # Equivalent approach: temporarily zero the e_c_values so no flip occurs,
    # then restore.  Simpler: just build the bias-point state directly.
    e_self = solver.potential.fe_field_from_polarization(
        float(np.mean(pin_states * fe_model.p_s_values)), v_diff=voltage_au
    )
    final_eval = {
        "residual": 0.0,
        "states": pin_states.copy(),
        "polarization": float(np.mean(pin_states * fe_model.p_s_values)),
        "guess_field": e_self,
        "self_consistent_field": e_self,
    }
    solver.potential.set_vdiff(voltage_au)
    state = solver._build_state_snapshot(
        voltage=voltage_v,
        voltage_au=voltage_au,
        previous_polarization=final_eval["polarization"],
        final_eval=final_eval,
        converged=True,
        used_fallback=False,
    )
    return state


# ---------------------------------------------------------------------------
# Test runner infra
# ---------------------------------------------------------------------------
class TestResult:
    def __init__(self, name: str):
        self.name = name
        self.status: str = "ERROR"   # PASS / FAIL / ERROR
        self.detail: str = ""

    def __str__(self):
        return f"[{self.status}] {self.name}: {self.detail}"


def _run(test_fn, name: str) -> TestResult:
    result = TestResult(name)
    print(f"\n========== {name} ==========")
    try:
        test_fn(result)
        if result.status == "ERROR":
            # test_fn finished without setting; mark as PASS by convention.
            result.status = "PASS"
            if not result.detail:
                result.detail = "ok"
    except AssertionError as e:
        result.status = "FAIL"
        result.detail = str(e).splitlines()[0] if str(e) else "assertion failed"
        traceback.print_exc()
    except Exception as e:
        result.status = "ERROR"
        result.detail = f"{type(e).__name__}: {e}"
        traceback.print_exc()
    print(f"---> {result.status}: {result.detail}")
    return result


# ---------------------------------------------------------------------------
# T1: Built-in flat-band asymmetry survives at P=0.
# ---------------------------------------------------------------------------
def test_t1_builtin_asymmetry_at_p_zero(result: TestResult):
    """At P_net = 0 the device should still be asymmetric (Ti vs Al WF)."""
    diode, potential, solver, fe_model = build_stack(il_nm=0.0, fe_nm=10.0,
                                                     num_domains=400, seed=0)
    force_polarization_zero(fe_model, seed=42)
    pinned_states = fe_model.state_values.copy()
    p_uc = AtomicUnits.convert_back_polarization(
        float(np.mean(pinned_states * fe_model.p_s_values))
    )
    print(f"  pinned net polarization: {p_uc:.4e} uC/cm^2 (target ~0)")

    evaluator = TransportEvaluator(diode, potential)
    state_pos = solve_at(solver, fe_model, pinned_states, +2.0)
    breakdown_pos = evaluator.evaluate(state_pos)

    state_neg = solve_at(solver, fe_model, pinned_states, -2.0)
    breakdown_neg = evaluator.evaluate(state_neg)

    j_pos = abs(breakdown_pos.total_a_per_m2)
    j_neg = abs(breakdown_neg.total_a_per_m2)
    print(f"  J(+2V)={breakdown_pos.total_a_per_m2:+.4e} A/m^2,  "
          f"J(-2V)={breakdown_neg.total_a_per_m2:+.4e} A/m^2")
    print(f"  |J(+2V)|={j_pos:.4e},  |J(-2V)|={j_neg:.4e}")

    # Avoid log(0).
    assert j_pos > 0 and j_neg > 0, (
        f"Both currents must be nonzero at P=0 to test asymmetry "
        f"(got j_pos={j_pos}, j_neg={j_neg})"
    )

    log_asym = abs(math.log10(j_pos) - math.log10(j_neg))
    print(f"  |log10(J+/J-)| = {log_asym:.4f}")

    # Analytic prediction: at +2V Al injects (barrier ~ W_Al - chi_AlScN);
    # at -2V Ti injects (barrier ~ W_Ti - chi_AlScN).  log10 ratio ~
    # (W_Ti - W_Al) / (kT * ln10) at room T:
    kT_eV = constants.k * 300.0 / constants.e
    delta_wf = 4.33 - 4.28  # eV
    analytic_log_ratio = delta_wf / (kT_eV * math.log(10))
    print(f"  analytic |log10 ratio| (TE-only, dWF=0.05 eV): {analytic_log_ratio:.4f}")

    # Spec: > 0.05 in log10 magnitude, and within a factor of 3 of analytic.
    assert log_asym > 0.05, (
        f"Expected polarization-free asymmetry > 0.05 decades, got {log_asym:.4f}"
    )
    assert (analytic_log_ratio / 3.0) <= log_asym <= (analytic_log_ratio * 3.0), (
        f"Asymmetry {log_asym:.4f} not within factor of 3 of analytic "
        f"{analytic_log_ratio:.4f}"
    )
    result.detail = (f"|log10 J+/J-|={log_asym:.3f}, analytic~{analytic_log_ratio:.3f}, "
                     f"P_net={p_uc:.2e} uC/cm^2")


# ---------------------------------------------------------------------------
# T2: Electroresistance: J(P_up) != J(P_down) at the same V.
# ---------------------------------------------------------------------------
def test_t2_electroresistance_signature(result: TestResult):
    diode, potential, solver, fe_model = build_stack(il_nm=0.0, fe_nm=10.0,
                                                     num_domains=400, seed=0)
    evaluator = TransportEvaluator(diode, potential)

    states_up = np.ones(fe_model.num_domains, dtype=np.int8)
    state_up = solve_at(solver, fe_model, states_up, +5.0)
    p_up_uc = AtomicUnits.convert_back_polarization(state_up.polarization_au)
    bd_up = evaluator.evaluate(state_up)

    states_dn = -np.ones(fe_model.num_domains, dtype=np.int8)
    state_dn = solve_at(solver, fe_model, states_dn, +5.0)
    p_dn_uc = AtomicUnits.convert_back_polarization(state_dn.polarization_au)
    bd_dn = evaluator.evaluate(state_dn)

    print(f"  P_up = {p_up_uc:+.3e} uC/cm^2,  P_down = {p_dn_uc:+.3e} uC/cm^2")
    print(f"  J(P_up,  +5V) = {bd_up.total_a_per_m2:+.4e} A/m^2")
    print(f"  J(P_down,+5V) = {bd_dn.total_a_per_m2:+.4e} A/m^2")
    print(f"  top barriers: up={bd_up.top_barrier_ev:.4f} eV, "
          f"down={bd_dn.top_barrier_ev:.4f} eV")
    print(f"  bot barriers: up={bd_up.bottom_barrier_ev:.4f} eV, "
          f"down={bd_dn.bottom_barrier_ev:.4f} eV")

    j_up = abs(bd_up.total_a_per_m2)
    j_dn = abs(bd_dn.total_a_per_m2)
    assert j_up > 0 and j_dn > 0, (
        f"Both J(P_up) and J(P_down) must be nonzero "
        f"(got j_up={j_up}, j_dn={j_dn})"
    )
    log_er = abs(math.log10(j_up) - math.log10(j_dn))
    print(f"  |log10(J_up/J_down)| = {log_er:.4f}")

    # Spec: > 0.3 decades.  Pre-C-COUPLE this fails because barriers ignore sigma_pol.
    assert log_er > 0.3, (
        f"Electroresistance |log10(J_up/J_down)| > 0.3 expected "
        f"(barrier sigma_pol coupling, Eq. 5/6); got {log_er:.4f}"
    )
    result.detail = f"|log10(J_up/J_down)|={log_er:.3f} at +5V"


# ---------------------------------------------------------------------------
# T3: Non-zero zero-bias current when P != 0 (depolarization-driven).
# ---------------------------------------------------------------------------
def test_t3_zero_bias_polarization_current(result: TestResult):
    diode, potential, solver, fe_model = build_stack(il_nm=0.0, fe_nm=10.0,
                                                     num_domains=400, seed=0)
    evaluator = TransportEvaluator(diode, potential)

    states_up = np.ones(fe_model.num_domains, dtype=np.int8)
    state0 = solve_at(solver, fe_model, states_up, 0.0)
    p_uc = AtomicUnits.convert_back_polarization(state0.polarization_au)
    e_fe_mv_cm = state0.fe_field_mv_cm

    bd0 = evaluator.evaluate(state0)
    j0 = abs(bd0.total_a_per_m2)
    print(f"  P_net = {p_uc:+.3e} uC/cm^2 (~+Pr)")
    print(f"  E_FE  = {e_fe_mv_cm:+.4e} MV/cm  (depolarizing field, V_app=0)")
    print(f"  Mechanism breakdown at V=0:")
    print(f"    thermionic = {bd0.thermionic_a_per_m2:+.4e}")
    print(f"    tunneling  = {bd0.tunneling_a_per_m2:+.4e}")
    print(f"    PF         = {bd0.poole_frenkel_a_per_m2:+.4e}")
    print(f"    TAT        = {bd0.trap_assisted_tunneling_a_per_m2:+.4e}")
    print(f"    SCLC       = {bd0.sclc_a_per_m2:+.4e}")
    print(f"    total      = {bd0.total_a_per_m2:+.4e}")
    print(f"  |J(V=0, P=+Pr)| = {j0:.4e} A/m^2")

    # Spec: |J(0)| > 1e-15 A/cm^2 == 1e-11 A/m^2.
    threshold_a_per_m2 = 1e-15 * 1e4
    assert j0 > threshold_a_per_m2, (
        f"|J(V=0, P=+Pr)| = {j0:.4e} A/m^2 must exceed {threshold_a_per_m2:.1e} A/m^2 "
        f"(depolarizing field tilts bands; abs/voltage_sign clamp must be removed)"
    )
    result.detail = (f"|J(V=0)|={j0:.3e} A/m^2, threshold {threshold_a_per_m2:.1e} A/m^2, "
                     f"E_FE={e_fe_mv_cm:.3e} MV/cm")


# ---------------------------------------------------------------------------
# T4: Image-force beta coefficient uses kappa_inf.
# ---------------------------------------------------------------------------
def test_t4_image_force_kappa_inf(result: TestResult):
    diode, potential, _solver, _fe_model = build_stack(il_nm=0.0, fe_nm=10.0,
                                                       num_domains=200, seed=0)
    evaluator = TransportEvaluator(diode, potential)

    alscn = Materials.alscn
    # Post-C-COUPLE this attribute exists; pre-C-COUPLE it is absent and the
    # test fails with AttributeError (which is the bug we are guarding).
    eps_inf = getattr(alscn, "eps_inf", None)
    assert eps_inf is not None, (
        "Materials.alscn.eps_inf is missing.  C-COUPLE must add an "
        "eps_inf field to Ferroelectric / Insulator (~ 4.6 for AlScN, "
        "~ 4.0 for HfOx) so the image-force lowering uses kappa_inf "
        "(polarization_barrier_coupling.md Eq. 15)."
    )
    print(f"  AlScN.eps_inf = {eps_inf}")

    # Probe the implementation at a representative oxide field.
    E_test_v_per_m = 1.0e9   # 10 MV/cm
    lowering_v = evaluator._schottky_lowering_ev(E_test_v_per_m, eps_inf)
    beta_impl = lowering_v / math.sqrt(E_test_v_per_m)

    # Analytic Schottky beta = sqrt(q / (4 pi eps0 kappa_inf)) [V * m^0.5].
    beta_analytic = math.sqrt(
        constants.e / (4.0 * math.pi * constants.epsilon_0 * eps_inf)
    )
    rel_err = abs(beta_impl - beta_analytic) / beta_analytic

    # If the bug were still present (using static kappa ~ 16-18) the ratio of
    # implementation beta / analytic-with-kappa_inf would be ~ sqrt(eps_inf/k_static)
    # ~ 0.5; we show that for context.
    beta_static_check = math.sqrt(
        constants.e / (4.0 * math.pi * constants.epsilon_0 * alscn.k)
    )
    print(f"  E_test = {E_test_v_per_m:.2e} V/m")
    print(f"  lowering(eps_inf) = {lowering_v:.4e} V")
    print(f"  beta_impl    = {beta_impl:.4e} V*m^0.5")
    print(f"  beta(kappa_inf={eps_inf}) = {beta_analytic:.4e} V*m^0.5")
    print(f"  beta(kappa_static={alscn.k}) = {beta_static_check:.4e} V*m^0.5  (bug ref)")
    print(f"  rel_err vs kappa_inf = {rel_err:.3e}")

    assert rel_err < 0.01, (
        f"Image-force beta off by {rel_err*100:.2f}% from analytic kappa_inf "
        f"prediction ({beta_analytic:.4e}); implementation gave {beta_impl:.4e}.  "
        f"Likely _schottky_lowering_ev still divides by static kappa."
    )
    result.detail = (f"beta_impl={beta_impl:.3e}, beta_inf={beta_analytic:.3e}, "
                     f"rel_err={rel_err:.2e}")


# ---------------------------------------------------------------------------
# T5: Per-mechanism breakdown sums to total.
# ---------------------------------------------------------------------------
def test_t5_breakdown_sums_to_total(result: TestResult):
    diode, potential, solver, fe_model = build_stack(il_nm=0.0, fe_nm=10.0,
                                                     num_domains=400, seed=0)
    evaluator = TransportEvaluator(diode, potential)

    states_up = np.ones(fe_model.num_domains, dtype=np.int8)
    voltages = [-3.0, -1.0, 0.0, +1.0, +3.0]

    max_rel_err = 0.0
    for V in voltages:
        state = solve_at(solver, fe_model, states_up, V)
        bd = evaluator.evaluate(state)
        mech_sum = (
            bd.thermionic_a_per_m2
            + bd.tunneling_a_per_m2
            + bd.poole_frenkel_a_per_m2
            + bd.trap_assisted_tunneling_a_per_m2
            + bd.sclc_a_per_m2
            + bd.background_leakage_a_per_m2
        )
        diff = bd.total_a_per_m2 - mech_sum
        denom = max(abs(bd.total_a_per_m2), 1e-30)
        rel_err = abs(diff) / denom
        max_rel_err = max(max_rel_err, rel_err)
        print(f"  V={V:+.1f}V: total={bd.total_a_per_m2:+.4e}, "
              f"mech_sum={mech_sum:+.4e}, "
              f"|delta|/|total|={rel_err:.2e}")

    print(f"  worst |delta|/|total| over sweep = {max_rel_err:.3e}")
    assert max_rel_err < 1e-9, (
        f"Mechanism breakdown does not sum to total (worst rel err {max_rel_err:.3e}). "
        f"Possible double-count between PF and TAT or missing mechanism."
    )
    result.detail = f"max rel err = {max_rel_err:.2e}"


# ---------------------------------------------------------------------------
# Entry point.
# ---------------------------------------------------------------------------
def main() -> int:
    tests = [
        ("T1_builtin_asymmetry_at_P_zero", test_t1_builtin_asymmetry_at_p_zero),
        ("T2_electroresistance_signature", test_t2_electroresistance_signature),
        ("T3_zero_bias_polarization_current", test_t3_zero_bias_polarization_current),
        ("T4_image_force_kappa_inf", test_t4_image_force_kappa_inf),
        ("T5_breakdown_sums_to_total", test_t5_breakdown_sums_to_total),
    ]

    results = []
    for name, fn in tests:
        results.append(_run(fn, name))

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for r in results:
        print(str(r))

    n_pass = sum(1 for r in results if r.status == "PASS")
    n_fail = sum(1 for r in results if r.status == "FAIL")
    n_err = sum(1 for r in results if r.status == "ERROR")
    print(f"\n{n_pass} PASS, {n_fail} FAIL, {n_err} ERROR  (total {len(results)})")

    # Return non-zero if any test did not pass, so CI / shells see a failure.
    return 0 if (n_fail == 0 and n_err == 0) else 1


# Optional pytest hooks: each test_t* function is also discoverable by pytest
# (it just ignores the `result` argument convention via a thin wrapper).
def _pytest_wrap(name, fn):
    def _inner():
        r = TestResult(name)
        fn(r)
    _inner.__name__ = "test_" + name.lower()
    return _inner


test_T1 = _pytest_wrap("T1", test_t1_builtin_asymmetry_at_p_zero)
test_T2 = _pytest_wrap("T2", test_t2_electroresistance_signature)
test_T3 = _pytest_wrap("T3", test_t3_zero_bias_polarization_current)
test_T4 = _pytest_wrap("T4", test_t4_image_force_kappa_inf)
test_T5 = _pytest_wrap("T5", test_t5_breakdown_sums_to_total)


if __name__ == "__main__":
    sys.exit(main())
