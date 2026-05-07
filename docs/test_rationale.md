# Stack-aware transport tests: rationale

This file accompanies `tests/test_stack_aware_transport.py`.  It explains, in
one paragraph each, what bug each test guards against and why the assertion is
phrased the way it is.  References below cite `polarization_barrier_coupling.md`
(the equation contract) and `transport_coupling_audit.md` (the punch list).

## Stack and conventions

All five tests instantiate the same nominal stack: Ti (top) / HfO_x (0 nm IL by
default) / AlScN (10 nm) / Al (bottom).  Voltage sign convention follows the
contract (`polarization_barrier_coupling.md` §1): `V_app > 0` means the top
electrode is at higher potential and electrons inject from the bottom (Al)
under forward bias.  `P > 0` means the polarization vector points toward the
bottom electrode.

The Preisach state is forced (not bisected) for each test.  Domains are pinned
via `Ferroelectric.set_states`, the self-consistent fields are solved
analytically by `Potential.fe_field_from_polarization`, and the snapshot is
built directly through `_build_state_snapshot` so that the requested
polarization configuration is preserved verbatim into the `BiasPointState`.

## T1 — Built-in flat-band asymmetry survives at P = 0

**Guards against:** a regression where polarization spuriously contributes to
the barrier at `P_net = 0`, or a regression that masks the electrode
work-function tilt entirely.  With `P = 0` the only source of asymmetry must
be electrode-only physics (`Eq. 5/6` reduce to zero), so `|log10 J(+V) − log10
J(−V)|` should reflect the work-function difference `W_Ti − W_Al ≈ 0.05–0.27
eV` (Ti = 4.33 eV, Al = 4.28 eV in `materials.py`).  The assertion is `> 0.05`
in log-decade and within a factor of 3 of the analytic thermionic prediction
`(W_Ti − W_Al) / (kT · ln 10) ≈ 0.84`.  A pure-tunneling-dominated regime can
also satisfy this test through the built-in field's effect on the WKB
integrand.

## T2 — Electroresistance: J(P_up) ≠ J(P_down) at the same V

**Guards against:** the dominant TER (tunneling electroresistance / FED ON-OFF)
bug.  `Eq. 5/6` of the contract require Tsymbal–Kohlstedt dipole steps
`Δφ_B = ∓ σ_pol · λ_TF · q / (ε_0 · κ_∞,M)` at the top and bottom Schottky
barriers.  Punch items 2, 3 and 14 of `transport_coupling_audit.md` flag that
the current barrier construction in `_barriers_ev` (`transport.py:96–119`)
ignores σ_pol entirely.  At `V = +5 V` the test compares J with all-up vs.
all-down domain configurations; the assertion is `|log10(J_up/J_down)| > 0.3`,
deliberately a low bar relative to the ~10–17 decades that a full TER signal
gives, so the test passes on first wiring of σ_pol into the barriers without
depending on the exact magnitude of the dipole shift.

## T3 — Non-zero zero-bias current when P ≠ 0

**Guards against:** the `abs()` / `_voltage_sign(V_app)` sign-stripping bug
(`transport.py:204, 223, 240, 268, 270, 277`, punch items 6, 7, 8, 11).  At
`V_app = 0` with `P = +Pr` the depolarizing field `E_FE = P t_IL /
(ε_IL t_FE + ε_FE t_IL)` from `Eq. 12` is nonzero and physically must drive a
finite leakage current.  Today every transport channel except TE multiplies by
`np.sign(V_app)` and returns 0 at `V_app = 0`.  Post-C-COUPLE the directional
prefactor must follow `np.sign(state.fe_field_au)` (or an equivalent supply-
function form), so the current at `V = 0, P = +Pr` becomes nonzero.  The
assertion floor `1 × 10⁻¹⁵ A/cm² = 1 × 10⁻¹¹ A/m²` is well above floating-point
noise but far below any DC measurement limit, so it triggers on any genuine
sign-preserving implementation while staying robust to model-mass artifacts.

## T4 — Image-force β coefficient uses κ_∞

**Guards against:** punch item 4 of `transport_coupling_audit.md` and `Eq. 15`
of `polarization_barrier_coupling.md`.  The Schottky / image-force lowering is
`Δφ_IF = sqrt(q³ E / (4 π ε_0 κ_∞))`, where the *optical* permittivity
`κ_∞ ≈ 4.6` for AlScN must be used because the image charge is induced by the
transiting carrier on a femtosecond timescale (faster than ionic relaxation).
The test calls `_schottky_lowering_ev` directly with a known field
`E = 1 × 10⁹ V/m` and the AlScN object's `eps_inf` field, then computes
`β_impl = lowering / sqrt(E)` and compares it to the analytic
`β_∞ = sqrt(q / (4 π ε_0 κ_∞))`.  The 1% relative-error gate catches the
"using `κ_static`" bug, which would shift β by `sqrt(κ_static / κ_∞)
≈ sqrt(18 / 4.6) ≈ 1.98` (~98 %).  Pre-C-COUPLE this test errors out at
`getattr(alscn, "eps_inf", None)` because the field has not yet been added to
`Ferroelectric` (punch item 15 in `transport_coupling_audit.md`, row "AlScN
optical k_∞" in `parameter_audit_20260506.md`).

## T5 — Per-mechanism breakdown sums to total

**Guards against:** double counting between mechanisms, especially PF and TAT
both claiming the same trap channel, or a missing mechanism whose current is
silently added into `total_a_per_m2` only.  This is a structural test on
`TransportEvaluator.evaluate`'s contract: `total_a_per_m2` must equal the sum
of `thermionic + tunneling + poole_frenkel + trap_assisted_tunneling + sclc +
background_leakage`.  The assertion floor `< 1e-9` relative error is loose
enough to absorb fp64 round-off and tight enough to catch any real arithmetic
discrepancy.  Today this test passes because the implementation literally
adds those six fields; the test guards future refactors (e.g. splitting PF FE
vs. PF IL, or adding an FN branch) from drifting the breakdown out of sync
with the total.

## Pre- vs. post-C-COUPLE expectations

| Test | Pre-C-COUPLE | Post-C-COUPLE | Diagnostic if failing post-merge |
|---|---|---|---|
| T1 | likely FAIL (currents underflow at small barriers / no built-in tilt in profile) | PASS | barriers blind to electrode WF, or tunneling profile not picking up `il_dv_bi` |
| T2 | FAIL | PASS | σ_pol not wired into `_barriers_ev` (Eq. 5/6), or `λ_TF` not multiplied in |
| T3 | FAIL | PASS | `abs()` / `_voltage_sign` still in PF/TAT/DT — current sign is V_app-locked |
| T4 | ERROR (no `eps_inf` attribute) | PASS | image-force still using static κ; check `_schottky_lowering_ev` callers |
| T5 | PASS | PASS | structural — should always pass; failing means double-count or missing mech |
