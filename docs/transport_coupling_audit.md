# Transport-to-Polarization Coupling Audit (C-COUPLE Wave-1)

**Stack:** Ti / HfO_x / AlScN / Al
**Scope:** Determine whether σ_pol (polarization-bound interface charge) and the self-consistent E_FE actually flow into the transport mechanisms (TE, DT, FN, PF, TAT, SCLC), and with what sign / κ / λ_TF.
**Method:** Read-only inspection of `src/transport_fed.py`, `src/transport_potential.py`, `src/transport_solver.py`, `src/transport.py`, `src/simulation_types.py`, `src/hysteresis_core.py`, `src/material_types.py`, `src/materials.py`. References below cite file:line and quote the relevant code verbatim.
**Companion contracts:** `docs/polarization_barrier_coupling.md` (Eq. 1–15), `docs/parameter_audit_20260506.md`.

---

## Q1. Where is σ_pol computed?

There is **no explicit `σ_pol` (polarization-bound interface charge) variable anywhere in the transport stack.** A `grep` for `sigma_pol`, `bound_charge`, `pol_charge`, `surface_charge` returns zero hits. The only polarization quantity propagated downstream is the volume-averaged FE polarization `polarization_au` and the *total* metal screening charge `sigma_s` (which lumps applied bias and bound charge together).

- `src/transport_fed.py:74`, `src/transport_fed.py:93–94`: the diode caches `self.fe_polarization = self.fe_model.avg_polarization()` — no per-interface bound charge is exposed.

  ```python
  self.fe_polarization = self.fe_model.avg_polarization()
  ...
  def get_polarization(self):
      self.fe_polarization = self.fe_model.avg_polarization()
      return self.fe_polarization
  ```

- `src/transport_potential.py:77–82`: `sigma_s` (the *metal* screening sheet) is computed from the `P · t / κ` voltage-divider expression. This is **not** σ_pol = ±P (Eq. 1, Eq. 2). It is the screening surface-charge density on the electrodes.

  ```python
  sigma_s = (
      fed.dl_polarization * fed.dl_thickness / fed.dl_k
      + fe_polarization * fed.fe_thickness / fed.fe_k
      + AtomicUnits.epsilon_0 * effective_v_diff
  ) / denominator
  ```

- `src/transport_solver.py:58–60, 77–79`: `screening_charge_au` is plumbed onto `BiasPointState`, but it is `sigma_s`, not σ_pol.

  ```python
  screening_charge_au = self.potential.screening_charge_from_polarization(
      polarization_au, v_diff=voltage_au
  )
  ```

**Verdict:** σ_pol,top = −P and σ_pol,bot = +P (Eq. 1, Eq. 2) are nowhere in the code. The only FE-polarization-derived charge that exists is `sigma_s` (the metal-side screening sheet), which is the wrong object to feed into a Tsymbal–Kohlstedt barrier shift.

---

## Q2. Does σ_pol enter the effective Schottky barrier (top / bottom)?

**No.** Trace from `state.polarization_au` into `_barriers_ev`:

- `src/transport.py:96–119`: barriers are built from `top_work_fxn − chi`, then *only* image-force-lowered by the layer field. Polarization never enters as a Δφ_B,top or Δφ_B,bot dipole shift.

  ```python
  top_base_barrier = AtomicUnits.hartree_to_ev(self.diode.top_work_fxn - top_affinity)
  bottom_base_barrier = AtomicUnits.hartree_to_ev(self.diode.bottom_work_fxn - self.diode.fe_chi)
  top_barrier = max(
      top_base_barrier
      - self.parameters.thermionic_top_schottky_scale
      * self._schottky_lowering_ev(
          top_field_v_per_m,
          self.diode.insulator_k if self.diode.insulator_thickness != 0 else self.diode.fe_k,
      ),
      0.0,
  )
  ```

- `src/transport.py:121–163` (`_thermionic_current`) calls `_barriers_ev(state)` and never reads `state.polarization_au`. The *only* polarization influence reaches TE indirectly through `state.il_field_au` / `state.fe_field_au` (which feed the image-force lowering).

There is **no Eq. 5 / Eq. 6 contribution** anywhere: no `λ_TF`, no `κ_∞,Ti`, no `κ_∞,Al`, no `±P · q · λ_TF / (ε0 · κ_∞)` term added to `top_base_barrier` or `bottom_base_barrier`. The Thomas–Fermi screening lengths do exist in `transport_fed.py:48–71` (used to set the screening exponential decay region in the band diagram), but they never multiply σ_pol into a barrier shift. The κ used for image-force lowering is the **static** `insulator_k` / `fe_k` (κ_static), not κ_∞ — see Q6.

**Verdict:** TE barrier is decoupled from σ_pol. Polarization affects TE only via the indirect, magnitude-only (`abs(...)`) image-force enhancement of the field.

---

## Q3. Does Poole–Frenkel use the self-consistent E_FE?

**Partially yes for the magnitude, but the field magnitude only — sign-stripped.** The PF kernel reads `state.fe_field_au` directly (so it is *not* a polarization-blind voltage divider), but it then takes `abs(...)`, so the polarization-induced flip of E_FE never alters the magnitude branch.

- `src/transport.py:214–242` (`_poole_frenkel_component`):

  ```python
  field_au = state.fe_field_au
  epsilon_r = self.diode.fe_k
  ...
  field_v_per_m = abs(self._field_to_v_per_m(field_au))
  ...
  lowering_ev = self._schottky_lowering_ev(field_v_per_m, epsilon_r)
  activation_ev = max(trap_depth_ev - lowering_ev, 0.0)
  ...
  current = (
      ...
      * field_v_per_m
      * np.exp(-activation_ev / max(thermal_ev, 1e-9))
      * self._voltage_sign(state.voltage_v)
  )
  ```

So the *magnitude* of E_FE used in PF includes the polarization contribution from Eq. 12 (because `state.fe_field_au` is built by the self-consistent solver at `transport_solver.py:80–83` from `layer_fields["total_e_field_fe"]`). But the directional sign is overridden by `self._voltage_sign(state.voltage_v)`, **not** the actual sign of E_FE. With a polarization that flips E_FE direction at zero bias (depolarizing field, Eq. 12 `P > 0, V_app = 0 ⇒ E_FE > 0`), this routine returns 0 instead of a nonzero leakage.

Additional issues:
- The Frenkel lowering uses `_schottky_lowering_ev`, which in turn (`transport.py:67–73`) divides by `epsilon_r = self.diode.fe_k` — **the static κ**, not κ_∞,FE; and it uses the `4π` Schottky-image prefactor, not the `π` Poole–Frenkel prefactor required by the spec (cf. polarization_barrier_coupling §7, "PF" row: `Δφ_PF = sqrt(q³ E_FE / (π ε0 κ_∞,FE))`).

**Verdict:** PF gets the `|E_FE|` magnitude correctly via the self-consistent state, but (a) zero-bias polarization-driven leakage is killed by `voltage_sign(V_app)`, and (b) the Frenkel lowering uses the wrong κ and the wrong prefactor.

---

## Q4. Does Trap-Assisted Tunneling use the same E_FE?

**Yes for the WKB barrier shape and the field-enhancement multiplier; sign is again clamped to V_app.**

- `src/transport.py:256–270` (`_trap_assisted_tunneling_current`):

  ```python
  profile, barrier_ev = self._tunneling_profile(state)
  trap_depth_ev = AtomicUnits.hartree_to_ev(self.diode.fe_trap_depth)
  reduced_barrier_ev = np.maximum(barrier_ev - trap_depth_ev, 0.0)
  transmission = self._wkb_transmission(
      x_au=profile["x_au"],
      barrier_ev=reduced_barrier_ev,
      effective_mass=profile["effective_mass"],
  )
  field_enhancement = 1.0 + abs(self._field_to_v_per_m(state.fe_field_au)) / 1e8
  current = self.parameters.trap_assisted_prefactor_a_per_m2 * transmission * field_enhancement
  return float(current * self._voltage_sign(state.voltage_v)), transmission
  ```

`_tunneling_profile` (`transport.py:181–188`) calls `self.potential.barrier_region_profile(fe_polarization=state.polarization_au, ...)` — so the WKB integrand of the TAT inherits whatever polarization-tilt is (or is not) baked into the barrier profile. See Q5 for what is actually baked in.

`field_enhancement` again uses `abs(state.fe_field_au)` and the directional output is multiplied by `_voltage_sign(state.voltage_v)` — same polarization-blind sign convention as PF.

**Verdict:** TAT shares the polarization-coupling status of `barrier_region_profile`. Magnitude of E_FE flows through; sign is clamped to V_app, killing zero-bias polarization-driven TAT.

---

## Q5. Does the WKB integrand use a P-tilted barrier φ(x)?

**Partially. The barrier picks up a uniform tilt from the *total* FE field (Eq. 12), but it does NOT include the explicit Eq. 8/9 σ_pol/ε_FE depolarizing gradient as a separate term, and it includes no Eq. 5/6 dipole step at the FE/electrode interfaces.** The implementation collapses the polarization-coupling into the linear voltage drop `fe_dv_electrostatic`, which already contains the depolarizing contribution because Eq. 12 is the ε_IL·V_app + P·t_IL solution.

- `src/transport_potential.py:88–92`:

  ```python
  fe_dv_electrostatic = (
      (sigma_s - fe_polarization) / (fed.fe_k * AtomicUnits.epsilon_0) * fed.fe_thickness
      ...
  )
  ```

  This is `(σ_s − P) · t_FE / ε_FE`, which algebraically equals `E_FE · t_FE` from Eq. 12 — so the *slope* of the FE band edge inside the FE region is correct.

- `src/transport_potential.py:181–189` builds the FE band edge as a single linear segment `x · fe_dv_electrostatic / t_FE` plus the constant `fe_v_barrier`:

  ```python
  return (
      state["v_top_interface"]
      + state["il_dv_electrostatic"]
      + pos / fed.fe_thickness * state["fe_dv_electrostatic"]
  )
  ```

- `src/transport_potential.py:212–225` (`barrier_potential`) returns a constant `top_fermi_e + fe_v_barrier` inside the FE region:

  ```python
  if x <= fe_end:
      return fed.bottom_fermi_e + self.fe_v_barrier
  ```

  No Eq. 5 (`+σ_pol,top · λ_TF · q / (ε0 κ_∞,Ti)`) step at z = b, no Eq. 6 (`−σ_pol,bot · λ_TF · q / (ε0 κ_∞,Al)`) step at z = c.

- `src/transport.py:181–212` (`_tunneling_current`) and `transport.py:404–423` (audit version) use this barrier directly; the WKB integrand (`_wkb_transmission`, `transport.py:165–179`) sees only this trapezoidal/affine profile.

**Verdict:** The WKB barrier captures the *total field tilt* (Eq. 12) inside the FE — including the depolarizing-field contribution as a hidden ingredient of `(σ_s − P)`. It does **not** apply the Tsymbal–Kohlstedt dipole steps at the electrodes (Eq. 5, Eq. 6), which is the dominant source of TER asymmetry. So TER asymmetry is missing from the WKB.

---

## Q6. Which κ is used in image-force lowering?

**Static κ (κ_static), not κ_∞ as required by Eq. 15.**

- `src/transport.py:67–73`:

  ```python
  def _schottky_lowering_ev(self, field_v_per_m, epsilon_r):
      if field_v_per_m <= 0 or epsilon_r <= 0:
          return 0.0
      lowering_volts = np.sqrt(
          constants.e * field_v_per_m / (4 * np.pi * constants.epsilon_0 * epsilon_r)
      )
      return float(lowering_volts)
  ```

- The caller (`transport.py:104–116`) passes the static dielectric `self.diode.insulator_k` or `self.diode.fe_k`:

  ```python
  top_barrier = max(
      top_base_barrier
      - self.parameters.thermionic_top_schottky_scale
      * self._schottky_lowering_ev(
          top_field_v_per_m,
          self.diode.insulator_k if self.diode.insulator_thickness != 0 else self.diode.fe_k,
      ),
  ```

- `src/material_types.py:23–35` and `materials.py:73–79, 158–164`: `Insulator.k` is the static dielectric (e.g. `hfo2.k = 18` static; baseline_hfo2.k = 16.64 static). There is no κ_∞ field on `Insulator` or `Ferroelectric`. This matches the parameter audit's flag (`docs/parameter_audit_20260506.md`: missing κ_∞ fields).

**Verdict:** Image-force lowering uses κ_static, violating Eq. 15. With κ_static ≈ 18 (HfO_x) and κ_∞ ≈ 4.5, current Δφ_IF is underestimated by ~sqrt(18/4.5) ≈ 2×.

---

## Q7. Which m* is used in the WKB integrand for the HfO_x region?

`m*_HfO_x = 0.11 m_e` (free-electron units), set on the `Insulator` material.

- `src/materials.py:73–79` (and `158–164` for the non-baseline variant):

  ```python
  baseline_hfo2 = Insulator(
      k=16.64,
      chi=AtomicUnits.ev_to_hartree(2.0),
      m_eff=0.11,
      name="$HfO_2$",
      breakdown_field=0
  )
  ```

- `src/transport_fed.py:37`: stored as `self.insulator_m_eff = insulator.m_eff`.

- `src/transport_fed.py:76–90` (`m_eff(x)` map): the WKB integrand picks up `insulator_m_eff` for `top_screen_region < x ≤ il_end`:

  ```python
  if x <= il_end:
      return self.insulator_m_eff
  ```

- The WKB integrand at `transport.py:175–177` consumes this mass directly:

  ```python
  mass_mid = 0.5 * (effective_mass[:-1] + effective_mass[1:]) * constants.m_e
  exponent_density = np.sqrt(2.0 * mass_mid * barrier_mid_j) / constants.hbar
  ```

**Verdict:** HfO_x m* = 0.11 m_e (consistent with Robertson 2006 / common literature value, 0.08–0.20 m_e). This is plausible for HfO₂ but the legacy `baseline_hfo2` is used by `hysteresis_core.resolve_insulator("hfox")` (line 51), not the newer `materials.hfo2` (which has identical `m_eff = 0.11` so the choice is moot for m*).

---

## Q8. Thought experiment: P_FE → −P_FE, what changes per mechanism?

If we manually negate `polarization_au` on a `BiasPointState` (and propagate the flip to `fe_field_au`, `il_field_au`, `dl_field_au`, `screening_charge_*` consistently — which is how `_evaluate_candidate` does it at `transport_solver.py:22–34`):

| Mechanism | Affected by P-flip? | Why |
|---|---|---|
| **TE (`_thermionic_current`)** | **Partial / weak NO** | Barriers (`_barriers_ev`, `transport.py:96–119`) only depend on `|il_field|`, `|fe_field|` (`abs(...)`) plus electrode WF/χ. No σ_pol shift. So flipping P only changes the magnitude of E_IL/E_FE through Eq. 12/13 (which can change because `|E_FE|` can be larger or smaller depending on P sign at fixed V_app), and that propagates only into image-force lowering. **No direct dipole-step from σ_pol.** |
| **DT (`_tunneling_current`)** | **YES (magnitude only)** | `barrier_region_profile` reads `state.polarization_au` (`transport.py:181–188` → `transport_potential.py:88–92`). Slope of FE band tilts with sign of P (via `(σ_s − P)`). However, the directional prefactor is `bias_drive = self._voltage_sign(state.voltage_v)` (`transport.py:204`), so under fixed V_app the *sign* of the current does not flip with P; only the WKB transmission magnitude does. **No Eq. 5/Eq. 6 dipole step**, so TER asymmetry is partial. |
| **FN** | (not implemented as a separate kernel) | The WKB tunneling kernel is a single mechanism; there is no distinct FN branch in `transport.py`. So this row is N/A. |
| **PF (`_poole_frenkel_current`)** | **YES (magnitude only); sign forced to V_app** | `field_v_per_m = abs(self._field_to_v_per_m(field_au))` (`transport.py:223`) picks up P through `state.fe_field_au`. So `|E_FE|` changes with P-flip and modifies the Frenkel lowering and the linear `field_v_per_m` factor. But sign of current is `_voltage_sign(state.voltage_v)`, decoupled from P. At V_app = 0, P-flip gives no PF current sign change (it stays zero either way). |
| **TAT (`_trap_assisted_tunneling_current`)** | **YES (magnitude only)** | Same `barrier_region_profile` and same `abs(state.fe_field_au)` field_enhancement (`transport.py:268`), same `_voltage_sign(V_app)` directional prefactor. P-flip changes WKB transmission magnitude (via barrier slope) and field_enhancement magnitude. Sign of current does not flip with P. |
| **SCLC (`_sclc_current`)** | **YES (magnitude only)** | `voltage_drop_v = abs(AtomicUnits.convert_back_volts(state.fe_field_au * self.diode.fe_thickness))` (`transport.py:277`) → P-flip changes the magnitude of `E_FE · t_FE` through Eq. 12. Directional prefactor is `np.sign(state.voltage_v)` (`transport.py:288`), so sign is V_app-locked. |

**Decoupling summary.** Flipping P does change the WKB transmission and the magnitudes of PF/TAT/SCLC currents through the self-consistent E_FE — so transport is not *entirely* decoupled from the Preisach state. But the dominant TER signal (asymmetric Δφ_B between top and bottom electrodes from Eq. 5/Eq. 6) is missing entirely, every mechanism's *sign* is locked to V_app rather than reflecting σ_pol-driven directionality, image-force lowering uses the wrong κ, and the WKB integrand never sees a dipole step at the FE/electrode interfaces.

---

## PUNCH LIST (gaps for Wave-2 C-COUPLE to fix)

1. **σ_pol is never computed.** Add per-interface bound-charge variables σ_pol,top = −P and σ_pol,bot = +P next to where `polarization_au` is built. `src/transport_solver.py:57–79` (extend `BiasPointState` build to include `sigma_pol_top_au` / `sigma_pol_bot_au`); applies **Eq. 1, Eq. 2**.

2. **No Tsymbal–Kohlstedt dipole shift on the top barrier.** Add `Δφ_B,top = +σ_pol,top · λ_TF,Ti · q / (ε0 · κ_∞,Ti)` to the barrier construction. `src/transport.py:96–119` (`_barriers_ev`); applies **Eq. 5**.

3. **No Tsymbal–Kohlstedt dipole shift on the bottom barrier.** Add `Δφ_B,bot = −σ_pol,bot · λ_TF,Al · q / (ε0 · κ_∞,Al)`. `src/transport.py:96–119` (`_barriers_ev`); applies **Eq. 6**.

4. **Image-force lowering uses κ_static instead of κ_∞.** Replace `self.diode.insulator_k` / `self.diode.fe_k` with κ_∞ fields (which need to be added to `material_types.Insulator` / `Ferroelectric` and populated in `materials.py`). `src/transport.py:67–73` (`_schottky_lowering_ev`) and callers at `transport.py:104–116`; applies **Eq. 15**.

5. **PF Frenkel lowering uses the Schottky `4π` form and κ_static.** Replace with `Δφ_PF = sqrt(q³ E_FE / (π ε0 κ_∞,FE))`. `src/transport.py:229` (the `_schottky_lowering_ev` reuse inside `_poole_frenkel_component`); applies **§7 PF row** of `polarization_barrier_coupling.md`.

6. **PF current sign is `voltage_sign(V_app)`, killing zero-bias depolarizing-field leakage.** Switch to `np.sign(state.fe_field_au)` (or to a vector-J form `J = q n μ E` that carries the sign of E_FE). `src/transport.py:240` (and the `abs(...)` at line 223); applies **Eq. 12** (E_FE retains sign even at V_app = 0).

7. **TAT field-enhancement and current sign are V_app-locked.** Same fix as Punch 6: drop `abs(...)` on `state.fe_field_au` and replace `_voltage_sign(state.voltage_v)` with `np.sign(state.fe_field_au)`. `src/transport.py:268, 270`; applies **Eq. 12**.

8. **DT current sign is V_app-locked.** WKB transmission already varies with P, but `bias_drive = self._voltage_sign(state.voltage_v)` clamps the directional prefactor. Replace with a sign derived from σ_pol-shifted barriers (or from a Tsu–Esaki-style integral of the supply functions on each side). `src/transport.py:204`; applies **Eq. 5, Eq. 6, Eq. 9**.

9. **WKB barrier inside FE has no dipole step at the metal interfaces.** `Potential.barrier_potential` returns a flat `bottom_fermi_e + fe_v_barrier` across the FE. Insert the `Δφ_B,top`/`Δφ_B,bot` step contributions at `z = b` and `z = c`. `src/transport_potential.py:212–225` (and the symmetric handling in `wf_potential` lines 227–248); applies **Eq. 5, Eq. 6**.

10. **WKB barrier inside FE does not separately encode the σ_pol/ε_FE depolarizing gradient (Eq. 8/9).** Although Eq. 12 provides the right slope through `(σ_s − P)`, the Eq. 9 form requires the explicit `ΔU_pol(x) = −(P/ε_FE) q x` term so that future barrier-shape modifications (e.g. interface charge layers, stripe domains) compose correctly. `src/transport_potential.py:88–92` and `159–209` (refactor `electrostatic_potential` so the P-only contribution is an additive segment); applies **Eq. 8, Eq. 9**.

11. **TE barriers strip directional information via `abs(self._field_to_v_per_m(...))`.** With σ_pol-driven asymmetry (Eq. 5/6), top and bottom barriers move in *signed* opposite directions; the magnitude-only field obscures the asymmetry that drives TER. `src/transport.py:97–100`; applies **Eq. 5, Eq. 6, Eq. 12, Eq. 13**.

12. **PF / Schottky use `self.diode.fe_trap_depth` fallback when `poole_frenkel_trap_depth_ev` is None, but this fallback uses the AlScN trap depth even when the PF region is set to `"il"` (HfO_x).** `src/transport.py:224–228`; cross-flagged by `docs/parameter_audit_20260506.md` (HfO_x trap-depth fallback bug); applies **§7 PF row** (PF lives in the *FE* bulk per the contract; if PF is forced into IL the kernel needs an HfO_x trap depth, not AlScN's).

13. **`Potential._electrostatic_state` lumps applied bias and bound charge into `sigma_s` so that `(σ_s − P)` is the only quantity exposed.** Expose `E_FE` and `E_IL` directly from the Eq. 12 / Eq. 13 closed forms so that downstream consumers can tell apart the depolarizing and applied-field contributions. `src/transport_potential.py:67–127`; applies **Eq. 10, Eq. 11, Eq. 12, Eq. 13**.

14. **No `λ_TF` for the dielectric stack accounted for.** `transport_fed.py:47–71` computes `top_screening_len` / `bottom_screening_len` from electrode `e_f` / `n0` (Thomas–Fermi in the *metal*), which is correct for Eq. 5/Eq. 6 — but the values are never multiplied into a barrier shift, only used to define the screening-decay region of `electrostatic_potential` (`transport_potential.py:170–175, 200–207`). Wire these existing `λ_TF` values into the new Eq. 5/Eq. 6 terms. `src/transport.py:96–119`; applies **Eq. 5, Eq. 6**.

15. **Material classes lack κ_∞ fields entirely.** Cross-referenced with `docs/parameter_audit_20260506.md`; required by Punches 2, 3, 4, 5. `src/material_types.py:22–53` (add `k_optical` to `Insulator` and `Ferroelectric`) and `src/materials.py:56–172` (populate κ_∞ for HfO_x ≈ 4.5, AlScN ≈ 4.7, Al/Ti Drude tail ~ 1–few); applies **§3, §6**.

---

**End of audit.** Wave-2 implementation must address all 15 punch-list items.
