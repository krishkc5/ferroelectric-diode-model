# Polarization–Barrier Coupling: Sign-Convention Contract

**Stack:** Ti (top) / HfO\_x (interfacial layer, "IL") / AlScN (ferroelectric, "FE") / Al (bottom).
**Status:** Wave-2 implementation contract. Reference equations by `(Eq. N)` below.
**Units:** SI throughout, except where noted as atomic units (a.u.) at the WKB code boundary.

---

## 1. Geometry and sign convention

Coordinate axis: `+z` runs **from the top electrode (Ti) into the device toward the bottom electrode (Al)**. The stack ordering along `+z` is:

```
   z = 0  : Ti (top metal)         [n_top  = +z hat at Ti/IL interface, pointing INTO IL]
   z = a  : Ti / HfOx interface
   z = b  : HfOx / AlScN interface
   z = c  : AlScN / Al interface   [n_bot  = -z hat at FE/Al interface, pointing INTO FE bulk from Al side... see §2]
```

Thicknesses: `t_IL = b − a` (HfO\_x), `t_FE = c − b` (AlScN). All in [m].

**Polarization sign.** `P > 0` means the polarization vector `\vec P = P \hat z` points **along `+z`**, i.e. from top toward bottom. Equivalently, the positive bound-charge face of the FE is the **bottom (FE/Al) interface** and the negative bound-charge face is the **top (FE/IL) interface**.

**Voltage sign.** `V_app > 0` means the **top electrode (Ti) is at higher potential than the bottom (Al)**. Under this convention electrons are injected from the **bottom (Al)** under forward bias `V_app > 0` and from the **top (Ti)** under reverse bias `V_app < 0`. The conventional electric field across the stack under `V_app > 0` (with `P = 0`) points along `+z` (from high to low potential).

ASCII band diagram (electron energy, `V_app = 0`, `P > 0`):

```
   E
   ^         _____                              ____
   |  Ti ---/     \----HfOx----\               /     \---- Al
   |        |  φB_top         |\             /|         φB_bot
   |        |                 | \___________/ |
   |        |                 |    AlScN(FE)  |
   |                          (negative bound  (positive bound
   |                           charge face)     charge face)
   +---------------------------------------------------> z
```

The `P > 0` bound charges **lower** the FE conduction band on the bottom side and **raise** it on the top side relative to the `P = 0` reference (derived in §3).

---

## 2. Polarization-bound surface charge

Bound surface charge density at a FE face is `\sigma_{pol} = \vec P \cdot \hat n_{out}`, where `\hat n_{out}` is the **outward** normal of the FE region (pointing **out of** the FE).

- **Top FE face (FE/IL interface, z = b):** `\hat n_{out,top} = -\hat z` (points back into IL, i.e. opposite to `+z`).

  ```
  σ_pol,top = P · (−ẑ) = −P            [C/m²]                         (Eq. 1)
  ```

- **Bottom FE face (FE/Al interface, z = c):** `\hat n_{out,bot} = +\hat z`.

  ```
  σ_pol,bot = P · (+ẑ) = +P            [C/m²]                         (Eq. 2)
  ```

Sign table:

| State | `σ_pol,top` (top FE face) | `σ_pol,bot` (bottom FE face) |
|---|---|---|
| `P > 0` (P along +z) | **negative**, `−|P|` | **positive**, `+|P|` |
| `P < 0` (P along −z) | **positive**, `+|P|` | **negative**, `−|P|` |

Charge neutrality of the FE slab is preserved: `σ_pol,top + σ_pol,bot = 0` (Eq. 3). [Tsymbal & Kohlstedt, *Science* **313**, 181 (2006), Eq. 1.]

---

## 3. Effective Schottky barrier shifts (Tsymbal–Kohlstedt screening)

The bound charge `σ_pol` at each FE/electrode side is screened by the metal over a Thomas–Fermi length `λ_TF` [m]. The screening dipole shifts the metal Fermi level relative to the dielectric band edges by an electrostatic potential drop

```
ΔV_screen = σ_pol · λ_TF / (ε0 · κ_∞_metal)         [V]               (Eq. 4)
```

where `ε0 = 8.8541878e−12 F/m` and `κ_∞_metal` is the **electronic (optical, high-frequency) relative permittivity** of the metal screening region. The electron Schottky barrier `φ_B` is shifted by `Δφ_B = q · ΔV_screen` in [J]. In [eV] the numerical conversion is simply `Δφ_B [eV] = ΔV_screen [V]` (i.e. drop the `q` factor when both sides are read in their natural units — `ΔV_screen` in volts, `Δφ_B` in electron-volts).

Sign derivation. The electron Schottky barrier is `φ_B = E_CB,diel − E_F,metal`. A bound charge `σ_pol` at the metal/dielectric interface induces screening electrons in the metal over a Thomas–Fermi length, producing an electrostatic potential drop `ΔV_screen = σ_pol · λ_TF / (ε0 · κ_∞,metal)` across the screening region. **Positive** bound charge against the metal raises the bulk-metal potential relative to the metal/dielectric interface — i.e. raises `E_F,metal` relative to `E_CB,diel` — and therefore **lowers** the electron Schottky barrier. The same physics applies at both faces: the dipole-shift convention is uniform with a leading minus sign:

```
Δφ_B,top = − (σ_pol,top) · λ_TF,top · q / (ε0 · κ_∞,Ti)        [J]   (Eq. 5)
Δφ_B,bot = − (σ_pol,bot) · λ_TF,bot · q / (ε0 · κ_∞,Al)        [J]   (Eq. 6)
```

with the convention `φ_B^{eff} = φ_B^{(0)} + Δφ_B`. Substituting Eq. 1–2 for the two polarization states:

| State | `Δφ_B,top` (with σ_pol,top = −P) | `Δφ_B,bot` (with σ_pol,bot = +P) | Physical effect |
|---|---|---|---|
| `P > 0` (down) | `+|P| λ_TF,top q / (ε0 κ_∞,Ti)` (raised)  | `−|P| λ_TF,bot q / (ε0 κ_∞,Al)` (lowered) | Top **raised**, bottom **lowered** — opposite directions. |
| `P < 0` (up)   | `−|P| λ_TF,top q / (ε0 κ_∞,Ti)` (lowered) | `+|P| λ_TF,bot q / (ε0 κ_∞,Al)` (raised)  | Top **lowered**, bottom **raised** — opposite directions. |

**The TER signature is the opposite-sign shift at the two interfaces.** The leading minus is uniform; the asymmetry comes from `σ_pol,top = −P` and `σ_pol,bot = +P` (Eqs. 1–2), not from the prefactor sign. Switching the polarity of `P` flips the shift at each interface independently, and the *difference* (top − bottom) is what tilts the barrier and drives the diode's ON/OFF ratio. Tsymbal–Kohlstedt 2006; Zhuravlev–Sabirianov–Jaswal–Tsymbal, *PRL* **94**, 246802 (2005); Pantel & Alexe, *PRB* **82**, 134105 (2010), Eqs. 4–6.

**Why `κ_∞`, not `κ_static`.** The screening charge in a metal is purely **electronic** (free-carrier + bound-electron polarization). Ionic lattice polarization does not participate in metallic screening. The relative permittivity that enters Thomas–Fermi screening of an interfacial dipole is therefore the **optical / high-frequency** value `κ_∞ ≈ n²` (Mehta, Silverman & Jacobs, *J. Appl. Phys.* **44**, 4490 (1973), §II). For Ti and Al, `κ_∞ ~ 1`–a few (use Drude tail values). Using `κ_static` here would over-screen by an order of magnitude and is a known sign-of-effect bug.

### 3.1 Saturating-screening correction at high `σ_pol` (Eq. 4′)

Eq. 4 is a **linear-response** result. It assumes the screening-charge density induced in the metal stays small compared to the metal's free-electron density `n_0`, so the Thomas–Fermi approximation `δn(r) = δn_TF` holds. For perovskite ferroelectrics (BaTiO₃: `P_r ~ 26 µC/cm²`, BiFeO₃: `~60 µC/cm²`) Eq. 4 predicts dipole shifts of `~1–3 V`, comparable to the metal/dielectric Fermi-level pinning floor `Δ_pin ~ 0.5–1 V` reported in BTO/BFO FTJ literature (Pantel–Alexe 2010 *PRB* **82**, 134105; Wen et al. *Nat. Mater.* **12**, 617 (2013); Gruverman et al. *Nano Lett.* **9**, 3539 (2009)).

For AlScN (`P_r ~ 113 µC/cm²` at Sc ≈ 0.3), Eq. 4 predicts `ΔV_screen ≈ 5.7 V` with `λ_TF = 0.45 Å` and `κ_∞,metal = 1`. At that magnitude, the implied screening sheet charge density `ρ_screen ≈ σ_pol / λ_TF ~ 2.5×10²⁹ m⁻³` exceeds the metal's own free-electron density (Ti, Al: `n_0 ~ 7×10²⁸ m⁻³`) — i.e. the linear-response assumption fails. Physically, what cuts the divergence is **Fermi-level pinning at the metal/dielectric interface**: once the dipole drives the metal Fermi level toward the dielectric band edge, interfacial trap states pin it, and additional bound charge no longer translates into proportional band shift.

A standard phenomenological fix that respects the linear limit at small `σ_pol` and saturates at the pinning floor at large `σ_pol` is:

```
ΔV_screen = ΔV_pin · tanh( ΔV_lin / ΔV_pin )         [V]               (Eq. 4′)
```

where `ΔV_lin = σ_pol · λ_TF / (ε0 · κ_∞,metal)` is the linear formula (Eq. 4) and `ΔV_pin` is the **Fermi-level-pinning floor** of the specific metal/dielectric interface, in [V]. Limits: `tanh(x) ≈ x` for `|x| ≪ 1` recovers Eq. 4 in the BTO/BFO regime; `tanh(x) → ±1` for `|x| ≫ 1` saturates at `±ΔV_pin` in the AlScN regime.

**Eq. 5/6 are unchanged** under Eq. 4′: the leading minus and the σ_pol sign convention propagate through `tanh` unchanged. Only the magnitude saturates.

**`ΔV_pin` value.** The pinning floor is interface-specific. For nitride/metal interfaces the canonical range is `0.5 – 1.0 V` (Casamento *APL* **120**, 152901 (2022); Tsymbal–Kohlstedt 2006 review §III). A defensible default is `ΔV_pin = 1.0 V` for both Ti/AlScN and Al/AlScN, with the caveat that the actual pinning depends on interface chemistry (e.g. interfacial Ti–O reduction layer can shift it). This single parameter replaces the implicit `max(barrier, 0)` clamp that was load-bearing in the linear-formula implementation.

**Domain of validity.** Eq. 4′ is a phenomenological interpolation between two limits. It correctly recovers (a) the linear-response Tsymbal–Kohlstedt regime at low `σ_pol` and (b) the pinned-Fermi-level regime at high `σ_pol`. It does **not** model the full nonlinear screening problem (which requires Lindhard or self-consistent DFT) and so should not be used to predict TER ratios more precisely than `~factor 2–3`. For our purposes — making the model physically self-consistent across electrode/IL choice and across BTO/BFO/AlScN polarization scales — it is sufficient.

---

## 4. Polarization-induced flat-band tilt across the FE

Inside the FE, the bound volume charge is zero (uniform `P`), but the bound *surface* charges produce a uniform "depolarizing" field if unscreened, and a residual tilt of the band edges if screening is finite. Treating the FE as a slab with surface charges `σ_pol,top = −P` at `x = 0` (top face, `z = b`) and `σ_pol,bot = +P` at `x = t_FE` (bottom face, `z = c`), and defining the local coordinate

```
x = z − b,    x ∈ [0, t_FE]                                            (Eq. 7)
```

the polarization-only contribution to the conduction-band edge (electron potential energy) inside the FE, in the limit of perfectly compensating electrodes, is

```
ΔU_pol(x) = − (σ_pol,bot / ε_FE) · q · x         [J]                  (Eq. 8)
          = − (P / ε_FE) · q · x                  for P > 0
```

with `ε_FE = ε0 · κ_FE` (use the **static** `κ_FE` here, since this is the bulk-FE response, not metal screening) [F/m]. With the sign of Eq. 8, `P > 0` produces a band edge that **decreases linearly with `x`** (drops from top to bottom), consistent with §1 (positive bound charge at the bottom face attracts electrons → lower electron potential energy on the bottom). This is the **gradient that enters the WKB tunneling integrand** as an additive term to the applied-field tilt:

```
U_FE(x) = U_FE^{(0)} − q E_FE x + ΔU_pol(x)       [J]                  (Eq. 9)
```

with `E_FE` from §5. Pantel & Alexe 2010, Eq. 7. Tsymbal & Kohlstedt 2006, Fig. 1.

---

## 5. Self-consistent fields with `σ_pol` present

Treat the IL and FE as a series capacitor stack with a sheet of bound charge `σ_pol,top = −P` at the IL/FE interface (the top FE face). Displacement-field continuity across that interface, including the bound-charge step, gives

```
D_IL − D_FE = σ_free,interface  (assumed 0 in the absence of interface traps)
ε_IL E_IL − ε_FE E_FE = σ_pol,top = −P        [C/m²]                  (Eq. 10)
```

(Sign: `D` jumps by the *free* surface charge; the bound charge appears as the difference between `D` and `ε0 E` is `P`, so writing the equation in terms of `ε_FE E_FE` already absorbs the bulk `P` of the FE — Eq. 10 is the correct boundary condition and reduces to the textbook form when `P → 0`.) [Sze & Ng, *Physics of Semiconductor Devices* (3rd ed.), Ch. 3, §3.2.2.]

Voltage-sum constraint (taking the top electrode at `+V_app`, bottom at 0, drops measured along `+z`):

```
V_app = E_IL · t_IL + E_FE · t_FE             [V]                     (Eq. 11)
```

Solving Eq. 10 and Eq. 11 simultaneously:

```
E_FE  = [ ε_IL · V_app + P · t_IL ] / [ ε_IL · t_FE + ε_FE · t_IL ]    [V/m]   (Eq. 12)
E_IL  = [ ε_FE · V_app − P · t_FE ] / [ ε_IL · t_FE + ε_FE · t_IL ]    [V/m]   (Eq. 13)
```

(Signs: `P > 0` adds to `E_FE` — positive `P` plus zero applied bias produces a residual depolarizing field of magnitude `P t_IL / (ε_IL t_FE + ε_FE t_IL)` along `+z` inside the FE, which is the standard incomplete-screening result. Conversely `P > 0` **subtracts** from `E_IL`. Confirm by plugging `V_app = 0, P > 0`: `E_FE > 0, E_IL < 0`, consistent with bound charges of `−P` at IL/FE pulling field lines from FE into IL.)

**`P → 0` limit.** Eq. 12–13 reduce to the polarization-free capacitive divider:

```
E_FE → ε_IL V_app / (ε_IL t_FE + ε_FE t_IL)
E_IL → ε_FE V_app / (ε_IL t_FE + ε_FE t_IL)                                    (Eq. 14)
```

i.e. `V_FE / V_IL = (ε_IL t_FE) / (ε_FE t_IL)` as required for series capacitors. [Sze & Ng Ch. 3, Eq. 51.]

---

## 6. Image-force lowering: which `κ`?

For an electron at distance `x` from a metal interface in a dielectric of relative permittivity `κ`, the image-charge potential energy is `−q² / (16 π ε0 κ x)` [J]. Combined with a uniform field `E` [V/m], the maximum-of-barrier lowering relative to the flat-band Schottky barrier is

```
Δφ_IF = sqrt( q³ E / (4 π ε0 κ_∞) )         [J]                       (Eq. 15)
```

(equivalently `sqrt(q E / (4 π ε0 κ_∞))` in [V] when divided by `q`).

**Why `κ_∞` (optical), not `κ_static`.** The image charge is induced by the **transiting carrier itself**. The carrier crosses the image-force region (`x_max ~ sqrt(q / (16 π ε0 κ E)) ~ 1 nm` at `E ~ 1 MV/cm`) in a transit time `τ_trans ~ x_max / v ~ 1 nm / 10⁶ m/s ~ 10⁻¹⁵ s = 1 fs`. This is far **shorter** than the lattice / ionic relaxation time (`~ 1/ω_TO ~ 10⁻¹³ s = 100 fs` for typical oxides). Therefore only the **electronic** response of the dielectric tracks the carrier, and the relevant permittivity is `κ_∞`, not `κ_static`. [Sze & Ng, *Physics of Semiconductor Devices* (3rd ed.), Ch. 3, §3.2.4 ("dynamic dielectric constant"); Mehta–Silverman–Jacobs 1973, §III.] For HfO\_x, `κ_∞ ≈ 4.0`; for AlScN, `κ_∞ ≈ 4.6` (authoritative values: `docs/parameter_audit_20260506.md`). Using `κ_static (~18–22 for HfO_x, ~16–18 for AlScN)` would *underestimate* `Δφ_IF` by `~ sqrt(κ_static / κ_∞) ≈ 2×` — a factor-of-2 sign-of-effect bug.

---

## 7. Equation-to-mechanism mapping (Wave-2 contract)

For each transport mechanism, the implementation must consume **exactly** the equations listed. The "injection side" is the electrode emitting electrons: top (Ti) for `V_app < 0`, bottom (Al) for `V_app > 0`.

| Mechanism | Consumes | Notes |
|---|---|---|
| **TE** (Thermionic Emission, Schottky) | `φ_B^{eff}` = `φ_B^{(0)} +` Eq. 5 (top inj.) or Eq. 6 (bot inj.); minus `Δφ_IF` from Eq. 15 with `κ_∞` of the **injection-side dielectric** (HfO\_x for top, AlScN for bot); `E` for image-force = `E_IL` (Eq. 13) for top, `E_FE` (Eq. 12) for bot. | Sze Ch. 3 Eq. 96. |
| **DT** (Direct Tunneling) | Eq. 5, 6 (both barrier heights); Eq. 9 for the barrier shape inside FE (and analog for IL); Eq. 12, 13 for tilts. WKB integrand uses `U(x) − E_F`. | Pantel–Alexe 2010 Eq. 8; Tsymbal–Kohlstedt 2006 Eq. 2. |
| **FN** (Fowler–Nordheim) | `φ_B^{eff}` from Eq. 5 (top inj.) or 6 (bot inj.); `E` = injection-side field, Eq. 12 or 13. No image-force in the standard FN form (already a triangular-barrier WKB). | Sze Ch. 8. |
| **PF** (Poole–Frenkel, bulk-trap-assisted) | `E_FE` from Eq. 12 only (PF lives in the FE bulk); trap depth `φ_T`; Frenkel lowering `Δφ_PF = sqrt(q³ E_FE / (π ε0 κ_∞,FE))` (note `π`, not `4π`, and `κ_∞,FE`). Polarization enters **only** through `E_FE`. | Sze Ch. 3 Eq. 137; Mehta–Silverman–Jacobs 1973 (PF in `κ_∞`). |
| **TAT** (Trap-Assisted Tunneling) | Eq. 5, 6 for endpoint barriers; Eq. 9, 12, 13 for the local potential at the trap site `x_T` inside FE; trap depth `φ_T` below FE CB. WKB legs use the same integrand as DT but split at `x_T`. | Standard two-step TAT, e.g. Pantel–Alexe 2010 §III.B. |

**Atomic-units boundary.** Equations 1–15 are in SI. The WKB / tunneling kernel in `src/` operates in a.u. (`ℏ = m_e = e = 1`, energies in Hartree, lengths in Bohr). Convert at the call boundary:
- Energy: 1 Ha = 4.3597447e−18 J = 27.2114 eV.
- Length: 1 Bohr = 5.29177e−11 m.
- Field: 1 Ha/(e·Bohr) = 5.14221e11 V/m.
The mass `m*` (FE/IL effective mass) enters as a dimensionless multiplier on the kinetic-energy term in a.u.

**Constants used (cite at first appearance in code).**
- `ε0 = 8.8541878128e−12 F/m` (CODATA 2018).
- `q = 1.602176634e−19 C` (exact, SI 2019).
- `κ_∞,HfOx ≈ 4.0` [Robertson, *Rep. Prog. Phys.* **69**, 327 (2006), Table 4]. Authoritative table is `docs/parameter_audit_20260506.md` row 14.
- `κ_∞,AlScN ≈ 4.6` [Ambacher et al., *J. Appl. Phys.* **130**, 045102 (2021), Table III; Deng et al., *Appl. Phys. Lett.* **102**, 112103 (2013)]. Authoritative table is `docs/parameter_audit_20260506.md` row 6.
- `κ_static,HfOx ≈ 18–22` (ALD amorphous; up to ~25 for crystallized monoclinic) [Robertson 2006]; `κ_static,AlScN ≈ 16–18` (Sc ≈ 0.3) [Fichtner 2019]. Do **not** use these in §3, §6 — listed only for comparison with `κ_∞`.
- `λ_TF,Ti ≈ 0.5–0.7 Å`, `λ_TF,Al ≈ 0.5 Å` [Ashcroft & Mermin Ch. 17; values used by Tsymbal–Kohlstedt 2006 Fig. 2].

---

**End of contract.** Wave-2 implementation must reference `(Eq. N)` numbers above when adding/modifying barrier and field expressions. Any deviation requires a follow-up memo.

---

## CHANGES (Wave-2 implementation map)

This section maps each contract equation to the file:line where the Wave-2
implementation lands. The equations themselves above are unchanged.

| Equation | What it is | Implemented at |
|---|---|---|
| Eq. 1 (`σ_pol,top = −P`) | Bound charge at top FE face | `src/transport_potential.py` `Potential.sigma_pol_top` (L151–160) |
| Eq. 2 (`σ_pol,bot = +P`) | Bound charge at bottom FE face | `src/transport_potential.py` `Potential.sigma_pol_bot` (L162–169) |
| Eq. 1, 2 propagation onto state | `σ_pol,top/bot` exposed on `BiasPointState` | `src/simulation_types.py` `BiasPointState.sigma_pol_top_au`, `sigma_pol_bot_au` (L37–38); populated in `src/transport_solver.py` `_build_state_snapshot` (~L120–127) |
| Eq. 4 (screening drop `ΔV`) | Folded into Eq. 5/6 helper | `src/transport.py` `_polarization_dipole_shift_ev` (L100–121) |
| Eq. 5 (`Δφ_B,top`) | Top barrier dipole shift, `−σ_pol,top·λ_TF·q/(ε0·κ_∞,Ti)` | `src/transport.py` `_barriers_ev` (L160–164 + L175–181); also audit_state mirror (~L460–474) |
| Eq. 6 (`Δφ_B,bot`) | Bottom barrier dipole shift, `−σ_pol,bot·λ_TF·q/(ε0·κ_∞,Al)` | `src/transport.py` `_barriers_ev` (L165–169 + L182–188); audit_state mirror (~L460–474) |
| Eq. 7 (FE local coordinate `x = z − b`) | `pos = x − il_end` inside `polarization_tilt_potential` | `src/transport_potential.py` `polarization_tilt_potential` (L274–300) |
| Eq. 8 (`ΔU_pol(x) = −(σ_pol,bot/ε_FE)·q·x`) | Polarization-only FE tilt | `src/transport_potential.py` `_electrostatic_state` `fe_dv_polarization` (L103–107); served via `polarization_tilt_potential` (L274–300) |
| Eq. 9 (`U_FE(x) = U^(0) − qE_FE_applied·x + ΔU_pol(x)`) | Composite WKB FE potential | `src/transport_potential.py` `total_potential` (L340–351) — sums `electrostatic_potential` (applied-only) + `polarization_tilt_potential` + `barrier_potential` + `wf_potential` |
| Eq. 10 (D-continuity with bound charge) | Implicit in displacement-divider `sigma_s`; the split into applied vs. polarization preserves it | `src/transport_potential.py` `_electrostatic_state` (L67–146) |
| Eq. 11 (voltage-sum constraint) | `denominator` and `sigma_s` solution | `src/transport_potential.py` `_electrostatic_state` (L70–81) |
| Eq. 12 (`E_FE` self-consistent) | `total_e_field_fe = (fe_dv_electrostatic_applied + fe_dv_polarization)/t_FE` | `src/transport_potential.py` `_electrostatic_state` (L125–127) |
| Eq. 13 (`E_IL` self-consistent) | `total_e_field_il = il_dv_electrostatic / t_IL` | `src/transport_potential.py` `_electrostatic_state` (L120–124) |
| Eq. 14 (`P → 0` divider, applied-only) | `fe_dv_electrostatic_applied = sigma_s/(ε_FE·ε0)·t_FE` | `src/transport_potential.py` `_electrostatic_state` (L98–102) |
| Eq. 15 (image-force `Δφ_IF` with κ_∞) | `_schottky_lowering_ev(field, κ_∞)`; callers pass `insulator_eps_inf`/`fe_eps_inf` | `src/transport.py` `_schottky_lowering_ev` (L67–80); call sites L172–174, L179, L186 |
| §7 PF row (`Δφ_PF` with `π` and κ_∞) | `_poole_frenkel_lowering_ev` (π, κ_∞,FE) | `src/transport.py` `_poole_frenkel_lowering_ev` (L82–98); used in `_poole_frenkel_component` (L314) |
| §7 PF row (PF in IL uses `insulator.trap_depth`) | `default_trap_depth_au = self.diode.insulator_trap_depth` for `region == "il"` | `src/transport.py` `_poole_frenkel_component` (L294–298); plumbing in `src/transport_fed.py` `insulator_trap_depth` (~L55–60) |
| Sign of E_FE through PF/TAT/SCLC/DT | `_field_sign(state.fe_field_au, ...)`; signed `field_v_per_m` in PF | `src/transport.py` `_tunneling_current` (L278), `_poole_frenkel_component` (L304–323), `_trap_assisted_tunneling_current` (L342–360), `_sclc_current` (L363–382) |

**Material parameter additions (parameter_audit_20260506.md punch list):**

| Punch | Field | Value | File:line |
|---|---|---|---|
| #1 | AlScN χ | 2.0 eV | `src/materials.py` `Materials.alscn` (L178); `Materials.baseline_alscn` (L74) |
| #2 | AlScN κ_∞ | 4.6 | `src/materials.py` `Materials.alscn` (L186); `Materials.baseline_alscn` (L80) |
| #3 | HfO_x κ_∞ | 4.0 | `src/materials.py` `Materials.hfo2` (L197); `Materials.baseline_hfo2` (L90) |
| #4 | HfO_x trap depth | 1.0 eV | `src/materials.py` `Materials.hfo2` (L201); `Materials.baseline_hfo2` (L92) |
| #5 | Al2O3 χ | 1.2 eV | `src/materials.py` `Materials.al2o3` (L161); `Materials.baseline_al2o3` (L59) |
| #6 | Al2O3 κ static | 9.0 | `src/materials.py` `Materials.al2o3` (L158) |
| #7 | Al2O3 κ_∞ | 3.1 | `src/materials.py` `Materials.al2o3` (L167); `Materials.baseline_al2o3` (L65) |
| #8 | Al2O3 trap depth | 1.4 eV | `src/materials.py` `Materials.al2o3` (L170); `Materials.baseline_al2o3` (L68) |
| eps_inf field | Insulator/Ferroelectric class | optional kwarg | `src/material_types.py` `Insulator` (L23), `Ferroelectric` (L40) |
| trap_depth on Insulator | Insulator class | optional kwarg | `src/material_types.py` `Insulator.__init__` (L23) |

**Notes.**
- Metal κ_∞ falls back to the existing electrode `k` field (kept at 1 for Ti/Al/Pd — Drude tail). To override, add an `eps_inf` attribute on `MetalElectrode` and `transport_fed.py` will pick it up via the `getattr` in `top_eps_inf`/`bottom_eps_inf` (L40–41).
- `screening_charge_au` (Eq. 4 metal-side total screening sheet) is unchanged; the new `sigma_pol_top_au`/`sigma_pol_bot_au` on `BiasPointState` are the FE-side bound-charge sheets used by Eq. 5/6.
