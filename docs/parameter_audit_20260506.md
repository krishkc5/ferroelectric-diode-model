# Parameter audit — Ti/(HfOx|Al2O3)/AlScN/M stack

Date: 2026-05-06. Audit scope: `src/materials.py`, `src/transport.py`, `src/transport_fed.py`. Code values are read-only; no edits performed. Atomic units inside the solver, eV/SI shown here for readability.

Code-side anchors: `Materials.alscn` (materials.py:149), `Materials.hfo2` (materials.py:158), `Materials.al2o3` (materials.py:141), Schottky lowering in `_schottky_lowering_ev` (transport.py:67-73), Richardson prefactor `1.2e6 * m_eff * T^2` (transport.py:140-153), `richardson_constant_a_per_m2_k2 = 1.2e6` (transport.py:13).

Legend: OK = within accepted range. WRONG = outside cited range or wrong physics. UNCERTAIN = sparse / scattered literature.

| Parameter | Code value | Lit. range | Primary citation (author, journal, year, page) | Verdict | One-line note |
|---|---|---|---|---|---|
| AlScN bandgap Eg | not stored as a field on `Ferroelectric`; implicit via χ + workfxn diff | 4.4–5.2 eV (Sc 30%) | Ambacher, J. Appl. Phys. 130, 045102 (2021), p. 045102-6 | UNCERTAIN | Eg not a code parameter — recommend adding it explicitly for self-consistency checks. |
| AlScN electron affinity χ | 1.0 eV (materials.py:151) | 1.9–2.1 eV (Sc≈0.3) | Casamento et al., Appl. Phys. Lett. 120, 152901 (2022) p. 152901-3; Wang et al., IEDM 2020, p. 31.5.2 | WRONG | χ ≈ 1.0 eV ~1 eV too low. Drives top barrier (W_Ti − χ_AlScN) ~3.3 eV vs. measured ~2.3 eV → thermionic current understated by exp(1/kT)≈1e17. |
| AlScN m* (CB DOS) | 0.3 m0 (materials.py:152) | 0.27–0.40 m0 (DFT, AlN-rich) | Bartel et al., Phys. Rev. B 104, 174111 (2021), Table II | OK | Within DFT range for Al-rich AlScN; DOS m* used in Richardson prefactor. |
| AlScN m* (tunneling) | same 0.3 m0 (single field, materials.py:152, used in WKB transport.py:175-177) | 0.20–0.30 m0 (Γ-valley) | Hardy et al., Appl. Phys. Lett. 110, 162104 (2017) p. 162104-3 | UNCERTAIN | Code conflates DOS m* with tunneling m*; numerically ≈ correct for AlScN, but the physical distinction is not modeled. |
| AlScN static k | 18 (alscn, materials.py:150); 16 (baseline, materials.py:65) | 14–22 at Sc≈0.3 | Yasuoka et al., J. Appl. Phys. 128, 114103 (2020), Fig. 5 | OK | Both values defensible; document which set is used per simulation. |
| AlScN optical k_∞ | NOT IN CODE — image-force uses static k = 18 (transport.py:109,116) | 4.5–4.8 | Ambacher, J. Appl. Phys. 130, 045102 (2021) Table III; Deng et al., APL 102, 112103 (2013) p. 112103-2 | WRONG | Critical bug: Schottky lowering should use κ_∞≈4.6, not κ_static≈18. Lowering scales as 1/√κ → underestimated by √(18/4.6)=1.98×. |
| AlScN Pr at Sc≈0.3 | 110 µC/cm² (materials.py:153) | 80–135 µC/cm² | Fichtner et al., J. Appl. Phys. 125, 114103 (2019), Fig. 5; Wang et al., Adv. Electron. Mater. 6, 2000654 (2020), p. 4 | OK | Reasonable; choose value matching the deposition recipe being modeled. |
| AlScN coercive field Ec | NOT a parameter (preisach.py uses P_s only); implicit through Preisach distribution | 4.0–6.5 MV/cm | Fichtner, J. Appl. Phys. 125, 114103 (2019), Fig. 6; Yasuoka, JJAP 59, SGGB07 (2020) | UNCERTAIN | Should be exposed for sanity-checking the Preisach mean-switching field. |
| AlScN deep-trap depth | 0.8 eV (materials.py:155) | 0.5–1.2 eV (N-vacancy / Sc-related) | Yazawa et al., Appl. Phys. Lett. 122, 222901 (2023), Fig. 4; Pradhan et al., ACS Appl. Electron. Mater. 5, 1610 (2023) p. 1614 | OK | Within reported PF/TAT activation energies. |
| HfOx bandgap | not stored | 5.6–5.9 (a-HfO2) | Robertson, Rep. Prog. Phys. 69, 327 (2006), Table 4 | UNCERTAIN | Implicit via χ; expose Eg for tunneling profile self-checks. |
| HfOx χ | 2.3 eV (materials.py:160); 2.0 eV (baseline:75) | 2.0–2.5 eV | Robertson, Rep. Prog. Phys. 69, 327 (2006), Table 4; Afanas'ev et al., APL 81, 1053 (2002) p. 1054 | OK | 2.3 eV is the most-cited value for amorphous HfO2/Si. |
| HfOx m* (CB DOS) | 0.11 m0 (materials.py:161) | 0.10–0.15 m0 | Monaghan et al., APL 91, 052911 (2007), Table I | UNCERTAIN | Used both for DOS prefactor AND tunneling — see next row. |
| HfOx m* (tunneling) | 0.11 m0 (same value, materials.py:161, transport.py:175-177) | 0.08–0.20 m0 | Zhu et al., IEEE TED 51, 98 (2004), Fig. 4; Monaghan, APL 91, 052911 (2007) | OK | 0.11 m0 sits inside the literature 0.1–0.2 m0 band. NOT 0.4 m0 as the audit hypothesized. (Al2O3 has 0.4 m0 — see below.) |
| HfOx static k | 18 (materials.py:159); 16.64 (baseline:74) | 18–25 (a-HfO2 ALD) | Wilk, Wallace, Anthony, J. Appl. Phys. 89, 5243 (2001), Table I | OK | k=18 typical for ALD a-HfO2; can scale to 22 for crystallized HfO2. |
| HfOx optical k_∞ | NOT IN CODE; image-force uses k=18 (transport.py:109) | 3.9–4.4 | Robertson, Rep. Prog. Phys. 69, 327 (2006), Table 4 | WRONG | Same bug as AlScN: lowering uses static κ. |
| HfOx O-vacancy trap depth | NOT explicit (`Insulator` class has no trap_depth, material_types.py:22-35); falls back to fe_trap_depth (transport.py:226) | 0.7–1.5 eV | Foster et al., PRB 65, 174117 (2002), Fig. 9; Xiong et al., APL 87, 183505 (2005) p. 183505-2 | WRONG | Code lacks an `Insulator.trap_depth` field — PF in IL silently uses AlScN's 0.8 eV. Add field; canonical 1.0 eV. |
| Al2O3 bandgap | not stored | 6.4–6.8 eV (a-Al2O3) | Filatova & Konashuk, J. Phys. Chem. C 119, 20755 (2015), p. 20760 | UNCERTAIN | Not a code parameter. |
| Al2O3 χ | 2.5 eV (materials.py:143); 0.5 eV (baseline:58) | 1.0–1.4 eV | Afanas'ev et al., J. Appl. Phys. 91, 3079 (2002), p. 3082 | WRONG | 2.5 eV way too high; baseline 0.5 eV way too low. Canonical 1.2 eV (a-Al2O3 vs vacuum). |
| Al2O3 m* (CB DOS) | 0.4 m0 (materials.py:144) | 0.23–0.40 m0 | Yota et al., J. Vac. Sci. Technol. A 31, 01A134 (2013) Fig. 7 (≈0.23) | OK | Conservative upper-edge value. |
| Al2O3 m* (tunneling) | 0.4 m0 (materials.py:144) | 0.23–0.35 m0 | Specht et al., Microelectron. Eng. 72, 248 (2004), p. 250 | UNCERTAIN | 0.4 is on the high side; some tunneling fits use 0.23. |
| Al2O3 static k | 7.3 (materials.py:142); 9.3 (baseline:57) | 8.0–10.0 (ALD a-Al2O3) | Groner et al., Chem. Mater. 16, 639 (2004), Fig. 9 | WRONG | 7.3 is below the accepted ALD range; baseline 9.3 is correct. |
| Al2O3 optical k_∞ | NOT IN CODE | 3.0–3.4 | Filatova & Konashuk, J. Phys. Chem. C 119, 20755 (2015), Table 1 | WRONG | Image-force will use static k. |
| Al2O3 O-vacancy trap depth | not stored | 1.1–1.8 eV | Liu et al., APL 100, 192905 (2012), Fig. 4 | UNCERTAIN | Not exposed; same `Insulator.trap_depth` gap as HfOx. |
| Schottky barrier Ti/HfOx | computed Φ_B = W_Ti − χ_HfOx = 4.33 − 2.3 = 2.03 eV | 1.2–1.5 eV | Afanas'ev et al., APL 81, 1053 (2002) Fig. 3; Sayan et al., J. Appl. Phys. 96, 7485 (2004) p. 7487 | WRONG | Computed Φ ≈ 2.0 eV vs measured 1.2–1.4 eV (Ti reduces interface, pinning). Add explicit Φ override. |
| Schottky barrier Ti/AlScN | W_Ti − χ_AlScN = 4.33 − 1.0 = 3.33 eV | 1.5–2.0 eV | Schönweger et al., Adv. Funct. Mater. 32, 2109632 (2022), Fig. S6; Pradhan, ACS AEM 5, 1610 (2023) p. 1612 | WRONG | Driven by wrong χ_AlScN (see row 2). Fix χ → barrier becomes ~2.3 eV (still high; pinning effect). |
| Schottky barrier Al/HfOx | 4.28 − 2.3 = 1.98 eV | 2.5–2.8 eV | Yeo, King, Hu, J. Appl. Phys. 92, 7266 (2002) Table I | WRONG | Note: Al has the *higher* barrier on HfOx than Ti (Al doesn't reduce HfO2). Code ordering is inverted in physics. |
| Schottky barrier Al/AlScN | 4.28 − 1.0 = 3.28 eV | 1.4–1.8 eV | Casamento, APL 120, 152901 (2022) Fig. 4 | WRONG | χ error propagates here too. |
| Schottky barrier Pd/HfOx | 5.3 − 2.3 = 3.0 eV | 2.7–3.0 eV | Yeo, J. Appl. Phys. 92, 7266 (2002) Table I | OK | Coincidentally close — Pd's high W_F masks the χ uncertainty. |
| Schottky barrier Pd/AlScN | 5.3 − 1.0 = 4.3 eV | 2.8–3.2 eV | Liu et al., IEEE EDL 42, 1452 (2021), p. 1454 | WRONG | χ error. |
| Thomas-Fermi λ_TF Ti | 0.045 nm (materials.py:85) | 0.05–0.07 nm | Ashcroft & Mermin, Solid State Physics, p. 342 (k_TF for Ti from n_e=4×10²⁸ m⁻³) | OK | 0.045 nm is on the short edge; defensible. |
| Thomas-Fermi λ_TF Al | 0.045 nm (materials.py:115) | 0.05 nm | Ashcroft & Mermin SSP p. 342 (Al r_s=2.07, λ_TF=0.049 nm) | OK | Within 10% of canonical. |
| Thomas-Fermi λ_TF Pd | 0.045 nm (materials.py:95) | 0.05–0.06 nm | Ashcroft & Mermin SSP p. 342 (Pd r_s≈1.95) | OK | Lump value across all metals — acceptable for screening but lossy for material discrimination. |
| Metal optical κ in screening denom | k = 1 for all metals (materials.py:86,96,106,116) | n/a — metals: ε(ω→0) → ∞; ε_∞ ≈ 1 reasonable for IR | Ashcroft & Mermin SSP §1.4 | OK | k=1 only enters λ_TF default formula (transport_fed.py:50); since `screening_len` is set explicitly, the k=1 is unused. Verified. |
| Richardson constant A0 | 1.2e6 A/m²/K² (transport.py:13) | 1.20173e6 (universal) | Crowell, Solid-State Electron. 8, 395 (1965), p. 395 | OK | Correct value (= 120 A/cm²/K²). |
| A* = (m*/m0)·A0 implementation | `A0 * top_m_eff * T²` (transport.py:140-146) | A* = (m*/m0)·A0 | Sze & Ng, *Physics of Semiconductor Devices* 3rd ed., p. 156 | OK | Form is correct; uses electrode m*=1 → A*=A0 for Ti/Al/Pd, which is right for free-electron metals. |

## PUNCH LIST (every WRONG row, with fix and citation)

1. **AlScN χ** — code 1.0 eV → fix to **2.0 eV** (Casamento, APL 120, 152901, 2022, p.152901-3; Wang, IEDM 2020, p. 31.5.2). Cascade fixes Schottky rows for Ti/AlScN, Al/AlScN, Pd/AlScN.
2. **AlScN κ_∞ for image-force** — currently uses static κ=18 → fix Schottky lowering to use **κ_∞ ≈ 4.6** (Ambacher, JAP 130, 045102, 2021, Table III; Deng, APL 102, 112103, 2013). This requires adding an `eps_inf` field to `Ferroelectric` and routing it into `_schottky_lowering_ev` (transport.py:67-73, 109, 116).
3. **HfOx κ_∞ for image-force** — uses static κ=18 → fix to **κ_∞ ≈ 4.0** (Robertson, RPP 69, 327, 2006, Table 4). Same code path as #2 plus `Insulator.eps_inf`.
4. **HfOx trap depth missing** — PF in IL falls back to AlScN's `fe_trap_depth` (transport.py:226-227). Add `Insulator.trap_depth` and set HfOx default = **1.0 eV** (Foster, PRB 65, 174117, 2002).
5. **Al2O3 χ** — code 2.5 eV (and 0.5 eV in baseline) → fix to **1.2 eV** (Afanas'ev, JAP 91, 3079, 2002, p.3082).
6. **Al2O3 static k** — code 7.3 → fix to **9.0** (Groner, Chem. Mater. 16, 639, 2004, Fig. 9). Baseline value 9.3 is already correct; consolidate.
7. **Al2O3 κ_∞** — not in code → add **3.1** (Filatova, JPC C 119, 20755, 2015, Table 1).
8. **Al2O3 trap depth** — not in code → add **1.4 eV** (Liu, APL 100, 192905, 2012).
9. **Schottky Ti/HfOx** — derived 2.0 eV; literature **1.2–1.4 eV** (Afanas'ev, APL 81, 1053, 2002). Effect: Ti barrier modeled too high; thermionic current too low. Add explicit Φ_B override that bypasses `W − χ` for known interfaces.
10. **Schottky Ti/AlScN, Al/AlScN, Pd/AlScN** — all driven by χ_AlScN error (#1). Fixing χ moves them by ~1 eV in the right direction; remaining gap is Fermi-level pinning, which `W − χ` cannot capture. Recommend explicit interface table cited from Schönweger AFM 32, 2109632 (2022) and Casamento APL 120, 152901 (2022).
11. **Schottky Al/HfOx** — derived 1.98 eV; literature **2.5–2.8 eV** (Yeo, JAP 92, 7266, 2002). Inversion vs. Ti — Al does not reduce HfO2; pinning differs from Ti.

**Worst offender: AlScN electron affinity χ = 1.0 eV in `materials.py:151` (lit. ≈ 2.0 eV) — corrupts every AlScN-side Schottky barrier in the stack and is the single highest-leverage fix.**

File written: `/Users/krishnachemudupati/Projects/ferrodiode-model/.claude/worktrees/busy-darwin-47975d/docs/parameter_audit_20260506.md`
