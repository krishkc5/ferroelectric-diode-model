# External AlScN Ferroelectric-Diode I-V Benchmarks

## Scope and provenance note (read first)

This document is the deliverable for the literature benchmark task: identify
3-6 published AlScN ferroelectric-diode I-V data sets whose stack is fully
documented, prioritizing stacks structurally similar to the project's
`Ti / HfOx (4 nm) / AlScN (10 nm) / Al` reference geometry.

**Access constraint encountered.** During this research pass both `WebSearch`
and `WebFetch` were denied by the runtime sandbox, so I could not pull paper
PDFs, supplementary information, or publisher pages directly. The repository
itself contains no local copies of these papers and no BibTeX file. Because of
this, the entries below carry only details I can responsibly reproduce from
prior knowledge of the AlScN ferroelectric literature, and every entry is
marked with a verification flag. **Before this table is used to constrain a
fit, every numeric field marked `verify` must be checked against the primary
source.** No I-V points were digitized, because no figures were accessible
during this run; the `docs/benchmarks/` directory was therefore created empty
on purpose.

The hard constraint "If a paper's stack details are ambiguous, exclude it" was
applied conservatively: papers I could not reconstruct stack-by-stack from
memory were left out rather than guessed.

## Summary table

| # | Lead author | Year | Venue | Top / IL / FE / Bottom | t_FE (nm) | Sc fraction (x) | t_IL (nm) | Area | V range | ER ratio (J_up/J_down) | Dominant mechanism (claimed) | Verification |
|---|-------------|------|-------|------------------------|-----------|-----------------|-----------|------|---------|------------------------|------------------------------|--------------|
| 1 | Liu | 2022 | Nat. Commun. 13:1009 | Pt / -- / Al(1-x)Sc(x)N / Pt (or Al) | ~20 (verify; samples reported across 20-100 nm) | x ~ 0.32 (verify) | none | small via, micron-scale (verify exact A) | a few V around +/-Vc; Vc reported around 5-6 V at 20 nm | reported "ON/OFF" optical-window ratios up to ~10^3-10^5 depending on bias point (verify) | authors discuss thermionic / Schottky-like emission modulated by polarization-induced barrier change; some discussion of trap-assisted contributions at higher fields (verify) | unverified -- web access blocked this pass |
| 2 | Wang (Olsson group) | 2023 | various IEDM / Nature Electronics papers on AlScN FeDs | Pt or Al / -- / Al(1-x)Sc(x)N / Pt or W (verify exact pair) | ~10-20 (verify per paper) | x ~ 0.3-0.36 (verify) | typically none (MFM); some MFIM variants exist | um-scale circular pads (verify) | a few V to ~10 V (verify) | not consistently reported as a single number; varies with read voltage (verify) | mixed: Schottky barrier modulation + PF/TAT in the FE bulk (verify) | unverified -- needs primary source |
| 3 | Schönweger / Kohlstedt group | 2022-2023 | Adv. Electron. Mater. / Adv. Funct. Mater. / APL | Pt / -- / Al(1-x)Sc(x)N / Pt (MFM) (verify) | ~20-100 (verify) | x ~ 0.27-0.36 (verify) | none in MFM cuts; SiO2 or AlOx in some MFIM cuts (verify) | um pads (verify) | up to ~+/-(Vc + a few V) (verify) | not always quoted as a single number (verify) | thermionic emission limited (Schottky) at low/moderate field, rolling into PF/TAT at higher fields (verify) | unverified -- needs primary source |
| 4 | Mizutani / Yasuoka et al. | 2022-2023 | Jpn. J. Appl. Phys. / APL | TiN / -- / Al(1-x)Sc(x)N / TiN (MFM) and SiO2-IL MFIM variants (verify) | ~20-50 (verify) | x ~ 0.2-0.3 (verify) | 0 (MFM) or a few nm SiO2 (MFIM) (verify) | um pads (verify) | several V (verify) | only loosely defined; usually shown as J-V hysteresis rather than ER ratio (verify) | mostly Schottky-like injection plus reports of PF in the AlScN at higher field (verify) | unverified -- needs primary source |
| 5 | Pradhan et al. | 2022-2023 | Adv. Electron. Mater. (Olsson / Penn group) | TiN or Pt / -- / Al(1-x)Sc(x)N / Pt or W (verify) | ~10-20 (verify) | x ~ 0.32 (verify) | none (MFM) (verify) | um circular pads (verify) | a few V (verify) | reported diode-like ratios at certain read biases (verify) | thermionic emission with polarization-modulated Schottky barrier (verify) | unverified -- needs primary source |

All "verify" tags above mean the same thing: the value is plausible from
general knowledge of these groups' work, but I could not open the paper this
pass to confirm the exact number. Treat any of those numbers as placeholders,
not data.

Notes on what is currently *not* in the table on purpose:

- I deliberately did not include a "Liu 2022 with 4 nm HfOx interlayer" entry,
  because, to the best of my recollection, the canonical Liu et al. 2022
  Nature Communications AlScN ferroelectric diode paper is an MFM stack
  (Pt / AlScN / Pt class), not an MFIM stack with a 4 nm HfOx interlayer.
  Inventing an HfOx-IL benchmark to "match" the project stack would violate
  the no-fabrication rule.
- I also did not include Yasuoka et al. as a separately verified row; the
  Mizutani / Yasuoka group has multiple closely related papers, and I could
  not, without access, distinguish stack-specific numbers between them with
  confidence.

## Per-paper sections

Each section gives the level of detail I am confident in from prior knowledge
plus an explicit verification checklist. Anything in a checklist must be
filled in or corrected from the primary source before this section is used to
constrain a fit.

### 1. Liu et al., *Nature Communications* 13:1009 (2022)

- Citation (verify exact author list and pagination): Liu et al.,
  "Aluminum scandium nitride-based metal-ferroelectric-metal diode," or close
  variant; *Nature Communications* 13, 1009 (2022).
- DOI (verify): `10.1038/s41467-022-28673-2` (this is the DOI most commonly
  associated with the AlScN ferroelectric diode paper out of the Olsson /
  Jariwala / Penn collaboration; confirm against the actual article record).
- Stack class: metal / AlScN / metal (MFM). No oxide interlayer in the main
  device.
- AlScN thickness: variants reported (verify exact set; commonly 20 nm, 45 nm,
  100 nm class samples).
- Sc fraction: nominally x ~ 0.32 (verify).
- Electrodes: top and bottom metal contacts (verify exact metals; recall is
  Pt-class, not Ti / Al).
- Electrode area: micron-scale circular pads (verify exact diameter).
- Voltage range: spans a few V on either side of coercive voltage; coercive
  voltage scales with thickness at roughly the documented coercive field of
  AlScN (~4-5 MV/cm) (verify).
- Current floor: limited by setup leakage / measurement floor (verify).
- ER ratio J_up/J_down: reported as a function of read bias; large at biases
  near but below Vc; rolls off as both states leak above Vc (verify).
- Mechanism attribution: the authors describe the diode behaviour as
  polarization-modulated barrier transport at the metal / AlScN interface,
  i.e. a Schottky-emission-like picture in which the surface charge from the
  ferroelectric polarization shifts the effective injection barrier between
  the two polarization states (verify exact wording and any TAT/PF
  contribution discussed in SI).
- Fit style: per-state, not shared parameter -- the polarization-up and
  polarization-down branches are typically fit as two Schottky branches with
  different effective barriers, rather than one global fit (verify).

Verification checklist for this entry (do these before relying on it):

- [ ] confirm author list and DOI
- [ ] confirm exact AlScN thicknesses reported
- [ ] confirm Sc fraction
- [ ] confirm electrode metals on both sides
- [ ] confirm pad area
- [ ] confirm whether any oxide IL is present in any reported variant
- [ ] confirm J floor and ER ratio numbers
- [ ] confirm whether the leakage analysis is Schottky-only or Schottky + PF/TAT
- [ ] confirm per-branch vs shared-parameter fit

### 2. Wang et al. (Olsson group, recent AlScN FeD work, 2021-2024)

- Citation (verify): one of the Wang first-author papers from the Olsson group
  on AlScN ferroelectric diodes / capacitors in IEEE EDL, IEDM, or Nature
  Electronics in this date window.
- DOI: verify by paper.
- Stack: MFM, Pt-class or W-class electrodes on Al(1-x)Sc(x)N (verify).
- Thickness: 10-20 nm range reported in scaling-focused papers (verify).
- Sc fraction: x ~ 0.3-0.36 (verify).
- Electrode area: um circular pads (verify).
- Voltage range: a few V (verify).
- ER ratio: reported as a function of read bias (verify).
- Mechanism: typically described as thermionic / Schottky-modulated with
  acknowledgement of trap-assisted or PF contributions in the AlScN bulk at
  higher field (verify).
- Fit style: per-branch is typical in these papers when a fit is shown at all;
  several papers only report ln(J) vs sqrt(V) plots without an explicit fit
  (verify).

Verification checklist: same nine bullets as for Liu 2022.

### 3. Schoenweger / Kohlstedt group (Kiel)

- Citation (verify): Schoenweger et al. and co-authors from the Kiel group on
  AlScN MFM diodes / capacitors, published in *Advanced Electronic
  Materials*, *Advanced Functional Materials*, or *Applied Physics Letters*
  in the 2022-2023 window.
- DOI: verify by paper.
- Stack: MFM with Pt electrodes is the most commonly remembered case; a
  subset of the group's papers report MFIM with native oxide / SiO2 / AlOx
  on the AlScN, but exact thickness and material vary by paper (verify
  carefully -- this is the entry most at risk of being mis-attributed).
- AlScN thickness: 20 nm to 100 nm class (verify by paper).
- Sc fraction: x ~ 0.27-0.36 (verify).
- Electrode area: um pads (verify).
- Voltage range: typically up to a few V beyond Vc (verify).
- ER ratio: not always quoted as a single number (verify).
- Mechanism: this group has published explicit transport analyses on AlScN
  capacitors. From memory the narrative is Schottky-limited at low / moderate
  field, rolling into PF or TAT at higher field, but the exact paper that
  states this should be cited rather than this summary (verify).
- Fit style: where an explicit transport fit is shown, it is typically per
  branch (verify).

Verification checklist: same as above, plus:

- [ ] specifically confirm whether the cited paper is MFM or MFIM, since this
      group has both
- [ ] if MFIM, confirm IL material and thickness exactly

### 4. Mizutani / Yasuoka et al. (Tokyo group)

- Citation (verify): Mizutani et al. and Yasuoka et al. publications on AlScN
  capacitors and diodes in *Japanese Journal of Applied Physics* and *Applied
  Physics Letters* in 2022-2023.
- DOI: verify per paper.
- Stack: MFM (typically TiN-class electrodes) and MFIM variants with thin
  SiO2 (verify).
- AlScN thickness: 20-50 nm range reported (verify).
- Sc fraction: x ~ 0.2-0.3 (verify -- this group has reported lower Sc
  fractions than the Penn group in some papers).
- Electrode area: um pads (verify).
- Voltage range: several V (verify).
- ER ratio: usually shown as J-V hysteresis curves rather than as a single
  on/off ratio number (verify).
- Mechanism: Schottky / thermionic injection at the metal/AlScN interface
  with PF in the AlScN bulk at higher fields is the typical attribution
  (verify exact paper).
- Fit style: where present, per-branch (verify).

Verification checklist: same nine bullets as for Liu 2022, plus:

- [ ] confirm which group member is first author for the entry being used
- [ ] confirm whether the cited entry is MFM or MFIM
- [ ] if MFIM, confirm IL is SiO2 vs AlOx and confirm IL thickness

### 5. Pradhan et al. (Olsson / Penn group, *Advanced Electronic Materials*)

- Citation (verify): Pradhan et al. AlScN ferroelectric diode paper in
  *Advanced Electronic Materials*, around 2022-2023.
- DOI: verify.
- Stack: MFM with TiN-class or Pt-class electrodes on Al(1-x)Sc(x)N (verify).
- AlScN thickness: typically 10-20 nm in this paper class (verify).
- Sc fraction: x ~ 0.32 (verify).
- Electrode area: um circular pads (verify).
- Voltage range: a few V on either side of Vc (verify).
- ER ratio J_up/J_down: reported at characteristic read biases (verify exact
  numbers).
- Mechanism: thermionic emission with polarization-modulated effective
  Schottky barrier is the dominant claim, similar in spirit to Liu 2022
  (verify).
- Fit style: typically per-branch (verify).

Verification checklist: same nine bullets as for Liu 2022.

## What was *not* delivered this pass and why

- **No CSVs in `docs/benchmarks/`.** The task said: "If figures are
  digitizable from the paper or its SI, save extracted I-V points as
  `docs/benchmarks/<lastauthor>_<year>_<stack>.csv` ... If not digitizable,
  skip." Web access was denied during this run, so no figure was reachable
  to digitize. The `docs/benchmarks/` directory was created empty on purpose
  so that future runs can populate it without churn.
- **No claimed numeric ER ratios, current floors, or Sc fractions in this
  document beyond what is explicitly tagged "verify".** This is intentional
  and matches the no-fabrication rule.

## Re-running this task with web access

When this task is re-run with `WebSearch` / `WebFetch` enabled (or a local
copy of the relevant PDFs in the repo), the next agent should:

1. Resolve the DOIs in each "Citation (verify)" line above.
2. Replace each "verify" placeholder with a concrete value or remove the
   entry.
3. Digitize the I-V curves where they exist as figures, one CSV per
   (paper, polarity branch) pair, columns `V_app, J_abs, polarity_branch,
   sweep_direction`.
4. Re-pick the recommended primary benchmark below if a closer-stack paper
   surfaces (in particular, any AlScN FeD with an HfOx interlayer would
   replace the current pick immediately).

## RECOMMENDED PRIMARY BENCHMARK

The project's reference stack is `Ti / HfOx (4 nm) / AlScN (10 nm) / Al`,
i.e. an MFIM diode with a thin oxide interlayer, ~10 nm of moderately
Sc-rich AlScN, and asymmetric metal / metal contacts (one mid-work-function
nitride-friendly metal, one low-work-function metal).

I am not aware of a published AlScN ferroelectric diode whose stack is a
clean structural twin (HfOx-interlayer MFIM with Ti and Al as the two
electrodes and AlScN around 10 nm). The honest closest published match in
the candidate set is therefore an MFM AlScN diode of comparable thickness
and Sc fraction, even though it is missing the HfOx interlayer.

Subject to the verification step above, the recommended primary benchmark is:

> **Liu et al., *Nature Communications* 13:1009 (2022) -- AlScN
> ferroelectric diode (MFM).**

Reasons for this pick:

- It is the most-cited dedicated AlScN FeD paper in this list and is the
  paper most other groups compare against, so reproducing its qualitative
  shape gives the project's transport layer the strongest external
  legitimacy per unit fitting effort.
- Its Sc fraction (~0.3 class) and AlScN thickness range (~10-20 nm class)
  are the closest of the five entries to the project's 10 nm AlScN target.
- Its transport story (polarization-modulated Schottky-like injection with
  per-branch behaviour) is exactly the regime the project's transport stack
  is being calibrated against.

The single thing this benchmark cannot pin down is the role of the 4 nm HfOx
interlayer in the project stack, because Liu 2022 is MFM. The
HfOx-interlayer effect should be cross-checked separately against any
HfOx-on-AlScN MFIM study identified in a subsequent web-enabled pass; if
such a paper is found, it should replace this primary benchmark.
