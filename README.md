# Ferroelectric Diode Preisach Model

This repository contains the active multidomain Preisach hysteresis workflow for
metal / interlayer / AlScN / metal stacks.

The current active model is:

- multidomain
- history-dependent
- electrostatically self-consistent
- quasistatic
- solved with a bisection-based fixed-point solver

The main output of the active workflow is a family of ferroelectric
polarization--voltage hysteresis loops across interlayer thickness, interlayer
material, electrode choice, and AlScN thickness.

## Provenance

The repository began from a mentor reference implementation.

Reference attribution:

- Zirun Han's reference implementation: [danny1078/fed_multidomain_presaich](https://github.com/danny1078/fed_multidomain_presaich)

The active root no longer mirrors that older tree exactly. It has been trimmed
and refocused around the current bisection-based multidomain hysteresis path.

## What Is Active

The active source code now lives under [`src`](src):

- [`src/atomicunits.py`](src/atomicunits.py)
- [`src/material_types.py`](src/material_types.py)
- [`src/materials.py`](src/materials.py)
- [`src/fed.py`](src/fed.py)
- [`src/potential.py`](src/potential.py)
- [`src/preisach.py`](src/preisach.py)
- [`src/self_consistent_solver.py`](src/self_consistent_solver.py)
- [`src/plot_il_hysteresis.py`](src/plot_il_hysteresis.py)
- [`src/plot_electrode_comparison.py`](src/plot_electrode_comparison.py)

The active documentation is:

- [`docs/ACTIVE_TREE_WALKTHROUGH.md`](docs/ACTIVE_TREE_WALKTHROUGH.md)
- [`docs/multidomain_preisach_hysteresis_report.tex`](docs/multidomain_preisach_hysteresis_report.tex)
- [`docs/multidomain_preisach_hysteresis_report.pdf`](docs/multidomain_preisach_hysteresis_report.pdf)
- [`docs/mifm_first_principles_report.tex`](docs/mifm_first_principles_report.tex)
- [`docs/mifm_first_principles_report.pdf`](docs/mifm_first_principles_report.pdf)

## Model Summary

The active model couples three pieces of physics:

1. A multidomain Preisach ferroelectric constitutive law.
2. One-dimensional electrostatics for the stack.
3. A self-consistent bisection solver that enforces agreement between domain
   switching and internal field.

At one voltage point, the model does the following:

1. Guess the internal ferroelectric field.
2. Update which domains would switch under that field.
3. Average those domains to obtain the net ferroelectric polarization.
4. Feed that polarization into the stack electrostatics.
5. Recompute the field implied by that polarization.
6. Adjust the field guess until the guessed field and the electrostatic field match.

The plotted quantity in the main hysteresis figures is the net ferroelectric
polarization of the AlScN layer.

## Key Equations

The active Preisach model samples a structural-disorder variable `c/a` for each
domain and maps it to local polarization and coercive field using

$$
P_s(c/a) = 333.33(c/a) - 400 \quad \mu\text{C}/\text{cm}^2,
$$

$$
E_c(c/a) = 3.16(c/a) - 1.1 \quad \text{MV}/\text{cm}.
$$

The net ferroelectric polarization is

$$
P_{\mathrm{FE}} = \frac{1}{N}\sum_{n=1}^{N} s_n P_{s,n}.
$$

The active electrostatic core computes screening charge as

$$
\sigma_s =
\frac{
\dfrac{P_d t_d}{k_d}
+
\dfrac{P_{\mathrm{FE}} t_f}{k_f}
+
\epsilon_0 V_a
}{
\dfrac{\lambda_t}{k_t}
+
\dfrac{\lambda_b}{k_b}
+
\dfrac{t_d}{k_d}
+
\dfrac{t_f}{k_f}
+
\dfrac{t_i}{k_i}
}.
$$

The internal ferroelectric field is then

$$
E_{\mathrm{FE}} =
\frac{\sigma_s - P_{\mathrm{FE}}}{k_f \epsilon_0}
+
\frac{\Delta V_{\mathrm{bi,FE}}}{t_f}.
$$

The active bisection solver searches for a field `E` satisfying the fixed-point
condition

$$
r(E) = E - E_{\mathrm{FE}}(P(E),V_a) = 0.
$$

## Active Material Set

The active figure workflow uses the baseline materials defined in
[`src/materials.py`](src/materials.py).

Important baseline values include:

- Ti work function: `4.33 eV`
- Al work function: `4.28 eV`
- Pd work function: `5.30 eV`
- baseline metal screening length: `0.1 bohr`
- `Test` screening length: `1.0 bohr`
- `Test2` screening length: `1.0 bohr`
- `Al2O3` dielectric constant: `9.3`
- `HfOx` dielectric constant: `16.64`
- `AlScN` dielectric constant: `16`

`Test2` is especially useful because it keeps Ti's work function while changing
screening length, which isolates screening effects from built-in work-function
effects.

## Figure Folders

The main generated figure folders are:

- [`45nm_figures`](45nm_figures)
- [`20nm_figures`](20nm_figures)
- [`10nm_figures`](10nm_figures)

These correspond to AlScN thickness.

Representative active figures include:

- [`45nm_figures/Ti_Al2O3_AlScN_Ti_il_hysteresis.png`](45nm_figures/Ti_Al2O3_AlScN_Ti_il_hysteresis.png)
- [`20nm_figures/Ti_Al2O3_AlScN_Ti_il_hysteresis_12V.png`](20nm_figures/Ti_Al2O3_AlScN_Ti_il_hysteresis_12V.png)
- [`10nm_figures/Ti_Al2O3_AlScN_Ti_il_hysteresis_6p5V.png`](10nm_figures/Ti_Al2O3_AlScN_Ti_il_hysteresis_6p5V.png)
- [`45nm_figures/Ti_HfOx_AlScN_Ti_il_hysteresis_20V.png`](45nm_figures/Ti_HfOx_AlScN_Ti_il_hysteresis_20V.png)
- [`20nm_figures/Ti_HfOx_AlScN_Ti_il_hysteresis_12V.png`](20nm_figures/Ti_HfOx_AlScN_Ti_il_hysteresis_12V.png)
- [`10nm_figures/Ti_HfOx_AlScN_Ti_il_hysteresis_6p5V.png`](10nm_figures/Ti_HfOx_AlScN_Ti_il_hysteresis_6p5V.png)
- [`45nm_figures/all_electrodes_0nm_overlay.png`](45nm_figures/all_electrodes_0nm_overlay.png)
- [`45nm_figures/all_electrodes_0nm_deltaP_summary.png`](45nm_figures/all_electrodes_0nm_deltaP_summary.png)

## Main Active Findings

The active figure set supports the following conclusions:

- Thicker interlayers suppress the hysteresis loop strongly.
- HfOx preserves much larger loops than Al2O3 at the same thickness because it
  is a better dielectric match to AlScN.
- Electrode work-function asymmetry can skew or shift the loop.
- Electrode screening length can matter strongly even when work function is held fixed.
- Thickness trends must always be interpreted together with the applied voltage
  window, because ease of driving and depolarization compete.

## Running the Main Hysteresis Driver

The main hysteresis script is
[`src/plot_il_hysteresis.py`](src/plot_il_hysteresis.py).

Example:

```bash
MPLCONFIGDIR=.cache/matplotlib \
UV_CACHE_DIR=.cache/uv \
PYTHONDONTWRITEBYTECODE=1 \
uv run --with matplotlib --with scipy python src/plot_il_hysteresis.py \
  --insulator al2o3 \
  --fe-thickness-nm 20 \
  --top-electrode ti \
  --bottom-electrode ti \
  --vmax 12 \
  --report-loop-summary
```

The electrode comparison script is
[`src/plot_electrode_comparison.py`](src/plot_electrode_comparison.py).

Example:

```bash
MPLCONFIGDIR=.cache/matplotlib \
UV_CACHE_DIR=.cache/uv \
PYTHONDONTWRITEBYTECODE=1 \
uv run --with matplotlib --with scipy python src/plot_electrode_comparison.py \
  --il-nm 0.0 \
  --output-dir 45nm_figures
```

## Other Important Folders

- [`old`](old): archived earlier code and figures
- [`single_domain`](single_domain): separate sandbox for the simpler single-domain study
- [`.tinytex`](.tinytex): repo-local LaTeX toolchain used to compile the reports

## Best Reading Order

If you are new to the project, the cleanest reading order is:

1. [`README.md`](README.md)
2. [`docs/ACTIVE_TREE_WALKTHROUGH.md`](docs/ACTIVE_TREE_WALKTHROUGH.md)
3. [`docs/multidomain_preisach_hysteresis_report.pdf`](docs/multidomain_preisach_hysteresis_report.pdf)
4. [`docs/mifm_first_principles_report.pdf`](docs/mifm_first_principles_report.pdf)

That path goes from orientation, to code understanding, to results, to formal physics.
