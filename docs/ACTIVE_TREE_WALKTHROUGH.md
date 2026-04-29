# Active Tree Walkthrough

This document explains the active source tree in plain English.

It is intentionally focused on the files that define, run, explain, and summarize the current multidomain Preisach hysteresis model. It does **not** document:

- archived material in `old/`
- the separate `single_domain/` sandbox
- generated LaTeX auxiliary files such as `.aux`, `.log`, `.out`, or `.toc`
- local tooling artifacts such as `.tinytex/`
- cached plotting directories such as `.mplcache/`

The goal is to make the active repository understandable to someone who is new both to the code and to ferroelectric hysteresis modeling.

The active Python code now lives under `src/`. The original multidomain
Preisach hysteresis workflow remains in its legacy files, while the newer DC
transport work lives in separate transport-specific modules so the older loop
generation path stays intact.

## The Repository at a Glance

The active code implements a **quasistatic multidomain ferroelectric hysteresis model**. It has four main layers:

1. **Unit conversion**
   The code works internally in atomic units, so it needs a reliable place to convert user-facing values such as volts, nanometers, polarization, and MV/cm.
2. **Material and stack definition**
   The code needs to know which electrodes, interlayer, and ferroelectric are present, along with their dielectric constants, work functions, screening lengths, and thicknesses.
3. **Physics solver**
   The code builds a multidomain Preisach ferroelectric, couples it to one-dimensional electrostatics, and solves the field-polarization self-consistency problem with a bisection algorithm.
4. **Figure generation and documentation**
   The code generates hysteresis figures, electrode-comparison figures, and the written reports that summarize what the model is doing and what it predicts.

## File-by-File Guide

## `README.md`

This is the project entry point for a human reader.

Its job is to explain:

- what the repository is for
- what kind of model is being used
- which files matter most
- what the physical picture is

Right now it still carries historical provenance to the original mentor codebase, but conceptually it should be read as the high-level orientation page for the active multidomain Preisach workflow.

## `src/atomicunits.py`

This file is the unit-conversion backbone of the model.

Why it exists:

- The electrostatics and polarization calculations are performed in atomic units.
- Experimental and user-facing quantities are almost never supplied in atomic units.
- The code therefore needs a single trustworthy place for all conversions.

What it contains:

- energy conversions such as eV ↔ Hartree
- length conversions such as nm ↔ bohr and m ↔ bohr
- voltage conversion between volts and atomic units
- polarization conversion between `uC/cm^2` and atomic units
- field conversion between `MV/cm` and atomic units
- density and current-density conversions

Why it matters physically:

- It keeps the rest of the code clean and compact.
- It reduces the chance of unit inconsistency.
- It makes the equations in the solver files much closer to their natural electrostatic form.

How to think about it:

If the rest of the model is the “physics engine,” `src/atomicunits.py` is the translation layer that lets lab-scale quantities and atomic-unit equations talk to one another safely.

## `src/material_types.py`

This file defines the basic material containers used throughout the repository.

It does not perform physics itself. Instead, it provides simple structured objects for:

- `MetalElectrode`
- `Insulator`
- `Ferroelectric`

Why it exists:

- The rest of the code should not have to pass around loose dictionaries or long lists of parameters.
- A `MetalElectrode` should always carry the same kinds of information.
- An `Insulator` should always carry the same kinds of information.

What each class represents:

- `MetalElectrode`
  Contains metal properties such as electron density, Fermi energy, effective mass, dielectric constant used in screening, work function, name, and screening length.
- `Insulator`
  Contains dielectric constant, electron affinity, effective mass, label, and breakdown field.
- `Ferroelectric`
  Contains remanent polarization, electron affinity, effective mass, dielectric constant, name, and trap depth.

Why it matters:

- It gives the rest of the repository a consistent language for “what a material is.”
- It makes it easy to swap one electrode or one interlayer material for another without rewriting the solver itself.

## `src/materials.py`

This file is the material library for the active tree.

It is where the specific Ti, Al, Pd, Test, Test2, Al2O3, HfO2, and AlScN objects are instantiated.

There are two families of material definitions inside this file:

- `baseline_*`
- non-baseline versions without the `baseline_` prefix

In the current active plotting workflow, the **baseline** materials are the ones that matter.

### What the baseline materials are used for

- `baseline_titanium_electrode`
- `baseline_aluminum_electrode`
- `baseline_palladium_electrode`
- `baseline_test_electrode`
- `baseline_test_electrode_2`
- `baseline_al2o3`
- `baseline_hfo2`
- `baseline_alscn`

These are the objects used by `src/plot_il_hysteresis.py` and by the electrode comparison script.

### Why this file matters

This file is where the most important physical constants for the active model are chosen:

- electrode work functions
- metal screening lengths
- interlayer dielectric constants
- ferroelectric dielectric constant
- the baseline remanent polarization stored for AlScN

### What the special test electrodes mean

- `Test`
  Changes both work function and screening length relative to Ti.
- `Test2`
  Keeps the Ti work function but increases the screening length.

`Test2` is especially useful because it isolates the effect of screening without changing the built-in work-function asymmetry.

### Practical interpretation

If you want to know why one figure differs from another because of a changed material stack, this is the first file to check.

## `src/fed.py`

This file defines the stack object: `FerroelectricDiode`.

What it does:

- stores the interlayer, ferroelectric, and dead-layer thicknesses
- stores the dielectric constants used in electrostatics
- stores the top and bottom work functions
- stores the top and bottom screening lengths
- keeps a reference to the active ferroelectric model

Why it matters:

- It is the object that says, “This is the specific capacitor-like stack we are solving right now.”
- It bridges the material library and the electrostatics code.

Important modeling choices embedded here:

- the dead-layer dielectric constant is taken as `k_FE / 2`
- dead-layer polarization is fixed to zero in the active path
- if a metal did not have an explicit screening length, the code could estimate one from electron density and Fermi energy

In practice, for the active figures, the screening length is already provided by the chosen electrode object, so the explicit values in `src/materials.py` dominate.

## `src/potential.py`

This is the electrostatic core of the active model.

It is one of the most important files in the repository.

### What problem it solves

Given:

- a stack geometry
- a current ferroelectric polarization
- an applied voltage

it computes:

- the electrode screening charge
- the ferroelectric internal field
- the built-in field contribution from electrode work-function difference

### Why it matters physically

This file is where the code implements the idea that:

> the ferroelectric does not feel the full externally applied voltage directly

Instead, the applied voltage is divided across:

- the top screening region
- the interlayer
- the ferroelectric
- the dead layer
- the bottom screening region

### Two key outputs

1. `screening_charge_from_polarization(...)`
   Computes the free screening charge that makes the series stack electrostatically self-consistent for a given polarization and voltage.

2. `fe_field_from_polarization(...)`
   Converts that screening charge into the field actually seen by the ferroelectric.

### Why the file is central to the IL-thickness study

Interlayer thickness and interlayer dielectric constant both enter directly into the electrostatic denominator. That is why changing the IL thickness or switching from Al2O3 to HfOx changes the loop so strongly.

### How to think about it

If `src/preisach.py` tells you how domains respond to field, `src/potential.py` tells you what field the stack actually creates.

## `src/preisach.py`

This file is the multidomain ferroelectric model.

It is where the hysteresis memory lives.

### What it does

For each ferroelectric domain, the file stores:

- a sampled structural ratio `c/a`
- a local saturation polarization value
- a local coercive field
- a current state of either `+1` or `-1`

### How the domain population is built

The code samples `c/a` from a normal distribution with:

- mean set by the active driver
- standard deviation set by the active driver

Then it maps the sampled `c/a` to:

- `P_s`
- `E_c`

using linear relations.

### What the active workflow uses

In the active hysteresis path, the code uses:

- `c/a mean = 1.54`
- `c/a std = 0.04`
- default map for `P_s(c/a)`
- default map for `E_c(c/a)`

This means the multidomain distribution is not arbitrary. It is the result of a structural distribution passed through constitutive maps.

### What the file does not do

It does not evolve the system dynamically in time.

It only defines:

- the domain population
- the allowed binary states
- the average polarization of the current state

The actual field-driven switching step is handled in the self-consistent solver.

## `src/self_consistent_solver.py`

This file contains the bisection-based self-consistent solver.

It is the numerical heart of the active model.

### What it solves

At each voltage point, the model needs to find a field `E` such that:

- if the domains are exposed to `E`, they switch into some new state
- that new state produces a polarization `P`
- that polarization, when passed into the electrostatics, produces the same field `E`

That is a fixed-point problem.

### How the solver approaches it

1. Start from the current domain configuration.
2. Guess an internal ferroelectric field.
3. Compute the candidate switched state at that guessed field.
4. Compute the candidate polarization.
5. Feed that polarization into `src/potential.py` to get the self-consistent electrostatic field.
6. Form the residual:
   `guessed field - self-consistent field`
7. Use bisection to search for the field where that residual becomes zero.

### Why bisection is used

The old relaxation-style approach could oscillate badly when polarization feedback was strong.

The bisection approach is more stable for a quasistatic fixed-point problem because it solves directly for the final self-consistent branch instead of stepping through a pseudo-time relaxation.

### Important physical limitation

This file is a **quasistatic equilibrium solver**, not a true switching-dynamics model.

That means:

- it is good for slow hysteresis-loop simulation
- it is not yet a model of transient switching kinetics

## `src/plot_il_hysteresis.py`

This is the main figure-generation script for the active repository.

If someone asks, “Which script produces the IL-thickness hysteresis curves we have been discussing?”, this is the answer.

### What it does

- parses command-line options
- chooses the top and bottom electrodes
- chooses the interlayer material
- chooses the ferroelectric thickness
- constructs the voltage sweep
- initializes the multidomain ferroelectric state
- runs the self-consistent solver for each IL thickness
- plots the final-cycle hysteresis curves
- optionally prints compact loop summaries

### Why it matters

This file is where modeling assumptions become actual experimental-style figures.

It controls:

- IL list
- voltage window
- number of domains
- random seed
- c/a spread
- plot title
- annotation block

### Why the annotation block is useful

The figure text box written by this script makes every output figure self-documenting. It records:

- the stack
- the thickness choices
- the sweep window
- the material parameters
- the Preisach settings
- the solver settings

That is extremely useful for group meetings and paper drafting because it reduces ambiguity about what each plot actually corresponds to.

## `src/plot_electrode_comparison.py`

This file is a targeted comparison script rather than a general hysteresis driver.

Its job is to answer a narrower question:

> how much does the hysteresis loop change when the electrodes are changed?

### What it does

For a fixed IL thickness, currently used at `0 nm`, it:

- runs every ordered top/bottom combination from the active electrode set
- generates an overlay plot
- generates a summary heatmap of mean `|Delta P|` relative to `Ti/Ti`

### Why it is valuable

It separates two different electrode effects:

- work-function asymmetry
- screening-length sensitivity

This script helped establish several important conclusions:

- Ti/Al changes the loop only slightly
- Ti/Pd and Pd/Ti create a much stronger asymmetry
- Test2 changes the loop strongly even though its work function matches Ti, which shows that screening length alone can matter a lot

### How to think about it

If `src/plot_il_hysteresis.py` is the main experiment script, `src/plot_electrode_comparison.py` is the diagnostic script for boundary-condition sensitivity.

## Transport-Specific Files

The transport milestone is intentionally isolated from the older hysteresis
files. The state dataclasses live in `src/simulation_types.py`, the ordered
sweep helpers live in `src/sweeps.py`, and the transport-side hysteresis
wrapper lives in `src/hysteresis_core.py`. The files
`src/transport_fed.py`, `src/transport_potential.py`, and
`src/transport_solver.py` are transport-aware counterparts to the older stack,
electrostatic, and self-consistent solver files. They add the barrier,
effective-mass, field-partition, and bias-point snapshot machinery that the
original hysteresis path does not need.

The actual transport equations live in `src/transport.py`. That file evaluates
thermionic emission, tunneling, Poole-Frenkel conduction, trap-assisted
tunneling, and SCLC at every bias point and sums them into a total current
without using bias-window weights. The main entry point for this transport path
is `src/plot_dciv_sweep.py`, which builds an ordered sweep such as
`0 -> -X -> 0 -> +X -> 0`, advances the multidomain Preisach state across that
trace, evaluates each current mechanism, and saves a figure containing both the
ferroelectric state evolution and the current breakdown.

## `docs/multidomain_preisach_hysteresis_report.tex`

This LaTeX file is the narrative report built directly around the active figure set.

Its job is to explain:

- what the current multidomain Preisach model is
- what assumptions it uses
- what the major generated figures show
- what trends are physically meaningful

It is intentionally more descriptive and results-oriented than the first-principles note.

Think of it as:

- the “what we ran and what we found” document

rather than:

- the “derive every equation from scratch” document

## `docs/mifm_first_principles_report.tex`

This LaTeX file is the formal physics note.

Its role is to explain the exact active model at the level of:

- variables
- equations
- constitutive assumptions
- self-consistency structure
- limitations

This is the document that should answer the question:

> what is the mathematical model we are actually solving?

It is the closest thing in the active tree to a theory document.

## `45nm_figures/`

This folder contains the generated figures for the `45 nm` AlScN case.

These include:

- IL-thickness sweeps
- electrode comparison overlays
- electrode-comparison summary heatmaps
- alternative bottom-electrode studies

This folder is important because it captures the thickest ferroelectric case in the active comparison set.

## `20nm_figures/`

This folder contains the generated figures for the `20 nm` AlScN case.

It is the intermediate thickness case and uses a narrower voltage window than the `45 nm` case so that the loops remain physically informative rather than trivially saturated.

## `10nm_figures/`

This folder contains the generated figures for the `10 nm` AlScN case.

This is the thinnest active ferroelectric case and therefore the one that shows the strongest competition between:

- stronger driving at fixed voltage
- stronger depolarization sensitivity

## How the Pieces Fit Together

The cleanest way to understand the active tree is as a pipeline.

### Step 1: define units and materials

- `src/atomicunits.py`
- `src/material_types.py`
- `src/materials.py`

These files define what the physical quantities are and how they are represented.

### Step 2: define the stack

- `src/fed.py`

This file says what geometry and boundary-condition problem is being solved.

### Step 3: define the ferroelectric state law

- `src/preisach.py`

This file defines how the multidomain ferroelectric can store memory.

### Step 4: define the electrostatics

- `src/potential.py`

This file converts polarization into internal field.

### Step 5: solve the coupled problem

- `src/self_consistent_solver.py`

This file enforces self-consistency between switching and electrostatics.

### Step 6: generate figures

- `src/plot_il_hysteresis.py`
- `src/plot_electrode_comparison.py`

These files turn the model into interpretable outputs.

### Step 7: explain the model and results

- `docs/multidomain_preisach_hysteresis_report.tex`
- `docs/mifm_first_principles_report.tex`

These documents translate code and figures into scientific language.

## The Most Important Big Picture

If someone only remembers one thing from this walkthrough, it should be this:

The active repository is not drawing hysteresis loops by hand.

It is solving a coupled problem in which:

- a multidomain ferroelectric stores switching history
- the stack electrostatics determine the internal field
- the internal field changes which domains switch
- and that changed polarization feeds back into the electrostatics

Everything in the active tree is organized around that loop.
