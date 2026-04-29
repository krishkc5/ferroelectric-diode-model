# Transport mechanism input audit

This note lists every direct input that the current transport implementation uses to generate a current value in `src/transport.py`. The goal is not to defend the model yet, but to make the dependency chain explicit so each quantity can be checked against the stack definition and the solved hysteresis state.

## Common state and stack inputs

Every transport mechanism depends on a solved `BiasPointState` from `src/transport_solver.py` and a stack object from `src/transport_fed.py`. The common state inputs are `voltage_v`, `voltage_au`, `polarization_au`, `polarization_uc_cm2`, `screening_charge_au`, `screening_charge_uc_cm2`, `fe_field_au`, `fe_field_mv_cm`, `il_field_au`, `il_field_mv_cm`, `dl_field_au`, `dl_field_mv_cm`, `fe_built_in_field_au`, `fe_built_in_field_mv_cm`, `il_built_in_field_au`, `il_built_in_field_mv_cm`, `dl_built_in_field_au`, and `dl_built_in_field_mv_cm`.

The common stack inputs are the top and bottom work functions, top and bottom Fermi energies, top and bottom carrier densities, top, bottom, insulator, and ferroelectric effective masses, top and bottom screening lengths, the insulator and ferroelectric dielectric constants, the insulator and ferroelectric electron affinities, the ferroelectric trap depth, and the insulator, ferroelectric, and dead-layer thicknesses.

The common model-parameter inputs are `temperature_k` and the derived thermal voltage `kT/q`.

## Thermionic emission

The thermionic current in `TransportEvaluator._thermionic_current` depends on the applied-voltage polarity, on the top and bottom interface fields through the Schottky lowering term, on the top and bottom work functions, on the insulator or ferroelectric electron affinity used at the injecting interface, on the insulator and ferroelectric dielectric constants used in the Schottky lowering, on the top and bottom electron effective masses, and on the temperature. The tunable model inputs are `richardson_constant_a_per_m2_k2`, `thermionic_prefactor_scale`, `thermionic_top_prefactor_scale`, `thermionic_bottom_prefactor_scale`, `thermionic_top_barrier_offset_ev`, `thermionic_bottom_barrier_offset_ev`, `thermionic_top_barrier_floor_ev`, `thermionic_bottom_barrier_floor_ev`, `thermionic_top_schottky_scale`, and `thermionic_bottom_schottky_scale`.

Written as an explicit dependency chain, the thermionic current uses the following scalar inputs at each bias point: `temperature_k`, `voltage_sign`, `top_field_v_per_m`, `bottom_field_v_per_m`, `top_work_function_ev`, `bottom_work_function_ev`, `top_affinity_ev`, `fe_affinity_ev`, `insulator_k`, `fe_k`, `top_m_eff`, `bottom_m_eff`, `richardson_constant_a_per_m2_k2`, `thermionic_prefactor_scale`, `thermionic_top_prefactor_scale`, `thermionic_bottom_prefactor_scale`, `thermionic_top_barrier_offset_ev`, `thermionic_bottom_barrier_offset_ev`, `thermionic_top_barrier_floor_ev`, `thermionic_bottom_barrier_floor_ev`, `thermionic_top_schottky_scale`, and `thermionic_bottom_schottky_scale`.

## Tunneling

The tunneling current in `TransportEvaluator._tunneling_current` depends on the full barrier profile returned by `Potential.barrier_region_profile`. That profile depends on `polarization_au`, `voltage_au`, the electrostatic partitioning across the metal screening regions, interlayer, ferroelectric, and dead layer, the barrier offsets from work function minus electron affinity, and the built-in work-function drop. The WKB transmission then depends on the sampled position array, the sampled barrier array, and the sampled effective-mass array. The carrier-flux factor depends on the top and bottom carrier densities, Fermi energies, and effective masses. The tunable model input is `tunneling_prefactor_scale`, and the numerical resolution input is `barrier_sampling_points`.

The direct scalar and vector inputs are `voltage_sign`, `barrier_sampling_points`, `x_au`, `barrier_ev(x)`, `effective_mass(x)`, `top_n0`, `bottom_n0`, `top_fermi_e`, `bottom_fermi_e`, `top_m_eff`, `bottom_m_eff`, and `tunneling_prefactor_scale`.

## Poole-Frenkel conduction

The Poole-Frenkel current in `TransportEvaluator._poole_frenkel_current` now supports three region modes: `fe`, `il`, and `both`. In all three cases it uses the magnitude of the solved local field, the common trap depth, the region dielectric constant in the Schottky-like lowering term, the temperature, and the mobility and trap-density prefactors. The tunable model inputs are `poole_frenkel_mobility_m2_per_v_s`, `poole_frenkel_trap_density_m3`, `poole_frenkel_prefactor_scale`, `poole_frenkel_fe_prefactor_scale`, `poole_frenkel_il_prefactor_scale`, and `poole_frenkel_field_region`.

The direct inputs are `voltage_sign`, `pf_field_v_per_m`, `trap_depth_ev`, the active region dielectric constant (`fe_k` or `insulator_k`), `temperature_k`, `poole_frenkel_mobility_m2_per_v_s`, `poole_frenkel_trap_density_m3`, `poole_frenkel_prefactor_scale`, and the region-specific PF scale. In `both` mode the total PF current is the sum of the FE-side and IL-side PF components.

## Trap-assisted tunneling

The trap-assisted tunneling current in `TransportEvaluator._trap_assisted_tunneling_current` uses the same barrier profile and effective-mass profile as the direct tunneling channel, but it reduces the barrier by the trap depth before recomputing the WKB transmission. It also uses a simple field-enhancement factor based on the magnitude of the ferroelectric field. The tunable model input is `trap_assisted_prefactor_a_per_m2`.

The direct inputs are `voltage_sign`, `x_au`, `barrier_ev(x)`, `effective_mass(x)`, `trap_depth_ev`, `pf_field_v_per_m`, `trap_assisted_prefactor_a_per_m2`, and the derived field-enhancement factor.

## Background leakage

The background leakage current in `TransportEvaluator._background_leakage_current` uses the total solved internal electrostatic drop across the insulator, ferroelectric, and dead-layer regions. It is meant to represent a small residual leakage path that remains even when the main PF and thermionic channels are weak. The tunable model input is `background_conductance_a_per_m2_v`.

The direct inputs are `il_field_au`, `fe_field_au`, `dl_field_au`, `insulator_thickness`, `fe_thickness`, `dl_thickness`, and `background_conductance_a_per_m2_v`.

## Space-charge-limited current

The SCLC current in `TransportEvaluator._sclc_current` uses the ferroelectric thickness, the ferroelectric dielectric constant, the ferroelectric field through the corresponding voltage drop across the ferroelectric, the sign of the applied voltage, and the model mobility and scale factor. The tunable model inputs are `sclc_mobility_m2_per_v_s` and `sclc_prefactor_scale`.

The direct inputs are `voltage_v`, `fe_field_au`, `fe_thickness`, `fe_k`, `sclc_mobility_m2_per_v_s`, and `sclc_prefactor_scale`.

## What to inspect during debugging

For thermionic emission, the first quantities to inspect are the base barriers, the Schottky lowerings, the effective barriers after the fitted offsets and barrier floors, and the resulting top and bottom injection currents. For tunneling, the barrier-profile minimum and maximum, the effective-mass range, and the WKB transmission are the most important diagnostics. For Poole-Frenkel conduction, the key quantities are the FE-side and IL-side fields, the corresponding lowerings and activation energies, and the relative FE/IL PF scales when `both` mode is active. For trap-assisted tunneling, the reduced barrier after subtracting the trap depth and the resulting transmission are the central checks. For background leakage, the key quantity is the net internal electrostatic drop across the active dielectric stack. For SCLC, the important quantities are the ferroelectric thickness, the inferred voltage drop across the ferroelectric, and the resulting `V^2/L^3` scaling.
