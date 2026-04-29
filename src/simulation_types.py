from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class BiasPointState:
    voltage_v: float
    voltage_au: float
    domain_states: np.ndarray
    polarization_au: float
    polarization_uc_cm2: float
    polarization_change_uc_cm2: float
    screening_charge_au: float
    screening_charge_uc_cm2: float
    screening_charge_c_m2: float
    fe_field_au: float
    fe_field_mv_cm: float
    il_field_au: float
    il_field_mv_cm: float
    dl_field_au: float
    dl_field_mv_cm: float
    fe_built_in_field_au: float
    fe_built_in_field_mv_cm: float
    il_built_in_field_au: float
    il_built_in_field_mv_cm: float
    dl_built_in_field_au: float
    dl_built_in_field_mv_cm: float
    residual_au: float
    residual_mv_cm: float
    converged: bool
    used_fallback: bool


@dataclass(frozen=True)
class TransportBreakdown:
    voltage_v: float
    temperature_k: float
    thermionic_a_per_m2: float
    tunneling_a_per_m2: float
    poole_frenkel_a_per_m2: float
    trap_assisted_tunneling_a_per_m2: float
    sclc_a_per_m2: float
    total_a_per_m2: float
    tunneling_transmission: float
    trap_tunneling_transmission: float
    top_barrier_ev: float
    bottom_barrier_ev: float
    background_leakage_a_per_m2: float = 0.0
