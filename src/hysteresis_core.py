from dataclasses import dataclass

import numpy as np

from atomicunits import AtomicUnits
from materials import Materials
from preisach import Ferroelectric as PreisachFerroelectric
from transport_fed import FerroelectricDiode
from transport_potential import Potential
from transport_solver import SelfConsistentSolver


ZERO_STRAIN_CA_MEAN = 1.54
DEFAULT_CA_STD = 0.04
DEFAULT_FE_THICKNESS_NM = 45.0
DEFAULT_NUM_DOMAINS = 15000
DEFAULT_MAX_ITER = 5000
DEFAULT_SOLVER_THRESHOLD_UC_CM2 = 0.5
DEFAULT_CYCLES = 2


@dataclass(frozen=True)
class HysteresisConfig:
    num_domains: int = DEFAULT_NUM_DOMAINS
    c_a_mean: float = ZERO_STRAIN_CA_MEAN
    c_a_std: float = DEFAULT_CA_STD
    seed: int = 0
    max_iter: int = DEFAULT_MAX_ITER
    threshold_uc_cm2: float = DEFAULT_SOLVER_THRESHOLD_UC_CM2
    top_electrode_name: str = "ti"
    bottom_electrode_name: str = "ti"
    insulator_name: str = "al2o3"
    fe_thickness_nm: float = DEFAULT_FE_THICKNESS_NM
    dead_layer_thickness_nm: float = 0.0


def resolve_metal_electrode(name: str):
    electrode_map = {
        "ti": Materials.baseline_titanium_electrode,
        "al": Materials.baseline_aluminum_electrode,
        "pd": Materials.baseline_palladium_electrode,
        "test": Materials.baseline_test_electrode,
        "test2": Materials.baseline_test_electrode_2,
    }
    return electrode_map[name]


def resolve_insulator(name: str):
    insulator_map = {
        "al2o3": Materials.baseline_al2o3,
        "hfox": Materials.baseline_hfo2,
    }
    return insulator_map[name]


def insulator_display_label(name: str) -> str:
    display_map = {
        "al2o3": "Al2O3",
        "hfox": "HfOx",
    }
    return display_map[name]


def initialize_domain_states(fe_model: PreisachFerroelectric, seed: int) -> None:
    rng = np.random.default_rng(seed)
    fe_model.set_states(rng.choice([-1, 1], size=fe_model.num_domains))


def create_ferroelectric_model(config: HysteresisConfig) -> PreisachFerroelectric:
    fe_model = PreisachFerroelectric(
        num_domains=config.num_domains,
        c_a_mean=config.c_a_mean,
        c_a_std=config.c_a_std,
        p_s_mean=None,
        p_s_std=None,
        e_c_mean=None,
        e_c_std=None,
        seed=config.seed,
    )
    if np.min(fe_model.p_s_values) <= 0:
        min_ps = AtomicUnits.convert_back_polarization(np.min(fe_model.p_s_values))
        raise ValueError(
            "The default Preisach Ps(c/a) map produced non-positive local saturation polarization "
            f"(minimum {min_ps:.3f} uC/cm^2)."
        )
    initialize_domain_states(fe_model, config.seed)
    return fe_model


class HysteresisSimulator:
    def __init__(self, il_nm: float, config: HysteresisConfig):
        self.il_nm = float(il_nm)
        self.config = config
        self.ferroelectric_model = create_ferroelectric_model(config)
        self.diode = FerroelectricDiode(
            insulator_thickness=AtomicUnits.nm_to_bohr(self.il_nm),
            fe_thickness=AtomicUnits.nm_to_bohr(config.fe_thickness_nm),
            dead_layer_thickness=AtomicUnits.nm_to_bohr(config.dead_layer_thickness_nm),
            top_electrode=resolve_metal_electrode(config.top_electrode_name),
            bottom_electrode=resolve_metal_electrode(config.bottom_electrode_name),
            insulator=resolve_insulator(config.insulator_name),
            ferroelectric=Materials.baseline_alscn,
            fe_model=self.ferroelectric_model,
        )
        self.potential = Potential(self.diode, 0.0)
        self.solver = SelfConsistentSolver(
            ferroelectric=self.ferroelectric_model,
            potential=self.potential,
            max_iter=config.max_iter,
            threshold=config.threshold_uc_cm2,
        )

    def simulate_trace(self, voltage_trace_v: np.ndarray):
        voltage_trace_v = np.asarray(voltage_trace_v, dtype=float)
        if voltage_trace_v.size == 0:
            return []
        self.potential.set_vdiff(AtomicUnits.convert_volts(voltage_trace_v[0]))
        return [self.solver.solve_with_state(float(voltage)) for voltage in voltage_trace_v]
