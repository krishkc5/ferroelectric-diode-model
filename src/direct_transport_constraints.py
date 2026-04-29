from dataclasses import dataclass

import numpy as np


GLOBAL_PARAMETER_ORDER = (
    "top_phi_ev",
    "top_gamma",
    "log10_top_prefactor_a",
    "bottom_phi_ev",
    "bottom_gamma",
    "log10_bottom_prefactor_a",
    "pf_trap_depth_ev",
    "pf_fe_beta",
    "log10_pf_fe_prefactor_a_per_v",
    "pf_il_beta",
    "log10_pf_il_prefactor_a_per_v",
    "log10_bg_cond",
    "log10_i_floor",
)


@dataclass(frozen=True)
class ConstraintPreset:
    name: str
    description: str
    initial_guess: np.ndarray
    lower_bounds: np.ndarray
    upper_bounds: np.ndarray


def _as_array(values):
    return np.asarray(values, dtype=float)


def broad_preset():
    return ConstraintPreset(
        name="broad",
        description="Loose bounds intended for exploratory fitting.",
        initial_guess=_as_array([1.0, 0.3, -8.0, 1.0, 0.2, -8.0, 0.8, 0.10, -10.0, 0.20, -10.0, -12.0, -14.0]),
        lower_bounds=_as_array([0.0, 0.0, -20.0, 0.0, 0.0, -20.0, 0.05, 0.0, -20.0, 0.0, -20.0, -20.0, -20.0]),
        upper_bounds=_as_array([4.0, 5.0, 4.0, 4.0, 5.0, 4.0, 2.0, 5.0, 4.0, 5.0, 4.0, 2.0, -6.0]),
    )


def physics_reasonable_preset():
    return ConstraintPreset(
        name="physics-reasonable",
        description=(
            "Moderately constrained ranges for Ti | HfOx | AlScN | Al direct transport fits. "
            "Use these as a starting point before tightening to literature."
        ),
        initial_guess=_as_array([0.8, 0.20, -8.0, 1.0, 0.20, -8.0, 0.7, 0.20, -9.0, 0.20, -10.0, -13.0, -14.0]),
        lower_bounds=_as_array([0.2, 0.0, -14.0, 0.2, 0.0, -14.0, 0.3, 0.0, -14.0, 0.0, -14.0, -18.0, -18.0]),
        upper_bounds=_as_array([2.0, 1.0, -4.0, 2.5, 1.0, -4.0, 1.2, 0.6, -6.0, 0.6, -6.0, -8.0, -10.0]),
    )


def literature_ready_preset():
    return ConstraintPreset(
        name="literature-ready",
        description=(
            "Tighter placeholder bounds for runs after literature values have been reviewed and updated."
        ),
        initial_guess=_as_array([0.8, 0.15, -8.0, 1.0, 0.15, -8.0, 0.7, 0.15, -9.0, 0.15, -10.0, -13.0, -14.0]),
        lower_bounds=_as_array([0.3, 0.02, -12.0, 0.4, 0.02, -12.0, 0.4, 0.02, -12.0, 0.02, -12.0, -16.0, -16.0]),
        upper_bounds=_as_array([1.5, 0.5, -5.0, 1.8, 0.5, -5.0, 1.0, 0.35, -7.0, 0.35, -7.0, -10.0, -11.0]),
    )


def get_constraint_preset(name: str):
    presets = {
        "broad": broad_preset(),
        "physics-reasonable": physics_reasonable_preset(),
        "literature-ready": literature_ready_preset(),
    }
    try:
        return presets[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown constraint preset '{name}'. Available: {', '.join(sorted(presets))}"
        ) from exc

