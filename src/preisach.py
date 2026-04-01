import numpy as np

from atomicunits import AtomicUnits


class Ferroelectric:
    def __init__(
        self,
        num_domains,
        c_a_mean,
        c_a_std,
        p_s_mean=None,
        p_s_std=None,
        e_c_mean=None,
        e_c_std=None,
        seed=0,
    ):
        np.random.seed(seed)
        c_a_ratios = np.random.normal(loc=c_a_mean, scale=c_a_std, size=num_domains)

        if p_s_mean is None and p_s_std is None:
            p_s_values = AtomicUnits.convert_polarization(333.33 * c_a_ratios - 400)
        else:
            slope = p_s_std / c_a_std
            intercept = p_s_mean - slope * c_a_mean
            p_s_values = AtomicUnits.convert_polarization(np.abs(slope * c_a_ratios + intercept))

        if e_c_mean is None and e_c_std is None:
            e_c_values = AtomicUnits.Mv_per_cm_to_atomic_units(3.16 * c_a_ratios - 1.1)
        else:
            slope = e_c_std / c_a_std
            intercept = e_c_mean - slope * c_a_mean
            e_c_values = AtomicUnits.Mv_per_cm_to_atomic_units(np.abs(slope * c_a_ratios + intercept))

        self.num_domains = num_domains
        self.c_a_ratios = c_a_ratios
        self.p_s_values = p_s_values
        self.e_c_values = e_c_values
        self.state_values = np.ones(num_domains, dtype=np.int8)

    def set_states(self, state_values):
        state_array = np.asarray(state_values, dtype=np.int8)
        if state_array.shape != (self.num_domains,):
            raise ValueError(
                f"Expected state array of shape {(self.num_domains,)}, got {state_array.shape}."
            )
        if not np.all(np.isin(state_array, (-1, 1))):
            raise ValueError("Ferroelectric domain states must be +/-1.")
        self.state_values = state_array.copy()

    def avg_polarization(self):
        return float(np.mean(self.state_values * self.p_s_values))
