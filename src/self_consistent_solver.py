import numpy as np

from atomicunits import AtomicUnits


class SelfConsistentSolver:
    def __init__(self, ferroelectric, potential, max_iter=5000, threshold=0.5):
        self.ferroelectric = ferroelectric
        self.potential = potential
        self.max_iter = max_iter
        self._polarization_tolerance = AtomicUnits.convert_polarization(threshold)

    def _candidate_states(self, base_states, e_field):
        states = base_states.copy()
        flip_down = (base_states == 1) & (e_field <= -self.ferroelectric.e_c_values)
        flip_up = (base_states == -1) & (e_field >= self.ferroelectric.e_c_values)
        states[flip_down] = -1
        states[flip_up] = 1
        return states

    def _evaluate_candidate(self, e_field, voltage_au, base_states):
        states = self._candidate_states(base_states, e_field)
        polarization = float(np.mean(states * self.ferroelectric.p_s_values))
        e_self_consistent = self.potential.fe_field_from_polarization(
            polarization, v_diff=voltage_au
        )
        return {
            "residual": e_field - e_self_consistent,
            "states": states,
            "polarization": polarization,
        }

    def _initial_bracket(self, voltage_au):
        p_max = float(np.mean(np.abs(self.ferroelectric.p_s_values)))
        extreme_fields = (
            self.potential.fe_field_from_polarization(-p_max, v_diff=voltage_au),
            self.potential.fe_field_from_polarization(p_max, v_diff=voltage_au),
        )
        max_ec = float(np.max(np.abs(self.ferroelectric.e_c_values)))
        margin = max_ec + max(abs(field) for field in extreme_fields) + AtomicUnits.Mv_per_cm_to_atomic_units(1.0)
        low = min(*extreme_fields, -max_ec) - margin
        high = max(*extreme_fields, max_ec) + margin
        return low, high

    def solve(self, voltage):
        voltage_au = AtomicUnits.convert_volts(voltage)
        base_states = self.ferroelectric.state_values.copy()
        previous_polarization = float(np.mean(base_states * self.ferroelectric.p_s_values))

        low, high = self._initial_bracket(voltage_au)
        low_eval = self._evaluate_candidate(low, voltage_au, base_states)
        high_eval = self._evaluate_candidate(high, voltage_au, base_states)

        for _ in range(24):
            if low_eval["residual"] <= 0 <= high_eval["residual"]:
                break
            span = high - low
            low -= span
            high += span
            low_eval = self._evaluate_candidate(low, voltage_au, base_states)
            high_eval = self._evaluate_candidate(high, voltage_au, base_states)
        else:
            final_eval = min((low_eval, high_eval), key=lambda candidate: abs(candidate["residual"]))
            self.ferroelectric.set_states(final_eval["states"])
            self.potential.set_vdiff(voltage_au)
            polarization_change = final_eval["polarization"] - previous_polarization
            return final_eval["polarization"], AtomicUnits.convert_back_polarization(polarization_change)

        best_eval = min((low_eval, high_eval), key=lambda candidate: abs(candidate["residual"]))
        final_eval = best_eval

        for _ in range(self.max_iter):
            if abs(high_eval["polarization"] - low_eval["polarization"]) <= self._polarization_tolerance:
                final_eval = min(
                    (low_eval, high_eval, best_eval), key=lambda candidate: abs(candidate["residual"])
                )
                break

            mid = 0.5 * (low + high)
            mid_eval = self._evaluate_candidate(mid, voltage_au, base_states)
            if abs(mid_eval["residual"]) < abs(best_eval["residual"]):
                best_eval = mid_eval

            if mid_eval["residual"] < 0:
                low, low_eval = mid, mid_eval
            else:
                high, high_eval = mid, mid_eval
        else:
            final_eval = min(
                (low_eval, high_eval, best_eval), key=lambda candidate: abs(candidate["residual"])
            )

        self.ferroelectric.set_states(final_eval["states"])
        self.potential.set_vdiff(voltage_au)
        polarization_change = final_eval["polarization"] - previous_polarization
        return final_eval["polarization"], AtomicUnits.convert_back_polarization(polarization_change)
