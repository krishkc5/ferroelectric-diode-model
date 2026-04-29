import numpy as np

from atomicunits import AtomicUnits
from simulation_types import BiasPointState


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
            "guess_field": e_field,
            "self_consistent_field": e_self_consistent,
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

    def _build_state_snapshot(
        self,
        voltage,
        voltage_au,
        previous_polarization,
        final_eval,
        converged,
        used_fallback,
    ):
        polarization_au = final_eval["polarization"]
        screening_charge_au = self.potential.screening_charge_from_polarization(
            polarization_au, v_diff=voltage_au
        )
        layer_fields = self.potential.layer_fields_from_polarization(
            polarization_au, v_diff=voltage_au
        )
        polarization_change_au = polarization_au - previous_polarization

        polarization_uc_cm2 = AtomicUnits.convert_back_polarization(polarization_au)
        polarization_change_uc_cm2 = AtomicUnits.convert_back_polarization(polarization_change_au)
        screening_charge_uc_cm2 = AtomicUnits.convert_back_polarization(screening_charge_au)

        return BiasPointState(
            voltage_v=float(voltage),
            voltage_au=voltage_au,
            domain_states=final_eval["states"].copy(),
            polarization_au=polarization_au,
            polarization_uc_cm2=polarization_uc_cm2,
            polarization_change_uc_cm2=polarization_change_uc_cm2,
            screening_charge_au=screening_charge_au,
            screening_charge_uc_cm2=screening_charge_uc_cm2,
            screening_charge_c_m2=AtomicUnits.uc_per_cm2_to_c_per_m2(screening_charge_uc_cm2),
            fe_field_au=layer_fields["total_e_field_fe"],
            fe_field_mv_cm=AtomicUnits.atomic_units_to_Mv_per_cm(layer_fields["total_e_field_fe"]),
            il_field_au=layer_fields["total_e_field_il"],
            il_field_mv_cm=AtomicUnits.atomic_units_to_Mv_per_cm(layer_fields["total_e_field_il"]),
            dl_field_au=layer_fields["total_e_field_dl"],
            dl_field_mv_cm=AtomicUnits.atomic_units_to_Mv_per_cm(layer_fields["total_e_field_dl"]),
            fe_built_in_field_au=(
                self.potential.fe_dv_bi / self.potential.fed.fe_thickness
                if self.potential.fed.fe_thickness != 0
                else 0.0
            ),
            fe_built_in_field_mv_cm=AtomicUnits.atomic_units_to_Mv_per_cm(
                self.potential.fe_dv_bi / self.potential.fed.fe_thickness
                if self.potential.fed.fe_thickness != 0
                else 0.0
            ),
            il_built_in_field_au=(
                self.potential.il_dv_bi / self.potential.fed.insulator_thickness
                if self.potential.fed.insulator_thickness != 0
                else 0.0
            ),
            il_built_in_field_mv_cm=AtomicUnits.atomic_units_to_Mv_per_cm(
                self.potential.il_dv_bi / self.potential.fed.insulator_thickness
                if self.potential.fed.insulator_thickness != 0
                else 0.0
            ),
            dl_built_in_field_au=(
                self.potential.dl_dv_bi / self.potential.fed.dl_thickness
                if self.potential.fed.dl_thickness != 0
                else 0.0
            ),
            dl_built_in_field_mv_cm=AtomicUnits.atomic_units_to_Mv_per_cm(
                self.potential.dl_dv_bi / self.potential.fed.dl_thickness
                if self.potential.fed.dl_thickness != 0
                else 0.0
            ),
            residual_au=final_eval["residual"],
            residual_mv_cm=AtomicUnits.atomic_units_to_Mv_per_cm(final_eval["residual"]),
            converged=converged,
            used_fallback=used_fallback,
        )

    def solve_with_state(self, voltage):
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
            return self._build_state_snapshot(
                voltage=voltage,
                voltage_au=voltage_au,
                previous_polarization=previous_polarization,
                final_eval=final_eval,
                converged=False,
                used_fallback=True,
            )

        best_eval = min((low_eval, high_eval), key=lambda candidate: abs(candidate["residual"]))
        final_eval = best_eval
        converged = False

        for _ in range(self.max_iter):
            if abs(high_eval["polarization"] - low_eval["polarization"]) <= self._polarization_tolerance:
                final_eval = min(
                    (low_eval, high_eval, best_eval), key=lambda candidate: abs(candidate["residual"])
                )
                converged = True
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
        return self._build_state_snapshot(
            voltage=voltage,
            voltage_au=voltage_au,
            previous_polarization=previous_polarization,
            final_eval=final_eval,
            converged=converged,
            used_fallback=not converged,
        )

    def solve(self, voltage):
        state = self.solve_with_state(voltage)
        return state.polarization_au, state.polarization_change_uc_cm2
