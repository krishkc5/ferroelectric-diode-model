import numpy as np

from atomicunits import AtomicUnits


class Potential:
    def __init__(self, fed, v_diff):
        self.fed = fed
        self.v_diff = v_diff
        self.sigma_s = 0.0

        self.il_dv_electrostatic = 0.0
        self.fe_dv_electrostatic = 0.0
        self.dl_dv_electrostatic = 0.0
        self.v_top_interface = 0.0

        self.total_e_field_il = 0.0
        self.total_e_field_fe = 0.0
        self.total_e_field_dl = 0.0

        self.il_v_barrier = fed.top_work_fxn - fed.insulator_chi
        self.fe_v_barrier = fed.bottom_work_fxn - fed.fe_chi
        self.dl_v_barrier = fed.bottom_work_fxn - fed.fe_chi

        d_wf = fed.bottom_work_fxn - fed.top_work_fxn
        self.il_dv_bi = self._built_in_drop(d_wf, "il")
        self.fe_dv_bi = self._built_in_drop(d_wf, "fe")
        self.dl_dv_bi = self._built_in_drop(d_wf, "dl")

        self.set_vdiff(v_diff)

    def _built_in_drop(self, work_function_delta, layer):
        fed = self.fed
        if layer == "il":
            target_thickness = fed.insulator_thickness
            target_k = fed.insulator_k
            denominator = (
                fed.insulator_thickness
                + fed.insulator_k / fed.fe_k * fed.fe_thickness
                + fed.insulator_k / fed.dl_k * fed.dl_thickness
            )
        elif layer == "fe":
            target_thickness = fed.fe_thickness
            target_k = fed.fe_k
            denominator = (
                fed.fe_thickness
                + fed.fe_k / fed.insulator_k * fed.insulator_thickness
                + fed.fe_k / fed.dl_k * fed.dl_thickness
            )
        elif layer == "dl":
            target_thickness = fed.dl_thickness
            target_k = fed.dl_k
            denominator = (
                fed.dl_thickness
                + fed.dl_k / fed.fe_k * fed.fe_thickness
                + fed.dl_k / fed.insulator_k * fed.insulator_thickness
            )
        else:
            raise ValueError(f"Unsupported built-in layer '{layer}'.")
        if target_thickness == 0 or target_k == 0 or denominator == 0:
            return 0.0
        return work_function_delta * target_thickness / denominator

    def _resolve_vdiff(self, v_diff):
        return self.v_diff if v_diff is None else v_diff

    def _electrostatic_state(self, fe_polarization, v_diff=None):
        fed = self.fed
        effective_v_diff = self._resolve_vdiff(v_diff)
        denominator = (
            fed.top_screening_len / fed.top_k
            + fed.bottom_screening_len / fed.bottom_k
            + fed.dl_thickness / fed.dl_k
            + fed.fe_thickness / fed.fe_k
            + fed.insulator_thickness / fed.insulator_k
        )
        sigma_s = (
            fed.dl_polarization * fed.dl_thickness / fed.dl_k
            + fe_polarization * fed.fe_thickness / fed.fe_k
            + AtomicUnits.epsilon_0 * effective_v_diff
        ) / denominator

        # Split the FE voltage drop into:
        #   (i) applied-field-only tilt: sigma_s / (eps_FE eps0) * t_FE
        #       This is Eq. 14 in the P -> 0 limit, with sigma_s reducing
        #       to the polarization-free divider charge.
        #   (ii) polarization-only tilt (Eq. 8 of polarization_barrier_coupling.md):
        #       Delta_U_pol(x) = -(sigma_pol_bot / eps_FE) q x = -P/(eps_FE eps0) q x
        # In the *electron PE* convention used by total_potential, the FE-region
        # electrostatic contribution is the applied-field part; ΔU_pol(x) is
        # carried by `polarization_tilt_potential` (added in `total_potential`).
        # The sum reproduces Eq. 12 exactly: combined slope = (sigma_s - P)/eps_FE.
        il_dv_electrostatic = (
            sigma_s / (fed.insulator_k * AtomicUnits.epsilon_0) * fed.insulator_thickness
            if fed.insulator_thickness != 0
            else 0.0
        )
        fe_dv_electrostatic_applied = (
            sigma_s / (fed.fe_k * AtomicUnits.epsilon_0) * fed.fe_thickness
            if fed.fe_thickness != 0
            else 0.0
        )
        fe_dv_polarization = (
            -fe_polarization / (fed.fe_k * AtomicUnits.epsilon_0) * fed.fe_thickness
            if fed.fe_thickness != 0
            else 0.0
        )
        # Backwards-compatible composite (used by callers that want the total drop).
        fe_dv_electrostatic = fe_dv_electrostatic_applied + fe_dv_polarization
        dl_dv_electrostatic = (
            (sigma_s - fed.dl_polarization) / (fed.dl_k * AtomicUnits.epsilon_0) * fed.dl_thickness
            if fed.dl_thickness != 0
            else 0.0
        )
        v_top_interface = (
            sigma_s * fed.top_screening_len / (AtomicUnits.epsilon_0 * fed.top_k)
            if fed.top_screening_len != 0
            else 0.0
        )
        total_e_field_il = (
            (il_dv_electrostatic + self.il_dv_bi) / fed.insulator_thickness
            if fed.insulator_thickness != 0
            else 0.0
        )
        total_e_field_fe = (
            (fe_dv_electrostatic + self.fe_dv_bi) / fed.fe_thickness if fed.fe_thickness != 0 else 0.0
        )
        total_e_field_dl = (
            (dl_dv_electrostatic + self.dl_dv_bi) / fed.dl_thickness if fed.dl_thickness != 0 else 0.0
        )
        return {
            "v_diff": effective_v_diff,
            "sigma_s": sigma_s,
            "il_dv_electrostatic": il_dv_electrostatic,
            "fe_dv_electrostatic": fe_dv_electrostatic,
            "fe_dv_electrostatic_applied": fe_dv_electrostatic_applied,
            "fe_dv_polarization": fe_dv_polarization,
            "dl_dv_electrostatic": dl_dv_electrostatic,
            "v_top_interface": v_top_interface,
            "total_e_field_il": total_e_field_il,
            "total_e_field_fe": total_e_field_fe,
            "total_e_field_dl": total_e_field_dl,
            "il_dv_bi": self.il_dv_bi,
            "fe_dv_bi": self.fe_dv_bi,
            "dl_dv_bi": self.dl_dv_bi,
        }

    def screening_charge_from_polarization(self, fe_polarization, v_diff=None):
        return self._electrostatic_state(fe_polarization, v_diff=v_diff)["sigma_s"]

    @staticmethod
    def sigma_pol_top(fe_polarization):
        """Bound surface charge density at the top FE face (FE/IL interface).

        Per polarization_barrier_coupling.md Eq. 1: outward normal at the
        top FE face is -z, so sigma_pol,top = P . (-z) = -P.
        Sign convention: P > 0 means polarization along +z (Ti -> Al).
        Units match `fe_polarization` (atomic units in this codebase).
        """
        return -fe_polarization

    @staticmethod
    def sigma_pol_bot(fe_polarization):
        """Bound surface charge density at the bottom FE face (FE/Al interface).

        Per polarization_barrier_coupling.md Eq. 2: outward normal at the
        bottom FE face is +z, so sigma_pol,bot = P . (+z) = +P.
        """
        return +fe_polarization

    def layer_fields_from_polarization(self, fe_polarization, v_diff=None):
        return self._electrostatic_state(fe_polarization, v_diff=v_diff)

    def fe_field_from_polarization(self, fe_polarization, v_diff=None):
        return self._electrostatic_state(fe_polarization, v_diff=v_diff)["total_e_field_fe"]

    def set_vdiff(self, v_diff):
        self.v_diff = v_diff
        state = self._electrostatic_state(self.fed.get_polarization(), v_diff=v_diff)
        self.sigma_s = state["sigma_s"]
        self.il_dv_electrostatic = state["il_dv_electrostatic"]
        self.fe_dv_electrostatic = state["fe_dv_electrostatic"]
        self.dl_dv_electrostatic = state["dl_dv_electrostatic"]
        self.v_top_interface = state["v_top_interface"]
        self.total_e_field_il = state["total_e_field_il"]
        self.total_e_field_fe = state["total_e_field_fe"]
        self.total_e_field_dl = state["total_e_field_dl"]

    def _boundaries(self):
        fed = self.fed
        top_screen_end = 5 * fed.top_screening_len
        il_end = top_screen_end + fed.insulator_thickness
        fe_end = il_end + fed.fe_thickness
        dl_end = fe_end + fed.dl_thickness
        bottom_screen_end = dl_end + 5 * fed.bottom_screening_len
        return top_screen_end, il_end, fe_end, dl_end, bottom_screen_end

    def electrostatic_potential(self, x, fe_polarization=None, v_diff=None):
        # Applied-field-driven contribution to the local potential. Inside the
        # FE region this uses only the `sigma_s / (eps_FE eps0)` slope (the
        # P -> 0 limit of Eq. 12, i.e. Eq. 14). The polarization-only tilt
        # ΔU_pol(x) = -(P/(eps_FE eps0)) q x (Eq. 8) is carried separately by
        # `polarization_tilt_potential`. Sum reproduces Eq. 12 exactly.
        fed = self.fed
        polarization = fed.get_polarization() if fe_polarization is None else fe_polarization
        state = self._electrostatic_state(polarization, v_diff=v_diff)
        top_screen_end, il_end, fe_end, dl_end, bottom_screen_end = self._boundaries()
        # Bottom-screen region must absorb the full V_app drop, but with the
        # FE polarization contribution split out we accumulate the applied-
        # field-only voltage drops up to the bottom electrode.
        v_through_dl_applied = (
            state["v_top_interface"]
            + state["il_dv_electrostatic"]
            + state["fe_dv_electrostatic_applied"]
            + state["dl_dv_electrostatic"]
        )

        if x <= 0:
            return 0.0
        if x <= top_screen_end:
            if fed.top_screening_len == 0:
                return 0.0
            return (
                state["sigma_s"]
                * fed.top_screening_len
                * np.exp(-abs(top_screen_end - x) / fed.top_screening_len)
                / (AtomicUnits.epsilon_0 * fed.top_k)
            )
        if x <= il_end:
            if fed.insulator_thickness == 0:
                return state["v_top_interface"]
            pos = x - top_screen_end
            return state["v_top_interface"] + pos / fed.insulator_thickness * state["il_dv_electrostatic"]
        if x <= fe_end:
            if fed.fe_thickness == 0:
                return state["v_top_interface"] + state["il_dv_electrostatic"]
            pos = x - il_end
            return (
                state["v_top_interface"]
                + state["il_dv_electrostatic"]
                + pos / fed.fe_thickness * state["fe_dv_electrostatic_applied"]
            )
        if x <= dl_end:
            if fed.dl_thickness == 0:
                return (
                    state["v_top_interface"]
                    + state["il_dv_electrostatic"]
                    + state["fe_dv_electrostatic_applied"]
                )
            pos = x - fe_end
            return (
                state["v_top_interface"]
                + state["il_dv_electrostatic"]
                + state["fe_dv_electrostatic_applied"]
                + pos / fed.dl_thickness * state["dl_dv_electrostatic"]
            )
        if x <= bottom_screen_end:
            # The total potential at the bottom electrode equals v_diff (=V_app),
            # which means the applied-only sum + bottom-screen drop must reach
            # v_diff *minus* the polarization tilt that we routed into
            # polarization_tilt_potential. Equivalently, we exit the FE/DL
            # stack at v_through_dl_applied and let the screening tail close.
            if fed.bottom_screening_len == 0:
                return v_through_dl_applied
            return (
                -state["sigma_s"]
                * fed.bottom_screening_len
                * np.exp(-abs(x - dl_end) / fed.bottom_screening_len)
                / (AtomicUnits.epsilon_0 * fed.bottom_k)
                + v_through_dl_applied
            )
        return v_through_dl_applied

    def polarization_tilt_potential(self, x, fe_polarization=None, v_diff=None):
        """Polarization-only tilt across the FE per Eq. 8 of the contract.

        Returns the (atomic-units) contribution to electron PE inside FE:
          ΔU_pol(x_local) = -(sigma_pol_bot / eps_FE) q x_local
                          = -(P / (eps_FE eps0)) q x_local
        with x_local measured from the top FE face. Outside FE this is
        constant (carries the cumulative tilt forward so total potential
        across the stack remains continuous and Eq. 12 is reproduced when
        summed with electrostatic_potential).
        """
        fed = self.fed
        polarization = fed.get_polarization() if fe_polarization is None else fe_polarization
        state = self._electrostatic_state(polarization, v_diff=v_diff)
        top_screen_end, il_end, fe_end, dl_end, bottom_screen_end = self._boundaries()
        fe_dv_polarization = state["fe_dv_polarization"]

        # No tilt before the FE region.
        if x <= il_end:
            return 0.0
        if x <= fe_end:
            if fed.fe_thickness == 0:
                return 0.0
            pos = x - il_end
            return pos / fed.fe_thickness * fe_dv_polarization
        # Past the FE: carry the full cumulative tilt (constant).
        return fe_dv_polarization

    def barrier_potential(self, x):
        fed = self.fed
        top_screen_end, il_end, fe_end, dl_end, _ = self._boundaries()
        no_il_stack = fed.insulator_thickness == 0

        if x <= top_screen_end:
            return 0.0
        if x <= il_end and not no_il_stack:
            return fed.top_fermi_e + self.il_v_barrier
        if x <= fe_end:
            return fed.bottom_fermi_e + self.fe_v_barrier
        if x <= dl_end:
            return fed.bottom_fermi_e + self.dl_v_barrier
        return fed.top_fermi_e - fed.bottom_fermi_e

    def wf_potential(self, x):
        fed = self.fed
        top_screen_end, il_end, fe_end, dl_end, _ = self._boundaries()

        if x <= top_screen_end:
            return 0.0
        if x <= il_end:
            if fed.insulator_thickness == 0:
                return 0.0
            pos = x - top_screen_end
            return pos / fed.insulator_thickness * self.il_dv_bi
        if x <= fe_end:
            if fed.fe_thickness == 0:
                return self.il_dv_bi
            pos = x - il_end
            return self.il_dv_bi + pos / fed.fe_thickness * self.fe_dv_bi
        if x <= dl_end:
            if fed.dl_thickness == 0:
                return self.il_dv_bi + self.fe_dv_bi
            pos = x - fe_end
            return self.il_dv_bi + self.fe_dv_bi + pos / fed.dl_thickness * self.dl_dv_bi
        return 0.0

    def total_potential(self, x, fe_polarization=None, v_diff=None):
        # Eq. 9 of polarization_barrier_coupling.md:
        #   U_FE(x) = U_FE^(0) - q E_FE_applied x + ΔU_pol(x)
        # The first term is `electrostatic_potential` (P -> 0 divider, Eq. 14).
        # ΔU_pol(x) is `polarization_tilt_potential` (Eq. 8). Their sum
        # reproduces Eq. 12 exactly while keeping the polarization tilt
        # explicit and composable.
        return (
            self.electrostatic_potential(x, fe_polarization=fe_polarization, v_diff=v_diff)
            + self.polarization_tilt_potential(x, fe_polarization=fe_polarization, v_diff=v_diff)
            + self.barrier_potential(x)
            + self.wf_potential(x)
        )

    def barrier_region_profile(self, fe_polarization, v_diff=None, num_points=1200):
        fed = self.fed
        start = 5 * fed.top_screening_len
        end = start + fed.barrier_thickness
        if end <= start:
            raise ValueError("Barrier profile requires a non-zero stack thickness.")

        x_au = np.linspace(start, end, num_points)
        sample_x = np.clip(
            x_au + max(1e-9, 1e-9 * max(1.0, end - start)),
            start + 1e-9,
            end - 1e-9,
        )
        total_potential_au = np.array(
            [self.total_potential(x, fe_polarization=fe_polarization, v_diff=v_diff) for x in x_au],
            dtype=float,
        )
        effective_mass = np.array([fed.m_eff(x) for x in sample_x], dtype=float)
        return {
            "x_au": x_au,
            "x_nm": AtomicUnits.bohr_to_nm(x_au),
            "total_potential_au": total_potential_au,
            "total_potential_ev": AtomicUnits.hartree_to_ev(total_potential_au),
            "effective_mass": effective_mass,
        }

    def full_profile(self, fe_polarization=None, v_diff=None, num_points=1600):
        fed = self.fed
        start = 0.0
        end = 5 * fed.top_screening_len + fed.barrier_thickness + 5 * fed.bottom_screening_len
        x_au = np.linspace(start, end, num_points)
        total_potential_au = np.array(
            [self.total_potential(x, fe_polarization=fe_polarization, v_diff=v_diff) for x in x_au],
            dtype=float,
        )
        return {
            "x_au": x_au,
            "x_nm": AtomicUnits.bohr_to_nm(x_au),
            "total_potential_au": total_potential_au,
            "total_potential_ev": AtomicUnits.hartree_to_ev(total_potential_au),
        }
