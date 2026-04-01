from atomicunits import AtomicUnits


class Potential:
    def __init__(self, fed, v_diff):
        self.fed = fed
        self.v_diff = v_diff
        self.total_e_field_fe = 0.0

        d_wf = fed.bottom_work_fxn - fed.top_work_fxn
        self.fe_dv_bi = d_wf * fed.fe_thickness / (
            fed.fe_thickness
            + fed.fe_k / fed.insulator_k * fed.insulator_thickness
            + fed.fe_k / fed.dl_k * fed.dl_thickness
        )

        self.set_vdiff(v_diff)

    def screening_charge_from_polarization(self, fe_polarization, v_diff=None):
        fed = self.fed
        effective_v_diff = self.v_diff if v_diff is None else v_diff
        return (
            fed.dl_polarization * fed.dl_thickness / fed.dl_k
            + fe_polarization * fed.fe_thickness / fed.fe_k
            + AtomicUnits.epsilon_0 * effective_v_diff
        ) / (
            fed.top_screening_len / fed.top_k
            + fed.bottom_screening_len / fed.bottom_k
            + fed.dl_thickness / fed.dl_k
            + fed.fe_thickness / fed.fe_k
            + fed.insulator_thickness / fed.insulator_k
        )

    def fe_field_from_polarization(self, fe_polarization, v_diff=None):
        sigma_s = self.screening_charge_from_polarization(fe_polarization, v_diff=v_diff)
        return (
            (sigma_s - fe_polarization) / (self.fed.fe_k * AtomicUnits.epsilon_0)
            + self.fe_dv_bi / self.fed.fe_thickness
        )

    def set_vdiff(self, v_diff):
        self.v_diff = v_diff
        fe_polarization = self.fed.get_polarization()
        self.total_e_field_fe = self.fe_field_from_polarization(fe_polarization, v_diff=v_diff)
