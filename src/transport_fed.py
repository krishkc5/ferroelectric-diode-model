import math

import scipy.constants as constants

from atomicunits import AtomicUnits


class FerroelectricDiode:
    def __init__(
        self,
        insulator_thickness,
        fe_thickness,
        dead_layer_thickness,
        top_electrode,
        bottom_electrode,
        insulator,
        ferroelectric,
        fe_model,
    ):
        self.insulator_thickness = insulator_thickness
        self.fe_thickness = fe_thickness
        self.dl_thickness = dead_layer_thickness
        self.barrier_thickness = insulator_thickness + fe_thickness + dead_layer_thickness
        self.insulator_k = insulator.k
        self.fe_k = ferroelectric.k
        self.dl_k = ferroelectric.k / 2
        self.top_k = top_electrode.k
        self.bottom_k = bottom_electrode.k
        self.insulator_chi = insulator.chi
        self.fe_chi = ferroelectric.chi
        self.top_work_fxn = top_electrode.w_f
        self.bottom_work_fxn = bottom_electrode.w_f
        self.top_fermi_e = top_electrode.e_f
        self.bottom_fermi_e = bottom_electrode.e_f
        self.top_n0 = top_electrode.n0
        self.bottom_n0 = bottom_electrode.n0
        self.insulator_m_eff = insulator.m_eff
        self.fe_m_eff = ferroelectric.m_eff
        self.top_m_eff = top_electrode.m_eff
        self.bottom_m_eff = bottom_electrode.m_eff
        self.fe_trap_depth = ferroelectric.trap_depth
        self.name = (
            f"{top_electrode.name}-{insulator.name}-{ferroelectric.name}-{bottom_electrode.name}"
        )
        self.fe_model = fe_model

        if top_electrode.screening_len is None:
            self.top_screening_len = AtomicUnits.m_to_bohr(
                math.sqrt(
                    top_electrode.k
                    * 2
                    * constants.epsilon_0
                    * AtomicUnits.hartree_to_joule(top_electrode.e_f)
                    / (3 * constants.e**2 * AtomicUnits.convert_back_density(top_electrode.n0))
                )
            )
        else:
            self.top_screening_len = top_electrode.screening_len

        if bottom_electrode.screening_len is None:
            self.bottom_screening_len = AtomicUnits.m_to_bohr(
                math.sqrt(
                    bottom_electrode.k
                    * 2
                    * constants.epsilon_0
                    * AtomicUnits.hartree_to_joule(bottom_electrode.e_f)
                    / (3 * constants.e**2 * AtomicUnits.convert_back_density(bottom_electrode.n0))
                )
            )
        else:
            self.bottom_screening_len = bottom_electrode.screening_len

        self.dl_polarization = 0.0
        self.fe_polarization = self.fe_model.avg_polarization()

    def m_eff(self, x):
        top_screen_region = 5 * self.top_screening_len
        il_end = top_screen_region + self.insulator_thickness
        fe_end = il_end + self.fe_thickness
        dl_end = fe_end + self.dl_thickness

        if x <= top_screen_region:
            return self.top_m_eff
        if x <= il_end:
            return self.insulator_m_eff
        if x <= fe_end:
            return self.fe_m_eff
        if x <= dl_end:
            return self.fe_m_eff
        return self.bottom_m_eff

    def get_polarization(self):
        self.fe_polarization = self.fe_model.avg_polarization()
        return self.fe_polarization
