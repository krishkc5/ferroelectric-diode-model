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
        self.insulator_k = insulator.k
        self.fe_k = ferroelectric.k
        self.dl_k = ferroelectric.k / 2
        self.top_k = top_electrode.k
        self.bottom_k = bottom_electrode.k
        self.top_work_fxn = top_electrode.w_f
        self.bottom_work_fxn = bottom_electrode.w_f
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

    def get_polarization(self):
        return self.fe_model.avg_polarization()
