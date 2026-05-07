from dataclasses import dataclass

import numpy as np
import scipy.constants as constants

from atomicunits import AtomicUnits
from simulation_types import BiasPointState, TransportBreakdown


@dataclass(frozen=True)
class TransportModelParameters:
    temperature_k: float = 300.0
    richardson_constant_a_per_m2_k2: float = 1.2e6
    thermionic_prefactor_scale: float = 1.0
    thermionic_top_prefactor_scale: float = 1.0
    thermionic_bottom_prefactor_scale: float = 1.0
    thermionic_barrier_offset_ev: float = 0.0
    thermionic_top_barrier_offset_ev: float = 0.0
    thermionic_bottom_barrier_offset_ev: float = 0.0
    thermionic_top_barrier_floor_ev: float = 0.0
    thermionic_bottom_barrier_floor_ev: float = 0.0
    thermionic_top_schottky_scale: float = 1.0
    thermionic_bottom_schottky_scale: float = 1.0
    # Eq. 4' (polarization_barrier_coupling.md §3.1): Fermi-level pinning floor for
    # the saturating-screening correction. ΔV_screen = ΔV_pin · tanh(ΔV_lin / ΔV_pin).
    # Default 1.0 V matches the canonical metal/nitride pinning range (Casamento APL
    # 120, 152901 (2022); Tsymbal-Kohlstedt 2006 review §III).
    fermi_level_pinning_top_v: float = 1.0
    fermi_level_pinning_bot_v: float = 1.0
    # Tung barrier-inhomogeneity sigma [eV] for thermionic emission.
    # A Gaussian distribution of barrier heights N(phi_bar, sigma^2) gives an
    # effective TE rate J = J_homog * exp(sigma^2 / (2 (kT)^2)). Real metal/oxide
    # interfaces have sigma_phi ~ 0.05-0.3 eV due to interface roughness, grain
    # variation, and dipole inhomogeneity. Default 0.0 = off (recover homogeneous
    # TE). Tung, Phys. Rev. B 45, 13509 (1992); Werner-Guettler, J. Appl. Phys.
    # 69, 1522 (1991).
    te_barrier_sigma_top_ev: float = 0.0
    te_barrier_sigma_bot_ev: float = 0.0
    tunneling_prefactor_scale: float = 1.0
    poole_frenkel_mobility_m2_per_v_s: float = 1.0e-10
    poole_frenkel_trap_density_m3: float = 5.0e23
    poole_frenkel_prefactor_scale: float = 1.0
    poole_frenkel_fe_prefactor_scale: float = 1.0
    poole_frenkel_il_prefactor_scale: float = 1.0
    poole_frenkel_trap_depth_ev: float | None = None
    poole_frenkel_field_region: str = "fe"
    trap_assisted_prefactor_a_per_m2: float = 5.0e8
    sclc_mobility_m2_per_v_s: float = 5.0e-10
    sclc_prefactor_scale: float = 1.0
    background_conductance_a_per_m2_v: float = 0.0
    barrier_sampling_points: int = 1500


class TransportEvaluator:
    def __init__(
        self,
        diode,
        potential,
        parameters: TransportModelParameters | None = None,
        enable_thermionic: bool = True,
        enable_tunneling: bool = True,
        enable_poole_frenkel: bool = True,
        enable_trap_assisted_tunneling: bool = True,
        enable_sclc: bool = True,
    ):
        self.diode = diode
        self.potential = potential
        self.parameters = parameters or TransportModelParameters()
        self.enable_thermionic = enable_thermionic
        self.enable_tunneling = enable_tunneling
        self.enable_poole_frenkel = enable_poole_frenkel
        self.enable_trap_assisted_tunneling = enable_trap_assisted_tunneling
        self.enable_sclc = enable_sclc

    @property
    def _thermal_voltage(self):
        return constants.k * self.parameters.temperature_k / constants.e

    def _field_to_v_per_m(self, field_au):
        return AtomicUnits.atomic_units_to_Mv_per_cm(field_au) * 1e8

    def _schottky_lowering_ev(self, field_v_per_m, epsilon_r):
        """Image-force (Schottky) barrier lowering. Eq. 15 of
        polarization_barrier_coupling.md: Delta_phi_IF = sqrt(q E / (4 pi eps0 kappa_inf)).

        Caller must pass kappa_infinity (optical permittivity), not kappa_static -
        the carrier-induced image charge follows electronic-only dielectric
        response (Mehta-Silverman-Jacobs 1973; Sze & Ng Ch. 3 sec 3.2.4).
        """
        if field_v_per_m <= 0 or epsilon_r <= 0:
            return 0.0
        lowering_volts = np.sqrt(
            constants.e * field_v_per_m / (4 * np.pi * constants.epsilon_0 * epsilon_r)
        )
        return float(lowering_volts)

    def _poole_frenkel_lowering_ev(self, field_v_per_m, epsilon_r):
        """Poole-Frenkel barrier lowering. polarization_barrier_coupling.md sec 7
        (PF row): Delta_phi_PF = sqrt(q^3 E_FE / (pi eps0 kappa_inf,FE)).

        Two differences from the Schottky form:
          (a) leading factor is pi, not 4 pi (Sze & Ng Ch. 3 Eq. 137);
          (b) kappa is the optical permittivity of the bulk dielectric
              (Mehta-Silverman-Jacobs 1973).
        Magnitude only - the carrier always sees a lower barrier on the
        side it is moving toward; sign-of-current is carried separately.
        """
        if field_v_per_m <= 0 or epsilon_r <= 0:
            return 0.0
        lowering_volts = np.sqrt(
            constants.e * field_v_per_m / (np.pi * constants.epsilon_0 * epsilon_r)
        )
        return float(lowering_volts)

    def _polarization_dipole_shift_ev(self, sigma_pol_au, lambda_tf_au, kappa_inf_metal, delta_v_pin_v):
        """Tsymbal-Kohlstedt screening shift of the electron Schottky barrier,
        with Fermi-level-pinning saturation (polarization_barrier_coupling.md Eq. 4').

        Linear T-K shift (Eq. 4):
            Delta_V_lin = - sigma_pol * lambda_TF / (eps0 * kappa_inf,metal)    [V]

        Saturating form (Eq. 4'), needed at AlScN-scale polarization where the linear
        formula leaves its domain of validity (BTO/BFO: ~1-3 V, fine; AlScN: ~5+ V,
        unphysical — see §3.1 of the contract):
            Delta_V_screen = Delta_V_pin * tanh(Delta_V_lin / Delta_V_pin)      [V]

        Limits: Delta_V_screen -> Delta_V_lin for |Delta_V_lin| << Delta_V_pin
        (BTO/BFO regime); Delta_V_screen -> ±Delta_V_pin for |Delta_V_lin| >>
        Delta_V_pin (AlScN regime, Fermi-level pinned).

        Per Eq. 5/6 the leading minus is uniform; asymmetry between top/bot interfaces
        comes from sigma_pol,top = -P, sigma_pol,bot = +P (Eq. 1, 2). Numerically
        Delta_phi_B [eV] = Delta_V_screen [V].

        Inputs in atomic units (sigma in atomic charge density, lambda in Bohr).
        kappa_inf_metal is the metal's *optical* permittivity ~1 (Drude tail).
        delta_v_pin_v is the Fermi-level pinning floor [V]; ~1.0 V for nitride/metal.
        """
        if kappa_inf_metal <= 0 or delta_v_pin_v <= 0:
            return 0.0
        sigma_si = AtomicUnits.convert_back_polarization(sigma_pol_au) / 1e2  # C/m^2
        lambda_si = AtomicUnits.bohr_to_m(lambda_tf_au)
        # Eq. 4 (linear): leading minus per Eq. 5/6 sign convention.
        delta_v_lin_v = -sigma_si * lambda_si / (constants.epsilon_0 * kappa_inf_metal)
        # Eq. 4' (saturating): tanh-cap at the pinning floor. Preserves sign of delta_v_lin.
        delta_v_screen_v = delta_v_pin_v * float(np.tanh(delta_v_lin_v / delta_v_pin_v))
        return delta_v_screen_v

    def _bias_drive(self, voltage_v):
        v_t = max(self._thermal_voltage, 1e-6)
        return float(np.tanh(np.clip(voltage_v / (2 * v_t), -50.0, 50.0)))

    def _voltage_sign(self, voltage_v):
        if abs(voltage_v) < 1e-15:
            return 0.0
        return float(np.sign(voltage_v))

    def _field_sign(self, field_au, fallback_voltage_v=0.0):
        if abs(field_au) > 1e-18:
            return float(np.sign(field_au))
        return self._voltage_sign(fallback_voltage_v)

    def _carrier_flux(self, carrier_density_au, fermi_energy_au, effective_mass_ratio):
        density_m3 = AtomicUnits.convert_back_density(carrier_density_au)
        fermi_energy_j = AtomicUnits.hartree_to_joule(fermi_energy_au)
        mass_kg = max(effective_mass_ratio, 1e-6) * constants.m_e
        fermi_velocity = np.sqrt(max(2 * fermi_energy_j / mass_kg, 0.0))
        return 0.25 * constants.e * density_m3 * fermi_velocity

    def _barriers_ev(self, state: BiasPointState):
        # Image-force lowering uses |E| of the injection-side dielectric
        # (Eq. 15: Delta_phi_IF magnitude is sign-independent of E).
        top_field_v_per_m = abs(
            self._field_to_v_per_m(state.il_field_au if self.diode.insulator_thickness != 0 else state.fe_field_au)
        )
        bottom_field_v_per_m = abs(self._field_to_v_per_m(state.fe_field_au))
        top_affinity = self.diode.insulator_chi if self.diode.insulator_thickness != 0 else self.diode.fe_chi
        top_base_barrier = AtomicUnits.hartree_to_ev(self.diode.top_work_fxn - top_affinity)
        bottom_base_barrier = AtomicUnits.hartree_to_ev(self.diode.bottom_work_fxn - self.diode.fe_chi)
        # Tsymbal-Kohlstedt polarization dipole shift (Eq. 5, Eq. 6). The
        # leading minus is uniform; asymmetry comes from sigma_pol,top = -P
        # (Eq. 1) vs sigma_pol,bot = +P (Eq. 2). Use the *metal*'s optical
        # permittivity for the screening dielectric (electronic-only response,
        # Mehta-Silverman-Jacobs 1973). For Ti/Al, kappa_inf,metal ~ 1 unless
        # materials.py supplies a better Drude-tail value.
        delta_phi_top_ev = self._polarization_dipole_shift_ev(
            state.sigma_pol_top_au,
            self.diode.top_screening_len,
            self.diode.top_eps_inf,
            self.parameters.fermi_level_pinning_top_v,
        )
        delta_phi_bot_ev = self._polarization_dipole_shift_ev(
            state.sigma_pol_bot_au,
            self.diode.bottom_screening_len,
            self.diode.bottom_eps_inf,
            self.parameters.fermi_level_pinning_bot_v,
        )
        # Image-force kappa is kappa_inf of the *injection-side dielectric*
        # (Eq. 15): top inj. = IL (HfO_x); bottom inj. = FE (AlScN).
        top_image_kappa = (
            self.diode.insulator_eps_inf if self.diode.insulator_thickness != 0 else self.diode.fe_eps_inf
        )
        top_barrier = max(
            top_base_barrier
            + delta_phi_top_ev
            - self.parameters.thermionic_top_schottky_scale
            * self._schottky_lowering_ev(top_field_v_per_m, top_image_kappa),
            0.0,
        )
        bottom_barrier = max(
            bottom_base_barrier
            + delta_phi_bot_ev
            - self.parameters.thermionic_bottom_schottky_scale
            * self._schottky_lowering_ev(bottom_field_v_per_m, self.diode.fe_eps_inf),
            0.0,
        )
        return top_barrier, bottom_barrier

    def _thermionic_current(self, state: BiasPointState):
        if not self.enable_thermionic:
            return 0.0, *self._barriers_ev(state)

        top_barrier_ev, bottom_barrier_ev = self._barriers_ev(state)
        top_offset_ev = (
            self.parameters.thermionic_barrier_offset_ev + self.parameters.thermionic_top_barrier_offset_ev
        )
        bottom_offset_ev = (
            self.parameters.thermionic_barrier_offset_ev + self.parameters.thermionic_bottom_barrier_offset_ev
        )
        top_barrier_ev = max(top_barrier_ev - top_offset_ev, self.parameters.thermionic_top_barrier_floor_ev, 0.0)
        bottom_barrier_ev = max(
            bottom_barrier_ev - bottom_offset_ev,
            self.parameters.thermionic_bottom_barrier_floor_ev,
            0.0,
        )
        temperature = self.parameters.temperature_k
        k_t_ev = constants.k * temperature / constants.e
        top_prefactor = (
            self.parameters.thermionic_prefactor_scale
            * self.parameters.thermionic_top_prefactor_scale
            * self.parameters.richardson_constant_a_per_m2_k2
            * self.diode.top_m_eff
            * temperature**2
        )
        bottom_prefactor = (
            self.parameters.thermionic_prefactor_scale
            * self.parameters.thermionic_bottom_prefactor_scale
            * self.parameters.richardson_constant_a_per_m2_k2
            * self.diode.bottom_m_eff
            * temperature**2
        )
        # Tung barrier-inhomogeneity enhancement (Tung 1992 PRB 45, 13509). A Gaussian
        # distribution of barrier heights N(phi_bar, sigma^2) over the device area
        # raises the integrated TE current by exp(sigma^2 / (2 (kT)^2)), which can be
        # many orders of magnitude. We cap the exponent at 80 (~exp(80) ~ 5.5e34)
        # to avoid floating-point overflow in unphysical fit excursions; sigma_phi
        # values larger than ~kT*sqrt(160) ~ 0.32 eV at 300K should be treated with
        # suspicion since they correspond to the high-end of the literature range.
        sigma_top_ev = max(self.parameters.te_barrier_sigma_top_ev, 0.0)
        sigma_bot_ev = max(self.parameters.te_barrier_sigma_bot_ev, 0.0)
        kt_safe = max(k_t_ev, 1e-9)
        top_inhom_log = min(0.5 * (sigma_top_ev / kt_safe) ** 2, 80.0)
        bot_inhom_log = min(0.5 * (sigma_bot_ev / kt_safe) ** 2, 80.0)
        top_injection = top_prefactor * np.exp(-top_barrier_ev / kt_safe + top_inhom_log)
        bottom_injection = bottom_prefactor * np.exp(-bottom_barrier_ev / kt_safe + bot_inhom_log)
        voltage_sign = self._voltage_sign(state.voltage_v)
        if voltage_sign > 0:
            current = top_injection
        elif voltage_sign < 0:
            current = -bottom_injection
        else:
            current = top_injection - bottom_injection
        return float(current), top_barrier_ev, bottom_barrier_ev

    def _wkb_transmission(self, x_au, barrier_ev, effective_mass):
        barrier_ev = np.maximum(np.asarray(barrier_ev, dtype=float), 0.0)
        if np.allclose(barrier_ev, 0.0):
            return 1.0

        x_m = AtomicUnits.bohr_to_m(np.asarray(x_au, dtype=float))
        dx = np.diff(x_m)
        if dx.size == 0:
            return 1.0

        barrier_mid_j = AtomicUnits.ev_to_joule(0.5 * (barrier_ev[:-1] + barrier_ev[1:]))
        mass_mid = 0.5 * (effective_mass[:-1] + effective_mass[1:]) * constants.m_e
        exponent_density = np.sqrt(2.0 * mass_mid * barrier_mid_j) / constants.hbar
        exponent = 2.0 * np.sum(exponent_density * dx)
        return float(np.exp(-np.clip(exponent, 0.0, 700.0)))

    def _tunneling_profile(self, state: BiasPointState):
        profile = self.potential.barrier_region_profile(
            fe_polarization=state.polarization_au,
            v_diff=state.voltage_au,
            num_points=self.parameters.barrier_sampling_points,
        )
        barrier_ev = np.maximum(profile["total_potential_ev"], 0.0)
        return profile, barrier_ev

    def _tunneling_current(self, state: BiasPointState):
        if not self.enable_tunneling:
            return 0.0, 0.0

        profile, barrier_ev = self._tunneling_profile(state)
        transmission = self._wkb_transmission(
            x_au=profile["x_au"],
            barrier_ev=barrier_ev,
            effective_mass=profile["effective_mass"],
        )
        flux_top = self._carrier_flux(self.diode.top_n0, self.diode.top_fermi_e, self.diode.top_m_eff)
        flux_bottom = self._carrier_flux(
            self.diode.bottom_n0, self.diode.bottom_fermi_e, self.diode.bottom_m_eff
        )
        # Direction follows the (signed) field across the FE: P-driven
        # depolarizing field at V_app = 0 produces nonzero tunneling current
        # of the correct sign (polarization_barrier_coupling.md Eq. 12).
        # When E_FE -> 0 we fall back to V_app to break degeneracy.
        bias_drive = self._field_sign(state.fe_field_au, fallback_voltage_v=state.voltage_v)
        current = (
            self.parameters.tunneling_prefactor_scale
            * 0.5
            * (flux_top + flux_bottom)
            * transmission
            * bias_drive
        )
        return float(current), transmission

    def _poole_frenkel_component(self, state: BiasPointState, region: str):
        # PF lives inside the bulk dielectric. polarization_barrier_coupling.md
        # sec 7 (PF row): Delta_phi_PF = sqrt(q^3 E / (pi eps0 kappa_inf)).
        # Note: pi (not 4 pi) in the denominator and kappa_inf (not kappa_static).
        # Trap depth: insulator-side PF must use insulator.trap_depth (HfO_x ~ 1.0 eV
        # / Al2O3 ~ 1.4 eV), not AlScN's fe_trap_depth (parameter audit punch #4, #8).
        if region == "il":
            field_au = state.il_field_au
            kappa_inf = self.diode.insulator_eps_inf
            region_scale = self.parameters.poole_frenkel_il_prefactor_scale
            default_trap_depth_au = self.diode.insulator_trap_depth
        else:
            field_au = state.fe_field_au
            kappa_inf = self.diode.fe_eps_inf
            region_scale = self.parameters.poole_frenkel_fe_prefactor_scale
            default_trap_depth_au = self.diode.fe_trap_depth
        # Carry the *signed* E_FE through the math: even at V_app = 0 with P != 0
        # there is a depolarizing field E_FE != 0 (Eq. 12) and PF must produce
        # nonzero current. Magnitude only enters the Frenkel sqrt.
        field_signed_v_per_m = self._field_to_v_per_m(field_au)
        field_v_per_m = abs(field_signed_v_per_m)
        trap_depth_ev = (
            self.parameters.poole_frenkel_trap_depth_ev
            if self.parameters.poole_frenkel_trap_depth_ev is not None
            else AtomicUnits.hartree_to_ev(default_trap_depth_au)
        )
        lowering_ev = self._poole_frenkel_lowering_ev(field_v_per_m, kappa_inf)
        activation_ev = max(trap_depth_ev - lowering_ev, 0.0)
        thermal_ev = constants.k * self.parameters.temperature_k / constants.e
        # J = q n mu E (drift) with field-dependent mobility e^{-phi_eff/kT}.
        # E carries its own sign (Eq. 12), so the current sign tracks E - not V_app.
        current = (
            self.parameters.poole_frenkel_prefactor_scale
            * region_scale
            * constants.e
            * self.parameters.poole_frenkel_trap_density_m3
            * self.parameters.poole_frenkel_mobility_m2_per_v_s
            * field_signed_v_per_m
            * np.exp(-activation_ev / max(thermal_ev, 1e-9))
        )
        return float(current), field_au, field_v_per_m, kappa_inf, lowering_ev, activation_ev

    def _poole_frenkel_current(self, state: BiasPointState):
        if not self.enable_poole_frenkel:
            return 0.0

        region = self.parameters.poole_frenkel_field_region
        if region == "both":
            fe_current, *_ = self._poole_frenkel_component(state, "fe")
            il_current, *_ = self._poole_frenkel_component(state, "il")
            return float(fe_current + il_current)
        current, *_ = self._poole_frenkel_component(state, region)
        return float(current)

    def _trap_assisted_tunneling_current(self, state: BiasPointState):
        if not self.enable_trap_assisted_tunneling:
            return 0.0, 0.0

        profile, barrier_ev = self._tunneling_profile(state)
        trap_depth_ev = AtomicUnits.hartree_to_ev(self.diode.fe_trap_depth)
        reduced_barrier_ev = np.maximum(barrier_ev - trap_depth_ev, 0.0)
        transmission = self._wkb_transmission(
            x_au=profile["x_au"],
            barrier_ev=reduced_barrier_ev,
            effective_mass=profile["effective_mass"],
        )
        # Field enhancement uses |E| (Schottky-like sqrt scaling); but the
        # current direction is set by the sign of E_FE, not V_app, so that
        # depolarizing-field-driven leakage at V_app = 0 (P != 0) is non-zero.
        field_signed_v_per_m = self._field_to_v_per_m(state.fe_field_au)
        field_enhancement = 1.0 + abs(field_signed_v_per_m) / 1e8
        direction = self._field_sign(state.fe_field_au, fallback_voltage_v=state.voltage_v)
        current = self.parameters.trap_assisted_prefactor_a_per_m2 * transmission * field_enhancement
        return float(current * direction), transmission

    def _sclc_current(self, state: BiasPointState):
        if not self.enable_sclc or self.diode.fe_thickness == 0:
            return 0.0

        thickness_m = AtomicUnits.bohr_to_m(self.diode.fe_thickness)
        # Mott-Gurney scales with V^2 (magnitude); current direction tracks the
        # sign of E_FE rather than V_app, so depolarizing-field leakage at zero
        # bias produces nonzero SCLC of the correct sign (Eq. 12).
        voltage_drop_v = abs(AtomicUnits.convert_back_volts(state.fe_field_au * self.diode.fe_thickness))
        direction = self._field_sign(state.fe_field_au, fallback_voltage_v=state.voltage_v)
        current = (
            self.parameters.sclc_prefactor_scale
            * 9.0
            / 8.0
            * constants.epsilon_0
            * self.diode.fe_k
            * self.parameters.sclc_mobility_m2_per_v_s
            * voltage_drop_v**2
            / max(thickness_m**3, 1e-30)
        )
        return float(direction * current)

    def _background_leakage_current(self, state: BiasPointState):
        conductance = self.parameters.background_conductance_a_per_m2_v
        if conductance <= 0:
            return 0.0, 0.0
        internal_drop_v = AtomicUnits.convert_back_volts(
            state.il_field_au * self.diode.insulator_thickness
            + state.fe_field_au * self.diode.fe_thickness
            + state.dl_field_au * self.diode.dl_thickness
        )
        current = conductance * internal_drop_v
        return float(current), float(internal_drop_v)

    def evaluate(self, state: BiasPointState):
        thermionic_a_per_m2, top_barrier_ev, bottom_barrier_ev = self._thermionic_current(state)
        tunneling_a_per_m2, tunneling_transmission = self._tunneling_current(state)
        poole_frenkel_a_per_m2 = self._poole_frenkel_current(state)
        tat_a_per_m2, trap_tunneling_transmission = self._trap_assisted_tunneling_current(state)
        sclc_a_per_m2 = self._sclc_current(state)
        background_leakage_a_per_m2, _ = self._background_leakage_current(state)
        total_a_per_m2 = (
            thermionic_a_per_m2
            + tunneling_a_per_m2
            + poole_frenkel_a_per_m2
            + tat_a_per_m2
            + sclc_a_per_m2
            + background_leakage_a_per_m2
        )
        return TransportBreakdown(
            voltage_v=state.voltage_v,
            temperature_k=self.parameters.temperature_k,
            thermionic_a_per_m2=float(thermionic_a_per_m2),
            tunneling_a_per_m2=float(tunneling_a_per_m2),
            poole_frenkel_a_per_m2=float(poole_frenkel_a_per_m2),
            trap_assisted_tunneling_a_per_m2=float(tat_a_per_m2),
            sclc_a_per_m2=float(sclc_a_per_m2),
            total_a_per_m2=float(total_a_per_m2),
            tunneling_transmission=float(tunneling_transmission),
            trap_tunneling_transmission=float(trap_tunneling_transmission),
            top_barrier_ev=float(top_barrier_ev),
            bottom_barrier_ev=float(bottom_barrier_ev),
            background_leakage_a_per_m2=float(background_leakage_a_per_m2),
        )

    def audit_state(self, state: BiasPointState):
        temperature = self.parameters.temperature_k
        thermal_voltage_v = self._thermal_voltage
        legacy_bias_drive = self._bias_drive(state.voltage_v)
        voltage_sign = self._voltage_sign(state.voltage_v)
        field_sign = self._field_sign(state.fe_field_au, state.voltage_v)

        top_field_v_per_m = abs(
            self._field_to_v_per_m(state.il_field_au if self.diode.insulator_thickness != 0 else state.fe_field_au)
        )
        bottom_field_v_per_m = abs(self._field_to_v_per_m(state.fe_field_au))
        top_affinity_ev = AtomicUnits.hartree_to_ev(
            self.diode.insulator_chi if self.diode.insulator_thickness != 0 else self.diode.fe_chi
        )
        fe_affinity_ev = AtomicUnits.hartree_to_ev(self.diode.fe_chi)
        top_work_function_ev = AtomicUnits.hartree_to_ev(self.diode.top_work_fxn)
        bottom_work_function_ev = AtomicUnits.hartree_to_ev(self.diode.bottom_work_fxn)
        top_fermi_ev = AtomicUnits.hartree_to_ev(self.diode.top_fermi_e)
        bottom_fermi_ev = AtomicUnits.hartree_to_ev(self.diode.bottom_fermi_e)

        top_base_barrier_ev = top_work_function_ev - top_affinity_ev
        bottom_base_barrier_ev = bottom_work_function_ev - fe_affinity_ev
        # Image-force: kappa_inf of the injection-side dielectric (Eq. 15).
        top_image_kappa = (
            self.diode.insulator_eps_inf if self.diode.insulator_thickness != 0 else self.diode.fe_eps_inf
        )
        top_schottky_lowering_ev = self._schottky_lowering_ev(top_field_v_per_m, top_image_kappa)
        bottom_schottky_lowering_ev = self._schottky_lowering_ev(bottom_field_v_per_m, self.diode.fe_eps_inf)
        top_schottky_effective_ev = self.parameters.thermionic_top_schottky_scale * top_schottky_lowering_ev
        bottom_schottky_effective_ev = (
            self.parameters.thermionic_bottom_schottky_scale * bottom_schottky_lowering_ev
        )
        # Tsymbal-Kohlstedt polarization dipole shift (Eq. 5, Eq. 6) with Eq. 4' saturation.
        delta_phi_top_ev = self._polarization_dipole_shift_ev(
            state.sigma_pol_top_au,
            self.diode.top_screening_len,
            self.diode.top_eps_inf,
            self.parameters.fermi_level_pinning_top_v,
        )
        delta_phi_bot_ev = self._polarization_dipole_shift_ev(
            state.sigma_pol_bot_au,
            self.diode.bottom_screening_len,
            self.diode.bottom_eps_inf,
            self.parameters.fermi_level_pinning_bot_v,
        )
        top_barrier_ev = max(
            top_base_barrier_ev + delta_phi_top_ev - top_schottky_effective_ev, 0.0
        )
        bottom_barrier_ev = max(
            bottom_base_barrier_ev + delta_phi_bot_ev - bottom_schottky_effective_ev, 0.0
        )
        top_offset_ev = (
            self.parameters.thermionic_barrier_offset_ev + self.parameters.thermionic_top_barrier_offset_ev
        )
        bottom_offset_ev = (
            self.parameters.thermionic_barrier_offset_ev + self.parameters.thermionic_bottom_barrier_offset_ev
        )
        top_barrier_effective_ev = max(top_barrier_ev - top_offset_ev, 0.0)
        bottom_barrier_effective_ev = max(bottom_barrier_ev - bottom_offset_ev, 0.0)
        top_barrier_effective_ev = max(top_barrier_effective_ev, self.parameters.thermionic_top_barrier_floor_ev)
        bottom_barrier_effective_ev = max(
            bottom_barrier_effective_ev, self.parameters.thermionic_bottom_barrier_floor_ev
        )

        k_t_ev = constants.k * temperature / constants.e
        top_prefactor_a_per_m2 = (
            self.parameters.thermionic_prefactor_scale
            * self.parameters.thermionic_top_prefactor_scale
            * self.parameters.richardson_constant_a_per_m2_k2
            * self.diode.top_m_eff
            * temperature**2
        )
        bottom_prefactor_a_per_m2 = (
            self.parameters.thermionic_prefactor_scale
            * self.parameters.thermionic_bottom_prefactor_scale
            * self.parameters.richardson_constant_a_per_m2_k2
            * self.diode.bottom_m_eff
            * temperature**2
        )
        top_injection_a_per_m2 = top_prefactor_a_per_m2 * np.exp(-top_barrier_effective_ev / max(k_t_ev, 1e-9))
        bottom_injection_a_per_m2 = bottom_prefactor_a_per_m2 * np.exp(
            -bottom_barrier_effective_ev / max(k_t_ev, 1e-9)
        )
        if voltage_sign > 0:
            thermionic_raw_a_per_m2 = top_injection_a_per_m2
        elif voltage_sign < 0:
            thermionic_raw_a_per_m2 = -bottom_injection_a_per_m2
        else:
            thermionic_raw_a_per_m2 = top_injection_a_per_m2 - bottom_injection_a_per_m2
        thermionic_a_per_m2 = thermionic_raw_a_per_m2 if self.enable_thermionic else 0.0

        profile, barrier_ev = self._tunneling_profile(state)
        tunneling_transmission = self._wkb_transmission(
            x_au=profile["x_au"],
            barrier_ev=barrier_ev,
            effective_mass=profile["effective_mass"],
        )
        flux_top_a_per_m2 = self._carrier_flux(self.diode.top_n0, self.diode.top_fermi_e, self.diode.top_m_eff)
        flux_bottom_a_per_m2 = self._carrier_flux(
            self.diode.bottom_n0, self.diode.bottom_fermi_e, self.diode.bottom_m_eff
        )
        # DT direction tracks signed E_FE (Eq. 12), not V_app (audit punch).
        tunneling_raw_a_per_m2 = (
            self.parameters.tunneling_prefactor_scale
            * 0.5
            * (flux_top_a_per_m2 + flux_bottom_a_per_m2)
            * tunneling_transmission
            * field_sign
        )
        tunneling_a_per_m2 = tunneling_raw_a_per_m2 if self.enable_tunneling else 0.0

        trap_depth_ev = AtomicUnits.hartree_to_ev(self.diode.fe_trap_depth)
        pf_region = self.parameters.poole_frenkel_field_region
        (
            pf_fe_raw_a_per_m2,
            pf_fe_field_au,
            pf_fe_field_v_per_m,
            pf_fe_k,
            pf_fe_lowering_ev,
            pf_fe_activation_ev,
        ) = self._poole_frenkel_component(state, "fe")
        (
            pf_il_raw_a_per_m2,
            pf_il_field_au,
            pf_il_field_v_per_m,
            pf_il_k,
            pf_il_lowering_ev,
            pf_il_activation_ev,
        ) = self._poole_frenkel_component(state, "il")
        if pf_region == "both":
            pf_field_au = pf_il_field_au
            pf_k = pf_il_k
            pf_field_v_per_m = pf_il_field_v_per_m
            pf_lowering_ev = min(pf_fe_lowering_ev, pf_il_lowering_ev)
            pf_activation_ev = min(pf_fe_activation_ev, pf_il_activation_ev)
            poole_frenkel_raw_a_per_m2 = pf_fe_raw_a_per_m2 + pf_il_raw_a_per_m2
        elif pf_region == "il":
            pf_field_au = pf_il_field_au
            pf_k = pf_il_k
            pf_field_v_per_m = pf_il_field_v_per_m
            pf_lowering_ev = pf_il_lowering_ev
            pf_activation_ev = pf_il_activation_ev
            poole_frenkel_raw_a_per_m2 = pf_il_raw_a_per_m2
        else:
            pf_field_au = pf_fe_field_au
            pf_k = pf_fe_k
            pf_field_v_per_m = pf_fe_field_v_per_m
            pf_lowering_ev = pf_fe_lowering_ev
            pf_activation_ev = pf_fe_activation_ev
            poole_frenkel_raw_a_per_m2 = pf_fe_raw_a_per_m2
        poole_frenkel_a_per_m2 = poole_frenkel_raw_a_per_m2 if self.enable_poole_frenkel else 0.0

        reduced_barrier_ev = np.maximum(barrier_ev - trap_depth_ev, 0.0)
        trap_tunneling_transmission = self._wkb_transmission(
            x_au=profile["x_au"],
            barrier_ev=reduced_barrier_ev,
            effective_mass=profile["effective_mass"],
        )
        field_enhancement = 1.0 + pf_field_v_per_m / 1e8
        # TAT direction tracks signed E_FE (Eq. 12), not V_app.
        trap_assisted_tunneling_raw_a_per_m2 = (
            self.parameters.trap_assisted_prefactor_a_per_m2
            * trap_tunneling_transmission
            * field_enhancement
            * field_sign
        )
        trap_assisted_tunneling_a_per_m2 = (
            trap_assisted_tunneling_raw_a_per_m2 if self.enable_trap_assisted_tunneling else 0.0
        )

        thickness_m = AtomicUnits.bohr_to_m(self.diode.fe_thickness)
        voltage_drop_v = abs(AtomicUnits.convert_back_volts(state.fe_field_au * self.diode.fe_thickness))
        # SCLC direction tracks signed E_FE (Eq. 12), not V_app.
        sclc_raw_a_per_m2 = (
            field_sign
            * self.parameters.sclc_prefactor_scale
            * 9.0
            / 8.0
            * constants.epsilon_0
            * self.diode.fe_k
            * self.parameters.sclc_mobility_m2_per_v_s
            * voltage_drop_v**2
            / max(thickness_m**3, 1e-30)
        )
        sclc_a_per_m2 = sclc_raw_a_per_m2 if self.enable_sclc else 0.0
        background_raw_a_per_m2, background_internal_drop_v = self._background_leakage_current(state)
        background_leakage_a_per_m2 = background_raw_a_per_m2

        total_a_per_m2 = (
            thermionic_a_per_m2
            + tunneling_a_per_m2
            + poole_frenkel_a_per_m2
            + trap_assisted_tunneling_a_per_m2
            + sclc_a_per_m2
            + background_leakage_a_per_m2
        )

        return {
            "voltage_v": float(state.voltage_v),
            "voltage_au": float(state.voltage_au),
            "temperature_k": float(temperature),
            "thermal_voltage_v": float(thermal_voltage_v),
            "legacy_bias_drive": float(legacy_bias_drive),
            "voltage_sign": float(voltage_sign),
            "field_sign": float(field_sign),
            "polarization_uc_cm2": float(state.polarization_uc_cm2),
            "polarization_au": float(state.polarization_au),
            "screening_charge_uc_cm2": float(state.screening_charge_uc_cm2),
            "screening_charge_c_m2": float(state.screening_charge_c_m2),
            "fe_field_mv_cm": float(state.fe_field_mv_cm),
            "il_field_mv_cm": float(state.il_field_mv_cm),
            "dl_field_mv_cm": float(state.dl_field_mv_cm),
            "fe_built_in_field_mv_cm": float(state.fe_built_in_field_mv_cm),
            "il_built_in_field_mv_cm": float(state.il_built_in_field_mv_cm),
            "dl_built_in_field_mv_cm": float(state.dl_built_in_field_mv_cm),
            "top_work_function_ev": float(top_work_function_ev),
            "bottom_work_function_ev": float(bottom_work_function_ev),
            "top_fermi_ev": float(top_fermi_ev),
            "bottom_fermi_ev": float(bottom_fermi_ev),
            "top_affinity_ev": float(top_affinity_ev),
            "fe_affinity_ev": float(fe_affinity_ev),
            "insulator_k": float(self.diode.insulator_k),
            "fe_k": float(self.diode.fe_k),
            "top_screening_len_nm": float(AtomicUnits.bohr_to_nm(self.diode.top_screening_len)),
            "bottom_screening_len_nm": float(AtomicUnits.bohr_to_nm(self.diode.bottom_screening_len)),
            "insulator_thickness_nm": float(AtomicUnits.bohr_to_nm(self.diode.insulator_thickness)),
            "fe_thickness_nm": float(AtomicUnits.bohr_to_nm(self.diode.fe_thickness)),
            "top_m_eff": float(self.diode.top_m_eff),
            "bottom_m_eff": float(self.diode.bottom_m_eff),
            "insulator_m_eff": float(self.diode.insulator_m_eff),
            "fe_m_eff": float(self.diode.fe_m_eff),
            "trap_depth_ev": float(trap_depth_ev),
            "pf_trap_depth_ev": float(
                self.parameters.poole_frenkel_trap_depth_ev
                if self.parameters.poole_frenkel_trap_depth_ev is not None
                else trap_depth_ev
            ),
            "top_field_v_per_m": float(top_field_v_per_m),
            "bottom_field_v_per_m": float(bottom_field_v_per_m),
            "top_base_barrier_ev": float(top_base_barrier_ev),
            "bottom_base_barrier_ev": float(bottom_base_barrier_ev),
            "top_schottky_lowering_ev": float(top_schottky_lowering_ev),
            "bottom_schottky_lowering_ev": float(bottom_schottky_lowering_ev),
            "thermionic_top_schottky_scale": float(self.parameters.thermionic_top_schottky_scale),
            "thermionic_bottom_schottky_scale": float(self.parameters.thermionic_bottom_schottky_scale),
            "top_schottky_effective_ev": float(top_schottky_effective_ev),
            "bottom_schottky_effective_ev": float(bottom_schottky_effective_ev),
            "top_barrier_ev": float(top_barrier_ev),
            "bottom_barrier_ev": float(bottom_barrier_ev),
            "thermionic_barrier_offset_ev": float(self.parameters.thermionic_barrier_offset_ev),
            "thermionic_top_barrier_offset_ev": float(self.parameters.thermionic_top_barrier_offset_ev),
            "thermionic_bottom_barrier_offset_ev": float(self.parameters.thermionic_bottom_barrier_offset_ev),
            "thermionic_top_barrier_floor_ev": float(self.parameters.thermionic_top_barrier_floor_ev),
            "thermionic_bottom_barrier_floor_ev": float(self.parameters.thermionic_bottom_barrier_floor_ev),
            "top_total_barrier_offset_ev": float(top_offset_ev),
            "bottom_total_barrier_offset_ev": float(bottom_offset_ev),
            "top_barrier_effective_ev": float(top_barrier_effective_ev),
            "bottom_barrier_effective_ev": float(bottom_barrier_effective_ev),
            "richardson_constant_a_per_m2_k2": float(self.parameters.richardson_constant_a_per_m2_k2),
            "thermionic_prefactor_scale": float(self.parameters.thermionic_prefactor_scale),
            "thermionic_top_prefactor_scale": float(self.parameters.thermionic_top_prefactor_scale),
            "thermionic_bottom_prefactor_scale": float(self.parameters.thermionic_bottom_prefactor_scale),
            "top_prefactor_a_per_m2": float(top_prefactor_a_per_m2),
            "bottom_prefactor_a_per_m2": float(bottom_prefactor_a_per_m2),
            "top_injection_a_per_m2": float(top_injection_a_per_m2),
            "bottom_injection_a_per_m2": float(bottom_injection_a_per_m2),
            "enable_thermionic": int(self.enable_thermionic),
            "thermionic_raw_a_per_m2": float(thermionic_raw_a_per_m2),
            "thermionic_a_per_m2": float(thermionic_a_per_m2),
            "barrier_sampling_points": int(self.parameters.barrier_sampling_points),
            "barrier_profile_min_ev": float(np.min(barrier_ev)),
            "barrier_profile_max_ev": float(np.max(barrier_ev)),
            "effective_mass_min": float(np.min(profile["effective_mass"])),
            "effective_mass_max": float(np.max(profile["effective_mass"])),
            "flux_top_a_per_m2": float(flux_top_a_per_m2),
            "flux_bottom_a_per_m2": float(flux_bottom_a_per_m2),
            "tunneling_prefactor_scale": float(self.parameters.tunneling_prefactor_scale),
            "enable_tunneling": int(self.enable_tunneling),
            "tunneling_transmission": float(tunneling_transmission),
            "tunneling_raw_a_per_m2": float(tunneling_raw_a_per_m2),
            "tunneling_a_per_m2": float(tunneling_a_per_m2),
            "pf_field_v_per_m": float(pf_field_v_per_m),
            "poole_frenkel_field_region": pf_region,
            "poole_frenkel_region_k": float(pf_k),
            "pf_lowering_ev": float(pf_lowering_ev),
            "pf_activation_ev": float(pf_activation_ev),
            "pf_fe_field_v_per_m": float(pf_fe_field_v_per_m),
            "pf_il_field_v_per_m": float(pf_il_field_v_per_m),
            "pf_fe_lowering_ev": float(pf_fe_lowering_ev),
            "pf_il_lowering_ev": float(pf_il_lowering_ev),
            "pf_fe_activation_ev": float(pf_fe_activation_ev),
            "pf_il_activation_ev": float(pf_il_activation_ev),
            "poole_frenkel_mobility_m2_per_v_s": float(self.parameters.poole_frenkel_mobility_m2_per_v_s),
            "poole_frenkel_trap_density_m3": float(self.parameters.poole_frenkel_trap_density_m3),
            "poole_frenkel_prefactor_scale": float(self.parameters.poole_frenkel_prefactor_scale),
            "poole_frenkel_fe_prefactor_scale": float(self.parameters.poole_frenkel_fe_prefactor_scale),
            "poole_frenkel_il_prefactor_scale": float(self.parameters.poole_frenkel_il_prefactor_scale),
            "enable_poole_frenkel": int(self.enable_poole_frenkel),
            "poole_frenkel_fe_raw_a_per_m2": float(pf_fe_raw_a_per_m2),
            "poole_frenkel_il_raw_a_per_m2": float(pf_il_raw_a_per_m2),
            "poole_frenkel_raw_a_per_m2": float(poole_frenkel_raw_a_per_m2),
            "poole_frenkel_a_per_m2": float(poole_frenkel_a_per_m2),
            "tat_reduced_barrier_min_ev": float(np.min(reduced_barrier_ev)),
            "tat_reduced_barrier_max_ev": float(np.max(reduced_barrier_ev)),
            "trap_assisted_prefactor_a_per_m2": float(self.parameters.trap_assisted_prefactor_a_per_m2),
            "field_enhancement": float(field_enhancement),
            "trap_tunneling_transmission": float(trap_tunneling_transmission),
            "enable_trap_assisted_tunneling": int(self.enable_trap_assisted_tunneling),
            "trap_assisted_tunneling_raw_a_per_m2": float(trap_assisted_tunneling_raw_a_per_m2),
            "trap_assisted_tunneling_a_per_m2": float(trap_assisted_tunneling_a_per_m2),
            "sclc_mobility_m2_per_v_s": float(self.parameters.sclc_mobility_m2_per_v_s),
            "sclc_prefactor_scale": float(self.parameters.sclc_prefactor_scale),
            "sclc_voltage_drop_v": float(voltage_drop_v),
            "sclc_thickness_m": float(thickness_m),
            "enable_sclc": int(self.enable_sclc),
            "sclc_raw_a_per_m2": float(sclc_raw_a_per_m2),
            "sclc_a_per_m2": float(sclc_a_per_m2),
            "background_conductance_a_per_m2_v": float(self.parameters.background_conductance_a_per_m2_v),
            "background_internal_drop_v": float(background_internal_drop_v),
            "background_raw_a_per_m2": float(background_raw_a_per_m2),
            "background_a_per_m2": float(background_leakage_a_per_m2),
            "total_a_per_m2": float(total_a_per_m2),
        }
