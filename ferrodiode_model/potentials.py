"""
potentials.py
Electrostatic + barrier potentials in a ferroelectric diode
"""

def thomas_fermi_length(kappa: float, Ef: float, rho0: float, eps0: float, e: float) -> float:
    """
    Compute Thomas-Fermi screening length (lambda).
    Args:
        kappa: dielectric constant of electrode
        Ef: Fermi energy [J]
        rho0: free electron density [m^-3]
        eps0: vacuum permittivity [F/m]
        e: elementary charge [C]
    Returns:
        lambda [m]
    """
    return ( (2 * eps0 * Ef) / (3 * (e**2) * rho0) * kappa )**0.5


def screening_charge(Pfe: float, Pdl: float, dfe: float, ddl: float, dox: float,
                     kfe: float, kdl: float, kox: float,
                     lam1: float, lam2: float, k1: float, k2: float,
                     Vext: float, eps0: float) -> float:
    """
    Solve for screening charge density (sigma_s).

    Args:
        Pfe: polarization in ferroelectric [C/m^2]
        Pdl: polarization in dead layer [C/m^2]
        dfe, ddl, dox: thickness of FE, DL, OX [m]
        kfe, kdl, kox: dielectric constants of FE, DL, OX
        lam1, lam2: screening lengths of electrodes [m]
        k1, k2: dielectric constants of electrodes
        Vext: external applied voltage [V]
        eps0: vacuum permittivity [F/m]

    Returns:
        sigma_s [C/m^2]
    """
    num = (Pdl * ddl / kdl) + (Pfe * dfe / kfe) + (eps0 * Vext)
    den = (lam1 / k1) + (lam2 / k2) + (ddl / kdl) + (dfe / kfe) + (dox / kox)
    return num / den


def depolarization_fields(sigma_s: float, Pfe: float, Pdl: float,
                          kfe: float, kdl: float, kox: float,
                          eps0: float) -> dict:
    """
    Compute depolarization fields in FE, DL, and OX layers.

    Args:
        sigma_s: screening charge density [C/m^2]
        Pfe: polarization in ferroelectric [C/m^2]
        Pdl: polarization in dead layer [C/m^2]
        kfe, kdl, kox: dielectric constants of FE, DL, OX
        eps0: vacuum permittivity [F/m]

    Returns:
        dict with keys:
            'Efe' [V/m] - field in ferroelectric
            'Edl' [V/m] - field in dead layer
            'Eox' [V/m] - field in oxide
    """
    Efe = (sigma_s - Pfe) / (kfe * eps0)
    Edl = (sigma_s - Pdl) / (kdl * eps0)
    Eox = sigma_s / (kox * eps0)
    return {"Efe": Efe, "Edl": Edl, "Eox": Eox}


def barrier_potential(phi1: float, phi2: float,
                      chi_fe: float, chi_ox: float) -> dict:
    """
    Compute barrier potentials at FE, DL, and OX interfaces.

    Args:
        phi1: work function of left electrode [eV]
        phi2: work function of right electrode [eV]
        chi_fe: electron affinity of FE [eV]
        chi_ox: electron affinity of oxide [eV]

    Returns:
        dict with keys:
            'Vb_fe' [eV] - barrier at FE interface
            'Vb_dl' [eV] - barrier at DL interface (same as FE)
            'Vb_ox' [eV] - barrier at oxide interface
    """
    Vb_fe = phi1 - chi_fe
    Vb_dl = phi1 - chi_fe
    Vb_ox = phi2 - chi_ox
    return {"Vb_fe": Vb_fe, "Vb_dl": Vb_dl, "Vb_ox": Vb_ox}


def total_potential(Vp: float, Vb: float) -> float:
    """
    Combine electrostatic and barrier potentials.

    Args:
        Vp: electrostatic potential [V]
        Vb: barrier potential [V or eV, must be consistent with Vp]
    Returns:
        total potential [same units as inputs]
    """
    return Vp + Vb
