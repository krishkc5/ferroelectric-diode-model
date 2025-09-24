"""
switching.py
Switching voltage analysis for a ferroelectric diode
"""

def switching_voltage(Pfe: float, Pdl: float,
                      dfe: float, ddl: float, dox: float,
                      kfe: float, kdl: float, kox: float,
                      lam1: float, lam2: float, k1: float, k2: float,
                      Ec: float, eps0: float) -> float:
    """
    Compute switching voltage Vsw.

    General case: includes electrode screening and dead layer.
    Special cases are recovered by parameter choices:
        - No screening: set lam1=lam2=0
        - No screening + no dead layer: set lam1=lam2=0, ddl=0, Pdl=0

    Args:
        Pfe, Pdl: polarizations [C/m^2]
        dfe, ddl, dox: thicknesses [m]
        kfe, kdl, kox: dielectric constants
        lam1, lam2: electrode screening lengths [m]
        k1, k2: dielectric constants of electrodes
        Ec: coercive field [V/m]
        eps0: vacuum permittivity [F/m]

    Returns:
        Vsw [V]
    """
    term1 = (Pfe * dfe) / (kfe * eps0)
    term2 = (Pdl * ddl) / (kdl * eps0)
    term3 = (dfe +
             (kfe / kdl) * ddl +
             (kfe / kox) * dox +
             (kfe / k1) * lam1 +
             (kfe / k2) * lam2)
    term4 = Ec - (Pfe / (kfe * eps0))
    return term1 + term2 + term3 * term4
