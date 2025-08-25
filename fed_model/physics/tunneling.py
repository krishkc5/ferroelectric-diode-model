import numpy as np
from scipy.integrate import quad
from scipy.constants import hbar, m_e, e

def wkb_transmission_coefficient(x, V, E):
    """Calculate tunneling probability using WKB approximation."""
    V_si = V * e  
    E_si = E * e
    
    mask = V_si > E_si
    if not np.any(mask):
        return 1.0
    
    dx = x[1] - x[0] if len(x) > 1 else 1e-10
    k_vals = np.sqrt(2 * m_e * np.abs(V_si - E_si) / hbar**2)
    
    integral = np.trapz(k_vals[mask], x[mask])
    T = np.exp(-2 * integral)
    
    return min(T, 1.0)

class TunnelingCalculator:
    def __init__(self):
        self.m_eff = m_e
    
    def set_effective_mass(self, m_eff_ratio):
        self.m_eff = m_eff_ratio * m_e
    
    def transmission_wkb(self, x, V, E):
        return wkb_transmission_coefficient(x, V, E)
    
    def current_density(self, x, V, voltage, temperature=300):
        kT = temperature * 8.617e-5  # eV
        E_max = 5 * kT
        energies = np.linspace(0, E_max, 100)
        
        def integrand(E):
            T = self.transmission_wkb(x, V, E)
            f1 = 1 / (1 + np.exp(E / kT))
            f2 = 1 / (1 + np.exp((E - voltage) / kT))
            return T * (f1 - f2)
        
        prefactor = e * m_e * kT / (2 * np.pi**2 * hbar**3)
        current_integral = quad(integrand, 0, E_max)[0]
        
        return prefactor * current_integral

tunnel_calc = TunnelingCalculator()
