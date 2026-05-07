class MetalElectrode:
    def __init__(self, n0, e_f, m_eff, k, w_f, name, screening_len):
        """
        :param n0: electron density
        :param e_f: fermi energy
        :param m_eff: effective mass of electron
        :param k: static (lattice) dielectric constant
        :param w_f: metal work function
        :param name: material name

        """
        self.n0 = n0
        self.e_f = e_f
        self.m_eff = m_eff
        self.k = k
        self.w_f = w_f
        self.name = name
        self.screening_len = screening_len



class Insulator:
    def __init__(self, k, chi, m_eff, name, breakdown_field, eps_inf=None, trap_depth=None):
        """
        :param k: static (low-frequency) dielectric constant
        :param chi: electron affinity
        :param m_eff: effective mass of electron
        :param name: material name
        :param breakdown_field: dielectric breakdown field of insulator
        :param eps_inf: optical / high-frequency relative permittivity (kappa_infinity);
            authoritative values in docs/parameter_audit_20260506.md.
        :param trap_depth: deep-trap depth (Hartree) for PF / TAT; see
            docs/parameter_audit_20260506.md.
        """
        self.k = k
        self.chi = chi
        self.m_eff = m_eff
        self.name = name
        self.breakdown_field = breakdown_field
        self.eps_inf = eps_inf
        self.trap_depth = trap_depth


class Ferroelectric:
    def __init__(self, p_r, chi, m_eff, k, name, trap_depth, eps_inf=None):
        """
        Describes a ferroelectric material
        :param p_r: remnant polarization
        :param chi: electron affinity
        :param m_eff: effective mass of electron
        :param k: static (low-frequency) dielectric constant
        :param name: material name
        :param trap_depth: deep-trap depth (Hartree) for PF / TAT
        :param eps_inf: optical / high-frequency relative permittivity (kappa_infinity);
            authoritative values in docs/parameter_audit_20260506.md.
        """
        self.k = k
        self.p_r = p_r
        self.chi = chi
        self.m_eff = m_eff
        self.name = name
        self.trap_depth = trap_depth
        self.eps_inf = eps_inf
