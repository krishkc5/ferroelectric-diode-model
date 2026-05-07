from material_types import *
from atomicunits import AtomicUnits


class Materials:
    baseline_titanium_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.33),
        screening_len=0.1,
        k=1,
        m_eff=1,
        name="Ti"
    )

    baseline_palladium_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(5.3),
        screening_len=0.1,
        k=1,
        m_eff=1,
        name="Pd"
    )

    baseline_aluminum_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.28),
        screening_len=0.1,
        k=1,
        m_eff=1,
        name="Al"
    )

    baseline_test_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.08),
        screening_len=1,
        k=1,
        m_eff=1,
        name="Test"
    )

    baseline_test_electrode_2 = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.33),
        screening_len=1,
        k=1,
        m_eff=1,
        name="Test2"
    )

    baseline_al2o3 = Insulator(
        k=9.3,
        # chi -> 1.2 eV (Afanas'ev, JAP 91, 3079, 2002). punch #5.
        chi=AtomicUnits.ev_to_hartree(1.2),
        m_eff=0.4,
        name="$Al_2O_3$",
        breakdown_field=0,
        # eps_inf = 3.1 (Filatova & Konashuk, J. Phys. Chem. C 119, 20755, 2015, Table 1).
        # docs/parameter_audit_20260506.md punch #7.
        eps_inf=3.1,
        # trap_depth = 1.4 eV (Liu et al., APL 100, 192905, 2012).
        # docs/parameter_audit_20260506.md punch #8.
        trap_depth=AtomicUnits.ev_to_hartree(1.4),
    )

    baseline_alscn = Ferroelectric(
        k=16,
        # 2.0 eV chi (Casamento APL 120, 152901, 2022). punch #1.
        chi=AtomicUnits.ev_to_hartree(2.0),
        m_eff=0.3,
        p_r=AtomicUnits.convert_polarization(110),
        name="AlScN",
        trap_depth=AtomicUnits.ev_to_hartree(0.8),
        # eps_inf = 4.6 (Ambacher 2021; Deng 2013). punch #2.
        eps_inf=4.6,
    )

    baseline_hfo2 = Insulator(
        k=16.64,
        chi=AtomicUnits.ev_to_hartree(2.0),
        m_eff=0.11,
        name="$HfO_2$",
        breakdown_field=0,
        # eps_inf = 4.0 (Robertson, RPP 69, 327, 2006, Table 4). punch #3.
        eps_inf=4.0,
        # trap_depth = 1.0 eV (Foster, PRB 65, 174117, 2002, Fig. 9). punch #4.
        trap_depth=AtomicUnits.ev_to_hartree(1.0),
    )

    titanium_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.33),
        screening_len=AtomicUnits.nm_to_bohr(0.045),
        k=1,
        m_eff=1,
        name="Ti"
    )

    palladium_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(5.3),
        screening_len=AtomicUnits.nm_to_bohr(0.045),
        k=1,
        m_eff=1,
        name="Pd"
    )

    chromium_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.5),
        screening_len=AtomicUnits.nm_to_bohr(0.045),
        k=1,
        m_eff=1,
        name="Cr"
    )

    aluminum_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.28),
        screening_len=AtomicUnits.nm_to_bohr(0.045),
        k=1,
        m_eff=1,
        name="Al"
    )

    test_electrode = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.08),
        screening_len=1,
        k=1,
        m_eff=1,
        name="Test"
    )

    test_electrode_2 = MetalElectrode(
        n0=AtomicUnits.convert_density(7e28),
        e_f=AtomicUnits.ev_to_hartree(4.354),
        w_f=AtomicUnits.ev_to_hartree(4.33),
        screening_len=1,
        k=1,
        m_eff=1,
        name="Test2"
    )

    al2o3 = Insulator(
        # k -> 9.0 (Groner et al., Chem. Mater. 16, 639, 2004, Fig. 9; ALD a-Al2O3
        # accepted range 8.0-10.0). docs/parameter_audit_20260506.md punch #6.
        k=9.0,
        # chi -> 1.2 eV (Afanas'ev et al., J. Appl. Phys. 91, 3079, 2002, p. 3082).
        # docs/parameter_audit_20260506.md punch #5.
        chi=AtomicUnits.ev_to_hartree(1.2),
        m_eff=0.4,
        name="$Al_2O_3$",
        breakdown_field=0,
        # eps_inf = 3.1 (Filatova & Konashuk, J. Phys. Chem. C 119, 20755, 2015, Table 1).
        # docs/parameter_audit_20260506.md punch #7.
        eps_inf=3.1,
        # trap_depth = 1.4 eV O-vacancy (Liu et al., APL 100, 192905, 2012, Fig. 4).
        # docs/parameter_audit_20260506.md punch #8.
        trap_depth=AtomicUnits.ev_to_hartree(1.4),
    )

    alscn = Ferroelectric(
        k=18,
        # chi -> 2.0 eV (Casamento et al., APL 120, 152901, 2022, p.152901-3;
        # Wang et al., IEDM 2020, p. 31.5.2).
        # docs/parameter_audit_20260506.md punch #1.
        chi=AtomicUnits.ev_to_hartree(2.0),
        m_eff=0.3,
        p_r=AtomicUnits.convert_polarization(110),
        name="AlScN",
        trap_depth=AtomicUnits.ev_to_hartree(0.8),
        # eps_inf = 4.6 (Ambacher, J. Appl. Phys. 130, 045102, 2021, Table III;
        # Deng et al., APL 102, 112103, 2013, p.112103-2).
        # docs/parameter_audit_20260506.md punch #2.
        eps_inf=4.6,
    )

    hfo2 = Insulator(
        k=18,
        chi=AtomicUnits.ev_to_hartree(2.3),
        m_eff=0.11,
        name="$HfO_2$",
        breakdown_field=0,
        # eps_inf = 4.0 (Robertson, Rep. Prog. Phys. 69, 327, 2006, Table 4).
        # docs/parameter_audit_20260506.md punch #3.
        eps_inf=4.0,
        # trap_depth = 1.0 eV O-vacancy canonical
        # (Foster et al., PRB 65, 174117, 2002, Fig. 9).
        # docs/parameter_audit_20260506.md punch #4.
        trap_depth=AtomicUnits.ev_to_hartree(1.0),
    )

    tio2 = Insulator(
        k=89.8,
        chi=AtomicUnits.ev_to_hartree(1.59),
        m_eff=5,
        name="$TiO_2$",
        breakdown_field=0
    )
