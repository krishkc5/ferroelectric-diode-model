import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from atomicunits import AtomicUnits
from fed import FerroelectricDiode
from materials import Materials
from potential import Potential
from preisach import Ferroelectric as PreisachFerroelectric
from self_consistent_solver import SelfConsistentSolver


ZERO_STRAIN_CA_MEAN = 1.54
DEFAULT_CA_STD = 0.04
DEFAULT_FE_THICKNESS_NM = 45.0
DEFAULT_NUM_DOMAINS = 15000
DEFAULT_MAX_ITER = 5000
DEFAULT_SOLVER_THRESHOLD_UC_CM2 = 0.5
DEFAULT_CYCLES = 2
DEFAULT_TITLE = "FE Hysteresis vs Interlayer Thickness"


def resolve_metal_electrode(name: str):
    electrode_map = {
        "ti": Materials.baseline_titanium_electrode,
        "al": Materials.baseline_aluminum_electrode,
        "pd": Materials.baseline_palladium_electrode,
        "test": Materials.baseline_test_electrode,
        "test2": Materials.baseline_test_electrode_2,
    }
    return electrode_map[name]


def resolve_insulator(name: str):
    insulator_map = {
        "al2o3": Materials.baseline_al2o3,
        "hfox": Materials.baseline_hfo2,
    }
    return insulator_map[name]


def insulator_display_label(name: str) -> str:
    display_map = {
        "al2o3": "Al2O3",
        "hfox": "HfOx",
    }
    return display_map[name]


def build_figure_dir(fe_thickness_nm: float) -> str:
    if float(fe_thickness_nm).is_integer():
        return f"{int(fe_thickness_nm)}nm_figures"
    return f"{str(fe_thickness_nm).replace('.', 'p')}nm_figures"


def build_stack_filename(args: argparse.Namespace) -> str:
    top_label = resolve_metal_electrode(args.top_electrode).name
    insulator_label = insulator_display_label(args.insulator)
    ferroelectric_label = "AlScN"
    bottom_label = resolve_metal_electrode(args.bottom_electrode).name
    return f"{top_label}_{insulator_label}_{ferroelectric_label}_{bottom_label}_il_hysteresis.png"


def build_voltage_sweep(vmax: float, step: float) -> np.ndarray:
    n_points = int(round((2.0 * vmax) / step)) + 1
    va_forward = np.linspace(-vmax, vmax, n_points)
    va_reverse = np.linspace(vmax, -vmax, n_points)
    return np.concatenate([va_forward, va_reverse])


def repeat_voltage_sweep(single_cycle: np.ndarray, cycles: int) -> np.ndarray:
    if cycles < 1:
        raise ValueError("The hysteresis sweep needs at least one cycle.")
    return np.tile(single_cycle, cycles)


def initialize_domain_states(fe_model: PreisachFerroelectric, seed: int) -> None:
    rng = np.random.default_rng(seed)
    state_values = rng.choice([-1, 1], size=fe_model.num_domains)
    fe_model.set_states(state_values)


def summarize_loop(il_nm: float, single_cycle: np.ndarray, p_curve: list[float]) -> None:
    midpoint = len(single_cycle) // 2
    forward_v = single_cycle[:midpoint]
    reverse_v = single_cycle[midpoint:]
    forward_p = np.asarray(p_curve[:midpoint], dtype=float)
    reverse_p = np.asarray(p_curve[midpoint:], dtype=float)
    zero_forward_index = int(np.argmin(np.abs(forward_v)))
    zero_reverse_index = int(np.argmin(np.abs(reverse_v)))
    print(
        f"[loop-summary] IL={il_nm:g}nm "
        f"Pmin={np.min(p_curve):.3f} Pmax={np.max(p_curve):.3f} "
        f"P(0V,fwd)={forward_p[zero_forward_index]:.3f} "
        f"P(0V,rev)={reverse_p[zero_reverse_index]:.3f}"
    )


def format_stack_text(args: argparse.Namespace) -> str:
    il_values = ", ".join(f"{il:g}" for il in args.il_nm)
    top_electrode = resolve_metal_electrode(args.top_electrode)
    bottom_electrode = resolve_metal_electrode(args.bottom_electrode)
    insulator = resolve_insulator(args.insulator)
    alscn = Materials.baseline_alscn

    top_wf_ev = AtomicUnits.hartree_to_ev(top_electrode.w_f)
    top_lambda_nm = AtomicUnits.bohr_to_nm(top_electrode.screening_len)
    bottom_wf_ev = AtomicUnits.hartree_to_ev(bottom_electrode.w_f)
    bottom_lambda_nm = AtomicUnits.bohr_to_nm(bottom_electrode.screening_len)
    insulator_chi_ev = AtomicUnits.hartree_to_ev(insulator.chi)
    alscn_chi_ev = AtomicUnits.hartree_to_ev(alscn.chi)
    alscn_pr_uc_cm2 = AtomicUnits.convert_back_polarization(alscn.p_r)
    alscn_trap_depth_ev = AtomicUnits.hartree_to_ev(alscn.trap_depth)
    nominal_top = (
        f"{args.top_thickness_nm:g} nm nominal"
        if args.top_thickness_nm is not None
        else "nominal thickness not specified"
    )
    nominal_bottom = (
        f"{args.bottom_thickness_nm:g} nm nominal"
        if args.bottom_thickness_nm is not None
        else "nominal thickness not specified"
    )
    insulator_label = insulator_display_label(args.insulator)

    return "\n".join([
        "Stack",
        f"{top_electrode.name} | {insulator_label} (t_IL sweep) | AlScN ({args.fe_thickness_nm:g} nm) | {bottom_electrode.name}",
        "",
        "Thickness / Bias",
        f"t_IL sweep = {il_values} nm",
        f"t_FE = {args.fe_thickness_nm:g} nm",
        "dead layer = 0 nm",
        f"top metal = {nominal_top}",
        f"bottom metal = {nominal_bottom}",
        "metal thickness is annotation only in the current model",
        f"V_a sweep = -{args.vmax:g} -> +{args.vmax:g} -> -{args.vmax:g} V",
        f"step = {args.vstep:g} V, cycles = {DEFAULT_CYCLES}",
        "",
        "Materials",
        f"{top_electrode.name}: WF = {top_wf_ev:.2f} eV, lambda = {top_lambda_nm:.4f} nm",
        f"{bottom_electrode.name}: WF = {bottom_wf_ev:.2f} eV, lambda = {bottom_lambda_nm:.4f} nm",
        f"{insulator_label}: k = {insulator.k:g}, chi = {insulator_chi_ev:.2f} eV, m* = {insulator.m_eff:g}",
        f"AlScN: k = {alscn.k:g}, chi = {alscn_chi_ev:.2f} eV, m* = {alscn.m_eff:g}",
        f"AlScN: Pr = {alscn_pr_uc_cm2:.0f} uC/cm^2, trap depth = {alscn_trap_depth_ev:.2f} eV",
        "",
        "Preisach / Solver",
        f"domains = {args.num_domains}, seed = {args.seed}",
        "initial state = random, same seed for all IL",
        f"c/a mean = {ZERO_STRAIN_CA_MEAN:g}, c/a std = {args.ca_std:g}",
        "Ps(c/a) = 333.33(c/a) - 400 uC/cm^2",
        "Ec(c/a) = 3.16(c/a) - 1.1 MV/cm",
        "solver = bisection, built-in WF term = on",
        f"tol = {DEFAULT_SOLVER_THRESHOLD_UC_CM2:g} uC/cm^2, max_iter = {DEFAULT_MAX_ITER}",
    ])


def format_assumptions_text(args: argparse.Namespace) -> str:
    return format_stack_text(args)


def run_hysteresis_for_il(il_nm: float, single_cycle: np.ndarray, args: argparse.Namespace) -> list[float]:
    fe_model = PreisachFerroelectric(
        num_domains=args.num_domains,
        c_a_mean=ZERO_STRAIN_CA_MEAN,
        c_a_std=args.ca_std,
        p_s_mean=None,
        p_s_std=None,
        e_c_mean=None,
        e_c_std=None,
        seed=args.seed,
    )
    if np.min(fe_model.p_s_values) <= 0:
        min_ps = AtomicUnits.convert_back_polarization(np.min(fe_model.p_s_values))
        raise ValueError(
            "The default Preisach Ps(c/a) map produced non-positive local saturation polarization "
            f"(minimum {min_ps:.3f} uC/cm^2)."
        )

    initialize_domain_states(fe_model=fe_model, seed=args.seed)

    diode = FerroelectricDiode(
        insulator_thickness=AtomicUnits.nm_to_bohr(il_nm),
        fe_thickness=AtomicUnits.nm_to_bohr(args.fe_thickness_nm),
        dead_layer_thickness=0.0,
        top_electrode=resolve_metal_electrode(args.top_electrode),
        bottom_electrode=resolve_metal_electrode(args.bottom_electrode),
        insulator=resolve_insulator(args.insulator),
        ferroelectric=Materials.baseline_alscn,
        fe_model=fe_model,
    )

    ec_vals = AtomicUnits.atomic_units_to_Mv_per_cm(fe_model.e_c_values)
    print(f"Median Ec (MV/cm), IL={il_nm:g}nm: {np.median(ec_vals):.6f}")

    full_sweep = repeat_voltage_sweep(single_cycle, DEFAULT_CYCLES)
    potential = Potential(diode, AtomicUnits.convert_volts(full_sweep[0]))
    solver = SelfConsistentSolver(
        ferroelectric=fe_model,
        potential=potential,
        max_iter=DEFAULT_MAX_ITER,
        threshold=DEFAULT_SOLVER_THRESHOLD_UC_CM2,
    )

    p_uc_cm2 = []
    for voltage in full_sweep:
        polarization_atomic, _ = solver.solve(voltage)
        p_uc_cm2.append(AtomicUnits.convert_back_polarization(polarization_atomic))

    final_cycle = p_uc_cm2[-len(single_cycle):]
    if args.report_loop_summary:
        summarize_loop(il_nm=il_nm, single_cycle=single_cycle, p_curve=final_cycle)
    return final_cycle


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate the raw multidomain Preisach hysteresis plot versus IL thickness."
    )
    parser.add_argument("--il-nm", nargs="+", type=float, default=[0.0, 0.1, 0.2, 0.5, 0.8, 1.0, 2.0, 3.0],
                        help="Interlayer thicknesses in nm.")
    parser.add_argument("--vmax", type=float, default=20.0,
                        help="Maximum sweep voltage magnitude.")
    parser.add_argument("--vstep", type=float, default=0.25,
                        help="Voltage step for the sweep.")
    parser.add_argument("--num-domains", type=int, default=DEFAULT_NUM_DOMAINS,
                        help="Number of Preisach domains in the structural distribution.")
    parser.add_argument("--ca-std", type=float, default=DEFAULT_CA_STD,
                        help="Std dev of the c/a distribution.")
    parser.add_argument("--seed", type=int, default=0,
                        help="Random seed for both domain generation and initial random states.")
    parser.add_argument("--top-electrode", choices=["ti", "al", "pd", "test", "test2"], default="ti",
                        help="Top electrode used in the stack annotation and electrostatics.")
    parser.add_argument("--bottom-electrode", choices=["ti", "al", "pd", "test", "test2"], default="ti",
                        help="Bottom electrode used in the stack annotation and electrostatics.")
    parser.add_argument("--insulator", choices=["al2o3", "hfox"], default="al2o3",
                        help="Interlayer material used in the stack annotation and electrostatics.")
    parser.add_argument("--fe-thickness-nm", type=float, default=DEFAULT_FE_THICKNESS_NM,
                        help="Ferroelectric thickness in nm.")
    parser.add_argument("--top-thickness-nm", type=float, default=None,
                        help="Nominal top metal thickness in nm for figure annotation only.")
    parser.add_argument("--bottom-thickness-nm", type=float, default=None,
                        help="Nominal bottom metal thickness in nm for figure annotation only.")
    parser.add_argument("--output", default=None,
                        help="Output image path.")
    parser.add_argument("--dpi", type=int, default=200, help="Figure DPI.")
    parser.add_argument("--title", default=DEFAULT_TITLE, help="Figure title.")
    parser.add_argument("--show-assumptions", action=argparse.BooleanOptionalAction, default=True,
                        help="Show the stack and simulation parameters in the figure margin.")
    parser.add_argument("--report-loop-summary", action="store_true",
                        help="Print a compact raw-loop summary for each IL thickness.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.output is None:
        args.output = f"{build_figure_dir(args.fe_thickness_nm)}/{build_stack_filename(args)}"
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    single_cycle = build_voltage_sweep(args.vmax, args.vstep)

    fig, ax = plt.subplots(figsize=(13.2, 7.0), dpi=args.dpi)
    sorted_il_values = sorted(args.il_nm)
    if len(sorted_il_values) == 1:
        color_positions = [0.5]
    else:
        color_positions = np.linspace(0.0, 1.0, len(sorted_il_values))
    color_map = {
        il_nm: plt.cm.rainbow(position)
        for il_nm, position in zip(sorted_il_values, color_positions)
    }

    for il_nm in sorted(args.il_nm, reverse=True):
        p_curve = run_hysteresis_for_il(il_nm, single_cycle, args)
        color = color_map.get(float(il_nm), None)
        ax.plot(single_cycle, p_curve, linewidth=2.8, alpha=1.0, color=color, label=f"{il_nm:g}nm IL", zorder=10 - il_nm)

    ax.set_title(args.title, pad=12)
    ax.set_xlabel("Applied Voltage (V)")
    ax.set_ylabel(r"FE Polarization ($\mu C/cm^2$)")
    x_padding = max(0.5, 0.08 * args.vmax)
    x_limit = args.vmax + x_padding
    if args.vmax <= 10:
        x_tick_step = 2
    elif args.vmax <= 20:
        x_tick_step = 5
    else:
        x_tick_step = 10
    ax.set_xlim(-x_limit, x_limit)
    ax.set_ylim(-120, 120)
    ax.set_xticks(np.arange(-int(args.vmax), int(args.vmax) + 1, x_tick_step))
    ax.set_yticks(np.arange(-100, 101, 25))
    handles, labels = ax.get_legend_handles_labels()
    order = np.argsort([float(label.replace("nm IL", "")) for label in labels])
    ax.legend(
        [handles[i] for i in order],
        [labels[i] for i in order],
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
        facecolor="white",
        edgecolor="#d9d9d9",
        framealpha=0.9,
    )
    if args.show_assumptions:
        fig.text(
            0.73,
            0.95,
            format_assumptions_text(args),
            ha="left",
            va="top",
            fontsize=8.4,
            family="monospace",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d9d9d9", alpha=0.95),
        )
    ax.tick_params(axis="both", which="major", labelsize=11)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    fig.tight_layout(rect=(0, 0, 0.69 if args.show_assumptions else 0.84, 1))
    fig.savefig(args.output, dpi=args.dpi)
    plt.close(fig)


if __name__ == "__main__":
    main()
