import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from hysteresis_core import (
    HysteresisConfig,
    HysteresisSimulator,
    insulator_display_label,
    resolve_metal_electrode,
)
from sweeps import build_waypoint_sweep
from transport import TransportEvaluator, TransportModelParameters


DEFAULT_WAYPOINTS = [-20.0, 20.0, -20.0]
SEGMENT_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]


def build_figure_dir(fe_thickness_nm: float):
    if float(fe_thickness_nm).is_integer():
        return f"{int(fe_thickness_nm)}nm_figures"
    return f"{str(fe_thickness_nm).replace('.', 'p')}nm_figures"


def build_output_path(args: argparse.Namespace):
    if args.output is not None:
        return Path(args.output)
    top_label = resolve_metal_electrode(args.top_electrode).name
    bottom_label = resolve_metal_electrode(args.bottom_electrode).name
    vmax = max(abs(float(value)) for value in args.waypoints)
    vmax_label = f"{vmax:g}".replace(".", "p")
    if args.il_nm == 0:
        filename = f"{top_label}_AlScN_{bottom_label}_dciv_{vmax_label}V.png"
    else:
        insulator_label = insulator_display_label(args.insulator)
        filename = f"{top_label}_{insulator_label}_AlScN_{bottom_label}_dciv_{vmax_label}V.png"
    return Path(build_figure_dir(args.fe_thickness_nm)) / filename


def a_per_m2_to_uA_per_um2(current_density):
    return np.asarray(current_density, dtype=float) * 1e-6


def build_segment_slices(waypoints, step):
    slices = []
    start_index = 0
    for index in range(len(waypoints) - 1):
        span = abs(float(waypoints[index + 1]) - float(waypoints[index]))
        n_points = int(round(span / step)) + 1
        segment_length = n_points if index == 0 else n_points - 1
        stop_index = start_index + segment_length
        slices.append(slice(start_index, stop_index))
        start_index = stop_index
    return slices


def default_csv_path(output_path: Path):
    return output_path.with_suffix(".csv")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a quasistatic DC I-V sweep using the active multidomain Preisach hysteresis engine."
    )
    parser.add_argument("--waypoints", nargs="+", type=float, default=DEFAULT_WAYPOINTS)
    parser.add_argument("--vstep", type=float, default=0.25)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--il-nm", type=float, default=1.0)
    parser.add_argument("--fe-thickness-nm", type=float, default=45.0)
    parser.add_argument("--top-electrode", choices=["ti", "al", "pd", "test", "test2"], default="ti")
    parser.add_argument("--bottom-electrode", choices=["ti", "al", "pd", "test", "test2"], default="ti")
    parser.add_argument("--insulator", choices=["al2o3", "hfox"], default="al2o3")
    parser.add_argument("--num-domains", type=int, default=15000)
    parser.add_argument("--ca-std", type=float, default=0.04)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--output", default=None)
    parser.add_argument("--csv-output", default=None)
    parser.add_argument(
        "--title",
        default="Quasistatic Bipolar DC I-V with Multidomain Preisach Hysteresis",
    )
    parser.add_argument("--report-summary", action="store_true")
    parser.add_argument("--thermionic", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--tunneling", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--poole-frenkel", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--trap-assisted-tunneling", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--sclc", action=argparse.BooleanOptionalAction, default=False)
    return parser.parse_args()


def main():
    args = parse_args()
    voltage_trace = build_waypoint_sweep(args.waypoints, args.vstep)
    config = HysteresisConfig(
        num_domains=args.num_domains,
        c_a_std=args.ca_std,
        seed=args.seed,
        top_electrode_name=args.top_electrode,
        bottom_electrode_name=args.bottom_electrode,
        insulator_name=args.insulator,
        fe_thickness_nm=args.fe_thickness_nm,
    )
    simulator = HysteresisSimulator(il_nm=args.il_nm, config=config)
    states = simulator.simulate_trace(voltage_trace)

    transport = TransportEvaluator(
        diode=simulator.diode,
        potential=simulator.potential,
        parameters=TransportModelParameters(temperature_k=args.temperature_k),
        enable_thermionic=args.thermionic,
        enable_tunneling=args.tunneling,
        enable_poole_frenkel=args.poole_frenkel,
        enable_trap_assisted_tunneling=args.trap_assisted_tunneling,
        enable_sclc=args.sclc,
    )
    breakdowns = [transport.evaluate(state) for state in states]

    voltage_v = np.asarray([state.voltage_v for state in states], dtype=float)
    polarization_uc_cm2 = np.asarray([state.polarization_uc_cm2 for state in states], dtype=float)
    total_uA_um2 = a_per_m2_to_uA_per_um2([entry.total_a_per_m2 for entry in breakdowns])
    thermionic_uA_um2 = a_per_m2_to_uA_per_um2([entry.thermionic_a_per_m2 for entry in breakdowns])
    tunneling_uA_um2 = a_per_m2_to_uA_per_um2([entry.tunneling_a_per_m2 for entry in breakdowns])
    pf_uA_um2 = a_per_m2_to_uA_per_um2([entry.poole_frenkel_a_per_m2 for entry in breakdowns])
    tat_uA_um2 = a_per_m2_to_uA_per_um2([entry.trap_assisted_tunneling_a_per_m2 for entry in breakdowns])
    sclc_uA_um2 = a_per_m2_to_uA_per_um2([entry.sclc_a_per_m2 for entry in breakdowns])

    output_path = build_output_path(args)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path = Path(args.csv_output) if args.csv_output is not None else default_csv_path(output_path)
    segment_slices = build_segment_slices(args.waypoints, args.vstep)

    fig, axes = plt.subplots(3, 1, figsize=(12.0, 11.5), dpi=args.dpi, sharex=True)
    fig.suptitle(args.title, fontsize=15, y=0.98)

    for segment_index, segment_slice in enumerate(segment_slices):
        color = SEGMENT_COLORS[segment_index % len(SEGMENT_COLORS)]
        label = (
            f"{args.waypoints[segment_index]:g} -> {args.waypoints[segment_index + 1]:g} V"
            if segment_index < len(args.waypoints) - 1
            else None
        )
        axes[0].plot(
            voltage_v[segment_slice],
            polarization_uc_cm2[segment_slice],
            color=color,
            linewidth=2.2,
            label=label,
        )
    axes[0].set_ylabel("P_FE (uC/cm^2)")
    axes[0].grid(alpha=0.25)
    axes[0].legend(loc="lower right", fontsize=9)

    for segment_index, segment_slice in enumerate(segment_slices):
        color = SEGMENT_COLORS[segment_index % len(SEGMENT_COLORS)]
        axes[1].plot(
            voltage_v[segment_slice],
            total_uA_um2[segment_slice],
            color=color,
            linewidth=2.2,
        )
    total_abs_uA_um2 = np.maximum(np.abs(total_uA_um2), 1e-20)
    for segment_index, segment_slice in enumerate(segment_slices):
        color = SEGMENT_COLORS[segment_index % len(SEGMENT_COLORS)]
        label = (
            f"{args.waypoints[segment_index]:g} -> {args.waypoints[segment_index + 1]:g} V"
            if segment_index < len(args.waypoints) - 1
            else None
        )
        axes[1].semilogy(
            voltage_v[segment_slice],
            total_abs_uA_um2[segment_slice],
            color=color,
            linewidth=2.4,
            label=label,
        )
    axes[1].set_ylabel("|J_total| (uA/um^2)")
    axes[1].grid(alpha=0.25)
    axes[1].legend(loc="lower right", fontsize=9)

    mechanism_series = [
        ("Thermionic", thermionic_uA_um2, "#d62728", 2.6, 1.0),
        ("Poole-Frenkel", pf_uA_um2, "#9467bd", 2.6, 1.0),
        ("Tunneling", tunneling_uA_um2, "#7f7f7f", 1.4, 0.5),
        ("TAT", tat_uA_um2, "#bdbdbd", 1.4, 0.5),
        ("SCLC", sclc_uA_um2, "#17becf", 1.4, 0.5),
    ]
    for label, values, color, linewidth, alpha in mechanism_series:
        axes[2].semilogy(
            voltage_v,
            np.maximum(np.abs(values), 1e-20),
            label=label,
            linewidth=linewidth,
            color=color,
            alpha=alpha,
        )
    axes[2].set_xlabel("Applied Voltage (V)")
    axes[2].set_ylabel("|J_i| (uA/um^2)")
    axes[2].grid(alpha=0.25)
    axes[2].legend(loc="upper left", ncol=2)

    if args.il_nm == 0:
        stack_line = (
            f"{resolve_metal_electrode(args.top_electrode).name} | "
            f"AlScN ({args.fe_thickness_nm:g} nm) | "
            f"{resolve_metal_electrode(args.bottom_electrode).name}"
        )
    else:
        stack_line = (
            f"{resolve_metal_electrode(args.top_electrode).name} | "
            f"{insulator_display_label(args.insulator)} ({args.il_nm:g} nm) | "
            f"AlScN ({args.fe_thickness_nm:g} nm) | "
            f"{resolve_metal_electrode(args.bottom_electrode).name}"
        )
    stack_text = "\n".join(
        [
            "Stack",
            stack_line,
            "",
            "Sweep",
            " -> ".join(f"{value:g} V" for value in args.waypoints),
            f"step = {args.vstep:g} V, T = {args.temperature_k:g} K",
            "",
            "Mechanisms",
            f"thermionic = {args.thermionic}",
            f"tunneling = {args.tunneling}",
            f"PF = {args.poole_frenkel}",
            f"TAT = {args.trap_assisted_tunneling}",
            f"SCLC = {args.sclc}",
        ]
    )
    fig.text(
        0.74,
        0.95,
        stack_text,
        ha="left",
        va="top",
        fontsize=8.8,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d9d9d9", alpha=0.95),
    )
    fig.tight_layout(rect=(0.0, 0.0, 0.72, 0.96))
    fig.savefig(output_path, dpi=args.dpi)
    plt.close(fig)

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "index",
                "segment_index",
                "voltage_v",
                "polarization_uc_cm2",
                "polarization_change_uc_cm2",
                "screening_charge_uc_cm2",
                "fe_field_mv_cm",
                "il_field_mv_cm",
                "residual_mv_cm",
                "converged",
                "used_fallback",
                "j_total_uA_per_um2",
                "j_thermionic_uA_per_um2",
                "j_tunneling_uA_per_um2",
                "j_pf_uA_per_um2",
                "j_tat_uA_per_um2",
                "j_sclc_uA_per_um2",
                "top_barrier_ev",
                "bottom_barrier_ev",
                "tunneling_transmission",
                "trap_tunneling_transmission",
            ]
        )
        for segment_index, segment_slice in enumerate(segment_slices):
            for local_index, global_index in enumerate(range(segment_slice.start, segment_slice.stop)):
                state = states[global_index]
                breakdown = breakdowns[global_index]
                writer.writerow(
                    [
                        global_index,
                        segment_index,
                        state.voltage_v,
                        state.polarization_uc_cm2,
                        state.polarization_change_uc_cm2,
                        state.screening_charge_uc_cm2,
                        state.fe_field_mv_cm,
                        state.il_field_mv_cm,
                        state.residual_mv_cm,
                        int(state.converged),
                        int(state.used_fallback),
                        a_per_m2_to_uA_per_um2(breakdown.total_a_per_m2).item(),
                        a_per_m2_to_uA_per_um2(breakdown.thermionic_a_per_m2).item(),
                        a_per_m2_to_uA_per_um2(breakdown.tunneling_a_per_m2).item(),
                        a_per_m2_to_uA_per_um2(breakdown.poole_frenkel_a_per_m2).item(),
                        a_per_m2_to_uA_per_um2(breakdown.trap_assisted_tunneling_a_per_m2).item(),
                        a_per_m2_to_uA_per_um2(breakdown.sclc_a_per_m2).item(),
                        breakdown.top_barrier_ev,
                        breakdown.bottom_barrier_ev,
                        breakdown.tunneling_transmission,
                        breakdown.trap_tunneling_transmission,
                    ]
                )

    if args.report_summary:
        print(
            f"[dciv-summary] Pmin={np.min(polarization_uc_cm2):.3f} "
            f"Pmax={np.max(polarization_uc_cm2):.3f} "
            f"J_total_max={np.max(np.abs(total_uA_um2)):.6e} uA/um^2"
        )
        print(
            f"[dciv-breakdown] thermionic={np.max(np.abs(thermionic_uA_um2)):.6e} "
            f"tunneling={np.max(np.abs(tunneling_uA_um2)):.6e} "
            f"pf={np.max(np.abs(pf_uA_um2)):.6e} "
            f"tat={np.max(np.abs(tat_uA_um2)):.6e} "
            f"sclc={np.max(np.abs(sclc_uA_um2)):.6e}"
        )
        print(f"[dciv-output] figure={output_path} csv={csv_path}")


if __name__ == "__main__":
    main()
