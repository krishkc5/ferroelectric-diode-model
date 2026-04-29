import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import least_squares

from direct_transport_constraints import (
    GLOBAL_PARAMETER_ORDER,
    get_constraint_preset,
)
from fit_dciv_experiment import (
    diameter_nm_to_area_m2,
    display_insulator_name,
    infer_segment_slices,
    load_experimental_data,
)


SEGMENT_COLORS = {
    "0→-V": "#1f77b4",
    "-V→0": "#9467bd",
    "0→+V": "#ff7f0e",
    "+V→0": "#2ca02c",
}
MECHANISM_NAMES = ("thermionic", "pf_fe", "pf_il", "background")


def thermal_ev(temperature_k: float) -> float:
    return 8.617333262145e-5 * temperature_k


def direct_transport_base_terms(
    voltage_v: np.ndarray,
    temperature_k: float,
    top_phi_ev: float,
    top_gamma_ev_per_sqrt_v: float,
    log10_top_prefactor_a: float,
    bottom_phi_ev: float,
    bottom_gamma_ev_per_sqrt_v: float,
    log10_bottom_prefactor_a: float,
    trap_depth_ev: float,
    pf_fe_beta_ev_per_sqrt_v: float,
    log10_pf_fe_prefactor_a_per_v: float,
    pf_il_beta_ev_per_sqrt_v: float,
    log10_pf_il_prefactor_a_per_v: float,
    log10_background_conductance_a_per_v: float,
    log10_current_floor_a: float,
):
    voltage_v = np.asarray(voltage_v, dtype=float)
    vmag = np.abs(voltage_v)
    sqrt_v = np.sqrt(np.maximum(vmag, 0.0))
    k_t_ev = max(thermal_ev(temperature_k), 1e-12)

    top_prefactor_a = 10.0 ** log10_top_prefactor_a
    bottom_prefactor_a = 10.0 ** log10_bottom_prefactor_a
    pf_fe_prefactor_a_per_v = 10.0 ** log10_pf_fe_prefactor_a_per_v
    pf_il_prefactor_a_per_v = 10.0 ** log10_pf_il_prefactor_a_per_v
    background_conductance_a_per_v = 10.0 ** log10_background_conductance_a_per_v
    current_floor_a = 10.0 ** log10_current_floor_a

    top_effective_barrier_ev = np.maximum(
        top_phi_ev - top_gamma_ev_per_sqrt_v * sqrt_v,
        0.0,
    )
    bottom_effective_barrier_ev = np.maximum(
        bottom_phi_ev - bottom_gamma_ev_per_sqrt_v * sqrt_v,
        0.0,
    )

    top_thermionic_a = top_prefactor_a * np.exp(-top_effective_barrier_ev / k_t_ev)
    bottom_thermionic_a = bottom_prefactor_a * np.exp(-bottom_effective_barrier_ev / k_t_ev)
    thermionic_a = np.where(voltage_v >= 0.0, top_thermionic_a, bottom_thermionic_a)

    pf_fe_activation_ev = np.maximum(
        trap_depth_ev - pf_fe_beta_ev_per_sqrt_v * sqrt_v,
        0.0,
    )
    pf_il_activation_ev = np.maximum(
        trap_depth_ev - pf_il_beta_ev_per_sqrt_v * sqrt_v,
        0.0,
    )
    pf_fe_a = pf_fe_prefactor_a_per_v * vmag * np.exp(-pf_fe_activation_ev / k_t_ev)
    pf_il_a = pf_il_prefactor_a_per_v * vmag * np.exp(-pf_il_activation_ev / k_t_ev)
    background_a = background_conductance_a_per_v * vmag
    current_floor_array_a = np.full_like(voltage_v, current_floor_a, dtype=float)

    return {
        "current_floor_a": current_floor_array_a,
        "thermionic_a": thermionic_a,
        "pf_fe_a": pf_fe_a,
        "pf_il_a": pf_il_a,
        "background_a": background_a,
        "top_effective_barrier_ev": top_effective_barrier_ev,
        "bottom_effective_barrier_ev": bottom_effective_barrier_ev,
        "pf_fe_activation_ev": pf_fe_activation_ev,
        "pf_il_activation_ev": pf_il_activation_ev,
    }


def build_output_path(args: argparse.Namespace):
    if args.output is not None:
        return Path(args.output)
    insulator_label = display_insulator_name(args.insulator)
    il_label = f"{args.il_nm:g}nm"
    return Path(
        f"{int(args.fe_thickness_nm)}nm_figures/"
        f"Ti_{insulator_label}_AlScN_Al_{il_label}_direct_transport_segmented_fit_5curves.png"
    )


def mean_abs_log10_error(measured_a: np.ndarray, model_a: np.ndarray):
    errors = [
        abs(math.log10(abs(model)) - math.log10(abs(measured)))
        for measured, model in zip(measured_a, model_a)
        if abs(measured) > 0 and abs(model) > 0
    ]
    return float(sum(errors) / len(errors)) if errors else float("nan")


def build_disjoint_segments(voltage_v: np.ndarray):
    coarse_segments = infer_segment_slices(voltage_v)
    if len(coarse_segments) < 3:
        return [(f"seg{i + 1}", segment) for i, segment in enumerate(coarse_segments)]

    first = coarse_segments[0]
    middle = slice(coarse_segments[1].start + 1, coarse_segments[1].stop)
    last = slice(coarse_segments[2].start + 1, coarse_segments[2].stop)

    segments = [("0→-V", first)]
    middle_voltage = voltage_v[middle]
    zero_matches = np.where(np.isclose(middle_voltage, 0.0, atol=1e-12))[0]

    if zero_matches.size > 0:
        zero_abs = middle.start + int(zero_matches[0])
        if middle.start < zero_abs + 1:
            segments.append(("-V→0", slice(middle.start, zero_abs + 1)))
        if zero_abs + 1 < middle.stop:
            segments.append(("0→+V", slice(zero_abs + 1, middle.stop)))
    else:
        negative_or_zero = np.where(middle_voltage <= 0.0)[0]
        if negative_or_zero.size > 0 and negative_or_zero[-1] + 1 < middle_voltage.size:
            split_abs = middle.start + int(negative_or_zero[-1] + 1)
            segments.append(("-V→0", slice(middle.start, split_abs)))
            segments.append(("0→+V", slice(split_abs, middle.stop)))
        else:
            segments.append(("-V→+V", middle))

    segments.append(("+V→0", last))
    return [(label, seg) for label, seg in segments if seg.stop > seg.start]


def build_segment_index_map(voltage_v: np.ndarray):
    segment_map = np.full(len(voltage_v), -1, dtype=int)
    segment_definitions = build_disjoint_segments(voltage_v)
    for segment_index, (_, segment_slice) in enumerate(segment_definitions):
        segment_map[segment_slice] = segment_index
    if np.any(segment_map < 0):
        raise ValueError("Failed to assign every voltage point to a fit segment.")
    return segment_definitions, segment_map


def normalized_segment_weights(raw_logs: np.ndarray):
    raw_logs = np.asarray(raw_logs, dtype=float)
    return 10.0 ** (raw_logs - np.mean(raw_logs))


def unpack_parameters(params, num_segments: int):
    global_params = [float(value) for value in params[:13]]
    offset = 13
    weight_params = {}
    for mechanism_name in MECHANISM_NAMES:
        mechanism_raw = np.asarray(params[offset : offset + num_segments], dtype=float)
        weight_params[mechanism_name] = normalized_segment_weights(mechanism_raw)
        offset += num_segments
    return global_params, weight_params


def build_weighted_model(base_terms, segment_index_map, weight_params):
    total_a = np.asarray(base_terms["current_floor_a"], dtype=float).copy()
    weighted_terms = {}
    for mechanism_name in MECHANISM_NAMES:
        mechanism_values = np.asarray(base_terms[f"{mechanism_name}_a"], dtype=float)
        segment_weights = weight_params[mechanism_name][segment_index_map]
        weighted = mechanism_values * segment_weights
        weighted_terms[mechanism_name] = weighted
        total_a += weighted
    return total_a, weighted_terms


def main():
    parser = argparse.ArgumentParser(
        description="Fit direct transport equations to workbook DC I-V traces using shared transport constants and segment-specific mechanism weights."
    )
    parser.add_argument("--data", default="dciv_experimental.xlsx")
    parser.add_argument("--device-diameter-nm", type=float, default=200.0)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--insulator", default="hfox")
    parser.add_argument("--il-nm", type=float, default=4.0)
    parser.add_argument("--fe-thickness-nm", type=float, default=10.0)
    parser.add_argument(
        "--constraint-preset",
        default="physics-reasonable",
        choices=("broad", "physics-reasonable", "literature-ready"),
    )
    parser.add_argument("--max-nfev", type=int, default=220)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    sweeps = load_experimental_data(Path(args.data))
    if len(sweeps) < 2:
        raise ValueError("Expected multiple sweeps in workbook for direct shared-parameter fitting.")

    segment_definitions, reference_segment_index_map = build_segment_index_map(
        sweeps[0].voltage_v
    )
    segment_labels = [label for label, _ in segment_definitions]
    num_segments = len(segment_definitions)
    for sweep in sweeps[1:]:
        _, segment_index_map = build_segment_index_map(sweep.voltage_v)
        if not np.array_equal(segment_index_map, reference_segment_index_map):
            raise ValueError("All sweeps must share the same voltage segmentation.")

    all_currents = np.concatenate([np.abs(sweep.current_a) for sweep in sweeps])
    current_scale_a = max(np.percentile(all_currents, 10), 1e-12)
    area_m2 = diameter_nm_to_area_m2(args.device_diameter_nm)
    constraint_preset = get_constraint_preset(args.constraint_preset)

    def residuals(params):
        global_params, weight_params = unpack_parameters(params, num_segments)
        (
            top_phi_ev,
            top_gamma,
            log10_top_pref,
            bottom_phi_ev,
            bottom_gamma,
            log10_bottom_pref,
            trap_depth_ev,
            pf_fe_beta,
            log10_pf_fe_pref,
            pf_il_beta,
            log10_pf_il_pref,
            log10_bg_cond,
            log10_i_floor,
        ) = global_params

        residual_parts = []
        for sweep in sweeps:
            _, segment_index_map = build_segment_index_map(sweep.voltage_v)
            base_terms = direct_transport_base_terms(
                voltage_v=sweep.voltage_v,
                temperature_k=args.temperature_k,
                top_phi_ev=top_phi_ev,
                top_gamma_ev_per_sqrt_v=top_gamma,
                log10_top_prefactor_a=log10_top_pref,
                bottom_phi_ev=bottom_phi_ev,
                bottom_gamma_ev_per_sqrt_v=bottom_gamma,
                log10_bottom_prefactor_a=log10_bottom_pref,
                trap_depth_ev=trap_depth_ev,
                pf_fe_beta_ev_per_sqrt_v=pf_fe_beta,
                log10_pf_fe_prefactor_a_per_v=log10_pf_fe_pref,
                pf_il_beta_ev_per_sqrt_v=pf_il_beta,
                log10_pf_il_prefactor_a_per_v=log10_pf_il_pref,
                log10_background_conductance_a_per_v=log10_bg_cond,
                log10_current_floor_a=log10_i_floor,
            )
            model_total_a, _ = build_weighted_model(
                base_terms=base_terms,
                segment_index_map=segment_index_map,
                weight_params=weight_params,
            )
            lhs = np.abs(model_total_a)
            rhs = np.abs(sweep.current_a)
            residual_parts.append(
                np.arcsinh(lhs / current_scale_a) - np.arcsinh(rhs / current_scale_a)
            )
        return np.concatenate(residual_parts)

    initial_params = np.concatenate(
        [
            constraint_preset.initial_guess.copy(),
            np.zeros(num_segments * len(MECHANISM_NAMES), dtype=float),
        ]
    )
    lower_bounds = np.concatenate(
        [
            constraint_preset.lower_bounds.copy(),
            np.full(num_segments * len(MECHANISM_NAMES), -6.0, dtype=float),
        ]
    )
    upper_bounds = np.concatenate(
        [
            constraint_preset.upper_bounds.copy(),
            np.full(num_segments * len(MECHANISM_NAMES), 6.0, dtype=float),
        ]
    )

    fit = least_squares(
        residuals,
        x0=initial_params,
        bounds=(lower_bounds, upper_bounds),
        max_nfev=args.max_nfev,
    )

    global_params, weight_params = unpack_parameters(fit.x, num_segments)
    (
        top_phi_ev,
        top_gamma,
        log10_top_pref,
        bottom_phi_ev,
        bottom_gamma,
        log10_bottom_pref,
        trap_depth_ev,
        pf_fe_beta,
        log10_pf_fe_pref,
        pf_il_beta,
        log10_pf_il_pref,
        log10_bg_cond,
        log10_i_floor,
    ) = global_params

    fitted_results = []
    for sweep in sweeps:
        segment_definitions, segment_index_map = build_segment_index_map(sweep.voltage_v)
        base_terms = direct_transport_base_terms(
            voltage_v=sweep.voltage_v,
            temperature_k=args.temperature_k,
            top_phi_ev=top_phi_ev,
            top_gamma_ev_per_sqrt_v=top_gamma,
            log10_top_prefactor_a=log10_top_pref,
            bottom_phi_ev=bottom_phi_ev,
            bottom_gamma_ev_per_sqrt_v=bottom_gamma,
            log10_bottom_prefactor_a=log10_bottom_pref,
            trap_depth_ev=trap_depth_ev,
            pf_fe_beta_ev_per_sqrt_v=pf_fe_beta,
            log10_pf_fe_prefactor_a_per_v=log10_pf_fe_pref,
            pf_il_beta_ev_per_sqrt_v=pf_il_beta,
            log10_pf_il_prefactor_a_per_v=log10_pf_il_pref,
            log10_background_conductance_a_per_v=log10_bg_cond,
            log10_current_floor_a=log10_i_floor,
        )
        model_total_a, weighted_terms = build_weighted_model(
            base_terms=base_terms,
            segment_index_map=segment_index_map,
            weight_params=weight_params,
        )
        fitted_results.append(
            (
                sweep,
                segment_definitions,
                segment_index_map,
                base_terms,
                weighted_terms,
                model_total_a,
            )
        )

    output_path = build_output_path(args)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_path.with_suffix(".csv")

    fig, axes = plt.subplots(3, 2, figsize=(14.5, 13.5), dpi=args.dpi, sharex=True, sharey=True)
    fig.suptitle(
        f"Segment-Weighted Direct Transport Fit Across {len(fitted_results)} Workbook Sweeps: "
        f"Ti | {display_insulator_name(args.insulator)} ({args.il_nm:g} nm) | "
        f"AlScN ({args.fe_thickness_nm:g} nm) | Al",
        fontsize=15,
        y=0.98,
    )
    plot_axes = list(axes.flat)
    sweep_errors = []
    for axis, (sweep, segment_defs, _, _, _, model_total_a) in zip(plot_axes, fitted_results):
        error = mean_abs_log10_error(sweep.current_a, model_total_a)
        sweep_errors.append(error)
        for segment_label, segment_slice in segment_defs:
            color = SEGMENT_COLORS.get(segment_label, "#444444")
            axis.semilogy(
                sweep.voltage_v[segment_slice],
                np.maximum(np.abs(sweep.current_a[segment_slice]), 1e-20),
                "o",
                markersize=3.2,
                color=color,
                alpha=0.72,
            )
            axis.semilogy(
                sweep.voltage_v[segment_slice],
                np.maximum(np.abs(model_total_a[segment_slice]), 1e-20),
                "-",
                linewidth=2.0,
                color=color,
            )
        axis.set_title(f"{sweep.name}\nmean |Δlog10 I| = {error:.3f}", fontsize=10.5)
        axis.grid(alpha=0.25)

    for axis in plot_axes[:5]:
        axis.set_ylabel("|I| (A)")
    plot_axes[3].set_xlabel("Applied Voltage (V)")
    plot_axes[4].set_xlabel("Applied Voltage (V)")

    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="None", color="#111111", markersize=5, label="Measured"),
        Line2D([0], [0], linestyle="-", color="#111111", linewidth=2, label="Fit"),
    ]
    legend_handles.extend(
        [
            Line2D([0], [0], linestyle="-", color=SEGMENT_COLORS[label], linewidth=2.2, label=label)
            for label in segment_labels
        ]
    )
    plot_axes[0].legend(handles=legend_handles, loc="upper left", fontsize=8.2, ncol=2)

    summary_axis = plot_axes[5]
    summary_axis.axis("off")
    fit_lines = [
        "Shared direct transport constants",
        f"num sweeps = {len(fitted_results)}",
        f"diameter = {args.device_diameter_nm:g} nm",
        f"area = {area_m2:.3e} m^2",
        f"T = {args.temperature_k:g} K",
        f"constraint preset = {constraint_preset.name}",
        f"top phi = {top_phi_ev:.4f} eV",
        f"top gamma = {top_gamma:.4f} eV / sqrt(V)",
        f"log10 top prefactor (A) = {log10_top_pref:.4f}",
        f"bottom phi = {bottom_phi_ev:.4f} eV",
        f"bottom gamma = {bottom_gamma:.4f} eV / sqrt(V)",
        f"log10 bottom prefactor (A) = {log10_bottom_pref:.4f}",
        f"PF trap depth = {trap_depth_ev:.4f} eV",
        f"PF FE beta = {pf_fe_beta:.4f} eV / sqrt(V)",
        f"log10 PF FE pref (A/V) = {log10_pf_fe_pref:.4f}",
        f"PF IL beta = {pf_il_beta:.4f} eV / sqrt(V)",
        f"log10 PF IL pref (A/V) = {log10_pf_il_pref:.4f}",
        f"log10 bg cond (A/V) = {log10_bg_cond:.4f}",
        f"log10 current floor (A) = {log10_i_floor:.4f}",
        "",
        "Geometric-mean-normalized segment weights",
    ]
    for segment_index, segment_label in enumerate(segment_labels):
        fit_lines.append(
            f"{segment_label}: TI={weight_params['thermionic'][segment_index]:.3f}, "
            f"PF_FE={weight_params['pf_fe'][segment_index]:.3f}, "
            f"PF_IL={weight_params['pf_il'][segment_index]:.3f}, "
            f"BG={weight_params['background'][segment_index]:.3f}"
        )
    fit_lines.extend(
        [
            "",
            f"overall mean |Δlog10 I| = {sum(sweep_errors) / len(sweep_errors):.3f}",
            f"cost = {fit.cost:.3e}",
            f"nfev = {fit.nfev}",
        ]
    )
    summary_axis.text(
        0.02,
        0.98,
        "\n".join(fit_lines),
        ha="left",
        va="top",
        fontsize=8.6,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d9d9d9", alpha=0.95),
        transform=summary_axis.transAxes,
    )

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
    fig.savefig(output_path, dpi=args.dpi)
    plt.close(fig)

    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "sweep_name",
                "segment_label",
                "voltage_v",
                "experimental_current_a",
                "model_total_a",
                "model_thermionic_weighted_a",
                "model_pf_fe_weighted_a",
                "model_pf_il_weighted_a",
                "model_background_weighted_a",
                "base_thermionic_a",
                "base_pf_fe_a",
                "base_pf_il_a",
                "base_background_a",
                "current_floor_a",
                "weight_thermionic",
                "weight_pf_fe",
                "weight_pf_il",
                "weight_background",
                "top_effective_barrier_ev",
                "bottom_effective_barrier_ev",
                "pf_fe_activation_ev",
                "pf_il_activation_ev",
            ]
        )
        for sweep, segment_defs, segment_index_map, base_terms, weighted_terms, model_total_a in fitted_results:
            segment_label_by_point = [segment_labels[index] for index in segment_index_map]
            for point_index, values in enumerate(
                zip(
                    sweep.voltage_v,
                    sweep.current_a,
                    model_total_a,
                    weighted_terms["thermionic"],
                    weighted_terms["pf_fe"],
                    weighted_terms["pf_il"],
                    weighted_terms["background"],
                    base_terms["thermionic_a"],
                    base_terms["pf_fe_a"],
                    base_terms["pf_il_a"],
                    base_terms["background_a"],
                    base_terms["current_floor_a"],
                    base_terms["top_effective_barrier_ev"],
                    base_terms["bottom_effective_barrier_ev"],
                    base_terms["pf_fe_activation_ev"],
                    base_terms["pf_il_activation_ev"],
                )
            ):
                segment_index = int(segment_index_map[point_index])
                writer.writerow(
                    [
                        sweep.name,
                        segment_label_by_point[point_index],
                        values[0],
                        values[1],
                        values[2],
                        values[3],
                        values[4],
                        values[5],
                        values[6],
                        values[7],
                        values[8],
                        values[9],
                        values[10],
                        weight_params["thermionic"][segment_index],
                        weight_params["pf_fe"][segment_index],
                        weight_params["pf_il"][segment_index],
                        weight_params["background"][segment_index],
                        values[11],
                        values[12],
                        values[13],
                        values[14],
                        values[15],
                    ]
                )

    summary_parts = [
        f"constraint_preset={constraint_preset.name}",
        f"top_phi_ev={top_phi_ev:.6f}",
        f"top_gamma={top_gamma:.6f}",
        f"log10_top_prefactor_a={log10_top_pref:.6f}",
        f"bottom_phi_ev={bottom_phi_ev:.6f}",
        f"bottom_gamma={bottom_gamma:.6f}",
        f"log10_bottom_prefactor_a={log10_bottom_pref:.6f}",
        f"pf_trap_depth_ev={trap_depth_ev:.6f}",
        f"pf_fe_beta={pf_fe_beta:.6f}",
        f"log10_pf_fe_prefactor_a_per_v={log10_pf_fe_pref:.6f}",
        f"pf_il_beta={pf_il_beta:.6f}",
        f"log10_pf_il_prefactor_a_per_v={log10_pf_il_pref:.6f}",
        f"log10_bg_cond={log10_bg_cond:.6f}",
        f"log10_i_floor={log10_i_floor:.6f}",
    ]
    for mechanism_name in MECHANISM_NAMES:
        for segment_index, segment_label in enumerate(segment_labels):
            safe_label = (
                segment_label.replace("→", "_to_")
                .replace("+", "plus")
                .replace("-", "minus")
            )
            summary_parts.append(
                f"{mechanism_name}_weight_{safe_label}={weight_params[mechanism_name][segment_index]:.6f}"
            )
    summary_parts.append(f"cost={fit.cost:.6e}")
    print("[segmented-direct-fit-summary] " + " ".join(summary_parts))
    print("[segmented-direct-fit-parameter-order] " + " ".join(GLOBAL_PARAMETER_ORDER))
    print(f"[segmented-direct-fit-output] figure={output_path} csv={csv_path}")


if __name__ == "__main__":
    main()
