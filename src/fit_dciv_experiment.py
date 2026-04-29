import argparse
import csv
import math
import zipfile
from dataclasses import dataclass
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

from hysteresis_core import HysteresisConfig, HysteresisSimulator
from sweeps import build_waypoint_sweep
from transport import TransportEvaluator, TransportModelParameters


@dataclass(frozen=True)
class ExperimentalSweep:
    name: str
    voltage_v: np.ndarray
    current_a: np.ndarray


def load_experimental_csv(path: Path):
    voltage = []
    current = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        fieldnames = {name.strip().lower(): name for name in reader.fieldnames or []}
        voltage_key = fieldnames.get("av", "AV")
        current_key = fieldnames.get("bi", "BI")
        for row in reader:
            voltage.append(float(row[voltage_key]))
            current.append(float(row[current_key]))
    return [ExperimentalSweep(name=path.stem, voltage_v=np.asarray(voltage, dtype=float), current_a=np.asarray(current, dtype=float))]


def _xlsx_cell_text(cell, ns):
    v = cell.find("a:v", ns)
    isel = cell.find("a:is", ns)
    if v is not None and v.text is not None:
        return v.text
    if isel is not None:
        return "".join(node.text or "" for node in isel.findall(".//a:t", ns))
    return ""


def load_experimental_xlsx(path: Path):
    ns = {
        "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
        "pr": "http://schemas.openxmlformats.org/package/2006/relationships",
    }
    sweeps = []
    with zipfile.ZipFile(path) as zf:
        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels.findall("pr:Relationship", ns)
        }
        for sheet in workbook.find("a:sheets", ns):
            sheet_name = sheet.attrib["name"]
            rel_id = sheet.attrib[f"{{{ns['r']}}}id"]
            target = rel_map[rel_id]
            sheet_path = target if target.startswith("xl/") else f"xl/{target}"
            root = ET.fromstring(zf.read(sheet_path))
            rows = root.findall(".//a:sheetData/a:row", ns)
            if not rows:
                continue
            header = [_xlsx_cell_text(cell, ns).strip().upper() for cell in rows[0].findall("a:c", ns)]
            if header[:2] != ["BI", "AV"]:
                continue
            voltage = []
            current = []
            for row in rows[1:]:
                cells = row.findall("a:c", ns)
                if len(cells) < 2:
                    continue
                bi = _xlsx_cell_text(cells[0], ns).strip()
                av = _xlsx_cell_text(cells[1], ns).strip()
                if not bi or not av:
                    continue
                current.append(float(bi))
                voltage.append(float(av))
            if voltage:
                sweeps.append(
                    ExperimentalSweep(
                        name=sheet_name,
                        voltage_v=np.asarray(voltage, dtype=float),
                        current_a=np.asarray(current, dtype=float),
                    )
                )
    return sweeps


def load_experimental_data(path: Path):
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return load_experimental_csv(path)
    if suffix == ".xlsx":
        return load_experimental_xlsx(path)
    raise ValueError(f"Unsupported experimental file type: {path.suffix}")


def infer_segment_slices(voltage_trace):
    voltage_trace = np.asarray(voltage_trace, dtype=float)
    if voltage_trace.size <= 1:
        return [slice(0, voltage_trace.size)]

    signs = []
    for index in range(1, len(voltage_trace)):
        dv = voltage_trace[index] - voltage_trace[index - 1]
        if abs(dv) < 1e-15:
            signs.append(0)
        else:
            signs.append(1 if dv > 0 else -1)

    slices = []
    start = 0
    current_sign = signs[0]
    for index in range(1, len(signs)):
        if signs[index] != 0 and current_sign != 0 and signs[index] != current_sign:
            slices.append(slice(start, index + 1))
            start = index
        if signs[index] != 0:
            current_sign = signs[index]
    slices.append(slice(start, len(voltage_trace)))
    return slices


def diameter_nm_to_area_m2(diameter_nm):
    radius_m = 0.5 * diameter_nm * 1e-9
    return math.pi * radius_m**2


def build_output_path(args: argparse.Namespace):
    if args.output is not None:
        return Path(args.output)
    insulator_label = display_insulator_name(args.insulator)
    il_label = f"{args.il_nm:g}nm" if args.il_nm > 0 else "0nm"
    pf_suffix_map = {"fe": "", "il": "_pfIL", "both": "_pfBoth"}
    pf_suffix = pf_suffix_map[args.pf_region]
    dataset_suffix = "_5curves" if Path(args.data).suffix.lower() == ".xlsx" else ""
    precond_suffix = (
        f"_precond{args.precondition_bipolar_cycles}"
        if args.precondition_bipolar_cycles > 0
        else ""
    )
    return Path(
        f"{int(args.fe_thickness_nm)}nm_figures/"
        f"Ti_{insulator_label}_AlScN_Al_{il_label}_experiment_fit{pf_suffix}{dataset_suffix}{precond_suffix}.png"
    )


def display_insulator_name(name: str):
    token = name.lower()
    if token == "hfox":
        return "HfOx"
    if token == "al2o3":
        return "Al2O3"
    return name


def model_currents_from_states(
    states,
    diode,
    potential,
    temperature_k,
    pf_region,
    top_barrier_offset_ev,
    bottom_barrier_offset_ev,
    top_barrier_floor_ev,
    top_schottky_scale,
    bottom_schottky_scale,
    top_thermionic_log10_scale,
    bottom_thermionic_log10_scale,
    pf_fe_log10_scale,
    pf_il_log10_scale,
    pf_trap_depth_ev,
    background_log10_conductance,
):
    evaluator = TransportEvaluator(
        diode=diode,
        potential=potential,
        parameters=TransportModelParameters(
            temperature_k=temperature_k,
            thermionic_top_barrier_offset_ev=top_barrier_offset_ev,
            thermionic_bottom_barrier_offset_ev=bottom_barrier_offset_ev,
            thermionic_top_barrier_floor_ev=top_barrier_floor_ev,
            thermionic_top_schottky_scale=top_schottky_scale,
            thermionic_bottom_schottky_scale=bottom_schottky_scale,
            thermionic_top_prefactor_scale=10 ** top_thermionic_log10_scale,
            thermionic_bottom_prefactor_scale=10 ** bottom_thermionic_log10_scale,
            poole_frenkel_fe_prefactor_scale=10 ** pf_fe_log10_scale,
            poole_frenkel_il_prefactor_scale=10 ** pf_il_log10_scale,
            poole_frenkel_trap_depth_ev=pf_trap_depth_ev,
            background_conductance_a_per_m2_v=10 ** background_log10_conductance,
            poole_frenkel_field_region=pf_region,
        ),
        enable_thermionic=True,
        enable_tunneling=False,
        enable_poole_frenkel=True,
        enable_trap_assisted_tunneling=False,
        enable_sclc=False,
    )
    return [evaluator.evaluate(state) for state in states]


def infer_voltage_step(voltage_trace_v):
    diffs = np.abs(np.diff(np.asarray(voltage_trace_v, dtype=float)))
    diffs = diffs[diffs > 1e-12]
    if diffs.size == 0:
        return 0.25
    return float(np.median(diffs))


def build_preconditioning_trace(experimental_voltage_v, cycles: int):
    if cycles <= 0:
        return np.asarray([], dtype=float)
    start_v = float(experimental_voltage_v[0])
    vmax = float(np.max(np.abs(experimental_voltage_v)))
    step = infer_voltage_step(experimental_voltage_v)
    single_cycle = build_waypoint_sweep([start_v, -vmax, vmax, -vmax, start_v], step)
    traces = [single_cycle]
    for _ in range(cycles - 1):
        traces.append(single_cycle[1:])
    return np.concatenate(traces)


def simulate_states_for_sweeps(sweeps, config, il_nm, precondition_cycles):
    states_by_sweep = []
    reference_simulator = None
    for sweep in sweeps:
        simulator = HysteresisSimulator(il_nm=il_nm, config=config)
        if precondition_cycles > 0:
            preconditioning_trace = build_preconditioning_trace(
                experimental_voltage_v=sweep.voltage_v,
                cycles=precondition_cycles,
            )
            simulator.simulate_trace(preconditioning_trace)
        states_by_sweep.append(simulator.simulate_trace(sweep.voltage_v))
        if reference_simulator is None:
            reference_simulator = simulator
    return reference_simulator, states_by_sweep


def main():
    parser = argparse.ArgumentParser(
        description="Fit the first-pass transport model to an experimental DC I-V sweep."
    )
    parser.add_argument("--data", default="dciv_experimental.csv")
    parser.add_argument("--device-diameter-nm", type=float, default=200.0)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--top-electrode", default="ti")
    parser.add_argument("--bottom-electrode", default="al")
    parser.add_argument("--insulator", default="hfox")
    parser.add_argument("--il-nm", type=float, default=4.0)
    parser.add_argument("--fe-thickness-nm", type=float, default=10.0)
    parser.add_argument("--num-domains", type=int, default=6000)
    parser.add_argument("--ca-std", type=float, default=0.04)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-nfev", type=int, default=80)
    parser.add_argument("--dpi", type=int, default=220)
    parser.add_argument("--output", default=None)
    parser.add_argument("--pf-region", choices=["fe", "il", "both"], default="fe")
    parser.add_argument(
        "--precondition-bipolar-cycles",
        type=int,
        default=0,
        help="Number of 0 -> -V -> +V -> -V -> 0 training cycles to run before the measured sweep.",
    )
    parser.add_argument(
        "--fit-signed-current",
        action="store_true",
        help="Fit signed current directly instead of current magnitude.",
    )
    args = parser.parse_args()

    experimental_sweeps = load_experimental_data(Path(args.data))
    if not experimental_sweeps:
        raise ValueError(f"No BI/AV experimental sweeps found in {args.data}")
    area_m2 = diameter_nm_to_area_m2(args.device_diameter_nm)

    config = HysteresisConfig(
        num_domains=args.num_domains,
        c_a_std=args.ca_std,
        seed=args.seed,
        top_electrode_name=args.top_electrode,
        bottom_electrode_name=args.bottom_electrode,
        insulator_name=args.insulator,
        fe_thickness_nm=args.fe_thickness_nm,
    )
    simulator, states_by_sweep = simulate_states_for_sweeps(
        sweeps=experimental_sweeps,
        config=config,
        il_nm=args.il_nm,
        precondition_cycles=args.precondition_bipolar_cycles,
    )

    all_currents = np.concatenate([np.abs(sweep.current_a) for sweep in experimental_sweeps])
    current_scale_a = max(np.percentile(all_currents, 10), 1e-12)

    def residuals(params):
        top_barrier_offset_ev = float(params[0])
        bottom_barrier_offset_ev = float(params[1])
        top_barrier_floor_ev = float(params[2])
        top_schottky_scale = float(params[3])
        bottom_schottky_scale = float(params[4])
        top_thermionic_log10_scale = float(params[5])
        bottom_thermionic_log10_scale = float(params[6])
        pf_fe_log10_scale = float(params[7])
        pf_il_log10_scale = float(params[8])
        pf_trap_depth_ev = float(params[9])
        background_log10_conductance = float(params[10])
        all_residuals = []
        for sweep, states in zip(experimental_sweeps, states_by_sweep):
            breakdowns = model_currents_from_states(
                states=states,
                diode=simulator.diode,
                potential=simulator.potential,
                temperature_k=args.temperature_k,
                pf_region=args.pf_region,
                top_barrier_offset_ev=top_barrier_offset_ev,
                bottom_barrier_offset_ev=bottom_barrier_offset_ev,
                top_barrier_floor_ev=top_barrier_floor_ev,
                top_schottky_scale=top_schottky_scale,
                bottom_schottky_scale=bottom_schottky_scale,
                top_thermionic_log10_scale=top_thermionic_log10_scale,
                bottom_thermionic_log10_scale=bottom_thermionic_log10_scale,
                pf_fe_log10_scale=pf_fe_log10_scale,
                pf_il_log10_scale=pf_il_log10_scale,
                pf_trap_depth_ev=pf_trap_depth_ev,
                background_log10_conductance=background_log10_conductance,
            )
            model_current_a = area_m2 * np.asarray([item.total_a_per_m2 for item in breakdowns], dtype=float)
            if args.fit_signed_current:
                lhs = model_current_a
                rhs = sweep.current_a
            else:
                lhs = np.abs(model_current_a)
                rhs = np.abs(sweep.current_a)
            all_residuals.append(np.arcsinh(lhs / current_scale_a) - np.arcsinh(rhs / current_scale_a))
        return np.concatenate(all_residuals)

    fit = least_squares(
        residuals,
        x0=np.array([1.5, 2.5, 0.0, 1.0, 1.0, 0.5, 0.5, 9.0, 9.0, 0.8, 1.0], dtype=float),
        bounds=(
            [0.0, 0.0, 0.0, 0.1, 0.1, -12.0, -12.0, -8.0, -8.0, 0.05, -12.0],
            [4.5, 4.5, 2.0, 10.0, 10.0, 24.0, 24.0, 20.0, 20.0, 2.0, 12.0],
        ),
        max_nfev=args.max_nfev,
    )

    fitted_top_barrier_offset_ev = float(fit.x[0])
    fitted_bottom_barrier_offset_ev = float(fit.x[1])
    fitted_top_barrier_floor_ev = float(fit.x[2])
    fitted_top_schottky_scale = float(fit.x[3])
    fitted_bottom_schottky_scale = float(fit.x[4])
    fitted_top_thermionic_log10_scale = float(fit.x[5])
    fitted_bottom_thermionic_log10_scale = float(fit.x[6])
    fitted_pf_fe_log10_scale = float(fit.x[7])
    fitted_pf_il_log10_scale = float(fit.x[8])
    fitted_pf_trap_depth_ev = float(fit.x[9])
    fitted_background_log10_conductance = float(fit.x[10])
    fitted_results = []
    for sweep, states in zip(experimental_sweeps, states_by_sweep):
        breakdowns = model_currents_from_states(
            states=states,
            diode=simulator.diode,
            potential=simulator.potential,
            temperature_k=args.temperature_k,
            pf_region=args.pf_region,
            top_barrier_offset_ev=fitted_top_barrier_offset_ev,
            bottom_barrier_offset_ev=fitted_bottom_barrier_offset_ev,
            top_barrier_floor_ev=fitted_top_barrier_floor_ev,
            top_schottky_scale=fitted_top_schottky_scale,
            bottom_schottky_scale=fitted_bottom_schottky_scale,
            top_thermionic_log10_scale=fitted_top_thermionic_log10_scale,
            bottom_thermionic_log10_scale=fitted_bottom_thermionic_log10_scale,
            pf_fe_log10_scale=fitted_pf_fe_log10_scale,
            pf_il_log10_scale=fitted_pf_il_log10_scale,
            pf_trap_depth_ev=fitted_pf_trap_depth_ev,
            background_log10_conductance=fitted_background_log10_conductance,
        )
        fitted_results.append(
            {
                "sweep": sweep,
                "states": states,
                "breakdowns": breakdowns,
                "model_total_a": area_m2 * np.asarray([item.total_a_per_m2 for item in breakdowns], dtype=float),
                "model_thermionic_a": area_m2
                * np.asarray([item.thermionic_a_per_m2 for item in breakdowns], dtype=float),
                "model_pf_a": area_m2
                * np.asarray([item.poole_frenkel_a_per_m2 for item in breakdowns], dtype=float),
                "model_background_a": area_m2
                * np.asarray([item.background_leakage_a_per_m2 for item in breakdowns], dtype=float),
                "polarization_uc_cm2": np.asarray([state.polarization_uc_cm2 for state in states], dtype=float),
            }
        )

    output_path = build_output_path(args)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_path.with_suffix(".csv")
    def mean_abs_log10_error(exp_current, model_current):
        errors = [
            abs(math.log10(abs(model)) - math.log10(abs(exp)))
            for exp, model in zip(exp_current, model_current)
            if abs(exp) > 0 and abs(model) > 0
        ]
        return float(sum(errors) / len(errors)) if errors else float("nan")

    overall_errors = [
        mean_abs_log10_error(result["sweep"].current_a, result["model_total_a"])
        for result in fitted_results
    ]
    overall_mean_error = float(sum(overall_errors) / len(overall_errors))

    fig, axes = plt.subplots(3, 2, figsize=(14.0, 13.0), dpi=args.dpi, sharex=True, sharey=True)
    fig.suptitle(
        f"Shared-Parameter Experimental DC I-V Fit Across {len(fitted_results)} Sweeps: "
        f"Ti | {display_insulator_name(args.insulator)} ({args.il_nm:g} nm) | "
        f"AlScN ({args.fe_thickness_nm:g} nm) | Al",
        fontsize=15,
        y=0.98,
    )
    plot_axes = list(axes.flat)
    for axis, result, sweep_error in zip(plot_axes, fitted_results, overall_errors):
        sweep = result["sweep"]
        axis.semilogy(
            sweep.voltage_v,
            np.maximum(np.abs(sweep.current_a), 1e-20),
            "o",
            markersize=3.5,
            color="#d62728",
            alpha=0.78,
            label="Measured",
        )
        axis.semilogy(
            sweep.voltage_v,
            np.maximum(np.abs(result["model_total_a"]), 1e-20),
            "-",
            linewidth=2.1,
            color="#111111",
            label="Model",
        )
        axis.set_title(f"{sweep.name}\nmean |Δlog10 I| = {sweep_error:.3f}", fontsize=10.5)
        axis.grid(alpha=0.25)
    for axis in plot_axes[:5]:
        axis.set_ylabel("|I| (A)")
    plot_axes[4].set_xlabel("Applied Voltage (V)")
    plot_axes[3].set_xlabel("Applied Voltage (V)")
    plot_axes[0].legend(loc="upper left", fontsize=8.5)

    summary_axis = plot_axes[5]
    summary_axis.axis("off")
    fit_text = "\n".join(
        [
            "Shared fit across workbook sweeps",
            f"num sweeps = {len(fitted_results)}",
            f"diameter = {args.device_diameter_nm:g} nm",
            f"area = {area_m2:.3e} m^2",
            f"top barrier offset = {fitted_top_barrier_offset_ev:.3f} eV",
            f"bottom barrier offset = {fitted_bottom_barrier_offset_ev:.3f} eV",
            f"top barrier floor = {fitted_top_barrier_floor_ev:.3f} eV",
            f"top Schottky scale = {fitted_top_schottky_scale:.3f}",
            f"bottom Schottky scale = {fitted_bottom_schottky_scale:.3f}",
            f"log10 top TI scale = {fitted_top_thermionic_log10_scale:.3f}",
            f"log10 bottom TI scale = {fitted_bottom_thermionic_log10_scale:.3f}",
            f"log10 PF FE scale = {fitted_pf_fe_log10_scale:.3f}",
            f"log10 PF IL scale = {fitted_pf_il_log10_scale:.3f}",
            f"PF trap depth = {fitted_pf_trap_depth_ev:.3f} eV",
            f"log10 bg cond = {fitted_background_log10_conductance:.3f}",
            f"PF region = {args.pf_region}",
            f"mode = {'signed I' if args.fit_signed_current else '|I|'}",
            f"precondition cycles = {args.precondition_bipolar_cycles}",
            f"overall mean |Δlog10 I| = {overall_mean_error:.3f}",
            f"cost = {fit.cost:.3e}",
            f"nfev = {fit.nfev}",
        ]
    )
    summary_axis.text(
        0.02,
        0.98,
        fit_text,
        ha="left",
        va="top",
        fontsize=9.5,
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
                "voltage_v",
                "experimental_current_a",
                "model_total_a",
                "model_pf_a",
                "model_thermionic_a",
                "model_background_a",
                "polarization_uc_cm2",
            ]
        )
        for result in fitted_results:
            sweep = result["sweep"]
            for voltage, current_exp, current_model, current_pf, current_ti, current_bg, polarization in zip(
                sweep.voltage_v,
                sweep.current_a,
                result["model_total_a"],
                result["model_pf_a"],
                result["model_thermionic_a"],
                result["model_background_a"],
                result["polarization_uc_cm2"],
            ):
                writer.writerow(
                    [sweep.name, voltage, current_exp, current_model, current_pf, current_ti, current_bg, polarization]
                )

    print(
        f"[fit-summary] top_barrier_offset_ev={fitted_top_barrier_offset_ev:.6f} "
        f"bottom_barrier_offset_ev={fitted_bottom_barrier_offset_ev:.6f} "
        f"top_barrier_floor_ev={fitted_top_barrier_floor_ev:.6f} "
        f"top_schottky_scale={fitted_top_schottky_scale:.6f} "
        f"bottom_schottky_scale={fitted_bottom_schottky_scale:.6f} "
        f"log10_top_ti_scale={fitted_top_thermionic_log10_scale:.6f} "
        f"log10_bottom_ti_scale={fitted_bottom_thermionic_log10_scale:.6f} "
        f"log10_pf_fe_scale={fitted_pf_fe_log10_scale:.6f} "
        f"log10_pf_il_scale={fitted_pf_il_log10_scale:.6f} "
        f"pf_trap_depth_ev={fitted_pf_trap_depth_ev:.6f} "
        f"log10_bg_cond={fitted_background_log10_conductance:.6f} cost={fit.cost:.6e}"
    )
    print(f"[fit-output] figure={output_path} csv={csv_path}")


if __name__ == "__main__":
    main()
