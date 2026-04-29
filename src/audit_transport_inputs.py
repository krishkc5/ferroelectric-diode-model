import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from hysteresis_core import HysteresisConfig, HysteresisSimulator
from transport import TransportEvaluator, TransportModelParameters


@dataclass(frozen=True)
class ExperimentalSweep:
    voltage_v: np.ndarray
    current_a: np.ndarray


def diameter_nm_to_area_m2(diameter_nm):
    radius_m = 0.5 * diameter_nm * 1e-9
    return math.pi * radius_m**2


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
    return ExperimentalSweep(voltage_v=np.asarray(voltage, dtype=float), current_a=np.asarray(current, dtype=float))


def display_insulator_name(name: str):
    token = name.lower()
    if token == "hfox":
        return "HfOx"
    if token == "al2o3":
        return "Al2O3"
    return name


def build_default_output_path(args):
    pf_suffix_map = {"fe": "", "il": "_pfIL", "both": "_pfBoth"}
    pf_suffix = pf_suffix_map[args.pf_region]
    return Path(
        f"{int(args.fe_thickness_nm)}nm_figures/"
        f"Ti_{display_insulator_name(args.insulator)}_AlScN_Al_{args.il_nm:g}nm_transport_audit{pf_suffix}.csv"
    )


def build_audit_rows(experimental: ExperimentalSweep, simulator, evaluator, area_m2):
    states = simulator.simulate_trace(experimental.voltage_v)
    rows = []
    for state, experimental_current_a in zip(states, experimental.current_a):
        audit = evaluator.audit_state(state)
        audit["experimental_current_a"] = float(experimental_current_a)
        audit["experimental_current_density_a_per_m2"] = float(experimental_current_a / area_m2)
        audit["experimental_current_density_uA_per_um2"] = float(experimental_current_a * 1e6 / (area_m2 * 1e12))
        audit["model_total_a"] = float(audit["total_a_per_m2"] * area_m2)
        audit["model_thermionic_a"] = float(audit["thermionic_a_per_m2"] * area_m2)
        audit["model_tunneling_a"] = float(audit["tunneling_a_per_m2"] * area_m2)
        audit["model_pf_a"] = float(audit["poole_frenkel_a_per_m2"] * area_m2)
        audit["model_tat_a"] = float(audit["trap_assisted_tunneling_a_per_m2"] * area_m2)
        audit["model_sclc_a"] = float(audit["sclc_a_per_m2"] * area_m2)
        audit["model_background_a"] = float(audit["background_a_per_m2"] * area_m2)
        rows.append(audit)
    return rows


def summarize_rows(rows):
    representative_targets = [0.0, -9.5, -5.0, 5.0, 9.5]
    summary = []
    for target in representative_targets:
        best = min(rows, key=lambda row: abs(row["voltage_v"] - target))
        summary.append(best)
    return summary


def main():
    parser = argparse.ArgumentParser(description="Dump all realized transport inputs for audit/debugging.")
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
    parser.add_argument("--top-barrier-offset-ev", type=float, default=1.548389)
    parser.add_argument("--bottom-barrier-offset-ev", type=float, default=0.0)
    parser.add_argument("--top-barrier-floor-ev", type=float, default=0.0)
    parser.add_argument("--top-schottky-scale", type=float, default=1.0)
    parser.add_argument("--bottom-schottky-scale", type=float, default=1.0)
    parser.add_argument("--log10-top-ti-scale", type=float, default=0.867826)
    parser.add_argument("--log10-bottom-ti-scale", type=float, default=0.867826)
    parser.add_argument("--log10-pf-fe-scale", type=float, default=9.415499)
    parser.add_argument("--log10-pf-il-scale", type=float, default=9.415499)
    parser.add_argument("--pf-trap-depth-ev", type=float, default=0.8)
    parser.add_argument("--log10-bg-cond", type=float, default=-12.0)
    parser.add_argument("--pf-region", choices=["fe", "il", "both"], default="fe")
    parser.add_argument("--output", default=None)
    args = parser.parse_args()

    experimental = load_experimental_csv(Path(args.data))
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
    simulator = HysteresisSimulator(il_nm=args.il_nm, config=config)
    evaluator = TransportEvaluator(
        diode=simulator.diode,
        potential=simulator.potential,
        parameters=TransportModelParameters(
            temperature_k=args.temperature_k,
            thermionic_top_barrier_offset_ev=args.top_barrier_offset_ev,
            thermionic_bottom_barrier_offset_ev=args.bottom_barrier_offset_ev,
            thermionic_top_barrier_floor_ev=args.top_barrier_floor_ev,
            thermionic_top_schottky_scale=args.top_schottky_scale,
            thermionic_bottom_schottky_scale=args.bottom_schottky_scale,
            thermionic_top_prefactor_scale=10 ** args.log10_top_ti_scale,
            thermionic_bottom_prefactor_scale=10 ** args.log10_bottom_ti_scale,
            poole_frenkel_fe_prefactor_scale=10 ** args.log10_pf_fe_scale,
            poole_frenkel_il_prefactor_scale=10 ** args.log10_pf_il_scale,
            poole_frenkel_trap_depth_ev=args.pf_trap_depth_ev,
            background_conductance_a_per_m2_v=10 ** args.log10_bg_cond,
            poole_frenkel_field_region=args.pf_region,
        ),
        enable_thermionic=True,
        enable_tunneling=False,
        enable_poole_frenkel=True,
        enable_trap_assisted_tunneling=False,
        enable_sclc=False,
    )

    rows = build_audit_rows(experimental, simulator, evaluator, area_m2)
    output_path = Path(args.output) if args.output is not None else build_default_output_path(args)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(f"[audit-output] csv={output_path}")
    print("[audit-representative-points]")
    for row in summarize_rows(rows):
        print(
            f"V={row['voltage_v']:+.2f} V | P={row['polarization_uc_cm2']:+.3f} uC/cm^2 | "
            f"Efe={row['fe_field_mv_cm']:+.3f} MV/cm | "
            f"Iexp={row['experimental_current_a']:+.3e} A | Imod={row['model_total_a']:+.3e} A | "
            f"TI={row['model_thermionic_a']:+.3e} A | PF={row['model_pf_a']:+.3e} A | "
            f"phi_top_eff={row['top_barrier_effective_ev']:.3f} eV | "
            f"phi_bot_eff={row['bottom_barrier_effective_ev']:.3f} eV | "
            f"PF_act={row['pf_activation_ev']:.3f} eV"
        )


if __name__ == "__main__":
    main()
