import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as path_effects
from matplotlib.lines import Line2D

from plot_il_hysteresis import build_voltage_sweep, run_hysteresis_for_il, resolve_metal_electrode


TOP_ELECTRODES = ["ti", "al", "pd", "test", "test2"]
BOTTOM_ELECTRODES = ["ti", "al", "pd", "test", "test2"]
TOP_COLORS = {
    "ti": "#1f77b4",
    "al": "#ff7f0e",
    "pd": "#2ca02c",
    "test": "#d62728",
    "test2": "#9467bd",
}
BOTTOM_STYLES = {
    "ti": ("-", "o"),
    "al": ("--", "s"),
    "pd": ("-.", "^"),
    "test": (":", "D"),
    "test2": ((0, (6, 2, 1.5, 2)), "x"),
}


def pretty(name: str) -> str:
    return resolve_metal_electrode(name).name


def build_args(top: str, bottom: str, il_nm: float) -> argparse.Namespace:
    return argparse.Namespace(
        il_nm=[il_nm],
        vmax=20.0,
        vstep=0.25,
        num_domains=15000,
        ca_std=0.04,
        seed=0,
        top_electrode=top,
        bottom_electrode=bottom,
        output=None,
        dpi=200,
        title="",
        show_assumptions=False,
        report_loop_summary=False,
    )


def compute_curves(il_nm: float):
    single_cycle = build_voltage_sweep(20.0, 0.25)
    curves = {}
    for top in TOP_ELECTRODES:
        for bottom in BOTTOM_ELECTRODES:
            args = build_args(top=top, bottom=bottom, il_nm=il_nm)
            curves[(top, bottom)] = np.asarray(
                run_hysteresis_for_il(il_nm, single_cycle, args), dtype=float
            )
    return single_cycle, curves


def add_encoding_legends(ax):
    top_handles = [
        Line2D([0], [0], color=TOP_COLORS[name], lw=2.8, label=pretty(name))
        for name in TOP_ELECTRODES
    ]
    bottom_handles = [
        Line2D(
            [0],
            [0],
            color="black",
            lw=2.5,
            linestyle=BOTTOM_STYLES[name][0],
            marker=BOTTOM_STYLES[name][1],
            markersize=5,
            label=pretty(name),
        )
        for name in BOTTOM_ELECTRODES
    ]
    legend_top = ax.legend(
        handles=top_handles,
        title="Top electrode (color)",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=True,
        facecolor="white",
        edgecolor="#d9d9d9",
        framealpha=0.95,
    )
    ax.add_artist(legend_top)
    ax.legend(
        handles=bottom_handles,
        title="Bottom electrode (style + marker)",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.48),
        frameon=True,
        facecolor="white",
        edgecolor="#d9d9d9",
        framealpha=0.95,
    )


def plot_overlay(single_cycle: np.ndarray, curves: dict[tuple[str, str], np.ndarray], il_nm: float, output: Path):
    reference_key = ("ti", "ti")
    fig, ax = plt.subplots(figsize=(13.8, 7.4), dpi=200)
    for top in TOP_ELECTRODES:
        for bottom in BOTTOM_ELECTRODES:
            key = (top, bottom)
            linestyle, marker = BOTTOM_STYLES[bottom]
            lw = 3.3 if key == reference_key else 1.8
            alpha = 1.0 if key == reference_key else 0.82
            zorder = 30 if key == reference_key else 10
            ax.plot(
                single_cycle,
                curves[key],
                color=TOP_COLORS[top],
                linestyle=linestyle,
                marker=marker,
                markevery=34,
                markersize=3.8,
                linewidth=lw,
                alpha=alpha,
                zorder=zorder,
            )

    ax.set_title(f"Electrode Comparison at {il_nm:g} nm IL")
    ax.set_xlabel("Applied Voltage (V)")
    ax.set_ylabel("FE Polarization (uC/cm^2)")
    ax.set_xlim(-21.6, 21.6)
    ax.set_ylim(-120, 120)
    ax.set_xticks(np.arange(-20, 21, 5))
    ax.set_yticks(np.arange(-100, 101, 25))
    ax.text(
        1.02,
        0.05,
        "Reference curve: Ti/Ti\n"
        "Test: lower WF + longer screening\n"
        "Test2: same WF as Ti, longer screening",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.5,
        family="monospace",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#d9d9d9", alpha=0.95),
    )
    add_encoding_legends(ax)
    fig.tight_layout(rect=(0.0, 0.0, 0.82, 1.0))
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def annotate_heatmap(ax, data: np.ndarray):
    max_value = float(np.max(data)) if np.max(data) > 0 else 1.0
    for row in range(data.shape[0]):
        for col in range(data.shape[1]):
            value = data[row, col]
            color = "white" if value > 0.55 * max_value else "black"
            stroke_color = "black" if color == "white" else "white"
            text = ax.text(
                col,
                row,
                f"{value:.1f}",
                ha="center",
                va="center",
                color=color,
                fontsize=11,
                fontweight="bold",
            )
            text.set_path_effects(
                [path_effects.withStroke(linewidth=1.25, foreground=stroke_color, alpha=0.9)]
            )


def plot_summary(curves: dict[tuple[str, str], np.ndarray], il_nm: float, output: Path):
    reference_curve = curves[("ti", "ti")]
    mean_abs = np.zeros((len(TOP_ELECTRODES), len(BOTTOM_ELECTRODES)))

    for row, top in enumerate(TOP_ELECTRODES):
        for col, bottom in enumerate(BOTTOM_ELECTRODES):
            delta = curves[(top, bottom)] - reference_curve
            mean_abs[row, col] = float(np.mean(np.abs(delta)))

    fig, ax = plt.subplots(figsize=(8.2, 6.8), dpi=200)
    image = ax.imshow(mean_abs, cmap="viridis", interpolation="nearest")
    annotate_heatmap(ax, mean_abs)
    ax.set_xticks(
        range(len(BOTTOM_ELECTRODES)),
        [pretty(name) for name in BOTTOM_ELECTRODES],
        rotation=35,
        ha="right",
        fontsize=11,
    )
    ax.set_yticks(range(len(TOP_ELECTRODES)), [pretty(name) for name in TOP_ELECTRODES], fontsize=11)
    ax.set_xlabel("Bottom electrode", fontsize=12)
    ax.set_ylabel("Top electrode", fontsize=12)
    ax.set_xticks(np.arange(-0.5, len(BOTTOM_ELECTRODES), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(TOP_ELECTRODES), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.6, alpha=0.75)
    ax.tick_params(which="minor", bottom=False, left=False)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="uC/cm^2")
    colorbar.ax.tick_params(labelsize=10)
    colorbar.set_label("uC/cm^2", size=11)

    fig.suptitle(f"Mean |Delta P| Relative to Ti/Ti at {il_nm:g} nm IL", y=0.98, fontsize=15)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Generate electrode-comparison figures at a fixed IL thickness.")
    parser.add_argument("--il-nm", type=float, default=0.0, help="IL thickness to compare.")
    parser.add_argument("--output-dir", default="45nm_figures", help="Directory for generated figures.")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    single_cycle, curves = compute_curves(il_nm=args.il_nm)
    plot_overlay(
        single_cycle=single_cycle,
        curves=curves,
        il_nm=args.il_nm,
        output=output_dir / f"all_electrodes_{args.il_nm:g}nm_overlay.png",
    )
    plot_summary(
        curves=curves,
        il_nm=args.il_nm,
        output=output_dir / f"all_electrodes_{args.il_nm:g}nm_deltaP_summary.png",
    )


if __name__ == "__main__":
    main()
