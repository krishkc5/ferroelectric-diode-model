import numpy as np


def _build_linear_segment(start: float, stop: float, step: float) -> np.ndarray:
    if step <= 0:
        raise ValueError("Voltage step must be positive.")
    span = abs(stop - start)
    n_points = int(round(span / step)) + 1
    return np.linspace(start, stop, n_points)


def build_waypoint_sweep(waypoints: list[float] | tuple[float, ...], step: float) -> np.ndarray:
    if len(waypoints) < 2:
        raise ValueError("At least two voltage waypoints are required.")

    segments = []
    for index in range(len(waypoints) - 1):
        segment = _build_linear_segment(waypoints[index], waypoints[index + 1], step)
        segments.append(segment if index == 0 else segment[1:])
    return np.concatenate(segments)


def build_bipolar_cycle(vmax: float, step: float) -> np.ndarray:
    return build_waypoint_sweep([-vmax, vmax, -vmax], step)


def repeat_voltage_trace(voltage_trace: np.ndarray, cycles: int) -> np.ndarray:
    if cycles < 1:
        raise ValueError("The voltage trace must be repeated at least once.")
    return np.tile(np.asarray(voltage_trace, dtype=float), cycles)
