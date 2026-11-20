#!/usr/bin/env python3
import argparse
import csv
import math
import re
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path
from xml.sax.saxutils import escape

from validate_acoustic import color, line_svg, read_nodal_tecplot, read_rows, write_rows

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
CASE = ROOT / "cylinder_aeolian"
OUTPUT = PROJECT / "validation"
DATA = OUTPUT / "data"
FIGURES = OUTPUT / "figures"
STROUHAL_REFERENCE = 0.183
MEAN_CD_REFERENCE = 1.32
LIFT_AMPLITUDE_REFERENCE = 0.52


def harmonic_amplitude(times, values, frequency):
    if len(times) < 3:
        return math.nan
    duration = times[-1] - times[0]
    real = 0.0
    imaginary = 0.0
    for index, (time, value) in enumerate(zip(times, values)):
        if index == 0:
            weight = (times[1] - time) / 2
        elif index == len(times) - 1:
            weight = (time - times[index - 1]) / 2
        else:
            weight = (times[index + 1] - times[index - 1]) / 2
        phase = 2.0 * math.pi * frequency * time
        real += weight * value * math.cos(phase)
        imaginary -= weight * value * math.sin(phase)
    return 2.0 * math.hypot(real, imaginary) / duration


def shedding_period(rows):
    start = max(20.0, 0.4 * float(rows[-1]["time"]))
    selected = [row for row in rows if float(row["time"]) >= start]
    times = [float(row["time"]) for row in selected]
    lift = [float(row["cl"]) for row in selected]
    half_window = 5
    times = times[half_window:-half_window]
    lift = [statistics.fmean(lift[index - half_window:index + half_window + 1])
            for index in range(half_window, len(lift) - half_window)]
    mean = statistics.fmean(lift)
    crossings = []
    for index in range(1, len(times)):
        before = lift[index - 1] - mean
        after = lift[index] - mean
        if before <= 0.0 < after:
            fraction = before / (before - after)
            crossings.append(times[index - 1] + fraction * (times[index] - times[index - 1]))
    if len(crossings) < 3:
        raise RuntimeError("Fewer than two complete shedding cycles were detected; extend the flow run")
    periods = [second - first for first, second in zip(crossings, crossings[1:])]
    period = statistics.median(periods[-min(5, len(periods)):])
    if not 1.0 < period < 20.0:
        raise RuntimeError(f"Implausible detected shedding period: {period}")
    return period


def prepare_source():
    rows = read_rows(CASE / "results" / "flow_forces.csv")
    if len(rows) < 100:
        raise RuntimeError("The flow force history is too short for acoustic source preparation")
    period = shedding_period(rows)
    end_time = float(rows[-1]["time"])
    start_time = end_time - 3.0 * period
    selected = [row for row in rows if float(row["time"]) >= start_time - 1e-10]
    mean_x = statistics.fmean(float(row["force_x"]) for row in selected)
    mean_y = statistics.fmean(float(row["force_y"]) for row in selected)
    source_rows = [{"time": float(row["time"]) - float(selected[0]["time"]),
                    "force_x": float(row["force_x"]) - mean_x,
                    "force_y": float(row["force_y"]) - mean_y} for row in selected]
    source_rows[0]["time"] = 0.0
    write_rows(CASE / "source" / "flow_dipole.csv", ("time", "force_x", "force_y"), source_rows)
    metadata = {"period": period, "frequency": 1.0 / period, "start_time": float(selected[0]["time"]),
                "end_time": end_time, "source_rows": source_rows}
    print(f"Detected shedding period {period:.8f}, St = {1.0 / period:.6f}; "
          f"source window {selected[0]['time']} to {end_time}")
    return metadata


def read_cell_zones(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    names = []
    nodes = {"X": [], "Y": []}
    cells = []
    values = defaultdict(list)
    while index < len(lines):
        line = lines[index].strip()
        if line.startswith("VARIABLES"):
            names = re.findall(r'"([^"]+)"', line)
            index += 1
            continue
        if not line.startswith("ZONE"):
            index += 1
            continue
        node_count = int(re.search(r"\bN=(\d+)", line).group(1))
        element_count = int(re.search(r"\bE=(\d+)", line).group(1))
        element_size = 3 if "FETRIANGLE" in line else 4
        index += 1
        coordinates = {}
        for name in names[:2]:
            coordinates[name] = [float(value) for value in lines[index:index + node_count]]
            index += node_count
        if not nodes["X"]:
            nodes = coordinates
        zone_values = {}
        for name in names[2:]:
            zone_values[name] = [float(value) for value in lines[index:index + element_count]]
            index += element_count
            values[name].extend(zone_values[name])
        for _ in range(element_count):
            cells.append(tuple(int(value) - 1 for value in lines[index].split()))
            if len(cells[-1]) != element_size:
                raise ValueError(f"Invalid connectivity in {path}")
            index += 1
    return nodes, cells, values


def scalar_contour_svg(path, nodes, cells, values, title, label, limits=None, marker=False, maximum_value=None):
    if limits:
        xmin, xmax, ymin, ymax = limits
    else:
        xmin, xmax = min(nodes["X"]), max(nodes["X"])
        ymin, ymax = min(nodes["Y"]), max(nodes["Y"])
    maximum = max(abs(min(values)), abs(max(values)), 1e-30)
    if maximum_value is not None:
        maximum = maximum_value
    width, height = 1000, 540
    left, top, plot_width, plot_height = 55, 55, 880, 390
    scale = min(plot_width / (xmax - xmin), plot_height / (ymax - ymin))
    drawn_width = scale * (xmax - xmin)
    drawn_height = scale * (ymax - ymin)
    offset_x = left + (plot_width - drawn_width) / 2
    offset_y = top + (plot_height - drawn_height) / 2
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="19">{escape(title)}</text>',
             f'<defs><clipPath id="plot"><rect x="{offset_x}" y="{offset_y}" width="{drawn_width}" height="{drawn_height}"/></clipPath></defs>',
             '<g clip-path="url(#plot)">']

    def point(node):
        return (offset_x + (nodes["X"][node] - xmin) * scale,
                offset_y + drawn_height - (nodes["Y"][node] - ymin) * scale)

    nodal = len(values) == len(nodes["X"])
    for cell_index, cell in enumerate(cells):
        cell_x = sum(nodes["X"][node] for node in cell) / len(cell)
        cell_y = sum(nodes["Y"][node] for node in cell) / len(cell)
        if not (xmin <= cell_x <= xmax and ymin <= cell_y <= ymax):
            continue
        coordinates = " ".join(f"{x:.3f},{y:.3f}" for x, y in map(point, cell))
        value = sum(values[node] for node in cell) / len(cell) if nodal else values[cell_index]
        band = round((value / maximum + 1.0) * 0.5 * 50) / 50
        parts.append(f'<polygon points="{coordinates}" fill="{color(band)}" stroke="none"/>')
    parts.append('</g>')
    parts.append(f'<rect x="{offset_x}" y="{offset_y}" width="{drawn_width}" height="{drawn_height}" fill="none" stroke="#333"/>')
    if marker:
        x = offset_x + (0.0 - xmin) * scale
        y = offset_y + drawn_height - (0.0 - ymin) * scale
        parts.append(f'<circle cx="{x}" cy="{y}" r="{0.5 * scale}" fill="none" stroke="#111" stroke-width="1.5"/>')
    for index in range(80):
        parts.append(f'<rect x="{300 + index * 5}" y="475" width="5" height="15" fill="{color(index / 79)}"/>')
    parts.extend(['<rect x="300" y="475" width="400" height="15" fill="none" stroke="#333"/>',
                  f'<text x="295" y="489" text-anchor="end" font-family="sans-serif" font-size="12">{-maximum:.3g}</text>',
                  '<text x="500" y="489" text-anchor="middle" font-family="sans-serif" font-size="12">0</text>',
                  f'<text x="705" y="489" text-anchor="start" font-family="sans-serif" font-size="12">{maximum:.3g}</text>',
                  f'<text x="500" y="515" text-anchor="middle" font-family="sans-serif" font-size="13">{escape(label)}</text>',
                  '</svg>'])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def write_report(metadata, metrics, cycle_rows):
    cycle_table = "\n".join(
        f"| {row['cycle']} | {row['start_time']:.4f}–{row['end_time']:.4f} | {row['mean_cd']:.6f} | {row['lift_amplitude']:.6f} |"
        for row in cycle_rows)
    st_error = 100.0 * (metrics["strouhal"] / STROUHAL_REFERENCE - 1.0)
    cd_error = 100.0 * (metrics["mean_cd"] / MEAN_CD_REFERENCE - 1.0)
    lift_error = 100.0 * (metrics["lift_amplitude"] / LIFT_AMPLITUDE_REFERENCE - 1.0)
    report = f"""# Cylinder vortex shedding and aeolian-tone propagation

## Result summary

The one-way hybrid workflow completed successfully. The transient incompressible cylinder flow produces a Kármán vortex street, and its mean-subtracted force history drives a compact Curle-type dipole in the convected acoustic solver.

- Measured Strouhal number: **{metrics['strouhal']:.6f}**, compared with **{STROUHAL_REFERENCE}** from Inoue and Hatakeyama ({st_error:+.2f}%).
- Mean drag coefficient: **{metrics['mean_cd']:.6f}**, compared with **{MEAN_CD_REFERENCE}** ({cd_error:+.2f}%).
- Lift half-range amplitude: **{metrics['lift_amplitude']:.6f}**, compared with **{LIFT_AMPLITUDE_REFERENCE}** ({lift_error:+.2f}%).
- Drag second-harmonic amplitude: **{metrics['drag_second_harmonic_amplitude']:.6f}**, compared with approximately **0.026** in the reference.
- Transverse/axial fundamental acoustic amplitude ratio at `r = 10 D`: **{metrics['acoustic_directivity_ratio']:.3f}**, confirming a dominant lift dipole.

The force comparison is quantitative. The acoustic calculation verifies propagation and dipole directivity, but it is not a compressible DNS or a grid-independent far-field pressure certification.

## Benchmark and discretization

The reference is O. Inoue and N. Hatakeyama, [Sound generation by a two-dimensional circular cylinder in a uniform flow](https://metodologiauntref.wordpress.com/wp-content/uploads/2015/09/inoue2002.pdf), *Journal of Fluid Mechanics* 471, 285–314, 2002, DOI `10.1017/S0022112002002124`.

The case uses `D = U = rho = 1`, `nu = 1/150`, `Re = 150`, `c = 5`, and `M = 0.2`. The near field uses the existing 2,304-cell body-fitted quadrilateral cylinder mesh and BDF2 with `dt = 0.02`. The acoustic domain is `[-20,20] x [-15,15]`, with 4,941 nodes and 7,200 mixed triangles/quadrilaterals, centered differencing with `dt = 0.02`, and maximum Courant number `0.24`. All acoustic exterior boundaries are absorbing.

The detected shedding period is **{metrics['shedding_period']:.8f}**. The acoustic source uses the last three detected periods, from flow time **{metadata['start_time']:.6f}** to **{metadata['end_time']:.6f}**, with the mean force removed and a one-time-unit cosine startup ramp.

## Cycle-by-cycle force development

| Cycle (oldest to latest) | Flow-time interval | Mean `Cd` | `Cl` half-range |
|---:|---|---:|---:|
{cycle_table}

These intervals show whether the selected source has reached a statistically stationary shedding state. The complete three-period harmonic fit gives a lift fundamental amplitude of **{metrics['lift_fundamental_amplitude']:.6f}**; this differs from the half-range when amplitude is still developing.

![Cylinder lift and drag histories](figures/cylinder_aeolian_forces.svg)

## Near-field vortex shedding

The vorticity plot uses a fixed `+/-2` color range to keep the wake structures visible despite much larger boundary-layer values at the cylinder wall. The underlying Tecplot data are not clipped.

![Cylinder wake vorticity](figures/cylinder_aeolian_vorticity.svg)

## Curle-type acoustic source

The flow output contains total, pressure, and viscous forces per unit span. For force `F` **of the fluid on the cylinder**, the force of the cylinder on the fluid is `-F`. With normalized compact Gaussian `g`, the acoustic right-hand side is

\\[
-\\nabla\\cdot\\left(-\\frac{{\\mathbf Fg}}{{\\rho}}\\right)
=\\nabla\\cdot\\left(\\frac{{\\mathbf Fg}}{{\\rho}}\\right).
\\]

The acoustic solver evaluates this dipole through the Gaussian gradient and linearly interpolates the measured force history. It does not use the incompressible pressure as radiating sound. Mean flow `U = (1,0)` is retained in the acoustic wave operator, so propagation includes the low-Mach Doppler effect.

## Acoustic pressure and directivity

The instantaneous contour shows alternating positive and negative pressure lobes generated by the oscillating transverse force, together with upstream/downstream asymmetry from convection.

![Flow-driven acoustic pressure](figures/cylinder_aeolian_pressure.svg)

At `r = 10 D`, the fundamental-pressure amplitudes are:

| Probe location | Fundamental amplitude |
|---|---:|
| `(0, +10)` | {metrics['acoustic_transverse_amplitude']:.8g} (mean with `(0,-10)`) |
| `(+/-10, 0)` | {metrics['acoustic_axial_amplitude']:.8g} (mean of both axial probes) |

![Measured and ideal lift-dipole directivity](figures/cylinder_aeolian_directivity.svg)

The directivity comparison normalizes the ideal transverse-dipole shape to the measured transverse amplitude. Axial probes also include weaker drag-dipole and near-field/numerical contributions, so the ideal curve is a shape reference rather than an exact error norm.

## Reproduction

```bash
make -j2
make cylinder-aeolian
```

To reuse the existing near-field result and rerun only acoustics:

```bash
python3 cases/run_cylinder_aeolian.py --acoustic-only
```

The machine-readable summary is `validation/cylinder_aeolian_metrics.csv`. Detailed force and probe data are under `validation/data/`; case-specific settings and source conventions are in `cases/cylinder_aeolian/README.md`.

## Limitations

- This is one-way Curle-type coupling, not two-way fluid–acoustic coupling or compressible DNS.
- The coarse near-field mesh and bounded flow domain affect force amplitudes even though the shedding frequency is accurate.
- The compact Gaussian dipole is a regularization of the cylinder loading; its width and the acoustic grid affect near-source pressure amplitudes.
- The acoustic domain is much smaller than the reference DNS radiation domain. First-order absorbing boundaries and low-order mixed-grid dispersion limit far-field comparisons.
- Acoustic pressure amplitudes are therefore reported for this implementation's source and domain, not asserted to match the reference DNS point-for-point.
"""
    (OUTPUT / "cylinder_aeolian.md").write_text(report, encoding="utf-8")


def analyze(metadata):
    period = metadata["period"]
    frequency = metadata["frequency"]
    force_rows = read_rows(CASE / "results" / "flow_forces.csv")
    selected = [row for row in force_rows if float(row["time"]) >= metadata["start_time"] - 1e-10]
    times = [float(row["time"]) for row in selected]
    cd = [float(row["cd"]) for row in selected]
    cl = [float(row["cl"]) for row in selected]
    mean_cd = statistics.fmean(cd)
    lift_amplitude = 0.5 * (max(cl) - min(cl))
    lift_fundamental = harmonic_amplitude(times, cl, frequency)
    drag_second = harmonic_amplitude(times, cd, 2.0 * frequency)
    cycle_rows = []
    for cycle in range(3):
        cycle_start = metadata["end_time"] - (3 - cycle) * period
        cycle_end = metadata["end_time"] - (2 - cycle) * period
        cycle_data = [row for row in selected if cycle_start - 1e-10 <= float(row["time"]) <= cycle_end + 1e-10]
        cycle_cd = [float(row["cd"]) for row in cycle_data]
        cycle_cl = [float(row["cl"]) for row in cycle_data]
        cycle_rows.append({"cycle": cycle + 1, "start_time": cycle_start, "end_time": cycle_end,
                           "mean_cd": statistics.fmean(cycle_cd),
                           "lift_amplitude": 0.5 * (max(cycle_cl) - min(cycle_cl))})
    write_rows(DATA / "cylinder_aeolian_cycles.csv",
               ("cycle", "start_time", "end_time", "mean_cd", "lift_amplitude"), cycle_rows)
    write_rows(DATA / "cylinder_aeolian_forces.csv", ("time", "cd", "cl"),
               [{key: row[key] for key in ("time", "cd", "cl")} for row in selected])
    line_svg(FIGURES / "cylinder_aeolian_forces.svg", [
        {"points": list(zip(times, cl)), "label": "Lift coefficient", "color": "#2166ac"},
        {"points": list(zip(times, cd)), "label": "Drag coefficient", "color": "#b2182b"},
    ], "Cylinder force history over three detected shedding periods", "time", "force coefficient", legend_left=True)

    probes = read_rows(CASE / "results" / "acoustic_probes.csv")
    grouped = defaultdict(list)
    for row in probes:
        grouped[int(row["probe"])].append(row)
    probe_summary = []
    amplitudes = {}
    for probe, rows in sorted(grouped.items()):
        rows = [row for row in rows if float(row["time"]) >= metadata["source_rows"][-1]["time"] - 2.0 * period]
        probe_times = [float(row["time"]) for row in rows]
        pressure = [float(row["p"]) for row in rows]
        amplitude = harmonic_amplitude(probe_times, pressure, frequency)
        x, y = float(rows[0]["x"]), float(rows[0]["y"])
        angle = math.degrees(math.atan2(y, x)) % 360
        amplitudes[(round(x, 6), round(y, 6))] = amplitude
        probe_summary.append({"probe": probe, "x": x, "y": y, "angle": angle,
                              "fundamental_amplitude": amplitude})
    write_rows(DATA / "cylinder_aeolian_probe_amplitudes.csv",
               ("probe", "x", "y", "angle", "fundamental_amplitude"), probe_summary)
    transverse = statistics.fmean([amplitudes[(0.0, 10.0)], amplitudes[(0.0, -10.0)]])
    axial = statistics.fmean([amplitudes[(10.0, 0.0)], amplitudes[(-10.0, 0.0)]])
    directivity = transverse / max(axial, 1e-30)
    angle_groups = defaultdict(list)
    for row in probe_summary:
        angle_groups[row["angle"]].append(row["fundamental_amplitude"])
    directivity_points = [(angle, statistics.fmean(values)) for angle, values in sorted(angle_groups.items())]
    line_svg(FIGURES / "cylinder_aeolian_directivity.svg", [
        {"points": directivity_points, "label": "Acoustic fundamental amplitude", "color": "#2166ac"},
        {"points": [(angle, transverse * math.sin(math.radians(angle)) ** 2) for angle in range(0, 361, 5)],
         "label": "Ideal transverse-dipole shape", "color": "#b2182b", "dash": "6 3"},
    ], "Acoustic directivity at r = 10 D", "angle (degrees)", "pressure amplitude")

    flow_step = len(read_rows(CASE / "results" / "flow_summary.csv"))
    nodes, cells, values = read_cell_zones(CASE / "results" / f"flow_{flow_step}.plt")
    scalar_contour_svg(FIGURES / "cylinder_aeolian_vorticity.svg", nodes, cells, values["vorticity"],
                       "Instantaneous cylinder-wake vorticity", "vorticity (color clipped at +/- 2)",
                       limits=(-2.0, 8.0, -3.0, 3.0), marker=True, maximum_value=2.0)
    acoustic_step = len(read_rows(CASE / "results" / "acoustic_summary.csv")) - 1
    acoustic = read_nodal_tecplot(CASE / "results" / f"acoustic_{acoustic_step}.plt")
    acoustic_nodes = {"X": acoustic["x"], "Y": acoustic["y"]}
    scalar_contour_svg(FIGURES / "cylinder_aeolian_pressure.svg", acoustic_nodes, acoustic["cells"], acoustic["p"],
                       f"Flow-driven acoustic pressure at t = {acoustic['time']:.6g}", "acoustic pressure", marker=True)

    metrics = {"shedding_period": period, "strouhal": frequency, "strouhal_reference": STROUHAL_REFERENCE,
               "mean_cd": mean_cd, "mean_cd_reference": MEAN_CD_REFERENCE,
               "lift_amplitude": lift_amplitude, "lift_amplitude_reference": LIFT_AMPLITUDE_REFERENCE,
               "lift_fundamental_amplitude": lift_fundamental, "drag_second_harmonic_amplitude": drag_second,
               "acoustic_transverse_amplitude": transverse, "acoustic_axial_amplitude": axial,
               "acoustic_directivity_ratio": directivity}
    write_rows(OUTPUT / "cylinder_aeolian_metrics.csv", tuple(metrics), [metrics])
    write_report(metadata, metrics, cycle_rows)
    failures = []
    if not 0.12 <= frequency <= 0.24:
        failures.append(f"Strouhal number {frequency:.4f} is outside the shedding sanity range")
    if not 0.15 <= lift_amplitude <= 0.8:
        failures.append(f"Lift amplitude {lift_amplitude:.4f} is outside the sanity range")
    if not 0.6 <= mean_cd <= 1.8:
        failures.append(f"Mean drag {mean_cd:.4f} is outside the sanity range")
    if directivity < 2.0:
        failures.append(f"Transverse/axial acoustic ratio {directivity:.4f} does not show a lift dipole")
    print(f"St = {frequency:.6f} (reference {STROUHAL_REFERENCE}), mean Cd = {mean_cd:.5f} "
          f"(reference {MEAN_CD_REFERENCE}), lift amplitude = {lift_amplitude:.5f} (reference {LIFT_AMPLITUDE_REFERENCE})")
    print(f"Lift fundamental = {lift_fundamental:.5f}, drag second harmonic = {drag_second:.5f}, "
          f"acoustic transverse/axial ratio = {directivity:.3f}")
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)
    print("Cylinder aeolian-tone workflow passed all flow and acoustic sanity checks")


def main():
    parser = argparse.ArgumentParser(description="Hybrid cylinder-flow/Curle-dipole acoustic workflow")
    parser.add_argument("--run", action="store_true", help="Run flow, source preparation, acoustics, and analysis")
    parser.add_argument("--acoustic-only", action="store_true", help="Prepare source and rerun only the acoustic solver")
    parser.add_argument("--prepare-only", action="store_true", help="Only prepare the acoustic force history")
    args = parser.parse_args()
    if args.run:
        if not (CASE / "mesh" / "acoustic.msh").exists():
            subprocess.run(["python3", str(ROOT / "generate_acoustic_meshes.py")], cwd=PROJECT, check=True)
        subprocess.run([str(PROJECT / "multiphysics"), str(CASE / "flow.in")], cwd=PROJECT, check=True)
    metadata = prepare_source()
    if args.prepare_only:
        return
    if args.run or args.acoustic_only:
        for old_output in (CASE / "results").glob("acoustic_*.plt"):
            old_output.unlink()
        subprocess.run([str(PROJECT / "multiphysics"), str(CASE / "acoustic.in")], cwd=PROJECT, check=True)
    analyze(metadata)


if __name__ == "__main__":
    main()
