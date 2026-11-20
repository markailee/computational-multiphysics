#!/usr/bin/env python3
import argparse
import csv
import math
import re
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
OUTPUT = PROJECT / "validation"
DATA = OUTPUT / "data"
FIGURES = OUTPUT / "figures"
PI = math.pi


def case(name, directory, input_name, stem, family, level=None, h=None):
    return {"name": name, "directory": ROOT / directory, "input": ROOT / directory / input_name,
            "stem": stem, "family": family, "level": level, "h": h}


CASES = [
    case("cavity_coarse", "acoustic_cavity", "coarse.in", "cavity_coarse", "cavity", "coarse", 0.1),
    case("cavity_medium", "acoustic_cavity", "case.in", "cavity", "cavity", "medium", 0.05),
    case("cavity_fine", "acoustic_cavity", "fine.in", "cavity_fine", "cavity", "fine", 0.025),
    case("cavity_triangles", "acoustic_cavity", "triangles.in", "cavity_triangles", "triangles"),
    case("cavity_quadrilaterals", "acoustic_cavity", "quadrilaterals.in", "cavity_quadrilaterals", "quadrilaterals"),
    case("cavity_newmark", "acoustic_cavity", "newmark.in", "cavity_newmark", "newmark"),
    case("plane_coarse", "acoustic_plane_wave", "coarse.in", "plane_wave_coarse", "plane", "coarse", 0.05),
    case("plane_medium", "acoustic_plane_wave", "case.in", "plane_wave", "plane", "medium", 0.025),
    case("plane_fine", "acoustic_plane_wave", "fine.in", "plane_wave_fine", "plane", "fine", 0.0125),
    case("convected_coarse", "acoustic_convected", "coarse.in", "convected_coarse", "convected", "coarse", 0.05),
    case("convected_medium", "acoustic_convected", "case.in", "convected", "convected", "medium", 0.025),
    case("convected_fine", "acoustic_convected", "fine.in", "convected_fine", "convected", "fine", 0.0125),
    case("absorbing", "acoustic_absorbing", "case.in", "absorbing", "absorbing"),
    case("rigid", "acoustic_absorbing", "rigid.in", "rigid", "rigid"),
]


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def result_path(spec, step=None):
    summary = spec["directory"] / "results" / f"{spec['stem']}_summary.csv"
    steps = len(read_rows(summary)) - 1
    if step is None:
        step = steps
    return spec["directory"] / "results" / f"{spec['stem']}_{step}.plt"


def read_nodal_tecplot(path):
    with path.open(encoding="utf-8") as stream:
        line = ""
        while line.strip() != "VARIABLES" and not line.startswith("VARIABLES"):
            line = stream.readline()
            if not line:
                raise ValueError(f"No variables in {path}")
        names = re.findall(r'"([^"]+)"', line)
        result = {"cells": [], "time": 0.0}
        for zone_line in stream:
            if not zone_line.startswith("ZONE"):
                continue
            node_count = int(re.search(r"\bN=(\d+)", zone_line).group(1))
            element_count = int(re.search(r"\bE=(\d+)", zone_line).group(1))
            result["time"] = float(re.search(r"SOLUTIONTIME=([0-9.eE+-]+)", zone_line).group(1))
            blocks = []
            for _ in names:
                blocks.append([float(stream.readline()) for _ in range(node_count)])
            if "x" not in result:
                for name, values in zip(names, blocks):
                    result[name.lower()] = values
                result["x"] = result.pop("x")
                result["y"] = result.pop("y")
            cell_size = 3 if "FETRIANGLE" in zone_line else 4
            for _ in range(element_count):
                connectivity = tuple(int(value) - 1 for value in stream.readline().split())
                if len(connectivity) != cell_size:
                    raise ValueError(f"Invalid connectivity in {path}")
                result["cells"].append(connectivity)
    required = ("x", "y", "p", "dual_area")
    if any(name not in result for name in required):
        raise ValueError(f"Missing nodal field in {path}")
    return result


def analytical(family, x, y, time):
    if family in ("cavity", "triangles", "quadrilaterals", "newmark"):
        omega = PI * math.sqrt(2.0)
        return math.cos(PI * x) * math.cos(PI * y) * math.cos(omega * time)
    if family in ("plane", "convected"):
        mach = 0.3 if family == "convected" else 0.0
        return math.sin(2.0 * PI * x - 2.0 * PI * (1.0 + mach) * time)
    if family in ("absorbing", "rigid"):
        return math.exp(-((x - 1.0 - time) ** 2) / (2.0 * 0.12 ** 2))
    raise ValueError(f"Unknown analytical family: {family}")


def error_metrics(result, family):
    exact = [analytical(family, x, y, result["time"]) for x, y in zip(result["x"], result["y"])]
    weights = result["dual_area"]
    numerator = sum(w * (p - reference) ** 2 for w, p, reference in zip(weights, result["p"], exact))
    denominator = sum(w * reference ** 2 for w, reference in zip(weights, exact))
    return {
        "l2_relative": math.sqrt(numerator / denominator),
        "linf_absolute": max(abs(p - reference) for p, reference in zip(result["p"], exact)),
        "exact": exact,
    }


def probe_rows(spec):
    return read_rows(spec["directory"] / "results" / f"{spec['stem']}_probe.csv")


def crossings(rows):
    times = [float(row["time"]) for row in rows]
    values = [float(row["p"]) for row in rows]
    result = []
    if abs(values[0]) < 1e-10:
        result.append(times[0])
    for i in range(1, len(values)):
        if values[i - 1] == 0.0 or values[i] == 0.0 or values[i - 1] * values[i] < 0.0:
            fraction = values[i - 1] / (values[i - 1] - values[i])
            crossing = times[i - 1] + fraction * (times[i] - times[i - 1])
            if not result or crossing - result[-1] > 1e-8:
                result.append(crossing)
    return result


def measured_frequency(rows):
    points = crossings(rows)
    if len(points) < 2:
        return math.nan
    return PI * (len(points) - 1) / (points[-1] - points[0])


def geometry_metrics(spec):
    rows = read_rows(spec["directory"] / "results" / f"{spec['stem']}_geometry.csv")
    return (max(float(row["area_error"]) for row in rows),
            max(float(row["gradient_sum"]) for row in rows))


def comparison_data(spec, result, exact):
    rows = [{"x": x, "y": y, "numerical": p, "analytical": reference, "error": p - reference,
             "dual_area": area} for x, y, p, reference, area in
            zip(result["x"], result["y"], result["p"], exact, result["dual_area"])]
    write_rows(DATA / f"{spec['name']}_comparison.csv",
               ("x", "y", "numerical", "analytical", "error", "dual_area"), rows)
    target_y = 0.5 if spec["family"] in ("cavity", "triangles", "quadrilaterals", "newmark", "absorbing", "rigid") else 0.25
    line = [{key: row[key] for key in ("x", "numerical", "analytical", "error")}
            for row in rows if abs(row["y"] - target_y) < 1e-9]
    line.sort(key=lambda row: row["x"])
    write_rows(DATA / f"{spec['name']}_line.csv", ("x", "numerical", "analytical", "error"), line)
    return line


def color(value):
    value = min(1.0, max(0.0, value))
    if value < 0.5:
        fraction = 2.0 * value
        rgb = (49 + 206 * fraction, 54 + 201 * fraction, 149 + 106 * fraction)
    else:
        fraction = 2.0 * (value - 0.5)
        rgb = (255 - 75 * fraction, 255 - 190 * fraction, 255 - 205 * fraction)
    return "#%02x%02x%02x" % tuple(round(component) for component in rgb)


def contour_svg(path, result, exact, title):
    triangles = []
    for cell in result["cells"]:
        if len(cell) == 3:
            triangles.append(cell)
        else:
            triangles.extend(((cell[0], cell[1], cell[2]), (cell[0], cell[2], cell[3])))
    xmin, xmax = min(result["x"]), max(result["x"])
    ymin, ymax = min(result["y"]), max(result["y"])
    vmax = max(abs(min(result["p"] + exact)), abs(max(result["p"] + exact)), 1e-14)
    width, height = 1100, 470
    panels = (("Numerical", result["p"]), ("Analytical", exact))
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="19">{escape(title)}</text>']
    for panel, (label, values) in enumerate(panels):
        left = 45 + panel * 530
        top = 55
        available_width, available_height = 480, 320
        scale = min(available_width / (xmax - xmin), available_height / (ymax - ymin))
        plot_width = scale * (xmax - xmin)
        plot_height = scale * (ymax - ymin)
        offset_x = left + (available_width - plot_width) / 2
        offset_y = top + (available_height - plot_height) / 2

        def point(node):
            return (offset_x + (result["x"][node] - xmin) * scale,
                    offset_y + plot_height - (result["y"][node] - ymin) * scale)

        for triangle in triangles:
            coordinates = " ".join(f"{x:.3f},{y:.3f}" for x, y in map(point, triangle))
            value = sum(values[node] for node in triangle) / 3
            band = round((value / vmax + 1.0) * 0.5 * 40) / 40
            parts.append(f'<polygon points="{coordinates}" fill="{color(band)}" stroke="none"/>')
        parts.append(f'<rect x="{offset_x:.3f}" y="{offset_y:.3f}" width="{plot_width:.3f}" height="{plot_height:.3f}" fill="none" stroke="#333"/>')
        parts.append(f'<text x="{offset_x + plot_width / 2:.3f}" y="{offset_y + plot_height + 28:.3f}" text-anchor="middle" font-family="sans-serif" font-size="16">{label}</text>')
        parts.append(f'<text x="{offset_x + plot_width / 2:.3f}" y="{offset_y + plot_height + 48:.3f}" text-anchor="middle" font-family="sans-serif" font-size="12">x</text>')
    for i in range(80):
        parts.append(f'<rect x="{350 + i * 5}" y="420" width="5" height="15" fill="{color(i / 79)}"/>')
    parts.extend(['<rect x="350" y="420" width="400" height="15" fill="none" stroke="#333"/>',
                  f'<text x="345" y="434" text-anchor="end" font-family="sans-serif" font-size="12">{-vmax:.3g}</text>',
                  '<text x="550" y="434" text-anchor="middle" font-family="sans-serif" font-size="12">0</text>',
                  f'<text x="755" y="434" text-anchor="start" font-family="sans-serif" font-size="12">{vmax:.3g}</text>',
                  '<text x="550" y="455" text-anchor="middle" font-family="sans-serif" font-size="12">acoustic pressure</text>',
                  '</svg>'])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def line_svg(path, series, title, xlabel, ylabel, width=900, height=500, legend_left=False):
    left, right, top, bottom = 80, 30, 55, 65
    all_points = [point for item in series for point in item["points"]]
    xmin, xmax = min(p[0] for p in all_points), max(p[0] for p in all_points)
    ymin, ymax = min(p[1] for p in all_points), max(p[1] for p in all_points)
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymax = ymin + 1.0
    padding = 0.08 * (ymax - ymin)
    ymin -= padding
    ymax += padding

    def point(x, y):
        return (left + (x - xmin) / (xmax - xmin) * (width - left - right),
                top + (ymax - y) / (ymax - ymin) * (height - top - bottom))

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="19">{escape(title)}</text>']
    for index in range(6):
        x_value = xmin + (xmax - xmin) * index / 5
        y_value = ymin + (ymax - ymin) * index / 5
        x, _ = point(x_value, ymin)
        _, y = point(xmin, y_value)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{height - bottom}" stroke="#e5e5e5"/>')
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}" stroke="#e5e5e5"/>')
        parts.append(f'<text x="{x:.2f}" y="{height - bottom + 20}" text-anchor="middle" font-family="sans-serif" font-size="11">{x_value:.3g}</text>')
        parts.append(f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="11">{y_value:.3g}</text>')
    parts.append(f'<rect x="{left}" y="{top}" width="{width - left - right}" height="{height - top - bottom}" fill="none" stroke="#333"/>')
    for item in series:
        points = " ".join(f"{x:.2f},{y:.2f}" for x, y in (point(*p) for p in item["points"]))
        dash = f' stroke-dasharray="{item["dash"]}"' if "dash" in item else ""
        parts.append(f'<polyline points="{points}" fill="none" stroke="{item["color"]}" stroke-width="1.8"{dash}/>')
    parts.extend([f'<text x="{width / 2}" y="{height - 15}" text-anchor="middle" font-family="sans-serif" font-size="14">{escape(xlabel)}</text>',
                  f'<text x="18" y="{height / 2}" text-anchor="middle" transform="rotate(-90 18 {height / 2})" font-family="sans-serif" font-size="14">{escape(ylabel)}</text>'])
    for index, item in enumerate(series):
        x = left + 20 if legend_left else width - 210
        y = top + 20 + 22 * index
        parts.append(f'<line x1="{x}" y1="{y}" x2="{x + 28}" y2="{y}" stroke="{item["color"]}" stroke-width="2"/>')
        parts.append(f'<text x="{x + 35}" y="{y + 4}" font-family="sans-serif" font-size="12">{escape(item["label"])}</text>')
    parts.append('</svg>')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def pulse_step(summary):
    return min(range(len(summary)), key=lambda index: abs(float(summary[index]["time"]) - 1.0))


def analyze(spec, failures):
    summary = read_rows(spec["directory"] / "results" / f"{spec['stem']}_summary.csv")
    step = pulse_step(summary) if spec["family"] in ("absorbing", "rigid") else None
    result = read_nodal_tecplot(result_path(spec, step))
    metrics = error_metrics(result, spec["family"])
    area_error, gradient_sum = geometry_metrics(spec)
    energy0 = float(summary[0]["energy"])
    energy_drift = (float(summary[-1]["energy"]) - energy0) / max(abs(energy0), 1e-14)
    row = {"case": spec["name"], "family": spec["family"], "level": spec["level"] or "",
           "h": spec["h"] if spec["h"] is not None else "", "nodes": len(result["x"]),
           "cells": len(result["cells"]), "time": result["time"], "l2_relative": metrics["l2_relative"],
           "linf_absolute": metrics["linf_absolute"], "observed_order": "",
           "frequency_relative_error": "", "phase_speed": "", "energy_drift": energy_drift,
           "max_area_error": area_error, "max_gradient_sum": gradient_sum,
           "reflection_ratio": "", "rigid_reflection_ratio": ""}
    if area_error > 1e-9 or gradient_sum > 1e-8:
        failures.append(f"{spec['name']}: geometry checks failed")
    if not math.isfinite(metrics["l2_relative"]):
        failures.append(f"{spec['name']}: non-finite error")
    if spec["family"] in ("cavity", "triangles", "quadrilaterals", "newmark"):
        frequency = measured_frequency(probe_rows(spec))
        exact_frequency = PI * math.sqrt(2.0)
        row["frequency_relative_error"] = (frequency - exact_frequency) / exact_frequency
        if metrics["l2_relative"] > 0.06 or abs(row["frequency_relative_error"]) > 0.02:
            failures.append(f"{spec['name']}: cavity accuracy gate failed")
        if abs(energy_drift) > 0.05:
            failures.append(f"{spec['name']}: energy drift gate failed")
    if spec["family"] in ("plane", "convected"):
        frequency = measured_frequency(probe_rows(spec))
        exact_speed = 1.3 if spec["family"] == "convected" else 1.0
        row["phase_speed"] = frequency / (2.0 * PI)
        row["frequency_relative_error"] = (frequency - 2.0 * PI * exact_speed) / (2.0 * PI * exact_speed)
        if metrics["l2_relative"] > 0.06 or abs(row["phase_speed"] - exact_speed) > 0.03:
            failures.append(f"{spec['name']}: plane-wave accuracy gate failed")
    if spec["name"] in ("cavity_medium", "plane_medium", "convected_medium"):
        line = comparison_data(spec, result, metrics["exact"])
        contour_svg(FIGURES / f"{spec['name']}_contour.svg", result, metrics["exact"],
                    f"{spec['name'].replace('_', ' ')} at t = {result['time']:.6g}")
        line_svg(FIGURES / f"{spec['name']}_line.svg", [
            {"points": [(float(item["x"]), float(item["numerical"])) for item in line], "label": "Numerical", "color": "#2166ac"},
            {"points": [(float(item["x"]), float(item["analytical"])) for item in line], "label": "Analytical", "color": "#b2182b", "dash": "6 3"},
        ], f"{spec['name'].replace('_', ' ')} centreline", "x", "pressure")
    return row


def absorption_analysis(rows_by_name, failures):
    histories = {}
    ratios = {}
    combined = []
    for name in ("absorbing", "rigid"):
        spec = next(item for item in CASES if item["name"] == name)
        history = [(float(row["time"]), float(row["p"])) for row in probe_rows(spec)]
        histories[name] = history
        incident = max(abs(p) for t, p in history if 2.0 <= t <= 2.9)
        reflected = max(abs(p) for t, p in history if t >= 3.1)
        ratios[name] = reflected / incident
        for t, p in history:
            combined.append({"time": t, "boundary": name, "pressure": p})
    write_rows(DATA / "acoustic_absorption_history.csv", ("time", "boundary", "pressure"), combined)
    line_svg(FIGURES / "acoustic_absorption_history.svg", [
        {"points": histories["absorbing"], "label": "Absorbing outlet", "color": "#2166ac"},
        {"points": histories["rigid"], "label": "Rigid outlet", "color": "#b2182b"},
    ], "Probe at x = 3.5: incident and reflected Gaussian pulse", "time", "pressure", legend_left=True)
    rows_by_name["absorbing"]["reflection_ratio"] = ratios["absorbing"]
    rows_by_name["absorbing"]["rigid_reflection_ratio"] = ratios["rigid"]
    rows_by_name["rigid"]["reflection_ratio"] = ratios["rigid"]
    if ratios["absorbing"] > 0.08 or ratios["rigid"] < 0.7:
        failures.append(f"absorption gate failed: absorbing={ratios['absorbing']:.4g}, rigid={ratios['rigid']:.4g}")

    spec = next(item for item in CASES if item["name"] == "absorbing")
    summary = read_rows(spec["directory"] / "results" / f"{spec['stem']}_summary.csv")
    result = read_nodal_tecplot(result_path(spec, pulse_step(summary)))
    metrics = error_metrics(result, "absorbing")
    comparison_data(spec, result, metrics["exact"])
    contour_svg(FIGURES / "acoustic_pulse_contour.svg", result, metrics["exact"],
                f"Directed Gaussian pulse at t = {result['time']:.6g}")
    rows = read_rows(DATA / "absorbing_line.csv")
    line_svg(FIGURES / "acoustic_pulse_line.svg", [
        {"points": [(float(row["x"]), float(row["numerical"])) for row in rows], "label": "Numerical", "color": "#2166ac"},
        {"points": [(float(row["x"]), float(row["analytical"])) for row in rows], "label": "Analytical", "color": "#b2182b", "dash": "6 3"},
    ], "Gaussian pulse centreline at t = 1", "x", "pressure")
    rows_by_name["absorbing"]["l2_relative"] = metrics["l2_relative"]
    rows_by_name["absorbing"]["linf_absolute"] = metrics["linf_absolute"]
    rows_by_name["absorbing"]["time"] = result["time"]
    if metrics["l2_relative"] > 0.04:
        failures.append(f"Gaussian pulse accuracy gate failed: {metrics['l2_relative']:.4g}")


def convergence_analysis(rows):
    for family in ("cavity", "plane", "convected"):
        group = sorted((row for row in rows if row["family"] == family), key=lambda row: row["h"], reverse=True)
        for fine, coarse in zip(group[1:], group):
            ratio = fine["l2_relative"] / coarse["l2_relative"]
            fine["observed_order"] = math.log(coarse["l2_relative"] / fine["l2_relative"]) / math.log(2.0) if ratio > 0 else math.nan
        if group[-1]["l2_relative"] >= group[0]["l2_relative"] / 2:
            raise RuntimeError(f"{family}: insufficient grid-error reduction")
        line_svg(FIGURES / f"acoustic_{family}_convergence.svg", [
            {"points": [(row["h"], row["l2_relative"]) for row in reversed(group)],
             "label": "Relative L2 error", "color": "#2166ac"},
        ], f"{family.replace('_', ' ')} grid convergence", "mesh spacing h", "relative L2 error")


def main():
    parser = argparse.ArgumentParser(description="Run and analyze acoustic analytical validations")
    parser.add_argument("--run", action="store_true", help="Generate meshes and run every acoustic case")
    args = parser.parse_args()
    failures = []
    if args.run:
        subprocess.run(["python3", str(ROOT / "generate_acoustic_meshes.py")], cwd=PROJECT, check=True)
        for spec in CASES:
            print(f"\n=== {spec['name']} ===", flush=True)
            subprocess.run([str(PROJECT / "multiphysics"), str(spec["input"])], cwd=PROJECT, check=True)
    rows = [analyze(spec, failures) for spec in CASES]
    rows_by_name = {row["case"]: row for row in rows}
    absorption_analysis(rows_by_name, failures)
    try:
        convergence_analysis(rows)
    except RuntimeError as error:
        failures.append(str(error))
    fieldnames = ("case", "family", "level", "h", "nodes", "cells", "time", "l2_relative", "linf_absolute",
                  "observed_order", "frequency_relative_error", "phase_speed", "energy_drift",
                  "max_area_error", "max_gradient_sum", "reflection_ratio", "rigid_reflection_ratio")
    write_rows(OUTPUT / "acoustic_metrics.csv", fieldnames, rows)
    convergence_rows = [{key: row[key] for key in ("case", "family", "level", "h", "l2_relative", "observed_order")}
                        for row in rows if row["level"]]
    write_rows(DATA / "acoustic_convergence.csv",
               ("case", "family", "level", "h", "l2_relative", "observed_order"), convergence_rows)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)
    print(f"Acoustic validation passed: {len(rows)} runs; metrics and figures written to {OUTPUT}")


if __name__ == "__main__":
    main()
