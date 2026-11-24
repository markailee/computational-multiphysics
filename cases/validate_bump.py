#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape

from validate import read_tecplot
from validate_acoustic import color, line_svg

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
REFERENCE = PROJECT / "reference" / "bump"
OUTPUT = PROJECT / "validation"
DATA = OUTPUT / "data"
FIGURES = OUTPUT / "figures"
R_GAS = 287.0
CP = 1005.0
GAMMA = CP / (CP - R_GAS)

CASES = [
    {"name": "subsonic", "directory": ROOT / "compressible_bump/subsonic", "stem": "bump_subsonic_final", "total_enthalpy": CP * 300.0 + 0.5 * 173.594354746921 ** 2},
    {"name": "transonic", "directory": ROOT / "compressible_bump/transonic", "stem": "bump_transonic_final", "total_enthalpy": CP * 289.814640852331 + 0.5 * 230.339754553282 ** 2},
    {"name": "supersonic", "directory": ROOT / "compressible_bump/supersonic", "stem": "bump_supersonic_final", "total_enthalpy": CP * 300.0 + 0.5 * 486.06419329138 ** 2},
]

FIELDS = ("case", "cells", "time_steps", "total_outer_iterations", "last_outer_iterations", "converged",
          "max_cell_imbalance", "net_boundary_flux", "total_enthalpy_closure", "relative_l2_total_enthalpy",
          "min_mach", "max_mach", "exit_upper_mach",
          "exit_lower_mach", "reference_upper_mach", "reference_lower_mach", "upper_error",
          "lower_error", "symmetry_error", "shock_x", "profile_rmse", "passed")


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames or tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def interpolate(points, x):
    for first, second in zip(points, points[1:]):
        if first[0] <= x <= second[0]:
            fraction = (x - first[0]) / (second[0] - first[0])
            return first[1] + fraction * (second[1] - first[1])
    return points[0][1] if x < points[0][0] else points[-1][1]


def mach(cell):
    if "Mach" in cell:
        return cell["Mach"]
    return math.hypot(cell["U_x"], cell["U_y"]) / math.sqrt(GAMMA * R_GAS * cell["T"])


def wall_profiles(cells):
    columns = {}
    for cell in cells:
        columns.setdefault(round(cell["x"], 12), []).append(cell)
    profiles = {"upper": [], "lower": []}
    for x, column in sorted(columns.items()):
        lower = min(column, key=lambda cell: cell["y"])
        upper = max(column, key=lambda cell: cell["y"])
        profiles["lower"].append((x, lower))
        profiles["upper"].append((x, upper))
    return profiles


def profile_rmse(profiles, reference):
    errors = []
    for wall, points in profiles.items():
        for x, cell in points:
            errors.append(mach(cell) - interpolate(reference[wall], x))
    return math.sqrt(sum(error * error for error in errors) / len(errors))


def shock_location(points):
    candidates = []
    for first, second in zip(points, points[1:]):
        if first[0] >= 1.4 and second[0] <= 2.2:
            slope = (mach(second[1]) - mach(first[1])) / (second[0] - first[0])
            candidates.append((slope, 0.5 * (first[0] + second[0])))
    return min(candidates)[1] if candidates else 0.0


def contour_svg(path, cells, title):
    width, height = 1050, 430
    left, top, plot_width, plot_height = 65, 55, 920, 275
    values = [mach(cell) for cell in cells]
    vmin, vmax = min(values), max(values)
    span = max(vmax - vmin, 1e-12)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="19">{escape(title)}</text>']
    for cell in cells:
        x = left + cell["x"] / 3.0 * plot_width
        y = top + (1.0 - cell["y"]) * plot_height
        shade = color((mach(cell) - vmin) / span)
        parts.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="2.3" fill="{shade}" stroke="none"/>')
    parts.append(f'<rect x="{left}" y="{top}" width="{plot_width}" height="{plot_height}" fill="none" stroke="#333"/>')
    parts.extend([f'<text x="{left}" y="{top + plot_height + 25}" font-family="sans-serif" font-size="12">0</text>',
                  f'<text x="{left + plot_width / 3}" y="{top + plot_height + 25}" font-family="sans-serif" font-size="12">1</text>',
                  f'<text x="{left + 2 * plot_width / 3}" y="{top + plot_height + 25}" font-family="sans-serif" font-size="12">2</text>',
                  f'<text x="{left + plot_width}" y="{top + plot_height + 25}" text-anchor="end" font-family="sans-serif" font-size="12">3</text>'])
    for index in range(80):
        parts.append(f'<rect x="{325 + index * 5}" y="365" width="5" height="15" fill="{color(index / 79)}"/>')
    parts.extend(['<rect x="325" y="365" width="400" height="15" fill="none" stroke="#333"/>',
                  f'<text x="320" y="379" text-anchor="end" font-family="sans-serif" font-size="12">{vmin:.3f}</text>',
                  f'<text x="730" y="379" font-family="sans-serif" font-size="12">{vmax:.3f}</text>',
                  '<text x="525" y="402" text-anchor="middle" font-family="sans-serif" font-size="12">Mach number</text>',
                  '</svg>'])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def analyze(specification, supersonic_reference):
    directory = specification["directory"]
    stem = specification["stem"]
    summaries = read_rows(directory / "results" / f"{stem}_summary.csv")
    result = directory / "results" / f"{stem}_{len(summaries)}.plt"
    if not result.exists():
        result = directory / "results" / f"{stem}.plt"
    cells = read_tecplot(result)
    profiles = wall_profiles(cells)
    name = specification["name"]
    row = {key: "" for key in FIELDS}
    inlet_enthalpy = specification["total_enthalpy"]
    row.update({"case": name, "cells": len(cells), "time_steps": len(summaries),
                "total_outer_iterations": int(summaries[-1]["iterations"]),
                "last_outer_iterations": int(summaries[-1]["iterations"]),
                "converged": int(summaries[-1]["converged"]),
                "max_cell_imbalance": float(summaries[-1]["max_cell_imbalance"]),
                "net_boundary_flux": abs(float(summaries[-1]["net_boundary_flux"])),
                "total_enthalpy_closure": max(
                    abs(cell["H0"] - (CP * cell["T"] + 0.5 * (cell["U_x"] ** 2 + cell["U_y"] ** 2))) /
                    max(abs(cell["H0"]), 1.0) for cell in cells),
                "relative_l2_total_enthalpy": math.sqrt(sum(
                    (cell["H0"] - inlet_enthalpy) ** 2 for cell in cells) / len(cells)) / inlet_enthalpy,
                "min_mach": min(mach(cell) for cell in cells), "max_mach": max(mach(cell) for cell in cells),
                "exit_upper_mach": mach(profiles["upper"][-1][1]),
                "exit_lower_mach": mach(profiles["lower"][-1][1]),
                "shock_x": shock_location(profiles["lower"]), "passed": 1})
    failures = []
    if not row["converged"]:
        failures.append("final time step is not converged")
    if row["max_cell_imbalance"] > 1e-4 or row["net_boundary_flux"] > 1e-4:
        failures.append("mass conservation exceeds 1e-4")
    if not math.isfinite(row["total_enthalpy_closure"]) or row["total_enthalpy_closure"] > 1e-9:
        failures.append("total-enthalpy/temperature closure is non-finite or exceeds 1e-9")
    if row["relative_l2_total_enthalpy"] > 0.03:
        failures.append("total-enthalpy relative L2 error exceeds 3%")
    for cell in cells:
        if not all(math.isfinite(cell[key]) for key in ("U_x", "U_y", "p", "T", "rho")) or \
                not math.isfinite(mach(cell)) or cell["p"] <= 0 or cell["T"] <= 0 or cell["rho"] <= 0:
            failures.append("non-finite or non-positive state")
            break

    if name == "subsonic":
        outlet_pressure_ratio = 101325.0 / 84606.375
        reference = math.sqrt((outlet_pressure_ratio ** ((1.4 - 1.0) / 1.4) - 1.0) / 0.2)
        row["reference_upper_mach"] = row["reference_lower_mach"] = reference
        row["upper_error"] = abs(row["exit_upper_mach"] - reference)
        row["lower_error"] = abs(row["exit_lower_mach"] - reference)
        lower = [(x, mach(cell)) for x, cell in profiles["lower"]]
        row["symmetry_error"] = max(abs(value - interpolate(lower, 3.0 - x)) for x, value in lower)
        if row["max_mach"] >= 1.0:
            failures.append("subsonic case contains Mach >= 1")
        if max(row["upper_error"], row["lower_error"]) > 0.08:
            failures.append("subsonic exit Mach error exceeds 0.08")
        if row["symmetry_error"] > 0.08:
            failures.append("subsonic fore/aft symmetry error exceeds 0.08")
    elif name == "transonic":
        row["reference_upper_mach"], row["reference_lower_mach"] = 0.670, 0.615
        row["upper_error"] = abs(row["exit_upper_mach"] - 0.670)
        row["lower_error"] = abs(row["exit_lower_mach"] - 0.615)
        if row["max_mach"] <= 1.0:
            failures.append("transonic case has no supersonic region")
        if not 1.4 <= row["shock_x"] <= 2.2:
            failures.append("transonic shock lies outside the reference region")
        if max(row["upper_error"], row["lower_error"]) > 0.10:
            failures.append("transonic exit Mach error exceeds 0.10")
    else:
        row["profile_rmse"] = profile_rmse(profiles, supersonic_reference)
        if row["max_mach"] < 1.4:
            failures.append("supersonic acceleration over the bump is absent")
        if row["min_mach"] < 0.5:
            failures.append("supersonic case has an excessively low Mach number")
        if row["exit_upper_mach"] < 1.0 or row["exit_lower_mach"] < 1.0:
            failures.append("supersonic outlet is not supersonic")
        if row["profile_rmse"] > 0.40:
            failures.append("supersonic wall-profile RMSE exceeds 0.40")
    if failures:
        row["passed"] = 0
    return row, cells, profiles, failures


def comparisons(rows, cells_by_name, profiles_by_name, supersonic_reference):
    for row in rows:
        name = row["case"]
        cells = cells_by_name[name]
        profiles = profiles_by_name[name]
        output_rows = []
        for wall in ("upper", "lower"):
            for x, cell in profiles[wall]:
                output_rows.append({"wall": wall, "x": x, "y": cell["y"], "mach": mach(cell),
                                    "p": cell["p"], "T": cell["T"], "rho": cell["rho"], "H0": cell["H0"],
                                    "U_x": cell["U_x"], "U_y": cell["U_y"]})
        write_rows(DATA / f"bump_{name}_walls.csv", output_rows)
        series = []
        for wall, shade in (("upper", "#2166ac"), ("lower", "#b2182b")):
            series.append({"points": [(x, mach(cell)) for x, cell in profiles[wall]],
                           "label": f"{wall} wall", "color": shade})
            if name == "supersonic":
                series.append({"points": supersonic_reference[wall], "label": f"{wall} reference (digitized)",
                               "color": shade, "dash": "6 3"})
        line_svg(FIGURES / f"bump_{name}_wall_mach.svg", series,
                 f"{name.capitalize()} circular-arc bump", "x", "Mach number")
        contour_svg(FIGURES / f"bump_{name}_mach_contour.svg", cells,
                    f"{name.capitalize()} circular-arc bump Mach field")


def main():
    parser = argparse.ArgumentParser(description="Validate subsonic, transonic, and supersonic bump cases")
    parser.add_argument("--run", action="store_true", help="Generate meshes and run all bump cases")
    args = parser.parse_args()
    if args.run:
        subprocess.run(["python3", str(ROOT / "generate_thermal_meshes.py")], cwd=PROJECT, check=True)
        for specification in CASES:
            print(f"\n=== bump {specification['name']} ===", flush=True)
            subprocess.run([str(PROJECT / "multiphysics"), str(specification["directory"])], cwd=PROJECT, check=True)
    supersonic_reference = {"upper": [], "lower": []}
    for item in read_rows(REFERENCE / "supersonic_mach_approx.csv"):
        supersonic_reference[item["wall"]].append((float(item["x"]), float(item["mach"])))
    failures = []
    rows = []
    cells_by_name = {}
    profiles_by_name = {}
    for specification in CASES:
        try:
            row, cells, profiles, case_failures = analyze(specification, supersonic_reference)
        except (OSError, ValueError, KeyError) as error:
            failures.append(f"{specification['name']}: {error}")
            continue
        rows.append(row)
        cells_by_name[specification["name"]] = cells
        profiles_by_name[specification["name"]] = profiles
        failures.extend(f"{specification['name']}: {failure}" for failure in case_failures)
        print(f"bump_{specification['name']}: iterations={row['total_outer_iterations']}, "
              f"exit=({row['exit_upper_mach']:.4f}, {row['exit_lower_mach']:.4f}), passed={row['passed']}")
    if len(rows) == len(CASES):
        comparisons(rows, cells_by_name, profiles_by_name, supersonic_reference)
    write_rows(OUTPUT / "bump_metrics.csv", rows, FIELDS)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)
    print(f"Bump validation passed: {len(rows)} regimes; results written to {OUTPUT}")


if __name__ == "__main__":
    main()
