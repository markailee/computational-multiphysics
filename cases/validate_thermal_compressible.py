#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
from pathlib import Path

from validate import read_tecplot
from validate_acoustic import line_svg

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
OUTPUT = PROJECT / "validation"
DATA = OUTPUT / "data"
FIGURES = OUTPUT / "figures"
R_GAS = 287.0
CP = 1005.0
P_IN = 100000.0
T_IN = 300.0
U_IN = 100.0
NU_REF = 2.243
U_REF = 16.178
V_REF = 19.617

CASES = [
    {"name": "boussinesq_coarse", "family": "cavity", "level": "coarse", "n": 16,
     "directory": ROOT / "boussinesq_cavity/coarse", "stem": "boussinesq_coarse"},
    {"name": "boussinesq_medium", "family": "cavity", "level": "medium", "n": 24,
     "directory": ROOT / "boussinesq_cavity/medium", "stem": "boussinesq_medium"},
    {"name": "boussinesq_fine", "family": "cavity", "level": "fine", "n": 32,
     "directory": ROOT / "boussinesq_cavity/fine", "stem": "boussinesq_fine"},
    {"name": "ideal_gas_conduction", "family": "conduction", "level": "", "n": 40, "h": 0.0025,
     "directory": ROOT / "compressible_conduction", "stem": "conduction"},
    {"name": "closed_heating", "family": "heating", "level": "", "n": 40, "h": 0.0025, "steps": 10,
     "directory": ROOT / "compressible_closed_heating", "stem": "closed_heating"},
    {"name": "nozzle_coarse", "family": "nozzle", "level": "coarse", "n": 40, "h": 0.1,
     "directory": ROOT / "compressible_nozzle/coarse", "stem": "nozzle_coarse"},
    {"name": "nozzle_medium", "family": "nozzle", "level": "medium", "n": 60, "h": 4.0 / 60.0,
     "directory": ROOT / "compressible_nozzle/medium", "stem": "nozzle_medium"},
    {"name": "nozzle_fine", "family": "nozzle", "level": "fine", "n": 80, "h": 0.05,
     "directory": ROOT / "compressible_nozzle/fine", "stem": "nozzle_fine"},
]

FIELDS = ("case", "family", "level", "h", "cells", "iterations", "converged", "max_cell_imbalance",
          "net_boundary_flux", "temperature_change", "density_change", "total_enthalpy_closure",
          "relative_l2_total_enthalpy", "relative_l2_temperature", "relative_l2_density",
          "relative_l2_pressure", "relative_l2_velocity", "observed_order",
          "nusselt", "nusselt_error_percent", "u_max", "u_max_error_percent", "v_max",
          "v_max_error_percent", "max_transverse_velocity", "passed")


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, rows, fieldnames=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames or tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def relative_l2(cells, name, exact, weight=None):
    numerator = 0.0
    denominator = 0.0
    for cell in cells:
        w = weight(cell) if weight else 1.0
        reference = exact(cell)
        numerator += w * (cell[name] - reference) ** 2
        denominator += w * reference ** 2
    return math.sqrt(numerator / denominator) if denominator > 0 else 0.0


def centerline(cells, n, vertical, name):
    groups = {}
    for cell in cells:
        if vertical and abs(cell["x"] - 0.5) <= 0.5 / n + 1e-10:
            groups.setdefault(round(cell["y"], 12), []).append(cell[name])
        if not vertical and abs(cell["y"] - 0.5) <= 0.5 / n + 1e-10:
            groups.setdefault(round(cell["x"], 12), []).append(cell[name])
    return sorted((position, sum(values) / len(values)) for position, values in groups.items())


def nozzle_exact(x):
    area = 0.2 - 0.025 * x
    lower, upper = U_IN, 500.0
    exponent = (CP - R_GAS) / R_GAS
    for _ in range(80):
        velocity = 0.5 * (lower + upper)
        temperature = T_IN + (U_IN ** 2 - velocity ** 2) / (2 * CP)
        residual = velocity * (temperature / T_IN) ** exponent - (0.2 / area) * U_IN
        if residual > 0:
            upper = velocity
        else:
            lower = velocity
    velocity = 0.5 * (lower + upper)
    temperature = T_IN + (U_IN ** 2 - velocity ** 2) / (2 * CP)
    pressure = P_IN * (temperature / T_IN) ** (CP / R_GAS)
    return {"U_x": velocity, "T": temperature, "p": pressure, "rho": pressure / (R_GAS * temperature)}


def profile(cells, name):
    groups = {}
    for cell in cells:
        groups.setdefault(round(cell["x"], 12), []).append(cell[name])
    return sorted((x, sum(values) / len(values)) for x, values in groups.items())


def analyze(specification):
    directory = specification["directory"]
    stem = specification["stem"]
    summaries = read_rows(directory / "results" / f"{stem}_summary.csv")
    summary = summaries[-1]
    suffix = f"_{specification['steps']}" if "steps" in specification else ""
    cells = read_tecplot(directory / "results" / f"{stem}{suffix}.plt")
    family = specification["family"]
    row = {key: "" for key in FIELDS}
    row.update({"case": specification["name"], "family": family, "level": specification["level"],
                "h": specification.get("h", 1.0 / specification["n"]), "cells": len(cells),
                "iterations": sum(int(item["iterations"]) for item in summaries),
                "converged": int(all(int(item["converged"]) == 1 for item in summaries)),
                "max_cell_imbalance": max(float(item["max_cell_imbalance"]) for item in summaries),
                "net_boundary_flux": max(abs(float(item["net_boundary_flux"])) for item in summaries),
                "temperature_change": float(summary["temperature_change"]),
                "density_change": float(summary["density_change"]),
                "max_transverse_velocity": max(abs(cell["U_y"]) for cell in cells), "passed": 1})
    failures = []
    if family != "cavity":
        row["total_enthalpy_closure"] = max(
            abs(cell["H0"] - (CP * cell["T"] + 0.5 * (cell["U_x"] ** 2 + cell["U_y"] ** 2))) /
            max(abs(cell["H0"]), 1.0) for cell in cells)
        if not math.isfinite(row["total_enthalpy_closure"]) or row["total_enthalpy_closure"] > 1e-9:
            failures.append("total-enthalpy/temperature closure is non-finite or exceeds 1e-9")
    if not row["converged"]:
        failures.append("solver convergence failed")
    if row["max_cell_imbalance"] > 1e-7 or row["net_boundary_flux"] > 1e-7:
        failures.append("mass conservation exceeds 1e-7")
    for cell in cells:
        if not all(math.isfinite(cell[key]) for key in ("U_x", "U_y", "p", "T", "rho")):
            failures.append("non-finite solution value")
            break
        if cell["T"] <= 0 or cell["rho"] <= 0 or cell["p"] <= 0 and family != "cavity":
            failures.append("non-positive compressible state")
            break

    if family == "cavity":
        n = specification["n"]
        hot_cells = [cell for cell in cells if abs(cell["x"] - 0.5 / n) < 1e-10]
        row["nusselt"] = sum((301.0 - cell["T"]) / cell["x"] for cell in hot_cells) / len(hot_cells)
        row["nusselt_error_percent"] = 100 * abs(row["nusselt"] - NU_REF) / NU_REF
        u_profile = centerline(cells, n, True, "U_x")
        v_profile = centerline(cells, n, False, "U_y")
        row["u_max"] = max(value for _, value in u_profile)
        row["v_max"] = max(value for _, value in v_profile)
        row["u_max_error_percent"] = 100 * abs(row["u_max"] - U_REF) / U_REF
        row["v_max_error_percent"] = 100 * abs(row["v_max"] - V_REF) / V_REF
        if row["nusselt_error_percent"] > 8:
            failures.append("hot-wall Nusselt error exceeds 8%")
        if row["u_max_error_percent"] > 15 or row["v_max_error_percent"] > 15:
            failures.append("mid-plane velocity maximum error exceeds 15%")
        row["_profiles"] = {"u": u_profile, "v": v_profile}
    elif family == "conduction":
        temperature = lambda cell: 300.0 + 1000.0 * cell["x"]
        density = lambda cell: 101325.0 / (R_GAS * temperature(cell))
        row["relative_l2_temperature"] = relative_l2(cells, "T", temperature)
        row["relative_l2_density"] = relative_l2(cells, "rho", density)
        row["relative_l2_pressure"] = relative_l2(cells, "p", lambda cell: 101325.0)
        row["relative_l2_velocity"] = max(math.hypot(cell["U_x"], cell["U_y"]) for cell in cells)
        if max(row[key] for key in ("relative_l2_temperature", "relative_l2_density",
                                    "relative_l2_pressure", "relative_l2_velocity")) > 1e-9:
            failures.append("zero-flow conduction/EOS error exceeds 1e-9")
    elif family == "heating":
        initial_density = 101325.0 / (R_GAS * 300.0)
        time = float(summary["time"])
        exact_temperature = 300.0 + 1000.0 * time / (initial_density * (CP - R_GAS))
        exact_pressure = initial_density * R_GAS * exact_temperature
        row["relative_l2_temperature"] = relative_l2(cells, "T", lambda cell: exact_temperature)
        row["relative_l2_density"] = relative_l2(cells, "rho", lambda cell: initial_density)
        row["relative_l2_pressure"] = relative_l2(cells, "p", lambda cell: exact_pressure)
        row["relative_l2_velocity"] = max(math.hypot(cell["U_x"], cell["U_y"]) for cell in cells)
        if max(row[key] for key in ("relative_l2_temperature", "relative_l2_density",
                                    "relative_l2_pressure", "relative_l2_velocity")) > 1e-8:
            failures.append("closed-heating analytical error exceeds 1e-8")
    else:
        inlet_enthalpy = CP * T_IN + 0.5 * U_IN ** 2
        row["relative_l2_total_enthalpy"] = relative_l2(
            cells, "H0", lambda cell: inlet_enthalpy, lambda cell: 0.2 - 0.025 * cell["x"])
        for name, output in (("T", "relative_l2_temperature"), ("rho", "relative_l2_density"),
                             ("p", "relative_l2_pressure"), ("U_x", "relative_l2_velocity")):
            row[output] = relative_l2(cells, name, lambda cell, key=name: nozzle_exact(cell["x"])[key],
                                      lambda cell: 0.2 - 0.025 * cell["x"])
        if row["relative_l2_total_enthalpy"] > 0.01:
            failures.append("nozzle total-enthalpy relative L2 error exceeds 1%")
        if row["relative_l2_pressure"] > 0.05 or row["relative_l2_density"] > 0.05:
            failures.append("nozzle pressure/density relative L2 error exceeds 5%")
        if row["relative_l2_temperature"] > 0.01 or row["relative_l2_velocity"] > 0.05:
            failures.append("nozzle temperature/velocity relative L2 error exceeds its tolerance")
    if failures:
        row["passed"] = 0
    return row, cells, failures


def grid_checks(rows, failures):
    for family in ("cavity", "nozzle"):
        group = sorted((row for row in rows if row["family"] == family), key=lambda row: row["h"], reverse=True)
        if len(group) != 3:
            continue
        if family == "cavity":
            if not group[-1]["nusselt_error_percent"] < group[0]["nusselt_error_percent"]:
                group[-1]["passed"] = 0
                failures.append("cavity Nusselt error did not decrease with refinement")
            if group[-1]["nusselt_error_percent"] > 5:
                group[-1]["passed"] = 0
                failures.append("fine-cavity Nusselt error exceeds 5%")
        else:
            for coarse, fine in zip(group, group[1:]):
                if fine["relative_l2_pressure"] > 0:
                    fine["observed_order"] = math.log(coarse["relative_l2_pressure"] / fine["relative_l2_pressure"]) / math.log(coarse["h"] / fine["h"])
            if not group[-1]["relative_l2_pressure"] < group[0]["relative_l2_pressure"]:
                group[-1]["passed"] = 0
                failures.append("nozzle pressure error did not decrease with refinement")


def comparisons(rows, cells_by_name):
    for row in rows:
        if row["family"] != "cavity":
            continue
        profiles = row.pop("_profiles")
        output_rows = []
        for component, values in profiles.items():
            for position, value in values:
                output_rows.append({"level": row["level"], "component": component,
                                    "position": position, "numerical": value})
        write_rows(DATA / f"boussinesq_{row['level']}_profiles.csv", output_rows)
    for specification in CASES:
        if specification["family"] == "nozzle":
            cells = cells_by_name[specification["name"]]
            output_rows = []
            for cell in cells:
                exact = nozzle_exact(cell["x"])
                output_rows.append({"x": cell["x"], "y": cell["y"],
                                    "p_numerical": cell["p"], "p_exact": exact["p"],
                                    "T_numerical": cell["T"], "T_exact": exact["T"],
                                    "rho_numerical": cell["rho"], "rho_exact": exact["rho"],
                                    "U_numerical": cell["U_x"], "U_exact": exact["U_x"],
                                    "H0_numerical": cell["H0"], "H0_exact": CP * T_IN + 0.5 * U_IN ** 2})
            write_rows(DATA / f"compressible_{specification['name']}.csv", output_rows)
    cells = cells_by_name["ideal_gas_conduction"]
    write_rows(DATA / "compressible_conduction.csv",
               [{"x": cell["x"], "T_numerical": cell["T"], "T_exact": 300 + 1000 * cell["x"],
                 "rho_numerical": cell["rho"], "rho_exact": 101325 / (R_GAS * (300 + 1000 * cell["x"]))}
                for cell in cells])

    cavity = sorted((row for row in rows if row["family"] == "cavity"), key=lambda row: row["h"])
    line_svg(FIGURES / "boussinesq_nusselt_grid.svg", [
        {"points": [(row["h"], row["nusselt"]) for row in cavity], "label": "Present", "color": "#2166ac"},
        {"points": [(row["h"], NU_REF) for row in cavity], "label": "de Vahl Davis", "color": "#b2182b", "dash": "6 3"},
    ], "Boussinesq cavity grid refinement", "mesh spacing h", "mean hot-wall Nusselt number")
    cells = cells_by_name["nozzle_fine"]
    numerical = {name: profile(cells, name) for name in ("p", "T", "rho", "U_x")}
    exact = {name: [(x, nozzle_exact(x)[name]) for x, _ in numerical["p"]] for name in numerical}
    normalizers = {"p": P_IN, "T": T_IN, "rho": P_IN / (R_GAS * T_IN), "U_x": U_IN}
    colors = {"p": "#2166ac", "T": "#b2182b", "rho": "#1b9e77", "U_x": "#7570b3"}
    series = []
    for name in ("p", "T", "rho", "U_x"):
        series.append({"points": [(x, value / normalizers[name]) for x, value in numerical[name]],
                       "label": f"{name} numerical", "color": colors[name]})
        series.append({"points": [(x, value / normalizers[name]) for x, value in exact[name]],
                       "label": f"{name} quasi-1D", "color": colors[name], "dash": "6 3"})
    line_svg(FIGURES / "compressible_nozzle_profiles.svg", series,
             "Subsonic ideal-gas nozzle", "x (m)", "normalized state", legend_left=True)


def main():
    parser = argparse.ArgumentParser(description="Validate Boussinesq and compressible SIMPLE models")
    parser.add_argument("--run", action="store_true", help="Generate meshes and run all cases")
    args = parser.parse_args()
    if args.run:
        subprocess.run(["python3", str(ROOT / "generate_thermal_meshes.py")], cwd=PROJECT, check=True)
        for specification in CASES:
            print(f"\n=== {specification['name']} ===", flush=True)
            subprocess.run([str(PROJECT / "multiphysics"), str(specification["directory"])], cwd=PROJECT, check=True)
    failures = []
    rows = []
    cells_by_name = {}
    for specification in CASES:
        try:
            row, cells, case_failures = analyze(specification)
        except (OSError, ValueError, KeyError) as error:
            failures.append(f"{specification['name']}: {error}")
            continue
        rows.append(row)
        cells_by_name[specification["name"]] = cells
        failures.extend(f"{specification['name']}: {failure}" for failure in case_failures)
        print(f"{specification['name']}: iterations={row['iterations']}, passed={row['passed']}")
    if len(rows) == len(CASES):
        grid_checks(rows, failures)
        comparisons(rows, cells_by_name)
    for row in rows:
        row.pop("_profiles", None)
    write_rows(OUTPUT / "thermal_compressible_metrics.csv", rows, FIELDS)
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)
    print(f"Thermal/compressible validation passed: {len(rows)} runs; results written to {OUTPUT}")


if __name__ == "__main__":
    main()
