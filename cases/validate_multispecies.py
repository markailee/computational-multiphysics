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
R_UNIVERSAL = 8.31446261815324


def spec(name, family, level, directory, stem, h=None, dt=None, steps=None):
    return {"name": name, "family": family, "level": level, "directory": ROOT / directory,
            "stem": stem, "h": h, "dt": dt, "steps": steps}


CASES = [
    spec("stefan_coarse", "stefan", "coarse", "species_stefan_tube/coarse", "stefan_coarse", 0.005),
    spec("stefan_medium", "stefan", "medium", "species_stefan_tube/medium", "stefan_medium", 0.0025),
    spec("stefan_fine", "stefan", "fine", "species_stefan_tube/fine", "stefan_fine", 0.00125),
    spec("correction_coarse", "correction", "coarse", "species_flux_correction/coarse", "correction_coarse", 0.005),
    spec("correction_medium", "correction", "medium", "species_flux_correction/medium", "correction_medium", 0.0025),
    spec("correction_fine", "correction", "fine", "species_flux_correction/fine", "correction_fine", 0.00125),
    spec("diffusion_coarse", "diffusion", "coarse", "species_diffusion/coarse", "diffusion_coarse", 0.002, 0.0002, 25),
    spec("diffusion_medium", "diffusion", "medium", "species_diffusion/medium", "diffusion_medium", 0.001, 0.00005, 100),
    spec("diffusion_fine", "diffusion", "fine", "species_diffusion/fine", "diffusion_fine", 0.0005, 0.0000125, 400),
    spec("air_properties", "air", "", "species_air_properties", "air_properties"),
]

MATERIALS = {
    "stefan": {"names": ("Vapor", "Carrier"), "molar_mass": (0.028, 0.028),
               "viscosity": (1.8e-5, 1.8e-5), "conductivity": (0.025, 0.025),
               "cp": (1000.0, 1000.0), "diffusivity": (2.0e-5, 2.0e-5),
               "temperature": 300.0, "pressure": R_UNIVERSAL * 300.0 / 0.028},
    "correction": {"names": ("A", "B"), "molar_mass": (0.028, 0.028),
                   "viscosity": (1.8e-5, 1.8e-5), "conductivity": (0.025, 0.025),
                   "cp": (1000.0, 1000.0), "diffusivity": (2.0e-5, 5.0e-6),
                   "temperature": 300.0, "pressure": R_UNIVERSAL * 300.0 / 0.028},
    "diffusion": {"names": ("A", "B", "C"), "molar_mass": (0.028, 0.028, 0.028),
                  "viscosity": (1.6e-5, 1.8e-5, 2.1e-5), "conductivity": (0.020, 0.025, 0.030),
                  "cp": (900.0, 1000.0, 1100.0), "diffusivity": (1.0e-3, 1.0e-3, 1.0e-3),
                  "temperature": 300.0, "pressure": R_UNIVERSAL * 300.0 / 0.028},
    "air": {"names": ("N2", "O2"), "molar_mass": (0.0280134, 0.0319988),
            "viscosity": (1.76e-5, 2.04e-5), "conductivity": (0.0259, 0.0266),
            "cp": (1040.0, 918.0), "diffusivity": (2.0e-5, 2.0e-5),
            "temperature": 300.0, "pressure": 101325.0},
}


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def analytical(family, name, x, time):
    if family == "stefan":
        velocity = 2.0e-5 / 0.1 * math.log(0.9 / 0.2)
        vapor = 1.0 - 0.2 * math.exp(velocity * x / 2.0e-5)
        return vapor if name == "Vapor" else 1.0 - vapor
    if family == "correction":
        diffusivity_a = 2.0e-5
        diffusivity_b = 5.0e-6
        y0 = 0.8
        yl = 0.1
        normalization = (diffusivity_a * (yl - y0) +
                         0.5 * (diffusivity_b - diffusivity_a) * (yl * yl - y0 * y0))
        lower, upper = yl, y0
        for _ in range(60):
            y = 0.5 * (lower + upper)
            position = (diffusivity_a * (y - y0) +
                        0.5 * (diffusivity_b - diffusivity_a) * (y * y - y0 * y0)) / normalization
            if position < x / 0.1:
                upper = y
            else:
                lower = y
        a = 0.5 * (lower + upper)
        return a if name == "A" else 1.0 - a
    if family == "diffusion":
        error_function = math.erf(x / (2.0 * math.sqrt(1.0e-3 * time)))
        return {"A": 0.6 - 0.4 * error_function,
                "B": 0.2 + 0.1 * error_function,
                "C": 0.2 + 0.3 * error_function}[name]
    raise ValueError(f"No analytical solution for {family}")


def reference_properties(material, fractions):
    count = len(fractions)
    molar = material["molar_mass"]
    moles = [fractions[k] / molar[k] for k in range(count)]
    mole_sum = sum(moles)
    x = [value / mole_sum for value in moles]
    mean_molar_mass = sum(x[k] * molar[k] for k in range(count))
    density = material["pressure"] * mean_molar_mass / (R_UNIVERSAL * material["temperature"])
    cp = sum(fractions[k] * material["cp"][k] for k in range(count))
    viscosity = 0.0
    for k in range(count):
        denominator = 0.0
        for l in range(count):
            ratio = material["viscosity"][k] / material["viscosity"][l]
            molar_ratio = molar[l] / molar[k]
            phi = ((1.0 + math.sqrt(ratio) * molar_ratio ** 0.25) ** 2 /
                   math.sqrt(8.0 * (1.0 + molar[k] / molar[l])))
            denominator += x[l] * phi
        viscosity += x[k] * material["viscosity"][k] / denominator
    arithmetic = sum(x[k] * material["conductivity"][k] for k in range(count))
    harmonic = sum(x[k] / material["conductivity"][k] for k in range(count))
    conductivity = 0.5 * (arithmetic + 1.0 / harmonic)
    return density, viscosity, conductivity, cp


def property_errors(cells, family):
    material = MATERIALS[family]
    errors = [0.0] * 5
    for cell in cells:
        fractions = [cell[f"Y_{name}"] for name in material["names"]]
        reference = reference_properties(material, fractions)
        numerical = [cell["rho_mix"], cell["mu_mix"], cell["lambda_mix"], cell["cp_mix"]]
        for index in range(4):
            errors[index] = max(errors[index], abs(numerical[index] - reference[index]) / abs(reference[index]))
        for k, name in enumerate(material["names"]):
            exact = material["diffusivity"][k]
            errors[4] = max(errors[4], abs(cell[f"D_{name}"] - exact) / exact)
    return errors


def result_path(specification):
    suffix = "" if specification["steps"] is None else f"_{specification['steps']}"
    return specification["directory"] / "results" / f"{specification['stem']}{suffix}.plt"


def profile(cells, name):
    groups = {}
    for cell in cells:
        groups.setdefault(round(cell["x"], 14), []).append(cell[name])
    return sorted((x, sum(values) / len(values)) for x, values in groups.items())


def analyze(specification):
    family = specification["family"]
    summary = read_rows(specification["directory"] / "results" / f"{specification['stem']}_summary.csv")
    cells = read_tecplot(result_path(specification))
    time = float(summary[-1]["time"])
    material = MATERIALS[family]
    numerator = 0.0
    denominator = 0.0
    linf = 0.0
    if family != "air":
        for cell in cells:
            for name in material["names"]:
                exact = analytical(family, name, cell["x"], time)
                error = cell[f"Y_{name}"] - exact
                numerator += error * error
                denominator += exact * exact
                linf = max(linf, abs(error))
        relative_l2 = math.sqrt(numerator / denominator)
    else:
        relative_l2 = ""
        linf = ""
    property_error = property_errors(cells, family)
    row = {
        "case": specification["name"], "family": family, "level": specification["level"],
        "h": specification["h"] if specification["h"] is not None else "",
        "cells": len(cells), "time": time,
        "dt": specification["dt"] if specification["dt"] is not None else "",
        "steps": specification["steps"] if specification["steps"] is not None else "",
        "iterations": sum(int(item["iterations"]) for item in summary),
        "relative_l2": relative_l2, "linf": linf, "observed_order": "",
        "composition_error": max(float(item["species_composition_error"]) for item in summary),
        "min_mass_fraction": min(float(item["minimum_mass_fraction"]) for item in summary),
        "max_mass_fraction": max(float(item["maximum_mass_fraction"]) for item in summary),
        "max_cell_imbalance": max(float(item["max_cell_imbalance"]) for item in summary),
        "net_boundary_flux": max(abs(float(item["net_boundary_flux"])) for item in summary),
        "density_error": property_error[0], "viscosity_error": property_error[1],
        "conductivity_error": property_error[2], "cp_error": property_error[3],
        "diffusivity_error": property_error[4],
        "converged": int(all(int(item["converged"]) == 1 for item in summary)), "passed": 1,
    }
    failures = []
    if not row["converged"]:
        failures.append("solver convergence failed")
    if row["composition_error"] > 1e-8:
        failures.append("mass-fraction closure exceeds 1e-8")
    if row["min_mass_fraction"] < -1e-9 or row["max_mass_fraction"] > 1.0 + 1e-9:
        failures.append("mass-fraction bounds exceed 1e-9")
    if row["max_cell_imbalance"] > 1e-8 or row["net_boundary_flux"] > 1e-8:
        failures.append("continuity conservation failed")
    if max(property_error) > 1e-8:
        failures.append("mixture-property relative error exceeds 1e-8")
    if family in ("stefan", "correction") and relative_l2 > 0.02:
        failures.append(f"{family} profile relative L2 error exceeds 2%")
    if family == "diffusion" and relative_l2 > 0.03:
        failures.append("transient-diffusion relative L2 error exceeds 3%")
    if failures:
        row["passed"] = 0
    return row, cells, failures


def convergence(rows, failures):
    for family in ("stefan", "correction", "diffusion"):
        group = sorted((row for row in rows if row["family"] == family),
                       key=lambda row: row["h"], reverse=True)
        for coarse, fine in zip(group, group[1:]):
            if fine["relative_l2"] > 0.0:
                fine["observed_order"] = math.log(coarse["relative_l2"] / fine["relative_l2"]) / math.log(2.0)
            if fine["relative_l2"] >= coarse["relative_l2"] / 1.5:
                fine["passed"] = 0
                failures.append(f"{family}: insufficient error reduction from {coarse['level']} to {fine['level']}")
        line_svg(FIGURES / f"multispecies_{family}_convergence.svg", [
            {"points": [(row["h"], row["relative_l2"]) for row in reversed(group)],
             "label": "Relative L2 error", "color": "#2166ac"},
        ], f"{family.replace('_', ' ')} grid convergence", "mesh spacing h", "relative L2 error")


def comparisons(rows_by_name, cells_by_name):
    for family, case_name in (("stefan", "stefan_fine"), ("correction", "correction_fine"),
                              ("diffusion", "diffusion_fine")):
        cells = cells_by_name[case_name]
        time = rows_by_name[case_name]["time"]
        rows = []
        for cell in cells:
            for name in MATERIALS[family]["names"]:
                exact = analytical(family, name, cell["x"], time)
                rows.append({"x": cell["x"], "y": cell["y"], "species": name,
                             "numerical": cell[f"Y_{name}"], "analytical": exact,
                             "error": cell[f"Y_{name}"] - exact})
        write_rows(DATA / f"multispecies_{family}_fine.csv",
                   ("x", "y", "species", "numerical", "analytical", "error"), rows)
        colors = ("#2166ac", "#1b9e77", "#b2182b")
        series = []
        for index, name in enumerate(MATERIALS[family]["names"]):
            numerical = profile(cells, f"Y_{name}")
            series.append({"points": numerical, "label": f"{name} numerical", "color": colors[index]})
            exact = [(x, analytical(family, name, x, time)) for x, _ in numerical]
            series.append({"points": exact, "label": f"{name} analytical", "color": colors[index], "dash": "6 3"})
        line_svg(FIGURES / f"multispecies_{family}_profiles.svg", series,
                 f"{family.replace('_', ' ')} species profiles", "x (m)", "mass fraction", legend_left=True)


def main():
    parser = argparse.ArgumentParser(description="Run and analyze multispecies verification cases")
    parser.add_argument("--run", action="store_true", help="Generate meshes and run every case")
    args = parser.parse_args()
    if args.run:
        subprocess.run(["python3", str(ROOT / "generate_multispecies_meshes.py")], cwd=PROJECT, check=True)
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
        error_text = f"{row['relative_l2']:.6e}" if row["relative_l2"] != "" else "property-only"
        print(f"{specification['name']}: L2={error_text}, Ysum={row['composition_error']:.3e}, "
              f"property max={max(row[key] for key in ('density_error', 'viscosity_error', 'conductivity_error', 'cp_error', 'diffusivity_error')):.3e}")
    if len(rows) == len(CASES):
        convergence(rows, failures)
        comparisons({row["case"]: row for row in rows}, cells_by_name)
    fieldnames = ("case", "family", "level", "h", "cells", "time", "dt", "steps", "iterations",
                  "relative_l2", "linf", "observed_order", "composition_error", "min_mass_fraction",
                  "max_mass_fraction", "max_cell_imbalance", "net_boundary_flux", "density_error",
                  "viscosity_error", "conductivity_error", "cp_error", "diffusivity_error", "converged", "passed")
    write_rows(OUTPUT / "multispecies_metrics.csv", fieldnames, rows)
    write_rows(DATA / "multispecies_convergence.csv",
               ("case", "family", "level", "h", "dt", "relative_l2", "observed_order"),
               [{key: row[key] for key in ("case", "family", "level", "h", "dt", "relative_l2", "observed_order")}
                for row in rows if row["level"]])
    if failures:
        for failure in failures:
            print(f"FAIL: {failure}")
        raise SystemExit(1)
    print(f"Multispecies validation passed: {len(rows)} runs; metrics and figures written to {OUTPUT}")


if __name__ == "__main__":
    main()
