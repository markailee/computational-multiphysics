#!/usr/bin/env python3
import csv
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def read_tecplot(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    names = []
    cells = []
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
        nodal = {}
        for name in names[:2]:
            nodal[name] = [float(value) for value in lines[index:index + node_count]]
            index += node_count
        values = {}
        for name in names[2:]:
            values[name] = [float(value) for value in lines[index:index + element_count]]
            index += element_count
        for element in range(element_count):
            connectivity = [int(value) - 1 for value in lines[index].split()]
            index += 1
            if len(connectivity) != element_size:
                raise ValueError(f"Invalid connectivity in {path}")
            cell = {name: values[name][element] for name in names[2:]}
            cell["x"] = sum(nodal["X"][node] for node in connectivity) / element_size
            cell["y"] = sum(nodal["Y"][node] for node in connectivity) / element_size
            cells.append(cell)
    if not cells:
        raise ValueError(f"No cells found in {path}")
    return cells


def summary(case):
    path = ROOT / case / "results" / f"{case}_summary.csv"
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))[-1]


def interpolate(cells, x, y, count=8):
    nearest = sorted(cells, key=lambda cell: (cell["x"] - x) ** 2 + (cell["y"] - y) ** 2)[:count]
    weights = [1.0 / max((cell["x"] - x) ** 2 + (cell["y"] - y) ** 2, 1e-12) for cell in nearest]
    return sum(weight * cell["U_x"] for weight, cell in zip(weights, nearest)) / sum(weights)


def reattachment(points):
    points = sorted(points, key=lambda cell: cell["x"])
    for first, second in zip(points, points[1:]):
        if first["U_x"] < 0 <= second["U_x"]:
            fraction = -first["U_x"] / (second["U_x"] - first["U_x"])
            return first["x"] + fraction * (second["x"] - first["x"])
    return None


def validate(case):
    cells = read_tecplot(ROOT / case / "results" / f"{case}.plt")
    result = summary(case)
    failures = []
    if int(result["converged"]) != 1:
        failures.append("SIMPLE did not converge")
    imbalance = float(result["max_cell_imbalance"])
    net_flux = abs(float(result["net_boundary_flux"]))
    if imbalance > 1e-8:
        failures.append(f"cell mass imbalance {imbalance:.3e} exceeds 1e-8")
    if net_flux > 1e-8:
        failures.append(f"net boundary mass flux {net_flux:.3e} exceeds 1e-8")
    for cell in cells:
        if not all(math.isfinite(cell[name]) for name in ("U_x", "U_y", "p")):
            failures.append("non-finite solution value")
            break

    details = ""
    if case == "cavity":
        velocity = interpolate(cells, 0.5, 0.5)
        reference = -0.20581
        error = abs(velocity - reference)
        details = f"u(0.5,0.5)={velocity:.5f}, Ghia Re=100 reference={reference:.5f}, |error|={error:.5f}"
        if error > 0.08:
            failures.append("cavity centreline velocity differs excessively from the reference")
    elif case == "backstep":
        lower_wall = [cell for cell in cells if cell["x"] > 0 and cell["y"] < 0.04]
        attachment = reattachment(lower_wall)
        minimum = min(cell["U_x"] for cell in lower_wall)
        details = f"minimum near-wall u={minimum:.5f}, estimated reattachment x={attachment}"
        if minimum >= 0 or attachment is None or not 0.1 < attachment < 8.0:
            failures.append("expected primary backstep recirculation was not found")
    elif case == "cylinder":
        wake = [cell for cell in cells if 0.5 < cell["x"] < 5.0 and abs(cell["y"]) < 0.15]
        minimum = min(cell["U_x"] for cell in wake)
        attachment = reattachment([cell for cell in wake if cell["x"] >= 0.5])
        details = f"minimum centreline wake u={minimum:.5f}, estimated wake reattachment x={attachment}"
        if minimum >= 0 or attachment is None:
            failures.append("expected steady cylinder wake recirculation was not found")
    print(f"{case}: iterations={result['iterations']}, max imbalance={imbalance:.3e}, net flux={net_flux:.3e}")
    print(f"  {details}")
    return failures


def main():
    failures = []
    for case in ("cavity", "backstep", "cylinder"):
        try:
            failures.extend(f"{case}: {failure}" for failure in validate(case))
        except (OSError, ValueError, KeyError) as error:
            failures.append(f"{case}: {error}")
    if failures:
        print("Validation failed:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1
    print("All conservation and benchmark checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
