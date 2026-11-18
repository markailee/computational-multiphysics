#!/usr/bin/env python3
import argparse
import bisect
import csv
import math
import subprocess
import sys
from pathlib import Path
from xml.sax.saxutils import escape

from validate import read_tecplot

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
OUTPUT = PROJECT / "validation"
DATA = OUTPUT / "data"
FIGURES = OUTPUT / "figures"
LEVELS = ("coarse", "medium", "fine")
MODELS = (("k_epsilon", "epsilon"), ("sst_k_omega", "omega"))


def cases():
    result = []
    for geometry in ("plate", "elbow"):
        for model, scale in MODELS:
            directory = ROOT / f"{geometry}_{model}"
            stem = f"{geometry}_{model}"
            for level in LEVELS:
                input_name = "case.in" if level == "medium" else f"{level}.in"
                result_name = f"{stem}.plt" if level == "medium" else f"{stem}_{level}.plt"
                result.append({
                    "geometry": geometry,
                    "model": model,
                    "scale": scale,
                    "level": level,
                    "name": f"{stem}_{level}",
                    "input": directory / input_name,
                    "result": directory / "results" / result_name,
                })
    return result


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_reference(path):
    points = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if not line or line.startswith("#") or line.lower().startswith(("variables", "zone")):
                continue
            parts = line.replace(",", " ").split()
            try:
                points.append((float(parts[0]), float(parts[1])))
            except (ValueError, IndexError):
                pass
    return sorted(points)


def interpolate(points, x, logarithmic=False):
    xs = [math.log(p[0]) if logarithmic else p[0] for p in points]
    ys = [p[1] for p in points]
    target = math.log(x) if logarithmic else x
    index = bisect.bisect_left(xs, target)
    if index == 0:
        return ys[0]
    if index >= len(points):
        return ys[-1]
    fraction = (target - xs[index - 1]) / (xs[index] - xs[index - 1])
    return ys[index - 1] + fraction * (ys[index] - ys[index - 1])


def error_metrics(pairs):
    if not pairs:
        return {"rmse_relative": math.nan, "max_relative": math.nan}
    rmse = math.sqrt(sum((numeric - reference) ** 2 for numeric, reference in pairs) / len(pairs))
    scale = math.sqrt(sum(reference ** 2 for _, reference in pairs) / len(pairs))
    nonzero = [(numeric, reference) for numeric, reference in pairs if abs(reference) > 1e-14]
    return {
        "rmse_relative": rmse / scale,
        "max_relative": max(abs(numeric - reference) / abs(reference)
                            for numeric, reference in nonzero) if nonzero else math.nan,
    }


def trapezoid(points, function):
    return sum((function(second) + function(first)) * 0.5 * (second[0] - first[0])
               for first, second in zip(points, points[1:]))


def plate_profiles(cells):
    groups = {}
    for cell in cells:
        groups.setdefault(round(cell["x"], 8), []).append(cell)
    return [(x, sorted(group, key=lambda cell: cell["y"])) for x, group in sorted(groups.items())]


def momentum_theta(profile, viscosity):
    freestream = profile[-1]["U_x"]
    points = [(0.0, 0.0)] + [(cell["y"], cell["U_x"] / freestream) for cell in profile] + [(1.0, 1.0)]
    theta = trapezoid(points, lambda point: point[1] * (1.0 - point[1]))
    return freestream * theta / viscosity


def nearest_profile(profiles, x):
    return min(profiles, key=lambda item: abs(item[0] - x))[1]


def sample_cells(cells, x, y):
    nearest = sorted(cells, key=lambda cell: (cell["x"] - x) ** 2 + (cell["y"] - y) ** 2)[:8]
    weights = [1.0 / max((cell["x"] - x) ** 2 + (cell["y"] - y) ** 2, 1e-16) for cell in nearest]
    total = sum(weights)
    return {name: sum(weight * cell[name] for weight, cell in zip(weights, nearest)) / total
            for name in ("U_x", "U_y", "p")}


def analyze_plate(spec, cells, wall):
    viscosity = 2.0e-7
    profiles = plate_profiles(cells)
    retheta = [(x, momentum_theta(profile, viscosity)) for x, profile in profiles]
    cf_reference = read_reference(PROJECT / "reference" / "flat_plate" / "cf_karman_schoenherr.dat")
    coles = read_reference(PROJECT / "reference" / "flat_plate" / "uplus_coles_retheta_10000.dat")
    plate_wall = sorted((row for row in wall if row["patch"] == "plate"), key=lambda row: float(row["x"]))
    cf_rows = []
    cf_pairs = []
    for row in plate_wall:
        x = float(row["x"])
        theta = interpolate(retheta, x)
        cf = float(row["tau_x"]) / 0.5
        reference = interpolate(cf_reference, theta)
        output = {
            "x": x,
            "retheta": theta,
            "cf": cf,
            "cf_reference": reference,
            "relative_error": (cf - reference) / reference,
            "y_plus": float(row["y_plus"]),
        }
        cf_rows.append(output)
        if 4000.0 <= theta <= 13000.0 and cf > 0.0:
            cf_pairs.append((cf, reference))
    metrics = error_metrics(cf_pairs)
    result = {
        "cf_rmse_relative": metrics["rmse_relative"],
        "cf_max_relative": metrics["max_relative"],
    }

    target_x = None
    for first, second in zip(retheta, retheta[1:]):
        if first[1] <= 10000.0 <= second[1]:
            fraction = (10000.0 - first[1]) / (second[1] - first[1])
            target_x = first[0] + fraction * (second[0] - first[0])
            break
    uplus_rows = []
    uplus_pairs = []
    if target_x is not None:
        profile = nearest_profile(profiles, target_x)
        friction = interpolate([(float(row["x"]), float(row["friction_velocity"])) for row in plate_wall], target_x)
        for cell in profile:
            y_plus = friction * cell["y"] / viscosity
            if y_plus < coles[0][0] or y_plus > coles[-1][0]:
                continue
            u_plus = cell["U_x"] / friction
            reference = interpolate(coles, y_plus, logarithmic=True)
            uplus_rows.append({"x": target_x, "y": cell["y"], "y_plus": y_plus,
                               "u_plus": u_plus, "u_plus_reference": reference})
            uplus_pairs.append((u_plus, reference))
    metrics = error_metrics(uplus_pairs)
    result["uplus_rmse_relative"] = metrics["rmse_relative"]
    result["uplus_max_relative"] = metrics["max_relative"]
    result["retheta_10000_x"] = target_x if target_x is not None else math.nan

    widths = []
    for index, row in enumerate(plate_wall):
        if len(plate_wall) == 1:
            widths.append(0.0)
        elif index == 0:
            widths.append(float(plate_wall[1]["x"]) - float(row["x"]))
        elif index == len(plate_wall) - 1:
            widths.append(float(row["x"]) - float(plate_wall[index - 1]["x"]))
        else:
            widths.append((float(plate_wall[index + 1]["x"]) - float(plate_wall[index - 1]["x"])) / 2)
    drag = sum(float(row["tau_x"]) * width for row, width in zip(plate_wall, widths))
    result["drag_coefficient"] = drag / 2.0
    yplus = [float(row["y_plus"]) for row in plate_wall if float(row["x"]) > 0.1]
    result["yplus_min"] = min(yplus)
    result["yplus_max"] = max(yplus)
    result["cap_cell_fraction"] = sum(cell["mu_t"] >= 1000.0 * viscosity * (1.0 - 1e-9) for cell in cells) / len(cells)

    write_rows(DATA / f"{spec['name']}_cf.csv",
               ("x", "retheta", "cf", "cf_reference", "relative_error", "y_plus"), cf_rows)
    write_rows(DATA / f"{spec['name']}_uplus.csv",
               ("x", "y", "y_plus", "u_plus", "u_plus_reference"), uplus_rows)
    return result, cf_rows, uplus_rows


def analyze_elbow(spec, cells, wall):
    reference = read_rows(PROJECT / "reference" / "elbow_enayet_1982" / "velocity_75deg_experimental.csv")
    angle = math.radians(75.0)
    tangent = (math.cos(angle), math.sin(angle))
    rows = []
    pairs = []
    for point in reference:
        eta = float(point["y_over_d"])
        radius = 2.3 + eta
        x = radius * math.sin(angle)
        y = 2.8 - radius * math.cos(angle)
        if eta <= 0.0 or eta >= 1.0:
            velocity = {"U_x": 0.0, "U_y": 0.0}
        else:
            velocity = sample_cells(cells, x, y)
        magnitude = math.hypot(velocity["U_x"], velocity["U_y"])
        streamwise = velocity["U_x"] * tangent[0] + velocity["U_y"] * tangent[1]
        expected = float(point["velocity_over_bulk"])
        rows.append({"y_over_d": eta, "x": x, "y": y, "velocity_magnitude": magnitude,
                     "streamwise_velocity": streamwise, "reference": expected,
                     "reference_uncertainty": float(point["velocity_uncertainty"]),
                     "error": magnitude - expected})
        pairs.append((magnitude, expected))
    metrics = error_metrics(pairs)
    result = {
        "elbow_rmse_relative": metrics["rmse_relative"],
        "elbow_max_relative": metrics["max_relative"],
    }

    inlet = []
    outlet = []
    for index in range(41):
        eta = index / 40
        inlet.append(sample_cells(cells, -4.5, -0.5 + eta)["p"])
        outlet.append(sample_cells(cells, 2.8, 2.8 + 9.5 + eta)["p"])
    pressure_drop = sum(inlet) / len(inlet) - sum(outlet) / len(outlet)
    result["pressure_drop"] = pressure_drop
    result["pressure_loss_coefficient"] = 2.0 * pressure_drop
    sampled = []
    for index in range(101):
        eta = index / 100
        radius = 2.3 + eta
        velocity = sample_cells(cells, radius * math.sin(angle), 2.8 - radius * math.cos(angle))
        sampled.append(velocity["U_x"] * tangent[0] + velocity["U_y"] * tangent[1])
    result["minimum_streamwise_velocity_75deg"] = min(sampled)
    yplus = [float(row["y_plus"]) for row in wall]
    result["yplus_min"] = min(yplus)
    result["yplus_max"] = max(yplus)
    viscosity = 1.0 / 43000.0
    result["cap_cell_fraction"] = sum(cell["mu_t"] >= 1000.0 * viscosity * (1.0 - 1e-9) for cell in cells) / len(cells)
    write_rows(DATA / f"{spec['name']}_75deg.csv",
               ("y_over_d", "x", "y", "velocity_magnitude", "streamwise_velocity", "reference",
                "reference_uncertainty", "error"), rows)
    return result, rows


def svg_plot(path, title, xlabel, ylabel, series, logarithmic=False):
    width, height = 860, 560
    left, right, top, bottom = 100, 30, 45, 75
    all_points = [point for item in series for point in item["points"]]
    if not all_points:
        return
    xs = [point[0] for point in all_points]
    ys = [point[1] for point in all_points]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    if logarithmic:
        xmin = max(xmin, 1e-12)
        xmin_log, xmax_log = math.log10(xmin), math.log10(xmax)
    else:
        xmin_log, xmax_log = xmin, xmax
    if abs(xmax_log - xmin_log) < 1e-14:
        xmax_log = xmin_log + 1.0
    if abs(ymax - ymin) < 1e-14:
        ymax = ymin + 1.0
    margin = 0.06 * (ymax - ymin)
    ymin -= margin
    ymax += margin

    def sx(value):
        value = math.log10(max(value, 1e-12)) if logarithmic else value
        return left + (value - xmin_log) / (xmax_log - xmin_log) * (width - left - right)

    def sy(value):
        return top + (ymax - value) / (ymax - ymin) * (height - top - bottom)

    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="{width / 2}" y="25" text-anchor="middle" font-family="sans-serif" font-size="20">{escape(title)}</text>']
    for index in range(6):
        fraction = index / 5
        y = top + fraction * (height - top - bottom)
        value = ymax - fraction * (ymax - ymin)
        lines.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}" stroke="#dddddd"/>')
        lines.append(f'<text x="{left - 10}" y="{y + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="12">{value:.3g}</text>')
        x = left + fraction * (width - left - right)
        lines.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{height - bottom}" stroke="#eeeeee"/>')
        if logarithmic:
            tick = 10 ** (xmin_log + fraction * (xmax_log - xmin_log))
        else:
            tick = xmin_log + fraction * (xmax_log - xmin_log)
        lines.append(f'<text x="{x:.2f}" y="{height - bottom + 22}" text-anchor="middle" font-family="sans-serif" font-size="12">{tick:.3g}</text>')
    lines.extend([
        f'<line x1="{left}" y1="{height - bottom}" x2="{width - right}" y2="{height - bottom}" stroke="black"/>',
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{height - bottom}" stroke="black"/>',
        f'<text x="{width / 2}" y="{height - 20}" text-anchor="middle" font-family="sans-serif" font-size="16">{escape(xlabel)}</text>',
        f'<text x="22" y="{height / 2}" text-anchor="middle" transform="rotate(-90 22 {height / 2})" font-family="sans-serif" font-size="16">{escape(ylabel)}</text>',
    ])
    for item in series:
        points = item["points"]
        color = item.get("color", "#0066aa")
        if item.get("line", True) and len(points) > 1:
            data = " ".join(f"{sx(x):.2f},{sy(y):.2f}" for x, y in points)
            lines.append(f'<polyline points="{data}" fill="none" stroke="{color}" stroke-width="2"/>')
        if item.get("markers", False):
            for x, y in points:
                lines.append(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="3.2" fill="{color}"/>')
    for index, item in enumerate(series):
        y = top + 18 * index + 10
        lines.append(f'<line x1="{width - 210}" y1="{y}" x2="{width - 185}" y2="{y}" stroke="{item.get("color", "#0066aa")}" stroke-width="2"/>')
        lines.append(f'<text x="{width - 178}" y="{y + 4}" font-family="sans-serif" font-size="12">{escape(item["name"])}</text>')
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="Run every case before post-processing")
    args = parser.parse_args()
    failures = []
    metrics_rows = []
    plot_data = {"plate": {}, "elbow": {}}
    for spec in cases():
        name = spec["name"]
        if args.run:
            print(f"Running {name}", flush=True)
            completed = subprocess.run([str(PROJECT / "multiphysics"), str(spec["input"])], cwd=PROJECT)
            if completed.returncode:
                failures.append(f"{name}: solver returned {completed.returncode}")
        result_path = spec["result"]
        summary_path = result_path.with_name(result_path.stem + "_summary.csv")
        wall_path = result_path.with_name(result_path.stem + "_wall.csv")
        if not result_path.exists() or not summary_path.exists() or not wall_path.exists():
            failures.append(f"{name}: missing result, summary, or wall output")
            continue
        summary = read_rows(summary_path)[-1]
        cells = read_tecplot(result_path)
        wall = read_rows(wall_path)
        row = {"case": name, "geometry": spec["geometry"], "model": spec["model"], "level": spec["level"],
               "iterations": int(summary["iterations"]), "converged": int(summary["converged"]),
               "max_cell_imbalance": float(summary["max_cell_imbalance"]),
               "net_boundary_flux": float(summary["net_boundary_flux"]),
               "turbulence_change": float(summary["turbulence_change"])}
        if not row["converged"]:
            failures.append(f"{name}: coupled solve did not converge")
        if row["max_cell_imbalance"] > 1e-8 or abs(row["net_boundary_flux"]) > 1e-8:
            failures.append(f"{name}: mass conservation check failed")
        required = ("U_x", "U_y", "p", "k", spec["scale"], "mu_t")
        for cell in cells:
            if not all(math.isfinite(cell[field]) for field in required):
                failures.append(f"{name}: non-finite solution value")
                break
        minimum_scale = min(cell[spec["scale"]] for cell in cells)
        minimum_k = min(cell["k"] for cell in cells)
        if minimum_scale <= 0.0 or minimum_k <= 0.0:
            failures.append(f"{name}: non-positive turbulence variable")
        row["minimum_k"] = minimum_k
        row["minimum_scale"] = minimum_scale
        row["maximum_mu_t"] = max(cell["mu_t"] for cell in cells)
        if spec["geometry"] == "plate":
            metrics, cf, uplus = analyze_plate(spec, cells, wall)
            plot_data["plate"].setdefault(spec["model"], {})[spec["level"]] = (cf, uplus)
            if metrics["cf_rmse_relative"] > 0.10:
                failures.append(f"{name}: Cf normalized RMSE exceeds 10%")
            if metrics["uplus_rmse_relative"] > 0.10:
                failures.append(f"{name}: u+ normalized RMSE exceeds 10%")
        else:
            metrics, profile = analyze_elbow(spec, cells, wall)
            plot_data["elbow"].setdefault(spec["model"], {})[spec["level"]] = profile
        row.update(metrics)
        metrics_rows.append(row)
        print(f"{name}: iterations={row['iterations']}, imbalance={row['max_cell_imbalance']:.3e}, "
              + (f"Cf RMSE={row['cf_rmse_relative']:.3%}, u+ RMSE={row['uplus_rmse_relative']:.3%}"
                 if spec["geometry"] == "plate" else f"75deg RMSE={row['elbow_rmse_relative']:.3%}"), flush=True)

    if metrics_rows:
        fields = ["case", "geometry", "model", "level"]
        fields += sorted({key for row in metrics_rows for key in row if key not in fields})
        write_rows(OUTPUT / "metrics.csv", fields, metrics_rows)

    reference_cf = read_reference(PROJECT / "reference" / "flat_plate" / "cf_karman_schoenherr.dat")
    reference_coles = read_reference(PROJECT / "reference" / "flat_plate" / "uplus_coles_retheta_10000.dat")
    cf_series = [{"name": "Karman-Schoenherr", "points": [point for point in reference_cf if 3500 <= point[0] <= 14000],
                  "color": "#111111", "markers": True}]
    uplus_series = [{"name": "Coles Re-theta=10000", "points": reference_coles,
                     "color": "#111111", "markers": False}]
    colors = {"k_epsilon": "#0066cc", "sst_k_omega": "#cc3300"}
    for model, levels in plot_data["plate"].items():
        if "medium" in levels:
            cf, uplus = levels["medium"]
            cf_series.append({"name": model, "points": [(float(row["retheta"]), float(row["cf"])) for row in cf],
                              "color": colors[model]})
            uplus_series.append({"name": model, "points": [(float(row["y_plus"]), float(row["u_plus"])) for row in uplus],
                                 "color": colors[model]})
    svg_plot(FIGURES / "plate_cf_comparison.svg", "Flat-plate skin friction",
             "Momentum-thickness Reynolds number", "Cf", cf_series)
    svg_plot(FIGURES / "plate_uplus_comparison.svg", "Flat-plate velocity at Re-theta=10000",
             "y+", "u+", uplus_series, logarithmic=True)

    elbow_reference = read_rows(PROJECT / "reference" / "elbow_enayet_1982" / "velocity_75deg_experimental.csv")
    elbow_series = [{"name": "Enayet experiment", "color": "#111111", "line": False, "markers": True,
                     "points": [(float(row["y_over_d"]), float(row["velocity_over_bulk"])) for row in elbow_reference]}]
    for model, levels in plot_data["elbow"].items():
        if "medium" in levels:
            elbow_series.append({"name": model, "color": colors[model],
                                 "points": [(float(row["y_over_d"]), float(row["velocity_magnitude"])) for row in levels["medium"]]})
    svg_plot(FIGURES / "elbow_75deg_comparison.svg", "90-degree elbow: 75-degree section",
             "y/d", "Velocity magnitude / bulk velocity", elbow_series)

    if failures:
        print("Turbulence validation failed:", file=sys.stderr)
        for failure in failures:
            print(f"  {failure}", file=sys.stderr)
        return 1
    print("All turbulence conservation and accuracy checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
