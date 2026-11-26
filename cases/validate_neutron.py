#!/usr/bin/env python3
import argparse
import csv
import math
import subprocess
from pathlib import Path
from xml.sax.saxutils import escape

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
OUTPUT = PROJECT / "validation"
PI = math.pi

CASES = [
    {"name": "attenuation_16", "directory": "neutron_attenuation", "input": "coarse.in", "stem": "attenuation_16", "kind": "attenuation", "h": 1 / 16},
    {"name": "attenuation_64", "directory": "neutron_attenuation", "input": "case.in", "stem": "attenuation_64", "kind": "attenuation", "h": 1 / 64},
    {"name": "attenuation_triangles", "directory": "neutron_attenuation", "input": "triangles.in", "stem": "attenuation_triangles", "kind": "attenuation", "h": 1 / 32},
    {"name": "sn_critical", "directory": "neutron_sn_critical", "input": "case.in", "stem": "sn_critical", "kind": "sn_critical"},
    {"name": "diffusion_16", "directory": "neutron_diffusion", "input": "coarse.in", "stem": "diffusion_16", "kind": "diffusion", "h": 1 / 16},
    {"name": "diffusion_32", "directory": "neutron_diffusion", "input": "case.in", "stem": "diffusion_32", "kind": "diffusion", "h": 1 / 32},
    {"name": "diffusion_two_group", "directory": "neutron_diffusion_two_group", "input": "case.in", "stem": "diffusion_two_group", "kind": "two_group", "h": 1 / 32},
    {"name": "diffusion_marshak", "directory": "neutron_diffusion_marshak", "input": "case.in", "stem": "diffusion_marshak", "kind": "marshak"},
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


def buckling(h):
    return 8.0 * math.sin(PI * h / 2.0) ** 2 / h**2 if h else 2.0 * PI**2


def one_group_reference(h):
    return 22.0 / (1.0 + buckling(h))


def two_group_reference(h):
    a11 = 1.2 + buckling(h)
    a22 = 2.0 + 0.5 * buckling(h)
    scattering = 0.2
    p1 = 0.8 / a11
    p2 = (0.2 + scattering * p1) / a22
    return 8.0 * p1 + 40.0 * p2, p1, p2


def normalized(values, target=1.0):
    return [value * target * len(values) / sum(values) for value in values]


def exact_fields(spec, cells):
    kind = spec["kind"]
    x = [float(row["x"]) for row in cells]
    y = [float(row["y"]) for row in cells]
    if kind == "attenuation":
        directions = read_rows(ROOT / spec["directory"] / "results" / f"{spec['stem']}_directions.csv")
        weight_sum = sum(float(row["weight"]) for row in directions)
        if abs(weight_sum - 4.0 * PI) > 1.0e-10:
            raise ValueError(f"{spec['name']}: angular weights do not sum to 4 pi")
        values = []
        for x_value in x:
            values.append(sum(float(row["weight"]) * math.exp(-x_value / float(row["omega_x"]))
                              for row in directions if float(row["omega_x"]) > 1.0e-14))
        return {"phi_total": values}
    if kind == "sn_critical":
        return {"phi_total": [1.0] * len(cells)}
    if kind in ("diffusion", "two_group"):
        shape = [math.sin(PI * x_value) * math.sin(PI * y_value) for x_value, y_value in zip(x, y)]
        if kind == "diffusion":
            return {"phi_total": normalized(shape)}
        _, p1, p2 = two_group_reference(spec["h"])
        total = normalized(shape)
        scale = 1.0 / (p1 + p2)
        return {"phi_g1": [scale * p1 * value for value in total],
                "phi_g2": [scale * p2 * value for value in total]}
    if kind == "marshak":
        b = -(math.sinh(1.0) + 0.5 * math.cosh(1.0)) / (math.cosh(1.0) + 0.5 * math.sinh(1.0))
        return {"phi_total": [math.cosh(x_value) + b * math.sinh(x_value) for x_value in x]}
    raise ValueError(f"Unknown case kind: {kind}")


def metrics(numerical, exact):
    numerator = sum((value - reference)**2 for value, reference in zip(numerical, exact))
    denominator = sum(reference**2 for reference in exact)
    return (math.sqrt(numerator / denominator) if denominator > 0.0 else math.sqrt(numerator),
            max(abs(value - reference) for value, reference in zip(numerical, exact)))


def analyze(spec, failures):
    directory = ROOT / spec["directory"] / "results"
    summary = read_rows(directory / f"{spec['stem']}_summary.csv")[0]
    cells = read_rows(directory / f"{spec['stem']}_cells.csv")
    exact = exact_fields(spec, cells)
    kind = spec["kind"]
    fields = ["phi_total"] if kind in ("attenuation", "sn_critical", "diffusion", "marshak") else ["phi_g1", "phi_g2"]
    numerical = {field: [float(row[field]) for row in cells] for field in fields}
    combined = [value for field in fields for value in numerical[field]]
    reference = [value for field in fields for value in exact[field]]
    l2, linf = metrics(combined, reference)
    k_eff = summary["k_eff"]
    k_reference = ""
    k_continuous = ""
    k_error = ""
    if kind == "sn_critical":
        k_reference = k_continuous = 1.2
    elif kind == "diffusion":
        k_reference = one_group_reference(spec["h"])
        k_continuous = one_group_reference(0.0)
    elif kind == "two_group":
        k_reference = two_group_reference(spec["h"])[0]
        k_continuous = two_group_reference(0.0)[0]
    if k_reference != "":
        k_error = abs(float(k_eff) - k_reference)
    row = {"case": spec["name"], "model": summary["model"], "k_eff": k_eff,
           "k_reference": k_reference, "k_continuous": k_continuous, "k_error": k_error,
           "l2_relative": l2, "linf_absolute": linf, "balance_relative": float(summary["balance_relative"]),
           "flux_min": float(summary["flux_min"]), "flux_max": float(summary["flux_max"]),
           "fission_production": float(summary["fission_production"]),
           "absorption": float(summary["absorption"]), "net_leakage": float(summary["net_leakage"]),
           "fission_power": float(summary["fission_power"]),
           "outer_iterations": int(summary["outer_iterations"])}
    balance = float(summary["balance_relative"])
    if balance > 1.0e-8:
        failures.append(f"{spec['name']}: particle balance {balance:.3e} exceeds 1e-8")
    if float(summary["flux_min"]) < -1.0e-8:
        failures.append(f"{spec['name']}: negative neutron flux {summary['flux_min']}")
    if k_error != "":
        tolerance = 2.0e-8 if kind == "sn_critical" else 5.0e-9
        if k_error > tolerance:
            failures.append(f"{spec['name']}: k_eff error {k_error:.3e} exceeds {tolerance:.1e}")
        if abs(float(k_eff) - float(summary["k_balance"])) > 1.0e-8:
            failures.append(f"{spec['name']}: eigenvalue and reaction-rate balance disagree")
    if kind == "attenuation":
        limit = 0.04 if spec["h"] >= 1 / 16 else 0.01 if spec["name"] == "attenuation_64" else 0.015
    else:
        limit = {"sn_critical": 1.0e-8, "diffusion": 1.0e-6, "two_group": 1.0e-6, "marshak": 0.005}[kind]
    if l2 > limit:
        failures.append(f"{spec['name']}: distribution L2 error {l2:.3e} exceeds {limit:.1e}")
    return row, {"spec": spec, "cells": cells, "numerical": numerical, "exact": exact, "summary": summary}


def line_svg(path, series, title, ylabel):
    width, height = 900, 500
    left, right, top, bottom = 85, 30, 55, 65
    points = [point for item in series for point in item["points"]]
    xmin, xmax = min(p[0] for p in points), max(p[0] for p in points)
    ymin, ymax = min(p[1] for p in points), max(p[1] for p in points)
    if xmax == xmin:
        xmax = xmin + 1.0
    if ymax == ymin:
        ymin, ymax = ymin - 0.5, ymax + 0.5
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
        xv = xmin + (xmax - xmin) * index / 5
        yv = ymin + (ymax - ymin) * index / 5
        x, _ = point(xv, ymin)
        _, y = point(xmin, yv)
        parts.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{height - bottom}" stroke="#e5e5e5"/>')
        parts.append(f'<line x1="{left}" y1="{y:.2f}" x2="{width - right}" y2="{y:.2f}" stroke="#e5e5e5"/>')
        parts.append(f'<text x="{x:.2f}" y="{height - bottom + 20}" text-anchor="middle" font-family="sans-serif" font-size="11">{xv:.3g}</text>')
        parts.append(f'<text x="{left - 8}" y="{y + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="11">{yv:.3g}</text>')
    parts.append(f'<rect x="{left}" y="{top}" width="{width - left - right}" height="{height - top - bottom}" fill="none" stroke="#333"/>')
    for item in series:
        data = " ".join(f"{x:.2f},{y:.2f}" for x, y in (point(*p) for p in item["points"]))
        dash = ' stroke-dasharray="6 4"' if item.get("dash") else ""
        parts.append(f'<polyline points="{data}" fill="none" stroke="{item["color"]}" stroke-width="1.8"{dash}/>')
    parts.extend([f'<text x="{width / 2}" y="{height - 15}" text-anchor="middle" font-family="sans-serif" font-size="14">x</text>',
                  f'<text x="18" y="{height / 2}" text-anchor="middle" transform="rotate(-90 18 {height / 2})" font-family="sans-serif" font-size="14">{escape(ylabel)}</text>'])
    for index, item in enumerate(series):
        x = width - 245
        y = top + 20 + 22 * index
        parts.append(f'<line x1="{x}" y1="{y}" x2="{x + 28}" y2="{y}" stroke="{item["color"]}" stroke-width="2"/>')
        parts.append(f'<text x="{x + 35}" y="{y + 4}" font-family="sans-serif" font-size="12">{escape(item["label"])}</text>')
    parts.append('</svg>')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def centreline(detail, field, exact=False):
    cells = detail["cells"]
    values = detail["exact" if exact else "numerical"][field]
    target_y = 0.5
    y_value = min({float(row["y"]) for row in cells}, key=lambda y: abs(y - target_y))
    pairs = [(float(row["x"]), value) for row, value in zip(cells, values) if abs(float(row["y"]) - y_value) < 1.0e-12]
    return sorted(pairs)


def make_figures(details):
    selected = {detail["spec"]["kind"]: detail for detail in details}
    selected["attenuation"] = next(detail for detail in details if detail["spec"]["name"] == "attenuation_64")
    selected["diffusion"] = next(detail for detail in details if detail["spec"]["name"] == "diffusion_32")
    styles = [("Numerical", "#1764ab", False), ("Analytical", "#c43c39", True)]
    for kind, filename, title, ylabel in [
        ("attenuation", "neutron_attenuation.svg", "2D SN absorbing-beam attenuation", "scalar flux"),
        ("diffusion", "neutron_diffusion.svg", "2D one-group diffusion eigenmode", "normalized scalar flux"),
        ("marshak", "neutron_diffusion_marshak.svg", "2D diffusion with Marshak vacuum boundary", "scalar flux"),
    ]:
        detail = selected[kind]
        series = [{"label": label, "color": color, "dash": dash,
                   "points": centreline(detail, "phi_total", exact)} for (label, color, dash), exact in zip(styles, (False, True))]
        line_svg(OUTPUT / "figures" / filename, series, title, ylabel)
    detail = selected["two_group"]
    series = []
    for field, label, color in [("phi_g1", "group 1", "#1764ab"), ("phi_g2", "group 2", "#2f8f4e")]:
        series.append({"label": label + " numerical", "color": color, "points": centreline(detail, field)})
        series.append({"label": label + " analytical", "color": color, "dash": True, "points": centreline(detail, field, True)})
    line_svg(OUTPUT / "figures" / "neutron_two_group.svg", series, "2D two-group diffusion eigenmode", "normalized scalar flux")


def write_report(rows, failures, one_order, attenuation_order):
    lines = ["# Neutron transport validation", "",
             "The cases exercise conservative SN upwind transport, exact axis reflection, multigroup diffusion, fixed-source and eigenvalue iteration, Dirichlet and Marshak boundaries, and triangle/quadrilateral meshes. Eigenvalue references use the matching discrete-grid buckling; continuous references are also listed.", "",
             "| Case | k_eff | Reference k | k error | Distribution L2 | Balance | Fission production | Absorption | Net leakage |",
             "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        k_eff = f"{float(row['k_eff']):.10f}" if row["k_eff"] else "—"
        reference = f"{float(row['k_reference']):.10f}" if row["k_reference"] != "" else "—"
        error = f"{float(row['k_error']):.3e}" if row["k_error"] != "" else "—"
        lines.append(f"| {row['case']} | {k_eff} | {reference} | {error} | {row['l2_relative']:.3e} | {row['balance_relative']:.3e} | {row['fission_production']:.9e} | {row['absorption']:.9e} | {row['net_leakage']:.9e} |")
    lines.extend(["", f"One-group continuous-eigenvalue refinement order: **{one_order:.4f}** (expected approximately 2).",
                  f"Attenuation L2 refinement order: **{attenuation_order:.4f}**.", "",
                  "The two-group reference solves $A_{11}\\phi_1=\\chi_1R$ and $A_{22}\\phi_2-\\Sigma_{s,1\\to2}\\phi_1=\\chi_2R$, with removal $\\Sigma_{r,g}=\\Sigma_{t,g}-\\Sigma_{s,g\\to g}$. The balance check independently evaluates $P_f/k=A+L$ for eigenvalue problems and $P_f+Q=A+L$ for fixed-source problems.", "",
                  "## Distribution plots", "",
                  "- [SN attenuation](figures/neutron_attenuation.svg)",
                  "- [One-group diffusion eigenmode](figures/neutron_diffusion.svg)",
                  "- [Two-group diffusion eigenmode](figures/neutron_two_group.svg)",
                  "- [Marshak fixed-source diffusion](figures/neutron_diffusion_marshak.svg)", ""])
    if failures:
        lines.extend(["## Failures", ""] + [f"- {failure}" for failure in failures])
    else:
        lines.append("**All neutron validation checks passed.**")
    (OUTPUT / "neutron_validation.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    args = parser.parse_args()
    if args.run:
        subprocess.run(["python3", str(ROOT / "generate_neutron_meshes.py")], cwd=PROJECT, check=True)
        for spec in CASES:
            case_file = ROOT / spec["directory"] / spec["input"]
            print(f"\n=== {spec['name']} ===", flush=True)
            subprocess.run([str(PROJECT / "multiphysics"), str(case_file)], cwd=PROJECT, check=True)
    failures = []
    details = [analyze(spec, failures) for spec in CASES]
    rows = [item[0] for item in details]
    details = [item[1] for item in details]
    fieldnames = list(rows[0])
    write_rows(OUTPUT / "neutron_metrics.csv", fieldnames, rows)
    by_name = {row["case"]: row for row in rows}
    coarse, fine = by_name["diffusion_16"], by_name["diffusion_32"]
    one_order = math.log(abs(float(coarse["k_eff"]) - float(coarse["k_continuous"])) /
                         abs(float(fine["k_eff"]) - float(fine["k_continuous"])), 2)
    attenuation_order = math.log(by_name["attenuation_16"]["l2_relative"] / by_name["attenuation_64"]["l2_relative"], 2)
    if not 1.7 <= one_order <= 2.3:
        failures.append(f"one-group eigenvalue convergence order {one_order:.3f} is outside [1.7, 2.3]")
    make_figures(details)
    write_report(rows, failures, one_order, attenuation_order)
    print(f"\nOne-group eigenvalue order: {one_order:.4f}; attenuation order: {attenuation_order:.4f}")
    if failures:
        for failure in failures:
            print("FAIL:", failure)
        return 1
    print("All neutron validation checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
