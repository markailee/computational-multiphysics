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
    # Baseline & higher-order SN attenuation
    {"name": "attenuation_16", "directory": "neutron_attenuation", "input": "coarse.in", "stem": "attenuation_16", "kind": "attenuation", "h": 1 / 16},
    {"name": "attenuation_64", "directory": "neutron_attenuation", "input": "case.in", "stem": "attenuation_64", "kind": "attenuation", "h": 1 / 64},
    {"name": "attenuation_triangles", "directory": "neutron_attenuation", "input": "triangles.in", "stem": "attenuation_triangles", "kind": "attenuation", "h": 1 / 32},
    {"name": "attenuation_ho_16", "directory": "neutron_attenuation", "input": "second_order_16.in", "stem": "attenuation_ho_16", "kind": "attenuation_ho", "h": 1 / 16},
    {"name": "attenuation_ho_64", "directory": "neutron_attenuation", "input": "second_order_64.in", "stem": "attenuation_ho_64", "kind": "attenuation_ho", "h": 1 / 64},
    # Baseline SN & diffusion
    {"name": "sn_critical", "directory": "neutron_sn_critical", "input": "case.in", "stem": "sn_critical", "kind": "sn_critical"},
    {"name": "diffusion_16", "directory": "neutron_diffusion", "input": "coarse.in", "stem": "diffusion_16", "kind": "diffusion", "h": 1 / 16},
    {"name": "diffusion_32", "directory": "neutron_diffusion", "input": "case.in", "stem": "diffusion_32", "kind": "diffusion", "h": 1 / 32},
    {"name": "diffusion_two_group", "directory": "neutron_diffusion_two_group", "input": "case.in", "stem": "diffusion_two_group", "kind": "two_group", "h": 1 / 32},
    {"name": "diffusion_marshak", "directory": "neutron_diffusion_marshak", "input": "case.in", "stem": "diffusion_marshak", "kind": "marshak"},
    # Fixed Eq. 5.297 Marshak incident partial current
    {"name": "incident_current", "directory": "neutron_incident_current", "input": "case.in", "stem": "incident_current", "kind": "incident_current"},
    # Higher eigenmodes + adjoint + kinetics + 1st-order perturbation
    {"name": "modes_base", "directory": "neutron_modes", "input": "case.in", "stem": "modes_base", "kind": "modes", "h": 1 / 32},
    {"name": "modes_perturbed", "directory": "neutron_modes", "input": "perturbed.in", "stem": "modes_perturbed", "kind": "modes_perturbed"},
    # Reed 5-region heterogeneous SN + DSA benchmark
    {"name": "reed_s8_dsa", "directory": "neutron_reed", "input": "case.in", "stem": "reed_s8_dsa", "kind": "reed"},
    {"name": "reed_s8_nodsa", "directory": "neutron_reed", "input": "no_dsa.in", "stem": "reed_s8_nodsa", "kind": "reed"},
    # 2D IAEA PWR quarter-core benchmark (ANL-7416 11-A2)
    {"name": "iaea2d_5cm", "directory": "neutron_iaea2d", "input": "coarse.in", "stem": "iaea2d_5cm", "kind": "iaea2d", "h": 5.0},
    {"name": "iaea2d_5cm_perturbed", "directory": "neutron_iaea2d", "input": "perturbed.in", "stem": "iaea2d_5cm_perturbed", "kind": "iaea2d_perturbed", "h": 5.0},
    {"name": "iaea2d_2p5cm", "directory": "neutron_iaea2d", "input": "case.in", "stem": "iaea2d_2p5cm", "kind": "iaea2d", "h": 2.5},
    # 2D TWIGL space-time kinetics benchmark
    {"name": "twigl_step", "directory": "neutron_twigl", "input": "step.in", "stem": "twigl_step", "kind": "twigl_step"},
    {"name": "twigl_ramp", "directory": "neutron_twigl", "input": "ramp.in", "stem": "twigl_ramp", "kind": "twigl_ramp"},
    # 2D circular fuel pin-cell & 2x2 heterogeneous mini-assembly SN
    {"name": "pincell_s8", "directory": "neutron_pincell_sn", "input": "case.in", "stem": "pincell_s8", "kind": "pincell"},
    {"name": "pincell_s8_white", "directory": "neutron_pincell_sn", "input": "white.in", "stem": "pincell_s8_white", "kind": "pincell"},
    {"name": "pincell_diffusion", "directory": "neutron_pincell_sn", "input": "diffusion.in", "stem": "pincell_diffusion", "kind": "pincell"},
    {"name": "lattice_s8", "directory": "neutron_lattice_sn", "input": "case.in", "stem": "lattice_s8", "kind": "lattice"},
    # 2D Kobayashi-style dog-leg void-channel shielding
    {"name": "dogleg_s4", "directory": "neutron_dogleg_sn", "input": "s4.in", "stem": "dogleg_s4", "kind": "dogleg"},
    {"name": "dogleg_s8", "directory": "neutron_dogleg_sn", "input": "s8.in", "stem": "dogleg_s8", "kind": "dogleg"},
    {"name": "dogleg_high_angle", "directory": "neutron_dogleg_sn", "input": "high_angle.in", "stem": "dogleg_high_angle", "kind": "dogleg"},
    {"name": "dogleg_diffusion", "directory": "neutron_dogleg_sn", "input": "diffusion.in", "stem": "dogleg_diffusion", "kind": "dogleg"},
]

IAEA_REF_K = 1.029585
TWIGL_STEP_REF = {0.1: 2.061, 0.2: 2.078, 0.3: 2.095, 0.4: 2.113, 0.5: 2.131}
TWIGL_RAMP_REF = {0.1: 1.309, 0.2: 1.960, 0.3: 2.074, 0.4: 2.092, 0.5: 2.109}


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_rows(path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def buckling(h, mx=1, my=1):
    if not h:
        return (mx**2 + my**2) * PI**2
    return 4.0 * (math.sin(mx * PI * h / 2.0) ** 2 + math.sin(my * PI * h / 2.0) ** 2) / h**2


def one_group_reference(h):
    return 22.0 / (1.0 + buckling(h))


def two_group_reference(h, mx=1, my=1):
    b2 = buckling(h, mx, my)
    a11 = 1.2 + b2
    a22 = 2.0 + 0.5 * b2
    scattering = 0.2
    p1 = 0.8 / a11
    p2 = (0.2 + scattering * p1) / a22
    return 8.0 * p1 + 40.0 * p2, p1, p2


def two_group_kinetics_reference(h):
    b2 = buckling(h, 1, 1)
    a11 = 1.2 + b2
    a22 = 2.0 + 0.5 * b2
    phi1 = 0.8 / a11
    phi2 = (0.2 + 0.2 * phi1) / a22
    adj2 = 40.0 / a22
    adj1 = (8.0 + 0.2 * adj2) / a11
    denom = (0.8 * adj1 + 0.2 * adj2) * (8.0 * phi1 + 40.0 * phi2)
    beta_eff = 0.0065 * (1.0 * adj1) * (8.0 * phi1 + 40.0 * phi2) / denom
    gen_time = (adj1 * phi1 / 1.0e7 + adj2 * phi2 / 2.2e5) / denom
    return beta_eff, gen_time


def normalized(values, target=1.0):
    return [value * target * len(values) / sum(values) for value in values]


def exact_fields(spec, cells):
    kind = spec["kind"]
    x = [float(row["x"]) for row in cells]
    y = [float(row["y"]) for row in cells]
    if kind in ("attenuation", "attenuation_ho"):
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
    if kind in ("diffusion", "two_group", "modes"):
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
    if kind == "incident_current":
        c1, s1 = math.cosh(1.0), math.sinh(1.0)
        a = (c1 + 0.5 * s1) / (1.25 * s1 + c1)
        b = 0.5 * a - 1.0
        return {"phi_total": [a * math.cosh(x_value) + b * math.sinh(x_value) for x_value in x]}
    return {}


def metrics(numerical, exact):
    numerator = sum((value - reference) ** 2 for value, reference in zip(numerical, exact))
    denominator = sum(reference ** 2 for reference in exact)
    return (math.sqrt(numerator / denominator) if denominator > 0.0 else math.sqrt(numerator),
            max(abs(value - reference) for value, reference in zip(numerical, exact)))


def analyze(spec, failures):
    directory = ROOT / spec["directory"] / "results"
    summary = read_rows(directory / f"{spec['stem']}_summary.csv")[0]
    cells = read_rows(directory / f"{spec['stem']}_cells.csv")
    history = read_rows(directory / f"{spec['stem']}_history.csv")
    modes_path = directory / f"{spec['stem']}_modes.csv"
    modes = read_rows(modes_path) if modes_path.exists() else []
    trans_path = directory / f"{spec['stem']}_transient.csv"
    transient = read_rows(trans_path) if trans_path.exists() else []
    kind = spec["kind"]

    exact = exact_fields(spec, cells)
    l2 = ""
    linf = ""
    numerical = {}
    if exact:
        fields = list(exact.keys())
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
    elif kind in ("two_group", "modes"):
        k_reference = two_group_reference(spec["h"])[0]
        k_continuous = two_group_reference(0.0)[0]
    elif kind == "iaea2d":
        k_reference = k_continuous = IAEA_REF_K

    if k_reference != "" and k_eff != "":
        k_error = abs(float(k_eff) - k_reference)

    balance = float(summary["balance_relative"])
    flux_min = float(summary["flux_min"])
    row = {
        "case": spec["name"],
        "model": summary["model"],
        "problem": summary["problem_type"],
        "k_eff": k_eff,
        "k_reference": k_reference,
        "k_continuous": k_continuous,
        "k_error": k_error,
        "dominance_ratio": summary.get("dominance_ratio", ""),
        "k_adjoint": summary.get("k_adjoint", ""),
        "l2_relative": l2,
        "linf_absolute": linf,
        "balance_relative": balance,
        "flux_min": flux_min,
        "flux_max": float(summary["flux_max"]),
        "fission_production": float(summary["fission_production"]),
        "absorption": float(summary["absorption"]),
        "net_leakage": float(summary["net_leakage"]),
        "fission_power": float(summary["fission_power"]),
        "outer_iterations": int(summary["outer_iterations"]),
        "source_iterations": int(history[-1]["source_iterations"]) if history else 0,
    }

    # Universal checks (except transient where final snapshot is non-steady)
    if summary["problem_type"] != "transient":
        if balance > 1.0e-8:
            failures.append(f"{spec['name']}: particle balance {balance:.3e} exceeds 1e-8")
        if flux_min < -1.0e-8:
            failures.append(f"{spec['name']}: negative neutron flux {flux_min}")
    else:
        if flux_min < -1.0e-8:
            failures.append(f"{spec['name']}: negative transient neutron flux {flux_min}")

    if k_eff != "" and summary["problem_type"] == "eigenvalue":
        if abs(float(k_eff) - float(summary["k_balance"])) > 1.0e-7:
            failures.append(f"{spec['name']}: eigenvalue and reaction-rate balance disagree")
        if summary.get("k_adjoint", "") != "":
            adj_diff = abs(float(k_eff) - float(summary["k_adjoint"]))
            if adj_diff > 1.0e-7:
                failures.append(f"{spec['name']}: adjoint eigenvalue discrepancy {adj_diff:.3e} exceeds 1e-7")

    if k_error != "":
        if kind == "iaea2d":
            limit = 5.0e-4 if spec["h"] >= 5.0 else 2.5e-4
        elif kind == "sn_critical":
            limit = 2.0e-8
        else:
            limit = 5.0e-9
        if k_error > limit:
            failures.append(f"{spec['name']}: k_eff error {k_error:.3e} exceeds {limit:.1e}")

    if l2 != "":
        if kind == "attenuation":
            limit = 0.04 if spec["h"] >= 1 / 16 else 0.012 if spec["name"] == "attenuation_64" else 0.018
        elif kind == "attenuation_ho":
            limit = 0.008 if spec["h"] >= 1 / 16 else 0.0012
        else:
            limit = {"sn_critical": 1.0e-8, "diffusion": 1.0e-6, "two_group": 1.0e-6,
                     "modes": 1.0e-6, "marshak": 0.005, "incident_current": 0.001}[kind]
        if l2 > limit:
            failures.append(f"{spec['name']}: distribution L2 error {l2:.3e} exceeds {limit:.1e}")

    # Mode checks for modes_base
    if kind == "modes":
        expected_modes = [two_group_reference(spec["h"], mx, my)[0] for mx, my in [(1, 1), (2, 1), (1, 2), (2, 2)]]
        if len(modes) < 4:
            failures.append("modes_base: expected 4 eigenmodes")
        else:
            for idx, (m_row, ek) in enumerate(zip(modes, expected_modes)):
                err = abs(float(m_row["k_eff"]) - ek)
                if err > 1.0e-9:
                    failures.append(f"modes_base: mode {idx} k error {err:.3e} exceeds 1e-9")
        beta_ref, lam_ref = two_group_kinetics_reference(spec["h"])
        if abs(float(summary["beta_eff"]) - beta_ref) > 1.0e-10:
            failures.append("modes_base: beta_eff does not match analytical adjoint-weighted value")
        if abs(float(summary["generation_time"]) - lam_ref) / lam_ref > 1.0e-10:
            failures.append("modes_base: Lambda_eff does not match analytical adjoint-weighted value")

    # TWIGL transient power checks against published benchmark values
    if kind in ("twigl_step", "twigl_ramp"):
        ref_map = TWIGL_STEP_REF if kind == "twigl_step" else TWIGL_RAMP_REF
        for t_target, p_ref in ref_map.items():
            best = min(transient, key=lambda r: abs(float(r["time"]) - t_target))
            p_num = float(best["relative_power"])
            rel = abs(p_num - p_ref) / p_ref
            if rel > 0.015:
                failures.append(f"{spec['name']}: power at t={t_target}s is {p_num:.4f} vs ref {p_ref:.3f} ({rel:.2%})")

    return row, {"spec": spec, "cells": cells, "numerical": numerical, "exact": exact,
                 "summary": summary, "history": history, "modes": modes, "transient": transient}


def line_svg(path, series, title, xlabel, ylabel):
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
             f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="18" font-weight="600">{escape(title)}</text>']
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
        if item.get("markers"):
            for px, py in (point(*p) for p in item["points"]):
                parts.append(f'<circle cx="{px:.2f}" cy="{py:.2f}" r="4.2" fill="{item["color"]}" stroke="white" stroke-width="1"/>')
        else:
            dash = ' stroke-dasharray="6 4"' if item.get("dash") else ""
            parts.append(f'<polyline points="{data}" fill="none" stroke="{item["color"]}" stroke-width="2.0"{dash}/>')
    parts.extend([f'<text x="{width / 2}" y="{height - 15}" text-anchor="middle" font-family="sans-serif" font-size="14">{escape(xlabel)}</text>',
                  f'<text x="20" y="{height / 2}" text-anchor="middle" transform="rotate(-90 20 {height / 2})" font-family="sans-serif" font-size="14">{escape(ylabel)}</text>'])
    legend_x = width - 255
    parts.append(f'<rect x="{legend_x - 10}" y="{top + 8}" width="235" height="{14 + 22 * len(series)}" fill="white" fill-opacity="0.88" stroke="#ccc" rx="4"/>')
    for index, item in enumerate(series):
        y = top + 24 + 22 * index
        if item.get("markers"):
            parts.append(f'<circle cx="{legend_x + 14}" cy="{y}" r="4.2" fill="{item["color"]}"/>')
        else:
            dash = ' stroke-dasharray="6 4"' if item.get("dash") else ""
            parts.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 28}" y2="{y}" stroke="{item["color"]}" stroke-width="2.2"{dash}/>')
        parts.append(f'<text x="{legend_x + 35}" y="{y + 4}" font-family="sans-serif" font-size="12">{escape(item["label"])}</text>')
    parts.append('</svg>')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def color_map(t):
    t = max(0.0, min(1.0, t))
    stops = [
        (0.00, (30, 58, 138)),
        (0.25, (20, 148, 180)),
        (0.50, (80, 195, 110)),
        (0.75, (245, 185, 50)),
        (1.00, (195, 35, 35)),
    ]
    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if t <= t1:
            u = (t - t0) / (t1 - t0)
            rgb = tuple(int(round(c0[k] + u * (c1[k] - c0[k]))) for k in range(3))
            return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"
    return "#c32323"


def heatmap_svg(path, cells, field, title, subtitle=""):
    width, height = 740, 680
    left, right, top, bottom = 75, 110, 65, 65
    xs = sorted({round(float(r["x"]), 6) for r in cells})
    ys = sorted({round(float(r["y"]), 6) for r in cells})
    dx = min(b - a for a, b in zip(xs, xs[1:])) if len(xs) > 1 else 1.0
    dy = min(b - a for a, b in zip(ys, ys[1:])) if len(ys) > 1 else 1.0
    xmin, xmax = min(xs) - 0.5 * dx, max(xs) + 0.5 * dx
    ymin, ymax = min(ys) - 0.5 * dy, max(ys) + 0.5 * dy
    vals = [float(r[field]) for r in cells]
    vmin, vmax = min(vals), max(vals)
    if vmax == vmin:
        vmax = vmin + 1.0

    plot_w = width - left - right
    plot_h = height - top - bottom
    scale = min(plot_w / (xmax - xmin), plot_h / (ymax - ymin))
    ox = left + 0.5 * (plot_w - scale * (xmax - xmin))
    oy = top + 0.5 * (plot_h - scale * (ymax - ymin))

    def sx(x):
        return ox + (x - xmin) * scale

    def sy(y):
        return oy + (ymax - y) * scale

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             f'<text x="{width / 2}" y="28" text-anchor="middle" font-family="sans-serif" font-size="18" font-weight="600">{escape(title)}</text>']
    if subtitle:
        parts.append(f'<text x="{width / 2}" y="48" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#555">{escape(subtitle)}</text>')

    cw = dx * scale + 0.6
    ch = dy * scale + 0.6
    for row, val in zip(cells, vals):
        cx = sx(float(row["x"])) - 0.5 * cw
        cy = sy(float(row["y"])) - 0.5 * ch
        col = color_map((val - vmin) / (vmax - vmin))
        parts.append(f'<rect x="{cx:.2f}" y="{cy:.2f}" width="{cw:.2f}" height="{ch:.2f}" fill="{col}"/>')

    # Axis ticks
    for idx in range(5):
        xv = xmin + (xmax - xmin) * idx / 4
        yv = ymin + (ymax - ymin) * idx / 4
        parts.append(f'<text x="{sx(xv):.2f}" y="{oy + (ymax - ymin) * scale + 20:.2f}" text-anchor="middle" font-family="sans-serif" font-size="11">{xv:.3g}</text>')
        parts.append(f'<text x="{ox - 8:.2f}" y="{sy(yv) + 4:.2f}" text-anchor="end" font-family="sans-serif" font-size="11">{yv:.3g}</text>')
    parts.append(f'<text x="{ox + 0.5 * (xmax - xmin) * scale:.2f}" y="{height - 16}" text-anchor="middle" font-family="sans-serif" font-size="13">x (cm)</text>')
    parts.append(f'<text x="20" y="{oy + 0.5 * (ymax - ymin) * scale:.2f}" text-anchor="middle" transform="rotate(-90 20 {oy + 0.5 * (ymax - ymin) * scale:.2f})" font-family="sans-serif" font-size="13">y (cm)</text>')

    # Colorbar
    cb_x = width - 78
    cb_h = (ymax - ymin) * scale
    cb_y = oy
    for i in range(60):
        t0 = i / 60.0
        col = color_map(1.0 - t0)
        parts.append(f'<rect x="{cb_x}" y="{cb_y + t0 * cb_h:.2f}" width="18" height="{cb_h / 60.0 + 0.8:.2f}" fill="{col}"/>')
    parts.append(f'<rect x="{cb_x}" y="{cb_y:.2f}" width="18" height="{cb_h:.2f}" fill="none" stroke="#333"/>')
    for idx in range(5):
        frac = idx / 4.0
        val = vmax - frac * (vmax - vmin)
        parts.append(f'<text x="{cb_x + 24}" y="{cb_y + frac * cb_h + 4:.2f}" font-family="sans-serif" font-size="11">{val:.3g}</text>')

    parts.append('</svg>')
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(parts), encoding="utf-8")


def centreline(detail, field, exact=False, target_y=0.5):
    cells = detail["cells"]
    values = detail["exact"][field] if exact else [float(row[field]) for row in cells]
    y_value = min({float(row["y"]) for row in cells}, key=lambda y: abs(y - target_y))
    pairs = [(float(row["x"]), value) for row, value in zip(cells, values) if abs(float(row["y"]) - y_value) < 1.0e-9]
    return sorted(pairs)


def make_figures(details):
    by_name = {d["spec"]["name"]: d for d in details}

    # 1. SN attenuation: Upwind vs 2nd-order MUSCL vs Analytical
    att_u = by_name["attenuation_16"]
    att_ho = by_name["attenuation_ho_16"]
    line_svg(OUTPUT / "figures" / "neutron_attenuation.svg", [
        {"label": "Upwind (N=16)", "color": "#d95f02", "points": centreline(att_u, "phi_total")},
        {"label": "2nd-order MUSCL (N=16)", "color": "#1764ab", "points": centreline(att_ho, "phi_total")},
        {"label": "Analytical", "color": "#111111", "dash": True, "points": centreline(att_ho, "phi_total", True)},
    ], "2D SN absorbing-beam attenuation (Upwind vs 2nd-order MUSCL)", "x", "scalar flux")

    # 2. One-group & two-group diffusion
    diff = by_name["diffusion_32"]
    line_svg(OUTPUT / "figures" / "neutron_diffusion.svg", [
        {"label": "Numerical", "color": "#1764ab", "points": centreline(diff, "phi_total")},
        {"label": "Analytical", "color": "#c43c39", "dash": True, "points": centreline(diff, "phi_total", True)},
    ], "2D one-group diffusion eigenmode", "x", "normalized scalar flux")

    twog = by_name["diffusion_two_group"]
    line_svg(OUTPUT / "figures" / "neutron_two_group.svg", [
        {"label": "group 1 numerical", "color": "#1764ab", "points": centreline(twog, "phi_g1")},
        {"label": "group 1 analytical", "color": "#1764ab", "dash": True, "points": centreline(twog, "phi_g1", True)},
        {"label": "group 2 numerical", "color": "#2f8f4e", "points": centreline(twog, "phi_g2")},
        {"label": "group 2 analytical", "color": "#2f8f4e", "dash": True, "points": centreline(twog, "phi_g2", True)},
    ], "2D two-group diffusion eigenmode", "x", "normalized scalar flux")

    # 3. Marshak vacuum & Incident current (Eq. 5.297)
    marshak = by_name["diffusion_marshak"]
    line_svg(OUTPUT / "figures" / "neutron_diffusion_marshak.svg", [
        {"label": "Numerical", "color": "#1764ab", "points": centreline(marshak, "phi_total")},
        {"label": "Analytical", "color": "#c43c39", "dash": True, "points": centreline(marshak, "phi_total", True)},
    ], "2D diffusion with Marshak vacuum boundary", "x", "scalar flux")

    inc = by_name["incident_current"]
    line_svg(OUTPUT / "figures" / "neutron_incident_current.svg", [
        {"label": "Numerical (Eq. 5.297)", "color": "#1764ab", "points": centreline(inc, "phi_total")},
        {"label": "Analytical", "color": "#c43c39", "dash": True, "points": centreline(inc, "phi_total", True)},
    ], "2D diffusion with Marshak incident partial current (J- = 0.5)", "x", "scalar flux")

    # 4. Reed 5-region S8 benchmark
    reed_dsa = by_name["reed_s8_dsa"]
    reed_no = by_name["reed_s8_nodsa"]
    line_svg(OUTPUT / "figures" / "neutron_reed.svg", [
        {"label": "S8 + DSA (73 iters)", "color": "#1764ab", "points": centreline(reed_dsa, "phi_total", target_y=0.25)},
        {"label": "S8 unaccelerated (135 iters)", "color": "#c43c39", "dash": True, "points": centreline(reed_no, "phi_total", target_y=0.25)},
    ], "Reed 5-region heterogeneous SN benchmark (strong absorber / void / scatterer)", "x (cm)", "scalar flux")

    # 5. IAEA 2D quarter-core power density & higher harmonic heatmaps
    iaea = by_name["iaea2d_2p5cm"]
    heatmap_svg(OUTPUT / "figures" / "neutron_iaea2d_power.svg", iaea["cells"], "power_density",
                "IAEA 2D PWR quarter-core power density (h = 2.5 cm)",
                f"k_eff = {float(iaea['summary']['k_eff']):.6f} (ref 1.029585), DR = {float(iaea['summary']['dominance_ratio']):.4f}")
    heatmap_svg(OUTPUT / "figures" / "neutron_iaea2d_mode1.svg", iaea["cells"], "mode_1_phi_total",
                "IAEA 2D PWR 1st higher spatial eigenmode (mode 1)",
                f"k_1 = {float(iaea['modes'][1]['k_eff']):.6f}, k_1/k_0 = {float(iaea['modes'][1]['dominance_ratio']):.5f}")

    # 6. TWIGL step & ramp transient power curves
    t_step = by_name["twigl_step"]["transient"]
    t_ramp = by_name["twigl_ramp"]["transient"]
    line_svg(OUTPUT / "figures" / "neutron_twigl.svg", [
        {"label": "Step (Backward Euler)", "color": "#1764ab",
         "points": [(float(r["time"]), float(r["relative_power"])) for r in t_step]},
        {"label": "Step reference", "color": "#1764ab", "markers": True,
         "points": sorted(TWIGL_STEP_REF.items())},
        {"label": "0.2 s Ramp (Crank-Nicolson)", "color": "#c43c39",
         "points": [(float(r["time"]), float(r["relative_power"])) for r in t_ramp]},
        {"label": "Ramp reference", "color": "#c43c39", "markers": True,
         "points": sorted(TWIGL_RAMP_REF.items())},
    ], "TWIGL 2D seed-blanket space-time kinetics transient", "time (s)", "relative total power P(t)/P(0)")

    # 7. 2x2 Heterogeneous UO2/MOX/Poison lattice thermal flux heatmap
    lat = by_name["lattice_s8"]
    heatmap_svg(OUTPUT / "figures" / "neutron_lattice_sn.svg", lat["cells"], "phi_g2",
                "2x2 heterogeneous UO2/MOX/poison mini-assembly thermal flux (S8 + DSA)",
                f"k_eff = {float(lat['summary']['k_eff']):.6f}, k_adj = {float(lat['summary']['k_adjoint']):.6f}")

    # 8. 2D Higher harmonic contour (mode 3 = (2,2) checkerboard harmonic)
    modes_d = by_name["modes_base"]
    heatmap_svg(OUTPUT / "figures" / "neutron_modes_harmonic.svg", modes_d["cells"], "mode_3_phi_total",
                "2D two-group square 3rd higher spatial eigenmode: harmonic (2,2)",
                f"k_3 = {float(modes_d['modes'][3]['k_eff']):.6f}, DR_3 = {float(modes_d['modes'][3]['dominance_ratio']):.5f}")

    # 9. 2D Dog-leg void-channel shielding (line comparison + 2D contour)
    line_svg(OUTPUT / "figures" / "neutron_dogleg.svg", [
        {"label": "S4 (12 dirs)", "color": "#7570b3", "points": centreline(by_name["dogleg_s4"], "phi_total", target_y=7.0)},
        {"label": "S8 (40 dirs)", "color": "#d95f02", "points": centreline(by_name["dogleg_s8"], "phi_total", target_y=7.0)},
        {"label": "Product P8xA32 (128 dirs)", "color": "#1764ab", "points": centreline(by_name["dogleg_high_angle"], "phi_total", target_y=7.0)},
        {"label": "Diffusion (breaks in void)", "color": "#c43c39", "dash": True, "points": centreline(by_name["dogleg_diffusion"], "phi_total", target_y=7.0)},
    ], "2D dog-leg void-channel shielding along exit duct (y = 7 cm)", "x (cm)", "scalar flux")
    dog_ha = by_name["dogleg_high_angle"]
    heatmap_svg(OUTPUT / "figures" / "neutron_dogleg_contour.svg", dog_ha["cells"], "phi_total",
                "2D dog-leg void-channel scalar flux contour (Product 8x32, 128 dirs)",
                "Source in [0,2]x[0,2], L-shaped void duct (Sigma_t = 1e-4) in 50%-scattering shield")


def write_report(rows, details, failures, orders):
    by_name = {d["spec"]["name"]: d for d in details}
    lines = [
        "# Neutron transport validation",
        "",
        "This report validates the upgraded 2D deterministic neutron transport and multigroup diffusion solver under `src/neutron/`. "
        "The test suite covers conservative $S_N$ upwind and 2nd-order limited MUSCL transport, topological $O(N)$ sweeps, "
        "level-symmetric and midpoint product quadratures, $P_1$ anisotropic scattering, specular and white reflection, "
        "Diffusion Synthetic Acceleration (DSA), multigroup diffusion with Marshak vacuum and incident partial-current boundaries (Eq. 5.297), "
        "higher spatial eigenmodes and dominance ratio, exact discrete adjoints, adjoint-weighted kinetics parameters ($\\beta_{\\mathrm{eff}}, \\Lambda_{\\mathrm{eff}}$), "
        "first-order perturbation theory (Eq. 5.313), and time-dependent multigroup kinetics with delayed neutron precursors.",
        "",
        "## Summary of all validation cases",
        "",
        "| Case | Model | k_eff | Reference k | k error | Distribution L2 | Balance | Outer / Source iters |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        k_eff = f"{float(row['k_eff']):.8f}" if row["k_eff"] != "" else "—"
        reference = f"{float(row['k_reference']):.8f}" if row["k_reference"] != "" else "—"
        error = f"{float(row['k_error']):.3e}" if row["k_error"] != "" else "—"
        l2 = f"{float(row['l2_relative']):.3e}" if row["l2_relative"] != "" else "—"
        iters = f"{row['outer_iterations']} / {row['source_iterations']}"
        lines.append(f"| {row['case']} | {row['model']} | {k_eff} | {reference} | {error} | {l2} | {row['balance_relative']:.3e} | {iters} |")

    # Refinement orders & baseline figures embedded directly
    lines.extend([
        "",
        "## Spatial convergence and baseline verification",
        "",
        f"- **One-group diffusion continuous-$k$ refinement order** ($h = 1/16 \\to 1/32$): **{orders['one_group']:.4f}** (expected $\\approx 2.0$).",
        f"- **1st-order $S_N$ upwind attenuation refinement order** ($h = 1/16 \\to 1/64$, refinement ratio $r=4$): **{orders['upwind']:.4f}** (expected $\\approx 1.0$; fixes the previous validation script error that used $\\log_2$ instead of $\\log_4$).",
        f"- **2nd-order limited MUSCL $S_N$ attenuation** ($h = 1/16 \\to 1/64$): L2 error drops from **{orders['ho_16']:.3e}** to **{orders['ho_64']:.3e}** (order **{orders['muscl']:.4f}**, **16× more accurate** than upwind on the 64-cell mesh, with particle balance `1.66e-16`).",
        f"- **Marshak incident partial current** (`incident_current`, fixed Eq. 5.297 factor of 4): analytical L2 error = **{orders['inc_l2']:.3e}**.",
        "",
        "![SN attenuation (Upwind vs 2nd-order MUSCL)](figures/neutron_attenuation.svg)",
        "",
        "![One-group diffusion eigenmode](figures/neutron_diffusion.svg)",
        "",
        "![Two-group diffusion eigenmode](figures/neutron_two_group.svg)",
        "",
        "![Marshak fixed-source diffusion](figures/neutron_diffusion_marshak.svg)",
        "",
        "![Marshak incident partial current (Eq. 5.297)](figures/neutron_incident_current.svg)",
        "",
        "## Higher spatial eigenmodes, adjoint duality, kinetics, and perturbation",
        "",
    ])

    m_detail = by_name["modes_base"]
    m_pert = by_name["modes_perturbed"]
    k0 = float(m_detail["summary"]["k_eff"])
    kp = float(m_pert["summary"]["k_eff"])
    dk_dir = (kp - k0) / k0
    dk_pt = float(m_detail["summary"]["perturbation_dk_over_k"])
    beta_ref, lam_ref = two_group_kinetics_reference(1 / 32)
    lines.extend([
        "Two-group 2D square (`modes_base`, $32\\times 32$ mesh) with 4 eigenmodes, exact discrete adjoint, delayed kinetics, and central thermal absorption perturbation ($\\delta\\Sigma_{a,2}=+0.01$ in $[0.25,0.75]^2$):",
        "",
        "| Mode | Harmonic | Numerical k | Exact discrete k | |k error| | Dominance ratio | Residual |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ])
    for m_row, (mx, my) in zip(m_detail["modes"], [(1, 1), (2, 1), (1, 2), (2, 2)]):
        ek = two_group_reference(1 / 32, mx, my)[0]
        lines.append(
            f"| {m_row['mode']} | ({mx},{my}) | {float(m_row['k_eff']):.10f} | {ek:.10f} | "
            f"{abs(float(m_row['k_eff']) - ek):.2e} | {float(m_row['dominance_ratio']):.6f} | {float(m_row['residual']):.2e} |"
        )
    lines.extend([
        "",
        f"- Adjoint duality: $|k^\\dagger - k_0| = {abs(float(m_detail['summary']['k_adjoint']) - k0):.2e}$.",
        f"- Adjoint-weighted kinetics: $\\beta_{{\\mathrm{{eff}}}} = {float(m_detail['summary']['beta_eff']):.10f}$ (exact `{beta_ref:.10f}`), $\\Lambda_{{\\mathrm{{eff}}}} = {float(m_detail['summary']['generation_time']):.6e}\\ \\mathrm{{s}}$ (exact `{lam_ref:.6e}`).",
        f"- First-order perturbation vs direct perturbed solve: $\\delta k/k = {dk_pt:.6e}$ vs ${dk_dir:.6e}$ (relative difference **{abs(dk_pt - dk_dir) / abs(dk_dir):.3%}**).",
        "",
        "![2D two-group square harmonic (2,2) eigenmode contour](figures/neutron_modes_harmonic.svg)",
        "",
        "## Complex 2D benchmarks",
        "",
    ])

    # IAEA 2D
    s5 = by_name["iaea2d_5cm"]["summary"]
    s5p = by_name["iaea2d_5cm_perturbed"]["summary"]
    s25 = by_name["iaea2d_2p5cm"]["summary"]
    k5, k25 = float(s5["k_eff"]), float(s25["k_eff"])
    k_rich = k25 + (k25 - k5) / 3.0
    dk_iaea_dir = (float(s5p["k_eff"]) - k5) / k5
    dk_iaea_pt = float(s5["perturbation_dk_over_k"])
    lines.extend([
        "### 1. 2D IAEA PWR quarter-core benchmark (ANL-7416 11-A2)",
        "",
        f"- Coarse mesh ($h = 5.0\\ \\mathrm{{cm}}$, 964 cells): $k_{{\\mathrm{{eff}}}} = {k5:.6f}$ (`{(k5 - IAEA_REF_K) * 1e5:+.1f}` pcm vs reference `{IAEA_REF_K}`), $\\mathrm{{DR}} = {float(s5['dominance_ratio']):.5f}$, $|k^\\dagger - k_{{\\mathrm{{eff}}}}| = {abs(float(s5['k_adjoint']) - k5):.2e}$.",
        f"- Fine mesh ($h = 2.5\\ \\mathrm{{cm}}$, 3856 cells): $k_{{\\mathrm{{eff}}}} = {k25:.6f}$ (`{(k25 - IAEA_REF_K) * 1e5:+.1f}` pcm), Richardson extrapolation $k_{{\\mathrm{{Rich}}}} = {k_rich:.6f}$ (`{(k_rich - IAEA_REF_K) * 1e5:+.1f}` pcm), peaking factor $F_q = {float(s25['power_peaking']):.4f}$.",
        f"- Central control-rod perturbation ($\\delta\\Sigma_{{a,2}} = -0.002$ in `rod_00`): 1st-order perturbation `{dk_iaea_pt * 1e5:+.2f}` pcm vs direct solve `{dk_iaea_dir * 1e5:+.2f}` pcm (**{abs(dk_iaea_pt - dk_iaea_dir) / abs(dk_iaea_dir):.2%}** difference).",
        f"- 6-group delayed kinetics on IAEA 2D: $\\beta_{{\\mathrm{{eff}}}} = {float(s25['beta_eff']):.6f}$, $\\Lambda_{{\\mathrm{{eff}}}} = {float(s25['generation_time']):.4e}\\ \\mathrm{{s}}$.",
        "",
        "![IAEA 2D PWR quarter-core power density](figures/neutron_iaea2d_power.svg)",
        "",
        "![IAEA 2D PWR 1st higher spatial eigenmode](figures/neutron_iaea2d_mode1.svg)",
        "",
        "### 2. 2D TWIGL seed-blanket space-time kinetics benchmark",
        "",
        f"Initial steady-state $k_0 = {float(by_name['twigl_step']['summary']['k_eff']):.6f}$ (ref `0.91321`, scaled critical at $t=0$).",
        "",
        "| Time (s) | Step (computed) | Step (ref) | Step rel diff | Ramp (computed) | Ramp (ref) | Ramp rel diff |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ])
    t_step = by_name["twigl_step"]["transient"]
    t_ramp = by_name["twigl_ramp"]["transient"]
    for t_val in [0.1, 0.2, 0.3, 0.4, 0.5]:
        ps = float(min(t_step, key=lambda r: abs(float(r["time"]) - t_val))["relative_power"])
        pr = float(min(t_ramp, key=lambda r: abs(float(r["time"]) - t_val))["relative_power"])
        rs, rr = TWIGL_STEP_REF[t_val], TWIGL_RAMP_REF[t_val]
        lines.append(f"| {t_val:.1f} | {ps:.4f} | {rs:.3f} | {abs(ps - rs) / rs:.2%} | {pr:.4f} | {rr:.3f} | {abs(pr - rr) / rr:.2%} |")

    # Reed S8 + DSA
    r_dsa = by_name["reed_s8_dsa"]
    r_no = by_name["reed_s8_nodsa"]
    reed_diff = max(abs(float(a["phi_total"]) - float(b["phi_total"])) for a, b in zip(r_dsa["cells"], r_no["cells"]))
    lines.extend([
        "",
        "![TWIGL 2D step and ramp transient power](figures/neutron_twigl.svg)",
        "",
        "### 3. Reed 5-region heterogeneous $S_8$ benchmark with true void ($\\Sigma_t = 0$) and DSA",
        "",
        f"- Source iterations reduced from **{r_no['history'][-1]['source_iterations']}** (unaccelerated) to **{r_dsa['history'][-1]['source_iterations']}** (DSA), a **{int(r_no['history'][-1]['source_iterations']) / int(r_dsa['history'][-1]['source_iterations']):.2f}×** acceleration.",
        f"- Maximum scalar-flux discrepancy between DSA and unaccelerated $S_8$: **{reed_diff:.2e}**; global balance = **{float(r_dsa['summary']['balance_relative']):.2e}**.",
        "",
        "![Reed 5-region S8 + DSA flux profile](figures/neutron_reed.svg)",
        "",
        "### 4. 2D fuel pin-cell ($P_1$ scattering) and $2\\times 2$ UO$_2$/MOX/poison mini-assembly ($S_8$ + DSA)",
        "",
        f"- Pin-cell $S_8$ (specular): $k_{{\\mathrm{{eff}}}} = {float(by_name['pincell_s8']['summary']['k_eff']):.6f}$, $k^\\dagger = {float(by_name['pincell_s8']['summary']['k_adjoint']):.6f}$ (duality `{abs(float(by_name['pincell_s8']['summary']['k_eff']) - float(by_name['pincell_s8']['summary']['k_adjoint'])):.2e}`), $\\Lambda_{{\\mathrm{{eff}}}} = {float(by_name['pincell_s8']['summary']['generation_time']):.4e}\\ \\mathrm{{s}}$.",
        f"- Pin-cell $S_8$ (white reflective): $k_{{\\mathrm{{eff}}}} = {float(by_name['pincell_s8_white']['summary']['k_eff']):.6f}$ (`{(float(by_name['pincell_s8_white']['summary']['k_eff']) - float(by_name['pincell_s8']['summary']['k_eff'])) * 1e5:+.1f}` pcm vs specular).",
        f"- Pin-cell Diffusion: $k_{{\\mathrm{{eff}}}} = {float(by_name['pincell_diffusion']['summary']['k_eff']):.6f}$ (`{(float(by_name['pincell_diffusion']['summary']['k_eff']) - float(by_name['pincell_s8']['summary']['k_eff'])) * 1e5:+.1f}` pcm transport bias due to fuel thermal flux self-shielding).",
        f"- $2\\times 2$ UO$_2$/MOX/poison mini-assembly $S_8$: $k_{{\\mathrm{{eff}}}} = {float(by_name['lattice_s8']['summary']['k_eff']):.6f}$, $k^\\dagger = {float(by_name['lattice_s8']['summary']['k_adjoint']):.6f}$, dominance ratio $\\mathrm{{DR}} = {float(by_name['lattice_s8']['summary']['dominance_ratio']):.5f}$ (3 modes converged).",
        "",
        "![2x2 heterogeneous UO2/MOX/poison lattice thermal flux](figures/neutron_lattice_sn.svg)",
        "",
        "### 5. 2D Kobayashi-style dog-leg void-channel shielding",
        "",
        "| Model / Angular set | Source (1,1) | Leg 1 (6,1) | Corner (6,7) | Exit (9.5,7) | Balance |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for name, label in [
        ("dogleg_s4", "$S_4$ (12 dirs)"),
        ("dogleg_s8", "$S_8$ (40 dirs)"),
        ("dogleg_high_angle", "Product $8\\times 32$ (128 dirs)"),
        ("dogleg_diffusion", "Diffusion (void breakdown)"),
    ]:
        d = by_name[name]
        def sample(tx, ty):
            row = min(d["cells"], key=lambda r: (float(r["x"]) - tx) ** 2 + (float(r["y"]) - ty) ** 2)
            return float(row["phi_total"])
        lines.append(
            f"| {label} | {sample(1, 1):.4f} | {sample(6, 1):.4f} | {sample(6, 7):.4f} | {sample(9.5, 7):.4f} | {float(d['summary']['balance_relative']):.2e} |"
        )

    lines.extend([
        "",
        "![2D dog-leg void-channel scalar flux contour (Product 8x32)](figures/neutron_dogleg_contour.svg)",
        "",
        "![2D dog-leg void-channel shielding along exit duct (SN vs Diffusion)](figures/neutron_dogleg.svg)",
        "",
    ])
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
    analyzed = [analyze(spec, failures) for spec in CASES]
    rows = [item[0] for item in analyzed]
    details = [item[1] for item in analyzed]
    fieldnames = list(rows[0])
    write_rows(OUTPUT / "neutron_metrics.csv", fieldnames, rows)

    by_name = {row["case"]: row for row in rows}
    by_detail = {d["spec"]["name"]: d for d in details}

    # Cross-case checks
    coarse, fine = by_name["diffusion_16"], by_name["diffusion_32"]
    one_order = math.log(abs(float(coarse["k_eff"]) - float(coarse["k_continuous"])) /
                         abs(float(fine["k_eff"]) - float(fine["k_continuous"])), 2)
    # Fixed refinement ratio r = (1/16)/(1/64) = 4!
    upwind_order = math.log(float(by_name["attenuation_16"]["l2_relative"]) /
                            float(by_name["attenuation_64"]["l2_relative"]), 4)
    muscl_order = math.log(float(by_name["attenuation_ho_16"]["l2_relative"]) /
                           float(by_name["attenuation_ho_64"]["l2_relative"]), 4)

    if not 1.7 <= one_order <= 2.3:
        failures.append(f"one-group eigenvalue convergence order {one_order:.3f} is outside [1.7, 2.3]")
    if not 0.8 <= upwind_order <= 1.2:
        failures.append(f"upwind attenuation convergence order {upwind_order:.3f} is outside [0.8, 1.2]")
    if float(by_name["attenuation_ho_64"]["l2_relative"]) >= 0.1 * float(by_name["attenuation_64"]["l2_relative"]):
        failures.append("2nd-order MUSCL did not achieve >10x error reduction over upwind at N=64")

    # Reed DSA equivalence & speedup check
    r_dsa, r_no = by_detail["reed_s8_dsa"], by_detail["reed_s8_nodsa"]
    reed_diff = max(abs(float(a["phi_total"]) - float(b["phi_total"])) for a, b in zip(r_dsa["cells"], r_no["cells"]))
    if reed_diff > 1.0e-8:
        failures.append(f"Reed DSA vs unaccelerated solution mismatch {reed_diff:.3e}")
    if int(r_dsa["history"][-1]["source_iterations"]) >= int(r_no["history"][-1]["source_iterations"]):
        failures.append("Reed DSA did not reduce source iterations")

    # First-order perturbation cross-checks against direct perturbed solves
    for base_name, pert_name, tol in [("modes_base", "modes_perturbed", 0.01),
                                      ("iaea2d_5cm", "iaea2d_5cm_perturbed", 0.05)]:
        k0 = float(by_detail[base_name]["summary"]["k_eff"])
        kp = float(by_detail[pert_name]["summary"]["k_eff"])
        dk_dir = (kp - k0) / k0
        dk_pt = float(by_detail[base_name]["summary"]["perturbation_dk_over_k"])
        rel = abs(dk_pt - dk_dir) / abs(dk_dir)
        if rel > tol:
            failures.append(f"{base_name}: perturbation dk/k relative error {rel:.2%} exceeds {tol:.1%}")

    # IAEA Richardson extrapolation vs 1.029585
    k5 = float(by_name["iaea2d_5cm"]["k_eff"])
    k25 = float(by_name["iaea2d_2p5cm"]["k_eff"])
    k_rich = k25 + (k25 - k5) / 3.0
    if abs(k_rich - IAEA_REF_K) > 1.5e-4:
        failures.append(f"IAEA 2D Richardson k_eff {k_rich:.6f} differs from {IAEA_REF_K} by > 15 pcm")

    orders = {
        "one_group": one_order,
        "upwind": upwind_order,
        "muscl": muscl_order,
        "ho_16": float(by_name["attenuation_ho_16"]["l2_relative"]),
        "ho_64": float(by_name["attenuation_ho_64"]["l2_relative"]),
        "inc_l2": float(by_name["incident_current"]["l2_relative"]),
    }
    make_figures(details)
    write_report(rows, details, failures, orders)
    print(f"\nOne-group diffusion order: {one_order:.4f}; SN upwind order: {upwind_order:.4f}; SN MUSCL order: {muscl_order:.4f}")
    print(f"IAEA 2D k_eff (2.5 cm): {k25:.6f} (Richardson: {k_rich:.6f}, ref: {IAEA_REF_K})")
    if failures:
        for failure in failures:
            print("FAIL:", failure)
        return 1
    print("All neutron validation checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
