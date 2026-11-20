#!/usr/bin/env python3
import math
from collections import defaultdict
from pathlib import Path

from generate_meshes import boundary_edges, write_msh

ROOT = Path(__file__).resolve().parent


def linspace(start, end, count):
    return [start + (end - start) * i / count for i in range(count + 1)]


def rectangle(path, length, height, nx, ny, composition="mixed", distort=True, xmin=0.0, ymin=0.0):
    if nx % 2 or ny % 2:
        raise ValueError("Acoustic validation grids require even node-cell counts")
    xs = linspace(xmin, xmin + length, nx)
    ys = linspace(ymin, ymin + height, ny)
    nodes = []
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            if distort and 0 < i < nx and 0 < j < ny:
                envelope = math.sin(math.pi * (x - xmin) / length) * math.sin(math.pi * (y - ymin) / height)
                centreline = math.sin(2.0 * math.pi * (y - ymin) / height)
                x += 0.12 * length / nx * envelope * centreline
                y += 0.10 * height / ny * envelope * centreline * math.sin(2.0 * math.pi * (x - xmin) / length)
            nodes.append((x, y))

    def node(i, j):
        return j * (nx + 1) + i

    cells = []
    for j in range(ny):
        for i in range(nx):
            a = node(i, j)
            b = node(i + 1, j)
            c = node(i + 1, j + 1)
            d = node(i, j + 1)
            triangular = composition == "triangles" or (composition == "mixed" and i >= nx // 2)
            if triangular:
                if (i + j) % 2:
                    cells.extend(((a, b, d), (b, c, d)))
                else:
                    cells.extend(((a, b, c), (a, c, d)))
            else:
                cells.append((a, b, c, d))

    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        y = (nodes[a][1] + nodes[b][1]) / 2
        if abs(x - xmin) < 1e-10:
            name = "left"
        elif abs(x - xmin - length) < 1e-10:
            name = "right"
        elif abs(y - ymin) < 1e-10:
            name = "bottom"
        else:
            name = "top"
        patches[name].append((a, b))
    write_msh(path, nodes, cells, patches)


def main():
    cavity = ROOT / "acoustic_cavity" / "mesh"
    for level, nx, ny in (("coarse", 20, 10), ("medium", 40, 20), ("fine", 80, 40)):
        rectangle(cavity / f"{level}_mixed.msh", 2.0, 1.0, nx, ny)
    rectangle(cavity / "medium_triangles.msh", 2.0, 1.0, 40, 20, "triangles")
    rectangle(cavity / "medium_quadrilaterals.msh", 2.0, 1.0, 40, 20, "quadrilaterals")

    for directory in ("acoustic_plane_wave", "acoustic_convected"):
        base = ROOT / directory / "mesh"
        for level, nx, ny in (("coarse", 40, 8), ("medium", 80, 16), ("fine", 160, 32)):
            rectangle(base / f"{level}_mixed.msh", 2.0, 0.5, nx, ny)

    pulse = ROOT / "acoustic_absorbing" / "mesh"
    rectangle(pulse / "mixed.msh", 4.0, 1.0, 80, 20)
    rectangle(pulse / "fine.msh", 4.0, 1.0, 160, 40)

    aeolian = ROOT / "cylinder_aeolian" / "mesh"
    rectangle(aeolian / "acoustic.msh", 40.0, 30.0, 80, 60, distort=False, xmin=-20.0, ymin=-15.0)


if __name__ == "__main__":
    main()
