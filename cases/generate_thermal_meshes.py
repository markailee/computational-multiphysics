#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path

from generate_meshes import boundary_edges, quad_regions, write_msh

ROOT = Path(__file__).resolve().parent


def cavity(level, n):
    grid = [[(i / n, j / n) for i in range(n + 1)] for j in range(n + 1)]
    nodes, cells = quad_regions([grid])
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        y = (nodes[a][1] + nodes[b][1]) / 2
        if abs(x) < 1e-12:
            name = "hot"
        elif abs(x - 1) < 1e-12:
            name = "cold"
        elif abs(y) < 1e-12:
            name = "bottom"
        else:
            name = "top"
        patches[name].append((a, b))
    write_msh(ROOT / "boussinesq_cavity" / "mesh" / f"cavity_{level}.msh", nodes, cells, patches)


def conduction():
    nx, ny = 40, 2
    grid = [[(0.1 * i / nx, 0.01 * j / ny) for i in range(nx + 1)] for j in range(ny + 1)]
    nodes, cells = quad_regions([grid])
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        name = "left" if abs(x) < 1e-12 else "right" if abs(x - 0.1) < 1e-12 else "walls"
        patches[name].append((a, b))
    write_msh(ROOT / "compressible_conduction" / "mesh" / "conduction.msh", nodes, cells, patches)


def nozzle(level, nx, ny):
    length = 4.0

    def half_height(x):
        return 0.1 - 0.0125 * x

    grid = [[(length * i / nx, half_height(length * i / nx) * (2 * j / ny - 1))
             for i in range(nx + 1)] for j in range(ny + 1)]
    nodes, cells = quad_regions([grid])
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        name = "inlet" if abs(x) < 1e-12 else "outlet" if abs(x - length) < 1e-12 else "walls"
        patches[name].append((a, b))
    write_msh(ROOT / "compressible_nozzle" / "mesh" / f"nozzle_{level}.msh", nodes, cells, patches)


if __name__ == "__main__":
    for level, n in (("coarse", 16), ("medium", 24), ("fine", 32)):
        cavity(level, n)
    conduction()
    for level, nx, ny in (("coarse", 40, 4), ("medium", 60, 6), ("fine", 80, 8)):
        nozzle(level, nx, ny)
