#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path

from generate_meshes import boundary_edges, linspace, quad_regions, triangulate_regions, write_msh

ROOT = Path(__file__).resolve().parent


def rectangle(case, name, lx, ly, nx, ny, triangles=False):
    regions = [(linspace(0.0, lx, nx), linspace(0.0, ly, ny))]
    nodes, cells = triangulate_regions(regions) if triangles else quad_regions([
        [[(x, y) for x in regions[0][0]] for y in regions[0][1]]
    ])
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        y = (nodes[a][1] + nodes[b][1]) / 2
        if abs(x) < 1e-12:
            patch = "left"
        elif abs(x - lx) < 1e-12:
            patch = "right"
        elif abs(y) < 1e-12:
            patch = "bottom"
        else:
            patch = "top"
        patches[patch].append((a, b))
    write_msh(ROOT / case / "mesh" / name, nodes, cells, patches)


def square(case, name, nx, ny, triangles=False):
    rectangle(case, name, 1.0, 1.0, nx, ny, triangles=triangles)


def iaea_2d_mesh(case, name, cell_size):
    """Generate the exact stepped 170 cm x 170 cm 2D IAEA PWR quarter-core mesh (ANL-7416 11-A2)."""
    # Quarter-core assembly bounds along x and y: 0, 10, 30, 50, 70, 90, 110, 130, 150, 170
    # Active blocks (i, j) where material != 0 (row j=0..8, col i=0..8):
    # Row 0..2: i=0..8 (x <= 170)
    # Row 3:    i=0..8 (x <= 170)  -- wait: check row 3 col 8: (50..70, 150..170) is reflector (4)!
    # Row 4..5: i=0..7 (x <= 150)
    # Row 6:    i=0..6 (x <= 130)
    # Row 7:    i=0..5 (x <= 110)
    # Row 8:    i=0..3 (x <= 70)
    edges = [0.0, 10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0]
    row_max_col = [9, 9, 9, 9, 8, 8, 7, 6, 4]  # number of active blocks in each row j

    # Build uniform coordinate arrays per block
    grids = []
    for j in range(9):
        y0, y1 = edges[j], edges[j + 1]
        ny = max(1, int(round((y1 - y0) / cell_size)))
        ys = linspace(y0, y1, ny)
        for i in range(row_max_col[j]):
            x0, x1 = edges[i], edges[i + 1]
            nx = max(1, int(round((x1 - x0) / cell_size)))
            xs = linspace(x0, x1, nx)
            grids.append([[(x, y) for x in xs] for y in ys])

    nodes, cells = quad_regions(grids)
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        y = (nodes[a][1] + nodes[b][1]) / 2
        if abs(x) < 1e-10:
            patch = "left"
        elif abs(y) < 1e-10:
            patch = "bottom"
        else:
            patch = "outer"
        patches[patch].append((a, b))
    write_msh(ROOT / case / "mesh" / name, nodes, cells, patches)


def twigl_mesh(case, name, h):
    """Generate the 80 cm x 80 cm TWIGL quarter-core mesh aligned with 0, 24, 56, 80 cm interfaces."""
    breaks = [0.0, 24.0, 56.0, 80.0]
    xs = [0.0]
    for a, b in zip(breaks[:-1], breaks[1:]):
        n = max(1, int(round((b - a) / h)))
        xs.extend(linspace(a, b, n)[1:])
    ys = list(xs)
    grid = [[(x, y) for x in xs] for y in ys]
    nodes, cells = quad_regions([grid])
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        y = (nodes[a][1] + nodes[b][1]) / 2
        if abs(x) < 1e-10:
            patch = "left"
        elif abs(y) < 1e-10:
            patch = "bottom"
        elif abs(x - 80.0) < 1e-10:
            patch = "right"
        else:
            patch = "top"
        patches[patch].append((a, b))
    write_msh(ROOT / case / "mesh" / name, nodes, cells, patches)


if __name__ == "__main__":
    # Baseline meshes
    square("neutron_attenuation", "quad_16.msh", 16, 4)
    square("neutron_attenuation", "quad_64.msh", 64, 8)
    square("neutron_attenuation", "triangles_32.msh", 32, 8, triangles=True)
    square("neutron_sn_critical", "quad_16.msh", 16, 16)
    square("neutron_diffusion", "quad_16.msh", 16, 16)
    square("neutron_diffusion", "quad_32.msh", 32, 32)
    square("neutron_diffusion_two_group", "quad_32.msh", 32, 32)
    square("neutron_diffusion_marshak", "quad_32.msh", 32, 8)

    # New complex benchmark meshes
    square("neutron_incident_current", "quad_32.msh", 32, 8)
    square("neutron_modes", "quad_32.msh", 32, 32)
    rectangle("neutron_reed", "reed_160.msh", 8.0, 1.0, 160, 2)
    iaea_2d_mesh("neutron_iaea2d", "iaea_5cm.msh", 5.0)
    iaea_2d_mesh("neutron_iaea2d", "iaea_2p5cm.msh", 2.5)
    twigl_mesh("neutron_twigl", "twigl_2cm.msh", 2.0)
    rectangle("neutron_pincell_sn", "pincell_24.msh", 1.26, 1.26, 24, 24)
    rectangle("neutron_lattice_sn", "lattice_32.msh", 2.52, 2.52, 32, 32)
    rectangle("neutron_dogleg_sn", "dogleg_40.msh", 10.0, 10.0, 40, 40)
