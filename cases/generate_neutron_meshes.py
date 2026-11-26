#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path

from generate_meshes import boundary_edges, linspace, quad_regions, triangulate_regions, write_msh

ROOT = Path(__file__).resolve().parent


def square(case, name, nx, ny, triangles=False):
    regions = [(linspace(0.0, 1.0, nx), linspace(0.0, 1.0, ny))]
    nodes, cells = triangulate_regions(regions) if triangles else quad_regions([
        [[(x, y) for x in regions[0][0]] for y in regions[0][1]]
    ])
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        y = (nodes[a][1] + nodes[b][1]) / 2
        if abs(x) < 1e-12:
            patch = "left"
        elif abs(x - 1.0) < 1e-12:
            patch = "right"
        elif abs(y) < 1e-12:
            patch = "bottom"
        else:
            patch = "top"
        patches[patch].append((a, b))
    write_msh(ROOT / case / "mesh" / name, nodes, cells, patches)


if __name__ == "__main__":
    square("neutron_attenuation", "quad_16.msh", 16, 4)
    square("neutron_attenuation", "quad_64.msh", 64, 8)
    square("neutron_attenuation", "triangles_32.msh", 32, 8, triangles=True)
    square("neutron_sn_critical", "quad_16.msh", 16, 16)
    square("neutron_diffusion", "quad_16.msh", 16, 16)
    square("neutron_diffusion", "quad_32.msh", 32, 32)
    square("neutron_diffusion_two_group", "quad_32.msh", 32, 32)
    square("neutron_diffusion_marshak", "quad_32.msh", 32, 8)
