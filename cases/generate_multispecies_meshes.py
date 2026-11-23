#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path

from generate_meshes import boundary_edges, quad_regions, write_msh

ROOT = Path(__file__).resolve().parent


def tube(case, mesh_name, nx, ny=1, length=0.1, width=0.004, family="stefan"):
    grid = [[(length * i / nx, width * j / ny) for i in range(nx + 1)] for j in range(ny + 1)]
    nodes, cells = quad_regions([grid])
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        if family == "air":
            name = "walls"
        elif abs(x) < 1e-12:
            name = "surface" if family == "stefan" else "inlet"
        elif family == "stefan" and abs(x - length) < 1e-12:
            name = "outlet"
        elif family == "stefan":
            name = "sides"
        else:
            name = "walls"
        patches[name].append((a, b))
    write_msh(ROOT / case / "mesh" / mesh_name, nodes, cells, patches)


if __name__ == "__main__":
    for level, nx in (("coarse", 20), ("medium", 40), ("fine", 80)):
        tube("species_stefan_tube", f"stefan_{level}.msh", nx, ny=2)
        tube("species_flux_correction", f"correction_{level}.msh", nx, ny=2)
    for level, nx in (("coarse", 50), ("medium", 100), ("fine", 200)):
        tube("species_diffusion", f"diffusion_{level}.msh", nx, family="diffusion")
    tube("species_air_properties", "air.msh", 4, length=0.01, width=0.002, family="air")
