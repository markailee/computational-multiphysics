#!/usr/bin/env python3
from collections import defaultdict
from math import cos, exp, pi, sin
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def boundary_edges(cells):
    edges = {}
    for cell in cells:
        for i, a in enumerate(cell):
            b = cell[(i + 1) % len(cell)]
            key = tuple(sorted((a, b)))
            if key in edges:
                edges[key][0] += 1
            else:
                edges[key] = [1, a, b]
    return [(a, b) for count, a, b in edges.values() if count == 1]


def write_msh(path, nodes, cells, patches):
    path.parent.mkdir(parents=True, exist_ok=True)
    patch_tags = {name: i + 1 for i, name in enumerate(patches)}
    fluid_tag = len(patch_tags) + 1
    with path.open("w", encoding="utf-8") as output:
        output.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        output.write(f"$PhysicalNames\n{len(patch_tags) + 1}\n")
        for name, tag in patch_tags.items():
            output.write(f'1 {tag} "{name}"\n')
        output.write(f'2 {fluid_tag} "fluid"\n$EndPhysicalNames\n')
        output.write(f"$Nodes\n{len(nodes)}\n")
        for i, (x, y) in enumerate(nodes, 1):
            output.write(f"{i} {x:.16g} {y:.16g} 0\n")
        output.write("$EndNodes\n")
        element_count = sum(len(faces) for faces in patches.values()) + len(cells)
        output.write(f"$Elements\n{element_count}\n")
        element = 1
        for name, faces in patches.items():
            for a, b in faces:
                output.write(f"{element} 1 2 {patch_tags[name]} {patch_tags[name]} {a + 1} {b + 1}\n")
                element += 1
        for cell in cells:
            element_type = 2 if len(cell) == 3 else 3
            connectivity = " ".join(str(node + 1) for node in cell)
            output.write(f"{element} {element_type} 2 {fluid_tag} {fluid_tag} {connectivity}\n")
            element += 1
        output.write("$EndElements\n")
    print(f"{path}: {len(nodes)} nodes, {len(cells)} cells")


def triangulate_regions(regions):
    nodes = []
    lookup = {}
    cells = []

    def node(x, y):
        key = (round(x, 12), round(y, 12))
        if key not in lookup:
            lookup[key] = len(nodes)
            nodes.append((x, y))
        return lookup[key]

    for xs, ys in regions:
        indices = [[node(x, y) for x in xs] for y in ys]
        for j in range(len(ys) - 1):
            for i in range(len(xs) - 1):
                a = indices[j][i]
                b = indices[j][i + 1]
                c = indices[j + 1][i + 1]
                d = indices[j + 1][i]
                if (i + j) % 2 == 0:
                    cells.extend(((a, b, c), (a, c, d)))
                else:
                    cells.extend(((a, b, d), (b, c, d)))
    return nodes, cells


def linspace(start, end, count):
    return [start + (end - start) * i / count for i in range(count + 1)]


def cavity():
    nodes, cells = triangulate_regions([(linspace(0.0, 1.0, 40), linspace(0.0, 1.0, 40))])
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        y = (nodes[a][1] + nodes[b][1]) / 2
        patches["top" if abs(y - 1.0) < 1e-10 else "walls"].append((a, b))
    write_msh(ROOT / "cavity" / "mesh" / "cavity.msh", nodes, cells, patches)


def backstep():
    regions = [
        (linspace(-1.0, 0.0, 12), linspace(0.5, 1.0, 12)),
        (linspace(0.0, 8.0, 96), linspace(0.0, 1.0, 24)),
    ]
    nodes, cells = triangulate_regions(regions)
    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        name = "inlet" if abs(x + 1.0) < 1e-10 else "outlet" if abs(x - 8.0) < 1e-10 else "walls"
        patches[name].append((a, b))
    write_msh(ROOT / "backstep" / "mesh" / "backstep.msh", nodes, cells, patches)


def cylinder():
    angular_count = 96
    radial_count = 24
    radius = 0.5
    xmin, xmax, ymin, ymax = -4.0, 16.0, -4.0, 4.0
    nodes = []
    for i in range(angular_count):
        theta = 2 * pi * i / angular_count
        dx, dy = cos(theta), sin(theta)
        candidates = []
        if abs(dx) > 1e-14:
            for x in (xmin, xmax):
                distance = (x - 0.0) / dx
                if distance > 0:
                    candidates.append(distance)
        if abs(dy) > 1e-14:
            for y in (ymin, ymax):
                distance = (y - 0.0) / dy
                if distance > 0:
                    candidates.append(distance)
        outer_radius = min(candidates)
        for j in range(radial_count + 1):
            fraction = (exp(2.0 * j / radial_count) - 1.0) / (exp(2.0) - 1.0)
            r = radius + (outer_radius - radius) * fraction
            nodes.append((r * dx, r * dy))

    def index(i, j):
        return (i % angular_count) * (radial_count + 1) + j

    cells = []
    for i in range(angular_count):
        for j in range(radial_count):
            cells.append((index(i, j), index(i, j + 1), index(i + 1, j + 1), index(i + 1, j)))

    patches = defaultdict(list)
    for a, b in boundary_edges(cells):
        x = (nodes[a][0] + nodes[b][0]) / 2
        ra = (nodes[a][0] ** 2 + nodes[a][1] ** 2) ** 0.5
        rb = (nodes[b][0] ** 2 + nodes[b][1] ** 2) ** 0.5
        if abs(ra - radius) < 1e-8 and abs(rb - radius) < 1e-8:
            name = "cylinder"
        elif abs(nodes[a][0] - xmin) < 1e-8 and abs(nodes[b][0] - xmin) < 1e-8:
            name = "inlet"
        elif abs(nodes[a][0] - xmax) < 1e-8 and abs(nodes[b][0] - xmax) < 1e-8:
            name = "outlet"
        else:
            name = "topBottom"
        patches[name].append((a, b))
    write_msh(ROOT / "cylinder" / "mesh" / "cylinder.msh", nodes, cells, patches)


if __name__ == "__main__":
    cavity()
    backstep()
    cylinder()
