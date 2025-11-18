# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 17:27:28 2025

@author: Micha

Final CHIMNEY mesh + XML creator
with inner hole, triangular mesh, and EXACT XML formatting.
"""

import math
from collections import defaultdict


# ---------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------
OUTER_W = 0.8
OUTER_H = 0.6

# hole (inner rectangle)
HOLE_X1, HOLE_X2 = 0.2, 0.6
HOLE_Y1, HOLE_Y2 = 0.2, 0.4

# BC values
T_INNER = 100.0
T_TOP = 30.0
H_CONV = 50.0
T_INF = 25.0


def generate_chimney_mesh(nx=4, ny=3):
    """
    Builds a structured rectangular grid,
    removes the internal 0.2×0.2 region,
    and splits each cell into 2 triangles.
    """

    dx = OUTER_W / nx
    dy = OUTER_H / ny

    n_nodes_x = nx + 1
    n_nodes_y = ny + 1

    # ------------------------
    # NODES
    # ------------------------
    nodes = []
    for j in range(n_nodes_y):
        for i in range(n_nodes_x):
            nodes.append((i * dx, j * dy))

    def nid(i, j):
        return j * n_nodes_x + i

    # ------------------------
    # TRIANGLES
    # ------------------------
    elems = []
    for j in range(ny):
        for i in range(nx):

            # center of cell
            xc = (i + 0.5) * dx
            yc = (j + 0.5) * dy

            # skip if inside hole
            if (HOLE_X1 < xc < HOLE_X2) and (HOLE_Y1 < yc < HOLE_Y2):
                continue

            # corner node indices
            n1 = nid(i, j)
            n2 = nid(i+1, j)
            n3 = nid(i+1, j+1)
            n4 = nid(i, j+1)

            # triangles with MATLAB orientation
            elems.append((n1, n2, n4))
            elems.append((n2, n3, n4))

    return nodes, elems


# ---------------------------------------------------------
# FIND BOUNDARY EDGES
# ---------------------------------------------------------
def find_boundary_edges(nodes, elems):
    edge_map = defaultdict(list)

    for e_idx, (n1, n2, n3) in enumerate(elems):
        # local edges (with local ID 1–3)
        edges = [
            (n1, n2, 1),
            (n2, n3, 2),
            (n3, n1, 3),
        ]
        for a, b, local in edges:
            key = (min(a, b), max(a, b))
            edge_map[key].append((e_idx, local))

    boundary_edges = []
    for (a, b), owners in edge_map.items():
        if len(owners) == 1:  # belongs to only one triangle
            e_idx, local = owners[0]
            boundary_edges.append((a, b, e_idx, local))

    return boundary_edges


# ---------------------------------------------------------
# CLASSIFY BOUNDARY
# ---------------------------------------------------------
def classify_edges(nodes, bnd):
    tol = 1e-10
    bottom = []
    right_ext = []
    top_ext = []
    left_ext = []
    inner = []

    for (a, b, e_idx, loc) in bnd:
        x1, y1 = nodes[a]
        x2, y2 = nodes[b]

        # BOTTOM y = 0
        if abs(y1) < tol and abs(y2) < tol:
            bottom.append((a, b, e_idx, loc))
        # RIGHT x = W
        elif abs(x1 - OUTER_W) < tol and abs(x2 - OUTER_W) < tol:
            right_ext.append((a, b, e_idx, loc))
        # TOP y = H
        elif abs(y1 - OUTER_H) < tol and abs(y2 - OUTER_H) < tol:
            top_ext.append((a, b, e_idx, loc))
        # LEFT x = 0  (insulated)
        elif abs(x1) < tol and abs(x2) < tol:
            left_ext.append((a, b, e_idx, loc))
        else:
            # INNER hole edges
            inner.append((a, b, e_idx, loc))

    return bottom, right_ext, top_ext, left_ext, inner


# ---------------------------------------------------------
# REMOVE UNUSED NODES + REINDEX
# ---------------------------------------------------------
def cleanup(nodes, elems, bottom, right, top, left, inner):
    used = set()
    for n1, n2, n3 in elems:
        used.update([n1, n2, n3])

    used = sorted(list(used))
    mapping = {old: new for new, old in enumerate(used)}

    new_nodes = [nodes[i] for i in used]

    new_elems = [(mapping[a], mapping[b], mapping[c]) for (a, b, c) in elems]

    def remap(edges):
        return [(mapping[a], mapping[b], e_idx, local) for (a, b, e_idx, local) in edges]

    return (
        new_nodes,
        new_elems,
        remap(bottom),
        remap(right),
        remap(top),
        remap(left),
        remap(inner)
    )


# ---------------------------------------------------------
# WRITE XML EXACTLY AS REQUIRED
# ---------------------------------------------------------
def write_xml(filename, nodes, elems, bottom, right, top, inner):
    with open(filename, "w", encoding="ISO-8859-1") as f:

        f.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        f.write('<SEMFE_spec>\n')
        f.write('  <Module type="heat conduction"/>\n\n')

        # MATERIAL
        f.write('  <Materials>\n')
        f.write('    <Material id="1" name="Brick">\n')
        f.write('      <conductivity>1.0</conductivity>\n')
        f.write('    </Material>\n')
        f.write('  </Materials>\n\n')

        # GEOMETRY
        f.write('  <Geometry>\n')
        f.write('    <Nodes>\n')
        for i, (x, y) in enumerate(nodes, start=1):
            f.write(f'      <node id="{i}" x="{x}" y="{y}" z="0.0"/>\n')
        f.write('    </Nodes>\n\n')

        f.write('    <Elements type="tri3" name="mesh">\n')
        for e_id, (a, b, c) in enumerate(elems, start=1):
            f.write(f'      <elem id="{e_id}">{a+1} {b+1} {c+1}</elem>\n')
        f.write('    </Elements>\n')
        f.write('  </Geometry>\n\n')

        # BOUNDARY CONDITIONS
        f.write('  <BoundaryConditions>\n')

        # ---- DIRICHLET ----
        f.write('    <Boundary>\n')

        # inner T=100
        inner_nodes = sorted({a for a,b,_,_ in inner} | {b for a,b,_,_ in inner})
        for n in inner_nodes:
            f.write(f'      <temperature node="{n+1}" value="{T_INNER}"/>\n')

        # top T=30
        top_nodes = sorted({a for a,b,_,_ in top} | {b for a,b,_,_ in top})
        for n in top_nodes:
            f.write(f'      <temperature node="{n+1}" value="{T_TOP}"/>\n')

        f.write('    </Boundary>\n\n')

        # ---- HEAT FLUX bottom = 0 ----
        f.write('    <HeatFlux>\n')
        for (a, b, e_idx, local) in bottom:
            f.write(f'      <flux elem="{e_idx+1}" edge="{local}" value="0.0"/>\n')
        f.write('    </HeatFlux>\n\n')

        # ---- CONVECTION at right ----
        f.write('    <Convection>\n')
        for (a, b, e_idx, local) in right:
            f.write(f'      <conv elem="{e_idx+1}" edge="{local}" h="{H_CONV}" Tinf="{T_INF}"/>\n')
        f.write('    </Convection>\n')

        f.write('  </BoundaryConditions>\n\n')

        f.write('  <Step name="step1" type="steady-state">\n')
        f.write('    <HeatSource>\n')
        f.write('    </HeatSource>\n')
        f.write('  </Step>\n')
        f.write('</SEMFE_spec>\n')

    print("XML saved:", filename)


# ---------------------------------------------------------
# MAIN
# ---------------------------------------------------------
if __name__ == "__main__":
    # Triangular mesh
    nodes, elems = generate_chimney_mesh(nx=20, ny=15)

    # boundaries
    bnd = find_boundary_edges(nodes, elems)
    bottom, right, top, left, inner = classify_edges(nodes, bnd)

    # cleanup & reindex
    nodes, elems, bottom, right, top, left, inner = cleanup(
        nodes, elems, bottom, right, top, left, inner
    )

    # WRITE XML
    write_xml("chimney.semfe", nodes, elems, bottom, right, top, inner)

    print("Nodes:", len(nodes))
    print("Elements:", len(elems))
