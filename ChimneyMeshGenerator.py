# -*- coding: utf-8 -*-
"""
Created on Sun Nov 16 17:27:28 2025

@author: Micha
"""

import xml.etree.ElementTree as ET

# ===========================
# Mesh parameters
# ===========================
Lx = 0.8
Ly = 0.6

# Mesh resolution
nx = 20   # number of divisions x
ny = 20   # number of divisions y

dx = Lx / nx
dy = Ly / ny

# Hole dimensions
hx1, hx2 = 0.2, 0.6
hy1, hy2 = 0.2, 0.4

# ===========================
# Generate nodes
# ===========================
nodes = []
node_id = {}
counter = 1

for j in range(ny+1):
    for i in range(nx+1):
        x = i * dx
        y = j * dy
        node_id[(i, j)] = counter
        nodes.append((counter, x, y, 0.0))
        counter += 1

# ===========================
# Generate triangular elements
# ===========================
elements = []
eid = 1

for j in range(ny):
    for i in range(nx):

        # Cell corners (bottom-left)
        x0 = i * dx
        y0 = j * dy
        x1 = (i+1) * dx
        y1 = (j+1) * dy

        # Skip if cell is fully inside the hole
        if (x0 >= hx1 and x1 <= hx2 and
            y0 >= hy1 and y1 <= hy2):
            continue

        # Node indices
        n1 = node_id[(i, j)]
        n2 = node_id[(i+1, j)]
        n3 = node_id[(i, j+1)]
        n4 = node_id[(i+1, j+1)]

        # Two triangles per cell
        elements.append((eid, n1, n2, n4)); eid += 1
        elements.append((eid, n1, n4, n3)); eid += 1

# ===========================
# Build XML structure
# ===========================
root = ET.Element("SEMFE_spec")
ET.SubElement(root, "Module", type="heat conduction")

materials = ET.SubElement(root, "Materials")
mat = ET.SubElement(materials, "Material", id="1", name="Steel")
ET.SubElement(mat, "conductivity").text = "1.5"

geom = ET.SubElement(root, "Geometry")

# ---- NODES ----
nodes_xml = ET.SubElement(geom, "Nodes")
for (nid, x, y, z) in nodes:
    ET.SubElement(nodes_xml, "node",
                  id=str(nid),
                  x=str(x), y=str(y), z=str(z))

# ---- ELEMENTS ----
elems_xml = ET.SubElement(geom, "Elements", type="tri3", name="mesh")
for (eid, a, b, c) in elements:
    el = ET.SubElement(elems_xml, "elem", id=str(eid))
    el.text = f"{a} {b} {c}"

# ===========================
# Boundary Conditions (INCLUDES HOLE = 100°C)
# ===========================
bc = ET.SubElement(root, "BoundaryConditions")

# ---- Dirichlet Boundary section ----
boundary = ET.SubElement(bc, "Boundary")

count_inner_bc = 0
count_outer_top = 0

for (nid, x, y, z) in nodes:

    # --- Inner hole at 100°C ---
    on_left   = abs(x - 0.2) < 1e-9 and (0.2 <= y <= 0.4)
    on_right  = abs(x - 0.6) < 1e-9 and (0.2 <= y <= 0.4)
    on_bottom = abs(y - 0.2) < 1e-9 and (0.2 <= x <= 0.6)
    on_top    = abs(y - 0.4) < 1e-9 and (0.2 <= x <= 0.6)

    if on_left or on_right or on_bottom or on_top:
        ET.SubElement(boundary, "temperature",
                      node=str(nid), value="100.0")
        count_inner_bc += 1
        continue  # avoid double counting

    # --- TOP outer boundary at 30°C (y = 0.6) ---
    if abs(y - 0.6) < 1e-9:
        ET.SubElement(boundary, "temperature",
                      node=str(nid), value="30.0")
        count_outer_top += 1

print("Inner hole BC nodes:", count_inner_bc)
print("Top boundary T=30°C nodes:", count_outer_top)


# ---- Convection at RIGHT WALL (x = 0.8) ----
convection = ET.SubElement(bc, "Convection")

for (nid, x, y, z) in nodes:
    if abs(x - 0.8) < 1e-9:     # right outer edge
        ET.SubElement(convection, "conv",
                      node=str(nid),
                      h="50.0",
                      Tinf="25.0")


# ---- Heat flux q = 0 at BOTTOM (y = 0) ----
heatflux = ET.SubElement(bc, "HeatFlux")

for (nid, x, y, z) in nodes:
    if abs(y - 0.0) < 1e-9:     # bottom edge
        ET.SubElement(heatflux, "flux",
                      node=str(nid),
                      value="0.0")

# ===========================
# Step section
# ===========================
step = ET.SubElement(root, "Step", name="step1", type="steady-state")
ET.SubElement(step, "HeatSource")

# ===========================
# WRITE XML FILE
# ===========================
tree = ET.ElementTree(root)
tree.write("chimney.semfe",
           encoding="ISO-8859-1",
           xml_declaration=True)

print("Mesh generated → chimney.semfe")
print("Total nodes:", len(nodes))
print("Total elements:", len(elements))
