#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The SEMFE Heat Transfer Solver
Computational Mechanics

Main Script
"""
import numpy as np
from PreProcessor import read_input_file
from Solver import assemble_global, apply_convection, apply_dirichlet
from Solver import apply_dirichlet_penalty1, apply_dirichlet_penalty2, apply_dirichlet_penalty3, apply_dirichlet_penalty4
from Solver import apply_heat_flux, solve_system
from PostProcessor import plot_mesh, plot_mesh_interactive, plot_temperature_field
from PostProcessor import export_temperature_csv


# Import model info
nodes, elems, materials, k, bcs = read_input_file('Part 2 - Chimney/chimney.semfe')

# Check Mesh Quality
plot_mesh_interactive(nodes, elems, show=True, filename='Part 2 - Chimney/interactive_mesh_chimney.html')

# Assemble global
K = assemble_global(nodes, elems, k=k)
f0 = np.zeros(nodes.shape[0])

# Apply BCs
# APPLY CONVECTION (Robin)
conv_bcs = bcs.get("convection", [])
Kc, fc = apply_convection(K, f0, nodes, elems, conv_bcs)
# APPLY HEAT FLUX (Neumann)
flux_bcs = bcs.get("heat_flux", [])
f2 = apply_heat_flux(fc, nodes, elems, flux_bcs)
# APPLY DIRICHLET (Strong)
temp_bcs = bcs.get("temperature", [])
bc_nodes = [node for node, val in temp_bcs]
bc_values = [val  for node, val in temp_bcs]

Kmod, fmod = apply_dirichlet(Kc, f2, bc_nodes, bc_values)

# Solve
u = solve_system(Kmod, fmod)

# Call it in main
plot_temperature_field(nodes, elems, u, filename='Part 2 - Chimney/temperature_field.png')
export_temperature_csv(nodes, u)
