# FEM Code for Heat Transfer on 2-D structures #

Finite Element Method (FEM) framework in Python for simulating steady-state heat transfer in 2D structures, built as a *triangular* mesh.
<!--The core implementation is based on a fork of the work by SotirisKak (original repository: https://github.com/SotirisKak/comp_mechII_Project1)-->

## Introduction ##
The code defines nodes and elements, assembles the global conductivity matrix, applies boundary conditions, and computes the temperature distribution across the domain. 
As **requirements**, we need a Python editor (such as Spyder).

The pipeline receives as input an XML file, which contains:
<!--- Conductivity constant
- Element connectivity
- Dirichlet, Newmann, and Robin BC-->

saved as .semfe

#### FEM Heat Transfer Theory ####


## Input file ##
For the first part of this project (Part 1 - Box), the input file was created manually. 

For a more complex triangular mesh, the file **ChimneyMeshGenerator.py** (located in the Part 2 – Chimney folder) generates a triangular mesh for a square domain.
This script was used for the second part of the project, so it includes parameters for creating a central hole - if required.

In the mesh generator script, as parameters, we set:
- The length and height of the square domain
- *plus the size of the hole*
- The number of divisions on the x and y axes
- The heat conductivity constant
- and the BC values.

## Main Script ##

## Outputs ##
