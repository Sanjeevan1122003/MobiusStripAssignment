# MobiusStripAssignment

# Short write-up

**Code Structure:**
The code is organised around a MobiusStrip class that encapsulates all relevant properties and methods for generating and analysing a Möbius strip surface. The constructor initialises parameters such as radius (R), width (w), and resolution (n), and computes the mesh grid coordinates.

The project is structured around a `MobiusStrip` class:
**Initialization (`__init__`)**: 
  Sets parameters — radius (`R`), width (`w`), and resolution (`n`). Generates mesh grids for parameters `u` and `v`.
  
**Mesh Computation (`compute_mesh`)**: 
  Calculates the 3D coordinates (`X`, `Y`, `Z`) of the Möbius strip surface.
  
**Surface Area Calculation (`compute_surface_area`)**: 
  Approximates surface area by calculating partial derivatives of the surface, computing the magnitude of their cross product (the local area element), and 
  numerically integrating over the parameter domain.
  
**Edge Length Calculation (`compute_edge_length`)**: 
  Calculates the length of one edge by sampling points along `v = w/2` and integrating the arc length.
  
**Visualisation (`plot`)**: 
  Plots the Möbius strip in 3D and displays computed parameters in a boxed text label at the top right of the figure.

---

**Surface Area Approximation:**

The Möbius strip surface is parametrised by two parameters, `u` (around the loop) and `v` (across the width). The surface area is approximated by:

1. Computing the partial derivatives of the surface with respect to `u` and `v`.
2. Taking the cross product of these derivative vectors at each mesh point to get the local surface normal vector.
3. Calculating the norm (magnitude) of the normal vector, which corresponds to the differential surface element area `dA`.
4. Numerically integrating these values over the parameter space using Simpson's rule to sum the total surface area.

---

**Challenges Faced:**

**Numerical Integration Accuracy:**
  Balancing the resolution of the mesh to achieve accurate numerical integration without excessive computation time.
  
**Correct Derivative Formulation:** 
  Deriving and implementing the parametric partial derivatives correctly was essential to compute valid surface normals.
  
**Text Box Positioning:** 
  Properly positioning the parameter text box inside the plot area, specifically at the top right corner, required using normalized figure coordinates and alignment options.\
  
**Edge Length Computation:** 
  Calculating the length of a nontrivial edge of the Möbius strip involved careful gradient calculations and integration along a parametric curve.

---
