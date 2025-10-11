# 🌀 Christoffel Symbol Calculator

A symbolic calculator built with **SymPy** and **Streamlit** to compute Christoffel symbols for arbitrary spacetime metrics in General Relativity.

## Features
- Compute Christoffel symbols from any metric tensor
- Predefined GR metrics: Schwarzschild, Kerr, FLRW, Minkowski, etc.
- Supports custom metric input
- Symbolic simplification and LaTeX output
- Easy web UI via Streamlit

## Example Metric
Schwarzschild 4D metric:
```python
[[-(1-2*M/r), 0, 0, 0],
 [0, 1/(1-2*M/r), 0, 0],
 [0, 0, r**2, 0],
 [0, 0, 0, r**2*sin(theta)**2]]
