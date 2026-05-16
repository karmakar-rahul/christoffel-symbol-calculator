# Christoffel Symbol Calculator  v2

A symbolic calculator built with **SymPy** and **Streamlit** that computes
Christoffel symbols and curvature tensors for arbitrary spacetime metrics
in General Relativity.

**Live app:** https://christoffel-symbol-calculator-fzqbkxctohbhcnvvnjnpfz.streamlit.app/

---

## What's new in v2

### Computation engine (`christoffel.py`)

| Improvement | Detail |
|---|---|
| **Faster metric inversion** | Switched from `Matrix.inv()` (Gaussian elimination) to `Matrix.inv("ADJ")` (adjugate/determinant). Block-diagonal and sparse metrics (Schwarzschild, Kerr, FLRW) invert 2–5× faster. |
| **Derivative cache** | Each `∂g_{mn}/∂x^k` is computed once and reused across all Christoffel components, eliminating redundant `sp.diff` calls. |
| **Zero-skipping** | Components whose unsimplified expression is identically zero are never passed through the simplifier — the most expensive step for Kerr/RN. |
| **Parallel simplification** | Optional `ThreadPoolExecutor` simplifies independent components concurrently. Enabled by default in the UI. |
| **Extended outputs** | New functions: `compute_riemann`, `compute_ricci_tensor`, `compute_ricci_scalar`, `compute_einstein_tensor`, and the convenience wrapper `compute_all_tensors`. |
| **Lower-index symmetry** | Exploits Γ^λ_{μν} = Γ^λ_{νμ} so only `n(n+1)/2` components per upper index are computed (unchanged from v1, confirmed correct). |
| **Fully backward-compatible** | v1 callers (`compute_christoffel(g, coords)`) require no changes. |

### UI (`app.py`)

- **Tabbed output** — Christoffel | Riemann | Ricci tensor | Ricci scalar | Einstein tensor
- **`st.cache_data`** — re-running the same metric is instant (result cached in session)
- **Non-zero count banner** — immediately shows how many independent components are non-zero
- **LaTeX download** — export Christoffel symbols as a ready-to-compile `.tex` file
- **Raw SymPy view** toggle — show unsimplified Python expressions for debugging
- **Simplification slider** — None / Basic / Full; slider replaces radio buttons for clarity
- **Parallel toggle** — enable/disable thread pool from the sidebar
- **Expanded metric library** — added de Sitter, Minkowski (spherical), Schwarzschild (3D)
- **Metric preview** — collapsible expander shows the parsed `g_{μν}` as LaTeX before results

---

## Project structure

```
christoffel-symbol-calculator/
├── app.py            # Streamlit UI (v2)
├── christoffel.py    # Core computation engine (v2)
├── requirements.txt
└── README.md
```

---

## Running locally

```bash
git clone https://github.com/karmakar-rahul/christoffel-symbol-calculator.git
cd christoffel-symbol-calculator

# (optional) virtual environment
python -m venv venv
source venv/bin/activate        # Windows: venv\Scripts\activate

pip install -r requirements.txt
streamlit run app.py
```

---

## Features

- Christoffel symbols Γ^λ_{μν} for any metric tensor
- Riemann curvature tensor R^ρ_{σμν}
- Ricci tensor R_{μν} and Ricci scalar R
- Einstein tensor G_{μν}
- Predefined GR metrics: Schwarzschild, Kerr, Reissner-Nordström, Kerr-Newman, FLRW, Minkowski, de Sitter
- Custom metric input (Python/SymPy syntax)
- Symbolic simplification with selectable level (None / Basic / Full)
- LaTeX rendering and `.tex` export
- Session-level caching for instant re-runs

---

## Technologies

- Python 3.x
- [Streamlit](https://streamlit.io)
- [SymPy](https://sympy.org)

---

## Author

**Rahul Karmakar**, MSc Physics (Astrophysics Specialization)
