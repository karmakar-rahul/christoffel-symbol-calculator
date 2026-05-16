"""
app.py  –  Christoffel Symbol Calculator  v2
Streamlit front-end.

v2 highlights
• st.cache_data on heavy computations – repeated runs for the same metric
  are instant (Streamlit caches across reruns within the same session).
• Tabbed output: Christoffel | Riemann | Ricci | Ricci Scalar | Einstein
• "Compute only Christoffel" vs "Compute all tensors" selector
• Non-zero filter with total count banner
• Download computed LaTeX as a .tex file
• Collapsible "raw SymPy" view alongside LaTeX for debugging
• Parallel simplification toggle (helpful for Full simplification level)
• Improved sidebar with metric preview and coordinate legend
• Dark-friendly CSS overrides for LaTeX panels

Author : Rahul Karmakar  (v1)
Updated: v2 improvements
"""

import streamlit as st
import sympy as sp

from christoffel import (
    compute_christoffel,
    compute_all_tensors,
    nonzero_christoffel,
    _resolve_coords,
    _safe_invert,
    compute_riemann,
    compute_ricci_tensor,
    compute_ricci_scalar,
    compute_einstein_tensor,
)
#  Page config
st.set_page_config(
    page_title="Christoffel Symbol Calculator v2",
    page_icon="🌀",
    layout="wide",
)

#  Custom CSS
st.markdown(
    """
    <style>
    /* ── General layout ── */
    .block-container { padding-top: 1.5rem; padding-bottom: 2rem; }

    /* ── LaTeX result card ── */
    .gamma-card {
        background: var(--secondary-background-color, #1e1e2e);
        border: 1px solid var(--primary-color, #5c5c8a);
        border-left: 4px solid var(--primary-color, #7c7cf0);
        border-radius: 8px;
        padding: 0.6rem 1rem;
        margin-bottom: 0.5rem;
        font-size: 1.05rem;
    }
    .gamma-label {
        font-family: monospace;
        font-size: 0.78rem;
        color: #aaaacc;
        margin-bottom: 0.2rem;
    }
    .zero-banner {
        background: #0e2a1f;
        border-left: 4px solid #2ecc71;
        border-radius: 6px;
        padding: 0.5rem 1rem;
        color: #2ecc71;
        font-weight: 600;
    }
    .count-banner {
        background: #1a1a3a;
        border-left: 4px solid #7c7cf0;
        border-radius: 6px;
        padding: 0.5rem 1rem;
        margin-bottom: 0.8rem;
        color: #c9c9ff;
    }
    /* ── Sidebar metric info ── */
    .metric-info {
        background: #12122a;
        border-radius: 8px;
        padding: 0.6rem 0.8rem;
        font-size: 0.82rem;
        color: #aaaacf;
        margin-bottom: 0.6rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)
#  Metric library
METRIC_LIBRARY: dict = {
    "Custom": {
        "metric": "[[1,0,0],[0,r**2,0],[0,0,r**2*sin(theta)**2]]",
        "coords": "r,theta,phi",
        "description": "Enter your own metric tensor",
        "dim": 3,
        "params": [],
    },
    "Schwarzschild (4D)": {
        "metric": (
            "[[-(1-2*M/r),0,0,0],"
            "[0,1/(1-2*M/r),0,0],"
            "[0,0,r**2,0],"
            "[0,0,0,r**2*sin(theta)**2]]"
        ),
        "coords": "t,r,theta,phi",
        "description": "Spherically symmetric vacuum black hole (M = mass)",
        "dim": 4,
        "params": ["M"],
    },
    "Schwarzschild (3D spatial)": {
        "metric": (
            "[[1/(1-2*M/r),0,0],"
            "[0,r**2,0],"
            "[0,0,r**2*sin(theta)**2]]"
        ),
        "coords": "r,theta,phi",
        "description": "Spatial hypersurface of Schwarzschild",
        "dim": 3,
        "params": ["M"],
    },
    "Reissner-Nordström": {
        "metric": (
            "[[-(1-2*M/r+Q**2/r**2),0,0,0],"
            "[0,1/(1-2*M/r+Q**2/r**2),0,0],"
            "[0,0,r**2,0],"
            "[0,0,0,r**2*sin(theta)**2]]"
        ),
        "coords": "t,r,theta,phi",
        "description": "Charged (non-rotating) black hole  (M = mass, Q = charge)",
        "dim": 4,
        "params": ["M", "Q"],
    },
    "Kerr": {
        "metric": (
            "[[-(1-2*M*r/(r**2+a**2*cos(theta)**2)),0,0,"
            "-2*M*r*a*sin(theta)**2/(r**2+a**2*cos(theta)**2)],"
            "[0,(r**2+a**2*cos(theta)**2)/(r**2-2*M*r+a**2),0,0],"
            "[0,0,r**2+a**2*cos(theta)**2,0],"
            "[-2*M*r*a*sin(theta)**2/(r**2+a**2*cos(theta)**2),0,0,"
            "(r**2+a**2+2*M*r*a**2*sin(theta)**2/(r**2+a**2*cos(theta)**2))*sin(theta)**2]]"
        ),
        "coords": "t,r,theta,phi",
        "description": "Rotating black hole  (M = mass, a = spin parameter)",
        "dim": 4,
        "params": ["M", "a"],
    },
    "Kerr-Newman": {
        "metric": (
            "[[-(1-2*M*r/(r**2+a**2*cos(theta)**2)+Q**2/(r**2+a**2*cos(theta)**2)),0,0,"
            "-2*M*r*a*sin(theta)**2/(r**2+a**2*cos(theta)**2)],"
            "[0,(r**2+a**2*cos(theta)**2)/(r**2-2*M*r+a**2+Q**2),0,0],"
            "[0,0,r**2+a**2*cos(theta)**2,0],"
            "[-2*M*r*a*sin(theta)**2/(r**2+a**2*cos(theta)**2),0,0,"
            "(r**2+a**2+((2*M*r-Q**2)*a**2*sin(theta)**2)/(r**2+a**2*cos(theta)**2))*sin(theta)**2]]"
        ),
        "coords": "t,r,theta,phi",
        "description": "Rotating charged black hole  (M, a, Q)",
        "dim": 4,
        "params": ["M", "a", "Q"],
    },
    "FLRW (flat)": {
        "metric": (
            "[[-1,0,0,0],[0,a**2,0,0],[0,0,a**2,0],[0,0,0,a**2]]"
        ),
        "coords": "t,r,theta,phi",
        "description": "Flat expanding universe  (a = a(t) scale factor)",
        "dim": 4,
        "params": ["a"],
    },
    "FLRW (spherical)": {
        "metric": (
            "[[-1,0,0,0],"
            "[0,a**2/(1-k*r**2),0,0],"
            "[0,0,a**2*r**2,0],"
            "[0,0,0,a**2*r**2*sin(theta)**2]]"
        ),
        "coords": "t,r,theta,phi",
        "description": "FLRW with curvature  (a = scale factor, k = curvature)",
        "dim": 4,
        "params": ["a", "k"],
    },
    "Minkowski (Cartesian)": {
        "metric": "[[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]",
        "coords": "t,x,y,z",
        "description": "Flat spacetime – all Christoffel symbols vanish",
        "dim": 4,
        "params": [],
    },
    "Minkowski (spherical)": {
        "metric": "[[−1,0,0,0],[0,1,0,0],[0,0,r**2,0],[0,0,0,r**2*sin(theta)**2]]",
        "coords": "t,r,theta,phi",
        "description": "Flat spacetime in spherical polar coordinates",
        "dim": 4,
        "params": [],
    },
    "de Sitter": {
        "metric": (
            "[[-(1-Lambda*r**2/3),0,0,0],"
            "[0,1/(1-Lambda*r**2/3),0,0],"
            "[0,0,r**2,0],"
            "[0,0,0,r**2*sin(theta)**2]]"
        ),
        "coords": "t,r,theta,phi",
        "description": "Cosmological constant spacetime  (Λ = Lambda)",
        "dim": 4,
        "params": ["Lambda"],
    },
}
#Sidebar
with st.sidebar:
    st.header("Settings")

    selected_metric = st.selectbox(
        "Preset metric:",
        list(METRIC_LIBRARY.keys()),
        help="Choose a predefined GR metric, or 'Custom' to enter your own.",
    )

    meta = METRIC_LIBRARY[selected_metric]
    st.markdown(
        f"<div class='metric-info'>"
        f"<b>Description:</b> {meta['description']}<br>"
        f"<b>Dimension:</b> {meta['dim']}D &nbsp;|&nbsp; "
        f"<b>Parameters:</b> {', '.join(meta['params']) if meta['params'] else 'none'}"
        f"</div>",
        unsafe_allow_html=True,
    )

    simplify_level = st.select_slider(
        "Simplification level:",
        options=["None", "Basic", "Full"],
        value="Basic",
        help=(
            "None → fastest (raw expressions)  |  "
            "Basic → cancel/expand balance  |  "
            "Full → sp.simplify (slow for Kerr/RN)"
        ),
    )

    parallel = st.toggle(
        "Parallel simplification",
        value=True,
        help="Use threads to simplify multiple components concurrently.",
    )

    compute_mode = st.radio(
        "Compute:",
        ["Christoffel symbols only", "All tensors (Riemann, Ricci, Einstein)"],
        index=0,
    )

    show_zero = st.toggle("Show zero components", value=False)
    show_raw = st.toggle("Show raw SymPy (debug)", value=False)

    st.divider()
    st.markdown(
        """
        **Symbol legend**
        - M — mass
        - a — spin per unit mass
        - Q — electric charge
        - k — spatial curvature (−1, 0, +1)
        - Λ (Lambda) — cosmological constant
        - a(t) — FLRW scale factor

        **Available functions:** `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`
        """
    )
    st.divider()
    st.caption("v2 · built with SymPy & Streamlit")

st.title("Christoffel Symbol Calculator")
st.write(
    "Compute Christoffel symbols and curvature tensors for any spacetime metric in General Relativity."
)

col_metric, col_coords = st.columns([3, 1])

with col_metric:
    metric_str = st.text_area(
        "Metric tensor  g_{μν}  (Python/SymPy list-of-lists):",
        value=meta["metric"],
        height=120,
        help="Use nested square brackets.  Example diagonal entry: -(1-2*M/r)",
    )

with col_coords:
    coord_input = st.text_input(
        "Coordinates (comma-separated):",
        value=meta["coords"],
        help="Must match the metric dimension.",
    )


#cached 
@st.cache_data(show_spinner=False)
def _cached_christoffel(metric_str: str, coords_str: str, level: str, par: bool):
    """Parse + compute Christoffel, return Gamma and coord list."""
    coords_raw = [c.strip() for c in coords_str.split(",")]
    coord_syms = [sp.Symbol(c) for c in coords_raw]
    allowed = {
        "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
        "exp": sp.exp, "log": sp.log, "sqrt": sp.sqrt,
        **{str(s): s for s in coord_syms},
        "M": sp.Symbol("M"), "a": sp.Symbol("a"),
        "Q": sp.Symbol("Q"), "k": sp.Symbol("k"),
        "Lambda": sp.Symbol("Lambda"), "pi": sp.pi,
    }
    g_list = eval(metric_str, {"__builtins__": {}}, allowed)  # noqa: S307
    g = sp.Matrix(g_list)
    Gamma = compute_christoffel(
        g, coord_syms, level, symmetrize=True, parallel=par
    )
    return Gamma, coord_syms, g


@st.cache_data(show_spinner=False)
def _cached_all_tensors(metric_str: str, coords_str: str, level: str, par: bool):
    """Parse + compute all tensors."""
    coords_raw = [c.strip() for c in coords_str.split(",")]
    coord_syms = [sp.Symbol(c) for c in coords_raw]
    allowed = {
        "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
        "exp": sp.exp, "log": sp.log, "sqrt": sp.sqrt,
        **{str(s): s for s in coord_syms},
        "M": sp.Symbol("M"), "a": sp.Symbol("a"),
        "Q": sp.Symbol("Q"), "k": sp.Symbol("k"),
        "Lambda": sp.Symbol("Lambda"), "pi": sp.pi,
    }
    g_list = eval(metric_str, {"__builtins__": {}}, allowed)  # noqa: S307
    g = sp.Matrix(g_list)
    result = compute_all_tensors(g, coord_syms, level, symmetrize=True, parallel=par)
    return result, coord_syms, g

if st.button("Compute", type="primary", use_container_width=False):
    with st.spinner("Computing … (large metrics like Kerr may take 10–30 s)"):
        try:
            if compute_mode.startswith("Christoffel"):
                Gamma, coord_syms, g = _cached_christoffel(
                    metric_str, coord_input, simplify_level, parallel
                )
                result = None
            else:
                result, coord_syms, g = _cached_all_tensors(
                    metric_str, coord_input, simplify_level, parallel
                )
                Gamma = result["Gamma"]
        except Exception as err:
            st.error(f"Error during computation:\n\n`{err}`")
            st.stop()

    dim = len(Gamma)
    with st.expander("Parsed metric tensor", expanded=False):
        st.latex(r"g_{\mu\nu} = " + sp.latex(g))

    # Tab layout 
    tabs_labels = ["Γ  Christoffel"]
    if result is not None:
        tabs_labels += ["ℝ  Riemann", "Ric  Ricci tensor", "R  Ricci scalar", "G  Einstein"]

    tabs = st.tabs(tabs_labels)

    #  Tab 0 : Christoffel 
    with tabs[0]:
        st.subheader("Christoffel Symbols  Γ^λ_{μν}")

        nz_list = nonzero_christoffel(Gamma, coord_syms)
        nz_count = len(nz_list)
        total_independent = dim * dim * (dim + 1) // 2
        zero_count = total_independent - nz_count

        st.markdown(
            f"<div class='count-banner'>"
            f"Non-zero independent components: <b>{nz_count}</b> / {total_independent}"
            f"</div>",
            unsafe_allow_html=True,
        )

        if nz_count == 0:
            st.markdown(
                "<div class='zero-banner'> All Christoffel symbols vanish — flat metric!</div>",
                unsafe_allow_html=True,
            )
        else:
            # Build LaTeX export string
            latex_lines = [
                r"\documentclass{article}",
                r"\usepackage{amsmath}",
                r"\begin{document}",
                r"\section*{Christoffel Symbols}",
            ]
            for label, val in nz_list:
                st.markdown(
                    f"<div class='gamma-card'><div class='gamma-label'>{label}</div>",
                    unsafe_allow_html=True,
                )
                st.latex(label + r" = " + sp.latex(val))
                if show_raw:
                    st.code(str(val), language="python")
                st.markdown("</div>", unsafe_allow_html=True)
                latex_lines.append(r"$$" + label + r" = " + sp.latex(val) + r"$$")

            latex_lines.append(r"\end{document}")
            st.download_button(
                "Download as LaTeX (.tex)",
                "\n".join(latex_lines),
                file_name="christoffel_symbols.tex",
                mime="text/plain",
            )

        if show_zero and zero_count > 0:
            st.markdown(f"**{zero_count} independent components are zero** (not shown)")

    # Riemann 
    if result is not None:
        with tabs[1]:
            st.subheader("Riemann Curvature Tensor  R^ρ_{σμν}")
            R_tensor = result["Riemann"]
            nz_riemann = []
            for rho in range(dim):
                for sigma in range(dim):
                    for mu in range(dim):
                        for nu in range(mu + 1, dim):
                            val = R_tensor[rho][sigma][mu][nu]
                            if val != sp.S.Zero:
                                c = [sp.latex(coord_syms[i]) for i in [rho, sigma, mu, nu]]
                                lbl = rf"R^{{{c[0]}}}{{{c[1]}}}{{{c[2]}}}{{{c[3]}}}"
                                nz_riemann.append((lbl, val))

            st.markdown(
                f"<div class='count-banner'>Non-zero independent components: <b>{len(nz_riemann)}</b></div>",
                unsafe_allow_html=True,
            )
            if not nz_riemann:
                st.markdown(
                    "<div class='zero-banner'> Riemann tensor vanishes — flat spacetime!</div>",
                    unsafe_allow_html=True,
                )
            for lbl, val in nz_riemann:
                st.markdown(f"<div class='gamma-card'><div class='gamma-label'>{lbl}</div>", unsafe_allow_html=True)
                st.latex(lbl + " = " + sp.latex(val))
                if show_raw:
                    st.code(str(val), language="python")
                st.markdown("</div>", unsafe_allow_html=True)

        # Ricci tensor
        with tabs[2]:
            st.subheader("Ricci Tensor  R_{μν}")
            Ric = result["Ricci"]
            nz_ricci = []
            for mu in range(dim):
                for nu in range(mu, dim):
                    val = Ric[mu][nu]
                    if val != sp.S.Zero:
                        cm, cn = sp.latex(coord_syms[mu]), sp.latex(coord_syms[nu])
                        lbl = rf"R_{{{cm}{cn}}}"
                        nz_ricci.append((lbl, val))

            st.markdown(
                f"<div class='count-banner'>Non-zero independent components: <b>{len(nz_ricci)}</b></div>",
                unsafe_allow_html=True,
            )
            if not nz_ricci:
                st.markdown(
                    "<div class='zero-banner'> Ricci tensor vanishes — Ricci flat metric!</div>",
                    unsafe_allow_html=True,
                )
            for lbl, val in nz_ricci:
                st.markdown(f"<div class='gamma-card'><div class='gamma-label'>{lbl}</div>", unsafe_allow_html=True)
                st.latex(lbl + " = " + sp.latex(val))
                if show_raw:
                    st.code(str(val), language="python")
                st.markdown("</div>", unsafe_allow_html=True)

        # Ricci scalar 
        with tabs[3]:
            st.subheader("Ricci Scalar  R")
            R_sc = result["R"]
            if R_sc == sp.S.Zero:
                st.markdown(
                    "<div class='zero-banner'> Ricci scalar R = 0</div>",
                    unsafe_allow_html=True,
                )
            else:
                st.latex(r"R = " + sp.latex(R_sc))
                if show_raw:
                    st.code(str(R_sc), language="python")

        # Einstein Tensor
        with tabs[4]:
            st.subheader("Einstein Tensor  G_{μν} = R_{μν} − ½ g_{μν} R")
            G = result["Einstein"]
            nz_G = []
            for mu in range(dim):
                for nu in range(mu, dim):
                    val = G[mu][nu]
                    if val != sp.S.Zero:
                        cm, cn = sp.latex(coord_syms[mu]), sp.latex(coord_syms[nu])
                        lbl = rf"G_{{{cm}{cn}}}"
                        nz_G.append((lbl, val))

            st.markdown(
                f"<div class='count-banner'>Non-zero independent components: <b>{len(nz_G)}</b></div>",
                unsafe_allow_html=True,
            )
            if not nz_G:
                st.markdown(
                    "<div class='zero-banner'>Einstein tensor vanishes.</div>",
                    unsafe_allow_html=True,
                )
            for lbl, val in nz_G:
                st.markdown(f"<div class='gamma-card'><div class='gamma-label'>{lbl}</div>", unsafe_allow_html=True)
                st.latex(lbl + " = " + sp.latex(val))
                if show_raw:
                    st.code(str(val), language="python")
                st.markdown("</div>", unsafe_allow_html=True)

else:
    st.info(
        "Select a metric (or enter a custom one), configure settings in the sidebar, then press **🚀 Compute**."
    )
