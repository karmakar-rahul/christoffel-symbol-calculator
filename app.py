import streamlit as st
import sympy as sp
from christoffel import compute_christoffel

st.set_page_config(page_title="Christoffel Symbol Calculator", page_icon="🌀")
st.title("Christoffel Symbol Calculator 🌀")
st.write("Calculate Christoffel symbols for any metric tensor in General Relativity")

# Sidebar with examples
with st.sidebar:
    st.header("Metric Library")
    st.markdown("""
    This calculator includes several important metrics from General Relativity:
    
    **Black Holes:**
    - Schwarzschild (non-rotating, uncharged)
    - Reissner-Nordström (charged)
    - Kerr (rotating)
    - Kerr-Newman (rotating + charged)
    
    **Cosmology:**
    - FLRW (expanding universe)
    - Minkowski (flat spacetime)
    
    **Parameters:**
    - M: mass
    - Q: electric charge
    - a: angular momentum per unit mass
    - k: spatial curvature (-1, 0, +1)
    
    **Functions:**
    Available: sin, cos, tan, exp, log, sqrt
    """)

# Main input area
st.subheader("Input Metric Tensor")

# Predefined metrics library
METRIC_LIBRARY = {
    "Custom": {
        "metric": "[[1,0,0],[0,r**2,0],[0,0,r**2*sin(theta)**2]]",
        "coords": "r,theta,phi",
        "description": "Enter your own metric"
    },
    "Schwarzschild (3D)": {
        "metric": "[[1/(1-2*M/r),0,0],[0,r**2,0],[0,0,r**2*sin(theta)**2]]",
        "coords": "r,theta,phi",
        "description": "Spherically symmetric vacuum solution (spatial part)"
    },
    "Schwarzschild (4D)": {
        "metric": "[[-(1-2*M/r),0,0,0],[0,1/(1-2*M/r),0,0],[0,0,r**2,0],[0,0,0,r**2*sin(theta)**2]]",
        "coords": "t,r,theta,phi",
        "description": "Spherically symmetric vacuum black hole"
    },
    "Reissner-Nordström": {
        "metric": "[[-(1-2*M/r+Q**2/r**2),0,0,0],[0,1/(1-2*M/r+Q**2/r**2),0,0],[0,0,r**2,0],[0,0,0,r**2*sin(theta)**2]]",
        "coords": "t,r,theta,phi",
        "description": "Charged black hole (M=mass, Q=charge)"
    },
    "Kerr": {
        "metric": "[[-(1-2*M*r/(r**2+a**2*cos(theta)**2)),0,0,-2*M*r*a*sin(theta)**2/(r**2+a**2*cos(theta)**2)],[0,(r**2+a**2*cos(theta)**2)/(r**2-2*M*r+a**2),0,0],[0,0,r**2+a**2*cos(theta)**2,0],[-2*M*r*a*sin(theta)**2/(r**2+a**2*cos(theta)**2),0,0,(r**2+a**2+2*M*r*a**2*sin(theta)**2/(r**2+a**2*cos(theta)**2))*sin(theta)**2]]",
        "coords": "t,r,theta,phi",
        "description": "Rotating black hole (M=mass, a=angular momentum per unit mass)"
    },
    "Kerr-Newman": {
        "metric": "[[-(1-2*M*r/(r**2+a**2*cos(theta)**2)+Q**2/(r**2+a**2*cos(theta)**2)),0,0,-2*M*r*a*sin(theta)**2/(r**2+a**2*cos(theta)**2)],[0,(r**2+a**2*cos(theta)**2)/(r**2-2*M*r+a**2+Q**2),0,0],[0,0,r**2+a**2*cos(theta)**2,0],[-2*M*r*a*sin(theta)**2/(r**2+a**2*cos(theta)**2),0,0,(r**2+a**2+((2*M*r - Q**2)*a**2*sin(theta)**2)/(r**2+a**2*cos(theta)**2))*sin(theta)**2]]",
        "coords": "t,r,theta,phi",
        "description": "Rotating charged black hole (M=mass, a=angular momentum, Q=charge)"
    },
    "FLRW (flat)": {
        "metric": "[[-1,0,0,0],[0,a**2,0,0],[0,0,a**2,0],[0,0,0,a**2]]",
        "coords": "t,r,theta,phi",
        "description": "Flat Friedmann-Lemaître-Robertson-Walker (a(t)=scale factor)"
    },
    "FLRW (spherical)": {
        "metric": "[[-1,0,0,0],[0,a**2/(1-k*r**2),0,0],[0,0,a**2*r**2,0],[0,0,0,a**2*r**2*sin(theta)**2]]",
        "coords": "t,r,theta,phi",
        "description": "FLRW with curvature (a(t)=scale factor, k=spatial curvature)"
    },
    "Minkowski (flat)": {
        "metric": "[[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]",
        "coords": "t,x,y,z",
        "description": "Flat spacetime (special relativity)"
    },
}

# Metric selector
selected_metric = st.selectbox(
    "Choose a metric:",
    options=list(METRIC_LIBRARY.keys()),
    help="Select from predefined metrics or choose Custom to enter your own"
)

st.info(f"ℹ️ {METRIC_LIBRARY[selected_metric]['description']}")

col1, col2 = st.columns([3, 1])

with col1:
    default_metric = METRIC_LIBRARY[selected_metric]["metric"]
    metric_str = st.text_area(
        "Metric tensor:",
        value=default_metric,
        height=100,
        help="Use square brackets for matrices. Diagonal metrics are faster to compute."
    )

with col2:
    st.write("**Coordinates:**")
    default_coords = METRIC_LIBRARY[selected_metric]["coords"]
    coord_input = st.text_input(
        "Order of coordinates",
        value=default_coords,
        help="Comma-separated list matching metric dimension"
    )

# Computation
if st.button("Compute Christoffel Symbols", type="primary"):
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    try:
        # Initialize symbols and functions
        status_text.text("Initializing symbols...")
        progress_bar.progress(5)
        r, theta, phi, t = sp.symbols('r theta phi t', real=True)
        x, y, z = sp.symbols('x y z', real=True)
        M, a_param, Lambda, Q, k = sp.symbols('M a Lambda Q k', real=True, positive=True)
        # Define scale factor a(t) for FLRW metrics; keep Kerr/Kerr-Newman 'a' as a constant parameter
        scale_factor = sp.Function('a')(t)
        
        # Use a(t) only for FLRW metrics; otherwise treat 'a' as constant spin parameter
        if "FLRW" in selected_metric:
            a_for_metric = scale_factor
        else:
            a_for_metric = a_param
        
        # Safe evaluation context with all necessary functions and symbols
        safe_dict = {
            "r": r, "theta": theta, "phi": phi, "t": t,
            "x": x, "y": y, "z": z,
            "M": M, "a": a_for_metric, "Lambda": Lambda, "Q": Q, "k": k,
            "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
            "exp": sp.exp, "log": sp.log, "ln": sp.log,
            "sqrt": sp.sqrt, "pi": sp.pi,
            "abs": sp.Abs, "sign": sp.sign
        }
        
        status_text.text("Parsing metric tensor...")
        progress_bar.progress(10)
        
        # Parse and validate metric with restricted eval environment
        try:
            metric = sp.Matrix(eval(metric_str, {"__builtins__": {}}, safe_dict))
            progress_bar.progress(20)
            status_text.text("Metric parsed successfully...")
        except SyntaxError:
            st.error("❌ Syntax Error: Please check your matrix format. Use [[row1], [row2], ...] format.")
            st.stop()
        except NameError as e:
            st.error(f"❌ Unknown symbol or function: {e}. Use only supported functions and coordinates.")
            st.stop()
        except Exception as e:
            st.error(f"❌ Error parsing metric: {e}")
            st.stop()

        status_text.text("Validating coordinate system...")
        progress_bar.progress(30)
        
        # Parse coordinate input
        coord_names = [c.strip() for c in coord_input.split(',')]
        if len(coord_names) != metric.shape[0]:
            st.error(f"❌ Number of coordinates ({len(coord_names)}) must match metric dimension ({metric.shape[0]})")
            st.stop()
        
        # Map coordinate names to symbols
        coord_map = {'t': t, 'r': r, 'theta': theta, 'phi': phi, 'x': x, 'y': y, 'z': z}
        try:
            coords_in_metric = [coord_map[name] for name in coord_names]
        except KeyError as e:
            st.error(f"❌ Unknown coordinate: {e}. Use: t, r, theta, phi, x, y, z")
            st.stop()
        
        # Validate metric properties
        if metric.shape[0] != metric.shape[1]:
            st.error("❌ Metric must be a square matrix!")
            st.stop()
        elif metric.shape[0] > 4:
            st.warning("⚠️ Large metrics may take considerable time to compute.")
        
        # Check if metric is symmetric
        if metric != metric.T:
            st.warning("⚠️ Warning: Metric should be symmetric. Using symmetrized version.")
            metric = (metric + metric.T) / 2
        
        # Display the metric
        st.success("✓ Metric tensor validation complete!")
        
        # Identify parameters from metric
        status_text.text("Analyzing metric parameters...")
        progress_bar.progress(40)
        free_syms = list(metric.free_symbols)
        params_in_metric = [s for s in free_syms if s not in coords_in_metric]
        
        st.write(f"**Coordinates:** {', '.join(coord_names)}")
        if params_in_metric:
            st.write(f"**Parameters:** {', '.join(str(s) for s in params_in_metric)}")
        
        st.write("#### Metric Tensor:")
        st.latex(r"g_{\mu\nu} = " + sp.latex(metric))
        
        # Compute Christoffel symbols with explicit coordinates
        status_text.text("Computing Christoffel symbols...")
        progress_bar.progress(50)
        
        def progress_callback(percentage, message):
            # Map the compute_christoffel progress (0-100) to our remaining progress (50-90)
            mapped_progress = 50 + int(percentage * 0.4)
            progress_bar.progress(mapped_progress)
            status_text.text(message)
        
        Gamma = compute_christoffel(metric, coord_symbols=coords_in_metric, progress_callback=progress_callback)
        
        status_text.text("Formatting results...")
        progress_bar.progress(90)
        
        # Display results
        st.write("### Non-zero Christoffel Symbols")
        st.write(r"Using the convention: $\Gamma^\lambda_{\mu\nu} = \frac{1}{2}g^{\lambda\sigma}\left(\frac{\partial g_{\sigma\mu}}{\partial x^\nu} + \frac{\partial g_{\sigma\nu}}{\partial x^\mu} - \frac{\partial g_{\mu\nu}}{\partial x^\sigma}\right)$")
        
        # Coordinate labels for LaTeX display
        coord_labels = []
        for name in coord_names:
            if name == 'theta':
                coord_labels.append(r'\theta')
            elif name == 'phi':
                coord_labels.append(r'\phi')
            else:
                coord_labels.append(name)
        
        # Count and display non-zero symbols
        non_zero_count = 0
        results = []
        displayed = set()  # Track which we've already shown (for symmetry)
        
        for l in range(len(Gamma)):
            for m in range(len(Gamma[l])):
                for n in range(len(Gamma[l][m])):
                    val = Gamma[l][m][n]
                    if val != 0:
                        # Use symmetry to avoid duplicates
                        key = (l, min(m, n), max(m, n))
                        if key not in displayed:
                            displayed.add(key)
                            non_zero_count += 1
                            results.append((l, m, n, val))
        
        if non_zero_count == 0:
            st.info("All Christoffel symbols are zero (flat spacetime).")
        else:
            st.info(f"Found {non_zero_count} non-zero independent symbols")
            # Display in a single container for better performance
            for l, m, n, val in results:
                st.latex(f"\\Gamma^{{{coord_labels[l]}}}_{{{coord_labels[m]} {coord_labels[n]}}} = {sp.latex(val)}")
        
        # Complete the progress bar
        status_text.text("Calculation complete!")
        progress_bar.progress(100)
        
        # Option to show all components (including zeros)
        st.write("---")
        if st.checkbox("Show all components (including zeros)"):
            st.write("#### Complete Christoffel Symbol Tensor:")
            for l in range(len(Gamma)):
                with st.expander(f"Γ^{coord_labels[l]}_{{μν}}"):
                    for m in range(len(Gamma[l])):
                        for n in range(len(Gamma[l][m])):
                            val = Gamma[l][m][n]
                            st.latex(f"\\Gamma^{{{coord_labels[l]}}}_{{{coord_labels[m]} {coord_labels[n]}}} = {sp.latex(val)}")

    except Exception as e:
        st.error(f"❌ Error: {e}")
        st.write("Please check your input and try again.")

# Footer
st.markdown("---")
st.caption("Built with Streamlit and SymPy | Christoffel symbols are computed symbolically")