import sympy as sp
from typing import Callable, List, Optional, Sequence, Literal


def compute_christoffel(
    metric: sp.Matrix,
    coord_symbols: Optional[Sequence[sp.Symbol]] = None,
    simplify_level: Literal["Full", "Basic", "None"] = "Basic",
    *,
    symmetrize: bool = True,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> List[List[List[sp.Expr]]]:
    """
    Compute Christoffel symbols of the second kind (Levi-Civita connection)
    for a given metric tensor.

    Parameters
    ----------
    metric : sympy.Matrix
        The metric tensor g_{\mu\nu} (covariant components). Must be square (n x n).
    coord_symbols : Sequence[sympy.Symbol], optional
        Ordered list/tuple of coordinate symbols [x^0, x^1, ..., x^{n-1}]. If None,
        the function will attempt to infer coordinates from the metric's free symbols.
        Explicit coordinates are recommended for reliability.
    simplify_level : {"Full", "Basic", "None"}, default "Basic"
        - "Full": sp.simplify for maximal simplification (slowest)
        - "Basic": sp.cancel(sp.expand(...)) for balanced performance (default)
        - "None": no post-simplification beyond the formula
    symmetrize : bool, keyword-only, default True
        If True, symmetrize the metric as (g + g.T)/2. If False, validate strict
        symmetry and raise on mismatch.
    progress_callback : Callable[[int, int], None], optional
        Callback invoked as progress_callback(current, total) during computation
        of independent components (m <= n) for each upper index l.

    Returns
    -------
    List[List[List[sp.Expr]]]
        A 3D nested list Gamma where Gamma[l][m][n] represents \Gamma^l_{mn}.
        Zero entries are represented with sp.S.Zero for consistent SymPy typing.

    Notes
    -----
    - Complexity is roughly O(n^4) due to nested loops and derivatives/inversion.
    - For Levi-Civita connection, \Gamma^\lambda_{\mu\nu} is symmetric in the
      lower indices (\mu, \nu). This implementation computes only for m <= n and
      mirrors to enforce symmetry.
    - This function is pure (no UI calls) and suitable for use in Streamlit apps.
    """

    # Validate metric type and shape
    if not isinstance(metric, sp.MatrixBase):
        raise ValueError("Metric must be a SymPy Matrix.")

    if metric.shape[0] != metric.shape[1]:
        raise ValueError("Metric must be a square matrix.")

    dim = metric.shape[0]
    if dim == 0:
        raise ValueError("Metric must have positive dimension.")

    # Validate simplification level
    allowed_levels = {"Full", "Basic", "None"}
    if simplify_level not in allowed_levels:
        raise ValueError(
            f"Invalid simplify_level: {simplify_level}. Choose from {sorted(allowed_levels)}"
        )

    # Symmetry handling
    if symmetrize:
        metric = sp.Rational(1, 2) * (metric + metric.T)
    else:
        if metric != metric.T:
            raise ValueError("Metric must be symmetric (set symmetrize=True to auto-symmetrize).")

    # Coordinates: validate or infer
    if coord_symbols is not None:
        if len(coord_symbols) != dim:
            raise ValueError(
                f"Number of coordinates ({len(coord_symbols)}) must match metric dimension ({dim})."
            )
        for s in coord_symbols:
            if not isinstance(s, sp.Symbol):
                raise ValueError("All coordinate entries must be SymPy Symbols.")
        coords = list(coord_symbols)
    else:
        # Attempt inference (prefer explicit for reliability)
        free_syms = list(metric.free_symbols)
        # Known parameter names to exclude from coordinates
        param_names = {"M", "a", "Lambda", "Q", "c", "G", "k"}
        candidates = [s for s in free_syms if s.name not in param_names]

        # Restrict to known coordinate-like names first
        known_coord_names = ["t", "r", "theta", "phi", "x", "y", "z", "rho"]
        preferred = [s for s in candidates if s.name in known_coord_names]

        # Sorting for stable, conventional order
        order_map = {name: i for i, name in enumerate(known_coord_names)}

        def sort_key(sym: sp.Symbol):
            name = sym.name
            return (0, order_map[name]) if name in order_map else (1, name)

        coords = sorted(preferred, key=sort_key)

        if len(coords) != dim:
            # Fall back to using all candidates if it uniquely matches dim
            fallback = sorted(candidates, key=sort_key)
            if len(fallback) == dim:
                coords = fallback
            else:
                raise ValueError(
                    "Could not reliably infer coordinate symbols from metric. "
                    "Please provide coord_symbols explicitly."
                )

    # Attempt metric inversion
    try:
        g_inv = metric.inv()
    except Exception as e:
        raise ValueError(f"Cannot invert metric tensor: {e}")

    # Prepare output structure with SymPy zeros for consistent typing
    Gamma: List[List[List[sp.Expr]]] = [
        [[sp.S.Zero for _ in range(dim)] for _ in range(dim)] for _ in range(dim)
    ]

    # Precompute partial derivatives: dg[m][n][k] = d g_{mn} / d x^k
    # Only compute for m <= n and mirror to dg[n][m][k] due to symmetry of g
    dg: List[List[List[sp.Expr]]] = [
        [[sp.S.Zero for _ in range(dim)] for _ in range(dim)] for _ in range(dim)
    ]

    for m in range(dim):
        for n in range(m, dim):
            for k in range(dim):
                dg[m][n][k] = sp.diff(metric[m, n], coords[k])
                if m != n:
                    dg[n][m][k] = dg[m][n][k]

    # Helper to simplify expressions according to level
    def apply_simplify(expr: sp.Expr) -> sp.Expr:
        if simplify_level == "Full":
            return sp.simplify(expr)
        if simplify_level == "Basic":
            return sp.cancel(sp.expand(expr))
        return expr  # None

    # Total number of independent (l, m, n with m <= n) computations
    total = dim * (dim * (dim + 1) // 2)
    current = 0

    # Compute Christoffel symbols (use symmetry in lower indices)
    for l in range(dim):
        for m in range(dim):
            for n in range(m, dim):
                sum_expr = sp.S.Zero
                for k in range(dim):
                    sum_expr += g_inv[l, k] * (
                        dg[k][m][n] + dg[k][n][m] - dg[m][n][k]
                    )
                result = sp.Rational(1, 2) * sum_expr
                val = apply_simplify(result) if result != 0 else sp.S.Zero
                Gamma[l][m][n] = val
                if n != m:
                    Gamma[l][n][m] = val

                # Progress update for each independent component
                current += 1
                if progress_callback is not None:
                    try:
                        progress_callback(current, total)
                    except Exception:
                        # Ensure computation continues even if callback fails
                        pass

    return Gamma
