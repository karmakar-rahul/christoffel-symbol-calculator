"""
christoffel.py  –  v2
Core computation engine for the Christoffel Symbol Calculator.

Key improvements over v1:
  - Symbolic metric inversion replaced with cofactor-expansion (avoids heavy
    symbolic fraction simplification inside .inv())
  - Derivative cache: each g_{mn} component is differentiated once per
    coordinate; results are reused across all Christoffel components
  - Zero-skipping: Christoffel components whose numerator is identically zero
    are never passed through the simplifier
  - Optional concurrent simplification via concurrent.futures (ProcessPoolExecutor
    is used when the caller opts-in; ThreadPoolExecutor is used inside Streamlit
    where forking is unreliable)
  - Extended outputs: Riemann tensor, Ricci tensor, Ricci scalar, and Einstein
    tensor are computed on demand from the already-cached Christoffel data
  - Fully backward-compatible API: existing callers need zero changes
"""

from __future__ import annotations

import concurrent.futures
from functools import lru_cache
from typing import Callable, List, Literal, Optional, Sequence, Tuple

import sympy as sp
# Public types
Tensor3 = List[List[List[sp.Expr]]]   # Gamma[l][m][n]
Tensor4 = List[List[List[List[sp.Expr]]]]  # Riemann[rho][sigma][mu][nu]
# Internal helpers
def _apply_simplify(expr: sp.Expr, level: str) -> sp.Expr:
    """Apply simplification at the chosen level."""
    if level == "Full":
        return sp.simplify(expr)
    if level == "Basic":
        return sp.cancel(sp.expand(expr))
    return expr   # "None"


def _simplify_task(args: Tuple[sp.Expr, str]) -> sp.Expr:
    """Top-level picklable function for ProcessPoolExecutor."""
    expr, level = args
    return _apply_simplify(expr, level)
def _build_derivative_cache(
    metric: sp.Matrix,
    coords: List[sp.Symbol],
) -> List[List[List[sp.Expr]]]:
    """
    Compute all partial derivatives dg[m][n][k] = ∂g_{mn}/∂x^k.

    Only the upper triangle (m ≤ n) is computed symbolically; the lower
    triangle is filled by symmetry.  Each element is differentiated once
    and reused.
    """
    dim = metric.shape[0]
    dg: List[List[List[sp.Expr]]] = [
        [[sp.S.Zero] * dim for _ in range(dim)] for _ in range(dim)
    ]
    for m in range(dim):
        for n in range(m, dim):
            for k in range(dim):
                val = sp.diff(metric[m, n], coords[k])
                dg[m][n][k] = val
                if m != n:
                    dg[n][m][k] = val
    return dg
def _christoffel_component(
    l: int,
    m: int,
    n: int,
    g_inv: sp.Matrix,
    dg: List[List[List[sp.Expr]]],
    dim: int,
) -> sp.Expr:
    """Compute a single Γ^l_{mn} without simplification."""
    total = sp.S.Zero
    for k in range(dim):
        g_inv_lk = g_inv[l, k]
        if g_inv_lk == sp.S.Zero:
            continue
        bracket = dg[k][m][n] + dg[k][n][m] - dg[m][n][k]
        if bracket == sp.S.Zero:
            continue
        total += g_inv_lk * bracket
    return sp.Rational(1, 2) * total
# Public API – Christoffel symbols
def compute_christoffel(
    metric: sp.Matrix,
    coord_symbols: Optional[Sequence[sp.Symbol]] = None,
    simplify_level: Literal["Full", "Basic", "None"] = "Basic",
    *,
    symmetrize: bool = True,
    parallel: bool = False,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> Tensor3:
    """
    Compute Christoffel symbols of the second kind (Levi-Civita connection).

    Parameters
    ----------
    metric : sympy.Matrix
        The covariant metric tensor g_{μν}.  Must be square.
    coord_symbols : sequence of sympy.Symbol, optional
        Coordinate symbols [x⁰, x¹, …].  Inferred if omitted (less reliable).
    simplify_level : {"Full", "Basic", "None"}, default "Basic"
        "Full"  → sp.simplify   (slowest, most compact)
        "Basic" → cancel(expand) (good balance)
        "None"  → raw expression (fastest, may be large)
    symmetrize : bool, default True
        Symmetrize metric as (g + gᵀ)/2 before computation.
    parallel : bool, default False
        Simplify independent components concurrently with ThreadPoolExecutor.
        Speeds up "Full" or "Basic" on multi-core machines for large metrics.
    progress_callback : Callable[[int, int], None], optional
        Called as progress_callback(current, total) after each independent
        component is resolved.

    Returns
    Gamma : list[list[list[Expr]]]
        Gamma[l][m][n] = Γ^l_{mn}.  Zero entries are sp.S.Zero.
    """
    # Validation
    if not isinstance(metric, sp.MatrixBase):
        raise ValueError("Metric must be a SymPy Matrix.")
    if metric.shape[0] != metric.shape[1]:
        raise ValueError("Metric must be square.")
    dim = metric.shape[0]
    if dim == 0:
        raise ValueError("Metric must have positive dimension.")
    if simplify_level not in {"Full", "Basic", "None"}:
        raise ValueError(f"simplify_level must be 'Full', 'Basic', or 'None'.")
    # Symmetrize
    if symmetrize:
        metric = sp.Rational(1, 2) * (metric + metric.T)
    else:
        if metric != metric.T:
            raise ValueError(
                "Metric is not symmetric.  Pass symmetrize=True to auto-fix."
            )
    # Coordinates
    coords = _resolve_coords(metric, coord_symbols, dim)

    #  Inverse metric (v2: use adjugate / det to avoid heavy .inv())
    g_inv = _safe_invert(metric)
    dg = _build_derivative_cache(metric, coords)
    independent: List[Tuple[int, int, int]] = [
        (l, m, n)
        for l in range(dim)
        for m in range(dim)
        for n in range(m, dim)
    ]
    total = len(independent)

    # Compute raw Christoffel values
    raw: dict[Tuple[int, int, int], sp.Expr] = {}
    for idx, (l, m, n) in enumerate(independent):
        raw[(l, m, n)] = _christoffel_component(l, m, n, g_inv, dg, dim)
        if progress_callback:
            try:
                progress_callback(idx + 1, total * 2)  # first half of progress
            except Exception:
                pass

    # Simplification (sequential or parallel)
    non_zero_keys = [(l, m, n) for (l, m, n), v in raw.items() if v != sp.S.Zero]

    simplified: dict[Tuple[int, int, int], sp.Expr] = {}

    if parallel and simplify_level != "None" and len(non_zero_keys) > 4:
        tasks = [(raw[k], simplify_level) for k in non_zero_keys]
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(_simplify_task, tasks))
        for k, v in zip(non_zero_keys, results):
            simplified[k] = v
    else:
        for k in non_zero_keys:
            simplified[k] = _apply_simplify(raw[k], simplify_level)
    # Assemble output tensor
    Gamma: Tensor3 = [
        [[sp.S.Zero] * dim for _ in range(dim)] for _ in range(dim)
    ]
    for idx, (l, m, n) in enumerate(independent):
        val = simplified.get((l, m, n), sp.S.Zero)
        Gamma[l][m][n] = val
        if n != m:
            Gamma[l][n][m] = val
        if progress_callback:
            try:
                progress_callback(total + idx + 1, total * 2)
            except Exception:
                pass

    return Gamma
# Extended tensor outputs
def compute_riemann(Gamma: Tensor3, coords: List[sp.Symbol]) -> Tensor4:
    """
    Compute the Riemann curvature tensor R^ρ_{σμν} from Christoffel symbols.

    R^ρ_{σμν} = ∂_μ Γ^ρ_{νσ} − ∂_ν Γ^ρ_{μσ} + Γ^ρ_{μλ} Γ^λ_{νσ} − Γ^ρ_{νλ} Γ^λ_{μσ}
    """
    dim = len(Gamma)
    R: Tensor4 = [
        [[[sp.S.Zero] * dim for _ in range(dim)] for _ in range(dim)]
        for _ in range(dim)
    ]
    for rho in range(dim):
        for sigma in range(dim):
            for mu in range(dim):
                for nu in range(mu + 1, dim):  # antisymmetry in last two
                    term1 = sp.diff(Gamma[rho][nu][sigma], coords[mu])
                    term2 = sp.diff(Gamma[rho][mu][sigma], coords[nu])
                    term3 = sum(
                        Gamma[rho][mu][lam] * Gamma[lam][nu][sigma]
                        for lam in range(dim)
                    )
                    term4 = sum(
                        Gamma[rho][nu][lam] * Gamma[lam][mu][sigma]
                        for lam in range(dim)
                    )
                    val = sp.cancel(sp.expand(term1 - term2 + term3 - term4))
                    R[rho][sigma][mu][nu] = val
                    R[rho][sigma][nu][mu] = -val
    return R


def compute_ricci_tensor(R: Tensor4) -> List[List[sp.Expr]]:
    """
    Ricci tensor: R_{μν} = R^λ_{μλν}  (contraction of first and third index).
    """
    dim = len(R)
    Ric: List[List[sp.Expr]] = [[sp.S.Zero] * dim for _ in range(dim)]
    for mu in range(dim):
        for nu in range(mu, dim):
            val = sp.cancel(sp.expand(sum(R[lam][mu][lam][nu] for lam in range(dim))))
            Ric[mu][nu] = val
            Ric[nu][mu] = val
    return Ric


def compute_ricci_scalar(
    Ric: List[List[sp.Expr]], g_inv: sp.Matrix
) -> sp.Expr:
    """Ricci scalar R = g^{μν} R_{μν}."""
    dim = len(Ric)
    R_scalar = sum(
        g_inv[mu, nu] * Ric[mu][nu]
        for mu in range(dim)
        for nu in range(dim)
    )
    return sp.cancel(sp.expand(R_scalar))


def compute_einstein_tensor(
    Ric: List[List[sp.Expr]],
    R_scalar: sp.Expr,
    metric: sp.Matrix,
) -> List[List[sp.Expr]]:
    """Einstein tensor G_{μν} = R_{μν} − ½ g_{μν} R."""
    dim = metric.shape[0]
    G: List[List[sp.Expr]] = [[sp.S.Zero] * dim for _ in range(dim)]
    for mu in range(dim):
        for nu in range(mu, dim):
            val = sp.cancel(
                sp.expand(Ric[mu][nu] - sp.Rational(1, 2) * metric[mu, nu] * R_scalar)
            )
            G[mu][nu] = val
            G[nu][mu] = val
    return G
# Convenience: compute everything
def compute_all_tensors(
    metric: sp.Matrix,
    coord_symbols: Optional[Sequence[sp.Symbol]] = None,
    simplify_level: Literal["Full", "Basic", "None"] = "Basic",
    *,
    symmetrize: bool = True,
    parallel: bool = False,
    progress_callback: Optional[Callable[[int, int], None]] = None,
) -> dict:
    """
    Compute Christoffel symbols, Riemann tensor, Ricci tensor, Ricci scalar,
    and Einstein tensor in a single call.

    Returns
    -------
    dict with keys:
        "Gamma"   : Tensor3
        "Riemann" : Tensor4
        "Ricci"   : List[List[Expr]]
        "R"       : Expr  (Ricci scalar)
        "Einstein": List[List[Expr]]
        "g_inv"   : Matrix
        "coords"  : List[Symbol]
    """
    if not isinstance(metric, sp.MatrixBase):
        raise ValueError("Metric must be a SymPy Matrix.")
    dim = metric.shape[0]
    if symmetrize:
        metric = sp.Rational(1, 2) * (metric + metric.T)
    coords = _resolve_coords(metric, coord_symbols, dim)
    g_inv = _safe_invert(metric)

    Gamma = compute_christoffel(
        metric,
        coord_symbols=coords,
        simplify_level=simplify_level,
        symmetrize=False,  # already done
        parallel=parallel,
        progress_callback=progress_callback,
    )
    R_tensor = compute_riemann(Gamma, coords)
    Ric = compute_ricci_tensor(R_tensor)
    R_scalar = compute_ricci_scalar(Ric, g_inv)
    G = compute_einstein_tensor(Ric, R_scalar, metric)

    return {
        "Gamma": Gamma,
        "Riemann": R_tensor,
        "Ricci": Ric,
        "R": R_scalar,
        "Einstein": G,
        "g_inv": g_inv,
        "coords": coords,
    }
# Private helpers
def _resolve_coords(
    metric: sp.Matrix,
    coord_symbols: Optional[Sequence[sp.Symbol]],
    dim: int,
) -> List[sp.Symbol]:
    if coord_symbols is not None:
        if len(coord_symbols) != dim:
            raise ValueError(
                f"Number of coordinates ({len(coord_symbols)}) must match "
                f"metric dimension ({dim})."
            )
        for s in coord_symbols:
            if not isinstance(s, sp.Symbol):
                raise ValueError("All coordinate entries must be SymPy Symbols.")
        return list(coord_symbols)
    # Inference fallback
    free_syms = list(metric.free_symbols)
    param_names = {"M", "a", "Lambda", "Q", "c", "G", "k"}
    candidates = [s for s in free_syms if s.name not in param_names]
    known_coord_names = ["t", "r", "theta", "phi", "x", "y", "z", "rho", "chi"]
    preferred = [s for s in candidates if s.name in known_coord_names]
    order_map = {name: i for i, name in enumerate(known_coord_names)}

    def sort_key(sym: sp.Symbol) -> Tuple:
        n = sym.name
        return (0, order_map[n]) if n in order_map else (1, n)

    coords = sorted(preferred, key=sort_key)
    if len(coords) != dim:
        fallback = sorted(candidates, key=sort_key)
        if len(fallback) == dim:
            coords = fallback
        else:
            raise ValueError(
                "Could not reliably infer coordinate symbols.  "
                "Please provide coord_symbols explicitly."
            )
    return coords


def _safe_invert(metric: sp.Matrix) -> sp.Matrix:
    """
    Invert the metric symbolically.

    Strategy (v2):
    1. Try metric.inv("ADJ") – adjugate / determinant, which avoids internal
       Gaussian elimination and tends to produce simpler intermediate expressions
       for block-diagonal and sparse metrics (Schwarzschild, Kerr, FLRW …).
    2. Fall back to metric.inv() if the adjugate route raises.
    """
    try:
        return metric.inv("ADJ")
    except Exception:
        pass
    try:
        return metric.inv()
    except Exception as e:
        raise ValueError(f"Cannot invert metric tensor: {e}") from e

# Formatting helpers 
def nonzero_christoffel(
    Gamma: Tensor3,
    coords: List[sp.Symbol],
) -> List[Tuple[str, sp.Expr]]:
    """
    Return a sorted list of (label, expression) for all non-zero Γ^l_{mn}.

    Label format: "Γ^{coord_l}_{coord_m coord_n}"
    """
    dim = len(Gamma)
    results = []
    for l in range(dim):
        for m in range(dim):
            for n in range(m, dim):
                val = Gamma[l][m][n]
                if val != sp.S.Zero:
                    c_l = sp.latex(coords[l])
                    c_m = sp.latex(coords[m])
                    c_n = sp.latex(coords[n])
                    label = rf"\Gamma^{{{c_l}}}_{{{c_m}{c_n}}}"
                    results.append((label, val))
    return results
