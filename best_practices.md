# 📘 Project Best Practices

## 1. Project Purpose
A Streamlit application for symbolic computation of Christoffel symbols (Γ^λ_{μν}) from a user-specified metric tensor in General Relativity. The app parses a metric matrix, validates coordinates, and uses SymPy to compute and render non-zero Christoffel symbols in LaTeX. Domain: differential geometry / GR tooling.

## 2. Project Structure
- app.py
  - Streamlit UI, metric presets (Metric Library), input parsing/validation, progress reporting, and LaTeX rendering.
  - Orchestrates computation by calling compute_christoffel from christoffel.py.
- christoffel.py
  - Core symbolic algorithm to compute Γ^λ_{μν} from a metric matrix.
  - Supports optional progress callbacks and adjustable simplification levels (Full/Basic/None).
- requirements.txt
  - Python dependencies (currently empty; see section 7 for recommendations).
- venv/
  - Local virtual environment (not committed typically).
- .qodo/
  - Workspace metadata (not part of runtime).

Conventions and roles:
- UI vs Core Separation: Keep Streamlit UI logic in app.py and symbolic math/computation in christoffel.py.
- Entry point: streamlit run app.py.
- Configuration: Currently none. Metric presets are embedded in app.py as a dictionary constant (METRIC_LIBRARY).

## 3. Test Strategy
Current state: No tests present. Recommended approach:
- Framework: pytest.
- Structure: tests/ directory mirroring module names (e.g., tests/test_christoffel.py, tests/test_app_parsing.py).
- Unit tests (priority):
  - Flat spacetime (Minkowski): All Γ^λ_{μν} = 0.
  - 2D polar coordinates metric: Validate known non-zero components.
  - Symmetry: Γ^λ_{μν} == Γ^λ_{νμ}.
  - Non-invertible metric: Raises a helpful error.
  - Coordinate handling: Explicit coord_symbols vs inferred (where applicable in christoffel.py).
  - Simplification modes: Results are mathematically equivalent across modes (where feasible).
- Integration tests:
  - End-to-end app input → output sanity using Streamlit scripting hooks or by factoring parsing into pure functions.
- Mocking guidelines:
  - For app logic, isolate computation and parsing into testable pure functions; avoid UI-level mocks when possible.
  - If needed, mock Streamlit components using lightweight stubs; avoid asserting presentation details, focus on computed values.
- Coverage expectations:
  - >90% for christoffel.py.
  - Key parsing/validation paths in app.py covered.

## 4. Code Style
- Language and typing:
  - Python 3.10+ recommended.
  - Add type hints for public functions in christoffel.py and any future utility modules.
- Formatting and linting:
  - Use Black (format), isort (imports), and Ruff/Flake8 (lint). Enforce pre-commit hooks.
- Naming conventions:
  - Modules: snake_case (app.py, christoffel.py).
  - Functions: snake_case (compute_christoffel).
  - Variables: snake_case; mathematical symbols may mirror domain notation (e.g., r, theta, phi) but avoid shadowing.
- Documentation:
  - Module and function docstrings with argument/return descriptions and domain notes (units, conventions).
  - Inline comments for non-obvious symbolic steps and performance choices.
- Error handling:
  - Validate input early (square, symmetric, invertible metric, coordinate count) and provide specific messages.
  - In app.py, prefer specific except clauses; avoid catching broad Exception unless reporting a user-friendly error and logging details.
- Security and parsing:
  - Current approach uses eval with __builtins__={} and a curated safe_dict. This is safer than raw eval but still evaluative.
  - Prefer using sympy.parsing.sympy_parser.parse_expr with standard transformations and a controlled local dict to avoid eval entirely.
  - Keep the whitelist of allowed functions/symbols explicit and minimal.

## 5. Common Patterns
- Separation of concerns: UI (Streamlit) vs core computation (SymPy) with a clean function boundary (compute_christoffel).
- Progress reporting: compute_christoffel accepts an optional progress_callback for granular UI feedback.
- Simplification strategy: Use Basic (sp.cancel(sp.expand(...))) by default to balance performance and readability; Full simplify is optional.
- Symmetry handling:
  - Metrics are symmetrized if the user input is not symmetric: metric = (metric + metric.T)/2.
  - Display deduplication: Only show independent Γ^λ_{μν} using symmetry in lower indices (μν).
- Coordinate mapping and labels:
  - Coordinates are validated against a fixed map {'t','r','theta','phi','x','y','z'}.
  - LaTeX uses special symbols for θ, φ for readability.

## 6. Do's and Don'ts
- Do
  - Validate metric shape and coordinate count before computation.
  - Warn on large dimensions; keep practical limit at ≤4 for responsiveness.
  - Keep the allowed function set small and documented (sin, cos, tan, exp, log, sqrt, etc.).
  - Preserve clear separation between UI and math core; keep compute_christoffel pure (no Streamlit calls).
  - Add unit tests for known metrics and properties (flatness, symmetry, invertibility errors).
  - Pin dependency versions and use a virtual environment.
  - Consider caching expensive symbolic operations if inputs repeat (Streamlit cache with careful keying by input strings).
- Don't
  - Don’t use unrestricted eval on user input.
  - Don’t silently swallow exceptions; surface clear actionable messages.
  - Don’t run heavy simplify() on every component by default; prefer Basic/None unless explicitly requested.
  - Don’t mix UI code into the computation module.
  - Don’t introduce unused imports (e.g., remove ast from app.py if not used).

## 7. Tools & Dependencies
- Key libraries
  - streamlit: UI framework for rapid interactive apps.
  - sympy: Symbolic mathematics for tensors and derivatives.
- Suggested requirements.txt (pin to tested versions):
  - streamlit>=1.30,<2
  - sympy>=1.12,<2
  - pytest>=7.0,<9 (for tests)
  - ruff>=0.5 (lint), black>=24.0 (format), isort>=5.13 (imports) [dev]
- Setup instructions
  - python -m venv venv
  - Windows: venv\Scripts\activate
  - pip install -r requirements.txt
  - Run: streamlit run app.py
  - Update requirements after changes: pip freeze > requirements.txt (or maintain by hand with pins)

## 8. Other Notes
- For LLM-generated changes:
  - Maintain the safe parsing approach (or move to SymPy parse_expr); never expand the eval scope beyond the curated whitelist.
  - Preserve the progress_callback pattern in compute_christoffel for UI responsiveness.
  - Keep display deduplication logic and LaTeX formatting conventions for θ, φ.
  - Ensure new metrics added to METRIC_LIBRARY include: metric (string matrix), coords (comma list), and a concise description.
  - Consider adding typing to compute_christoffel and factoring out small helpers for parsing to improve testability.
  - Be mindful of algorithmic complexity: nested loops are O(n^4) work overall due to derivatives/inversion; avoid >4D by default.
