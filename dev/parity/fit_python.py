"""Fit the parity battery with pycircstat2 (the oracle side).

Every case is fit with CLRegression's gam backend (hea), family=vmlss,
method="REML" -- Newton-REML, the same optimizer the R side pins (vM is
Tier 1 with l1-l4 on both sides). Results go to out/<case>_py.json.
"""

import inspect
import json
import pathlib

import numpy as np
import polars as pl

from pycircstat2.distributions import CircularLL, vmlss
from pycircstat2.regression import CLRegression

# ---------------------------------------------------------------------------
# TEMPORARY SHIM (remove once pycircstat2 is fixed upstream): hea's engine
# calls family.postproc(y, prior_weights=..., fitted=..., linear_predictors=
# ..., offset=..., intercept=...) but CircularLL.postproc(self, y, fitted)
# accepts only the first two -- every gam-backend fit crashes (pycircstat2's
# own test_katojones_family_dispatch fails the same way, 2026-06-12).
# Upstream one-liner in pycircstat2/distributions.py:
#   def postproc(self, y, fitted) -> dict:
# becomes
#   def postproc(self, y, fitted=None, **kwargs) -> dict:
# The shim is conditional on the narrow signature, so it self-disables once
# the upstream fix lands.
_pp = CircularLL.postproc
if not any(p.kind is inspect.Parameter.VAR_KEYWORD
           for p in inspect.signature(_pp).parameters.values()):
    def _postproc(self, y, fitted=None, **kwargs):
        return _pp(self, y, fitted)
    CircularLL.postproc = _postproc
# ---------------------------------------------------------------------------

HERE = pathlib.Path(__file__).parent
OUT = HERE / "out"
OUT.mkdir(exist_ok=True)

CASES = {
    "lin": {
        "formulas": ["y ~ x", "~ x"],
        "knots": None,
        "var": "x",
        "grid": np.linspace(0.0, 1.0, 101),
    },
    "smooth": {
        "formulas": ["y ~ s(x, k=10)", "~ s(x, k=10)"],
        "knots": None,
        "var": "x",
        "grid": np.linspace(0.0, 1.0, 101),
    },
    "cyclic": {
        "formulas": ["y ~ s(phi, bs='cc', k=10)", "~ s(phi, bs='cc', k=10)"],
        "knots": {"phi": [-np.pi, np.pi]},
        "var": "phi",
        "grid": np.linspace(-np.pi, np.pi, 101),
    },
    "small": {
        "formulas": ["y ~ s(x, k=8)", "~ s(x, k=8)"],
        "knots": None,
        "var": "x",
        "grid": np.linspace(0.0, 1.0, 101),
    },
}

for name, spec in CASES.items():
    df = pl.read_csv(HERE / "data" / f"{name}.csv")
    m = CLRegression(
        spec["formulas"], data=df, family=vmlss, method="REML",
        knots=spec["knots"],
    )
    g = m.gam_fit
    res = m.result
    grid_df = pl.DataFrame({spec["var"]: spec["grid"]})
    pp = m.predict_params(grid_df)

    edf_by = res.get("edf_by_smooth") or {}
    out = {
        "case": name,
        "formulas": spec["formulas"],
        "coef_names": list(g.bhat.columns),
        "coef": [float(v) for v in g.bhat.row(0)],
        "lpi": [[int(i) for i in np.asarray(ix)] for ix in g.lpi],
        "sp": [float(s) for s in np.atleast_1d(getattr(g, "sp", np.zeros(0)))],
        "edf_total": float(res["edf_total"]),
        "edf_by_smooth": {k: float(v) for k, v in edf_by.items()},
        "loglik": float(res["log_likelihood"]),
        "reml": float(res["reml"]),
        "grid": [float(v) for v in spec["grid"]],
        "mu_grid": [float(v) for v in pp["mu"]],
        "kappa_grid": [float(v) for v in pp["kappa"]],
    }
    path = OUT / f"{name}_py.json"
    path.write_text(json.dumps(out, indent=1))
    print(f"{name}: loglik={out['loglik']:.6f} edf={out['edf_total']:.4f} "
          f"sp={['%.5g' % s for s in out['sp']]} -> {path.name}")
