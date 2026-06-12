"""Compare R (circlss) and Python (pycircstat2) parity fits.

Reads out/<case>_py.json and out/<case>_r.json, reports the discrepancy on
every compared quantity, and exits nonzero if any pinned tolerance is
breached. With --freeze, copies the data CSVs and Python JSONs into
tests/testthat/fixtures/ so the testthat parity test is hermetic.

What "exactly the same result" means here (plan section 5.2): both sides
run Newton-REML on identical data; bitwise equality across two languages,
two BLASes and two iterative solvers is not the bar -- agreement orders of
magnitude below any statistically meaningful difference is, with the
thresholds pinned below and treated as regressions thereafter.

Two cross-engine conventions are normalized before comparing:

- REML scale: hea's REML_criterion is -2*log REML (deviance scale); mgcv's
  gcv.ubre is -log REML. R is compared against Python/2 (they then agree to
  ~1e-9).
- Thin-plate coefficient signs: the tp basis is built from an eigen
  decomposition whose column signs are arbitrary (mgcv does not guarantee
  them across LAPACK builds either). Coefficients are compared
  sign-aligned, min(|a-b|, |a+b|) per coefficient; the fitted function is
  invariant. Parametric and cyclic (cc) bases are deterministic, so their
  coefficients also match without alignment -- machine-tight tolerances
  apply (TIGHT class below).
"""

import json
import math
import pathlib
import shutil
import sys

HERE = pathlib.Path(__file__).parent
CASES = ["lin", "smooth", "cyclic", "small"]

# tolerance classes, pinned at ~10x the noise observed at v0.0.1
# (lin/cyclic observed 1e-14..1e-9; smooth/small tp cases observed up to
# ~3e-5 on loglik/sp -- REML-optimum flatness transmits the optimizer
# stopping tolerance into those, while the criterion itself and the fitted
# curves stay far tighter)
TIGHT = {  # deterministic bases: parametric terms, cc cyclic splines
    "coef": 1e-9, "sp_log": 1e-7, "edf_total": 1e-7, "edf_smooth": 1e-7,
    "loglik": 1e-8, "reml": 1e-7, "mu_grid": 1e-9, "kappa_grid_rel": 1e-9,
}
EIGEN = {  # thin-plate bases: eigen-decomposed, float path differs
    "coef": 1e-4, "sp_log": 3e-4, "edf_total": 3e-4, "edf_smooth": 3e-4,
    "loglik": 3e-4, "reml": 1e-6, "mu_grid": 3e-5, "kappa_grid_rel": 3e-5,
}
TOL = {"lin": TIGHT, "smooth": EIGEN, "cyclic": TIGHT, "small": EIGEN}


def wrap_diff(a, b):
    d = a - b
    return abs(math.atan2(math.sin(d), math.cos(d)))


def compare(case):
    py = json.loads((HERE / "out" / f"{case}_py.json").read_text())
    rr = json.loads((HERE / "out" / f"{case}_r.json").read_text())
    tol = TOL[case]
    rows = []
    ok = True

    def check(name, value, t):
        nonlocal ok
        good = value <= t
        ok = ok and good
        rows.append((name, value, t, "ok" if good else "FAIL"))

    if not rr.get("converged", False):
        ok = False
        rows.append(("converged", 0.0, 1.0, "FAIL"))

    assert len(py["coef"]) == len(rr["coef"]), "coefficient count differs"
    assert py["lpi"] == rr["lpi"], "linear predictor index layout differs"
    check("coef",
          max(min(abs(a - b), abs(a + b))
              for a, b in zip(py["coef"], rr["coef"])),
          tol["coef"])

    assert len(py["sp"]) == len(rr["sp"]), "smoothing parameter count differs"
    if py["sp"]:
        check("sp_log", max(abs(math.log(a / b))
                            for a, b in zip(rr["sp"], py["sp"])),
              tol["sp_log"])
    check("edf_total", abs(py["edf_total"] - rr["edf_total"]),
          tol["edf_total"])
    for lab, v in py.get("edf_by_smooth", {}).items():
        if lab in rr.get("edf_by_smooth", {}):
            check(f"edf[{lab}]", abs(v - rr["edf_by_smooth"][lab]),
                  tol["edf_smooth"])
    check("loglik", abs(py["loglik"] - rr["loglik"]), tol["loglik"])
    check("reml", abs(rr["reml"] - py["reml"] / 2.0), tol["reml"])
    check("mu_grid", max(wrap_diff(a, b)
                         for a, b in zip(py["mu_grid"], rr["mu_grid"])),
          tol["mu_grid"])
    check("kappa_grid_rel",
          max(abs(a - b) / abs(b)
              for a, b in zip(py["kappa_grid"], rr["kappa_grid"])),
          tol["kappa_grid_rel"])

    cls = "TIGHT" if tol is TIGHT else "EIGEN"
    print(f"\n== {case} ({cls}) ==")
    for name, value, t, status in rows:
        print(f"  {name:<18} {value:>12.3e}  (tol {t:.0e})  {status}")
    return ok


def freeze():
    fx = HERE.parent.parent / "tests" / "testthat" / "fixtures"
    fx.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        shutil.copy(HERE / "data" / f"{case}.csv", fx / f"{case}.csv")
        shutil.copy(HERE / "out" / f"{case}_py.json", fx / f"{case}_py.json")
    print(f"\nfroze {len(CASES)} cases into {fx}")


if __name__ == "__main__":
    results = {case: compare(case) for case in CASES}
    print("\n" + "-" * 40)
    for case, ok in results.items():
        print(f"{case:<10} {'PASS' if ok else 'FAIL'}")
    if not all(results.values()):
        sys.exit(1)
    if "--freeze" in sys.argv:
        freeze()
