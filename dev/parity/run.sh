#!/bin/sh
# Side-by-side parity runner (plan section 5.2): generate shared data once,
# fit the same models in pycircstat2 (Newton-REML via hea) and in circlss
# (Newton-REML via mgcv), compare within pinned tolerances, and freeze the
# fixtures the hermetic testthat parity test asserts against.
#
# Usage: ./run.sh [path-to-pycircstat2] (defaults to ../../../pycircstat2)
set -e
cd "$(dirname "$0")"

PYCS2="${1:-$(cd ../../../pycircstat2 && pwd)}"
PY="$PYCS2/.venv/bin/python"
[ -x "$PY" ] || { echo "no venv python at $PY" >&2; exit 1; }

echo "== generate shared data (fixed seed) =="
"$PY" gen_data.py

echo "== fit: pycircstat2 (oracle) =="
"$PY" fit_python.py

echo "== fit: circlss =="
Rscript fit_r.R .

echo "== compare + freeze fixtures =="
"$PY" compare.py --freeze
