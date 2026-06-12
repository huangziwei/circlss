"""Generate the v0.0.1 parity battery (plan section 5.2).

Data are generated ONCE here (fixed seed) and exchanged by file: R and
Python both read the same CSV, never re-simulate. Each case writes
data/<case>.csv with the response in radians (numpy vonmises branch,
(-pi, pi]) and the covariate.
"""

import csv
import pathlib

import numpy as np

OUT = pathlib.Path(__file__).parent / "data"
OUT.mkdir(exist_ok=True)

rng = np.random.default_rng(20260612)


def write(name, cols):
    path = OUT / f"{name}.csv"
    keys = list(cols)
    n = len(cols[keys[0]])
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for i in range(n):
            w.writerow([format(float(cols[k][i]), ".17g") for k in keys])
    print(f"wrote {path} (n={n})")


# A. linear covariate in both linear predictors (no smoothing: coefficient
#    parity is exact-model parity)
n = 256
x = rng.uniform(0.0, 1.0, n)
mu = 2.0 * np.arctan(-0.3 + 1.1 * x)
kappa = np.exp(0.9 + 0.6 * x)
write("lin", {"y": rng.vonmises(mu, kappa), "x": x})

# B. smooth covariate in both LPs, default thin-plate basis
n = 256
x = rng.uniform(0.0, 1.0, n)
mu = 2.0 * np.arctan(1.4 * np.sin(2.0 * np.pi * x))
kappa = np.exp(1.0 + 0.8 * np.cos(2.0 * np.pi * x))
write("smooth", {"y": rng.vonmises(mu, kappa), "x": x})

# C. circular covariate through a cyclic smooth, knots pinned to (-pi, pi)
n = 300
phi = rng.uniform(-np.pi, np.pi, n)
mu = 2.0 * np.arctan(0.9 * np.sin(phi))
kappa = np.exp(1.1 + 0.5 * np.cos(phi))
write("cyclic", {"y": rng.vonmises(mu, kappa), "phi": phi})

# D. small-n variant of B
n = 80
x = rng.uniform(0.0, 1.0, n)
mu = 2.0 * np.arctan(1.4 * np.sin(2.0 * np.pi * x))
kappa = np.exp(1.0 + 0.8 * np.cos(2.0 * np.pi * x))
write("small", {"y": rng.vonmises(mu, kappa), "x": x})

# ---- projected normal battery (v0.0.2) -- appended AFTER the vmlss cases
# so the rng stream above is unchanged and the frozen vmlss fixtures stay
# byte-identical. theta = atan2(mu2 + z2, mu1 + z1), z ~ N(0, 1).


def rpn(mu1, mu2):
    return np.arctan2(mu2 + rng.standard_normal(mu2.shape),
                      mu1 + rng.standard_normal(mu1.shape))


# E. linear covariate in both Cartesian mean components
n = 256
x = rng.uniform(0.0, 1.0, n)
write("pn_lin", {"y": rpn(0.8 + 1.2 * x, -0.4 + 1.5 * x), "x": x})

# F. smooth covariate, thin-plate basis
n = 256
x = rng.uniform(0.0, 1.0, n)
write("pn_smooth", {"y": rpn(1.8 * np.sin(2.0 * np.pi * x) + 0.5,
                             1.2 * np.cos(2.0 * np.pi * x) - 0.3),
                    "x": x})

# G. the winding case: mean direction = phi exactly (winding number 1),
#    constant concentration -- representable by pnlss, not by vmlss
n = 300
phi = rng.uniform(-np.pi, np.pi, n)
write("pn_cyclic", {"y": rpn(2.0 * np.cos(phi), 2.0 * np.sin(phi)),
                    "phi": phi})

# H. small-n variant of F
n = 80
x = rng.uniform(0.0, 1.0, n)
write("pn_small", {"y": rpn(1.8 * np.sin(2.0 * np.pi * x) + 0.5,
                            1.2 * np.cos(2.0 * np.pi * x) - 0.3),
                   "x": x})
