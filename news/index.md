# Changelog

## circlss 0.0.3

- [`wclss()`](https://huangziwei.github.io/circlss/reference/wclss.md):
  wrapped Cauchy family — the heavy-tailed counterpart of the von Mises
  (sharply peaked, fat circular tails; the robust choice under angular
  outliers). Tan-half location link, logit-linked mean resultant length
  `rho`; full fourth-order derivatives (Newton REML and EFS). The
  log-density uses the cancellation-free denominator
  `(1-rho)^2 + 4*rho*sin^2(d/2)`, exact as `rho -> 1`.
- The differential battery grows to twelve cases; all green.
- The parity instrument now runs both engines at tightened convergence
  (`conv.tol = 1e-11`; production defaults are untouched). This
  collapsed what had looked like cross-engine noise classes by up to
  five orders of magnitude — the former “EIGEN-class” spread (~3e-4) was
  the engines’ default 1e-6 stopping tolerance, not implementation
  disagreement. All twelve cases now agree to \<= 1.4e-7 on every
  quantity (typically 1e-9..1e-12), under a single pinned tolerance of
  2e-6. The remaining floors are genuine: thin-plate coefficient signs
  (eigen-arbitrary, compared sign-aligned) and ~1e-8 special-function
  ulps (scipy log_ndtr vs R pnorm) amplified at flat-gradient optima.

## circlss 0.0.2

- [`pnlss()`](https://huangziwei.github.io/circlss/reference/pnlss.md):
  projected normal family — the response is the angle of a bivariate
  normal with mean `(mu1, mu2)` and identity covariance; both Cartesian
  mean components carry identity links and their own linear predictors.
  No tan-half branch cut: the fitted mean direction can cross any angle
  and can *wind* around the circle, making `pnlss` the natural family
  for circular–circular regression with rotation-type association
  (`theta ~ phi`). Full fourth-order derivatives (Newton REML and EFS).
- The differential battery against pycircstat2 doubles to eight cases:
  the four `pnlss` cases (parametric, thin plate, cyclic-with-winding,
  small-n) agree to the same tolerance classes as `vmlss` — near machine
  precision on deterministic bases (the winding cyclic case to ~1e-10),
  optimizer-stopping tolerance on thin-plate ones. The TIGHT class’s
  coefficient/curve thresholds were recalibrated from 1e-9 to 1e-7 to
  cover the flatter parametric projected-normal likelihood
  (log-likelihood still agrees to 1e-13 there).

## circlss 0.0.1

Initial release.

- [`vmlss()`](https://huangziwei.github.io/circlss/reference/vmlss.md):
  von Mises location-scale general family for
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) — a Fisher–Lee
  tan-half link for the mean direction and a log link for the
  concentration, each parameter with its own linear predictor and
  smooths. Log-likelihood derivatives to fourth order, so full Newton
  `method = "REML"` works (as does `optimizer = "efs"`).
- Differentially tested against
  [pycircstat2](https://github.com/circstat/pycircstat2): identical
  models on identical data (Newton-REML both sides) agree to machine
  precision for parametric and cyclic-spline fits, and to optimizer
  stopping tolerance (~3e-5, fitted curves within 2e-6 radians) for
  thin-plate fits. Tolerances are pinned in the test suite.
- The derivative algebra is transcribed from pycircstat2’s
  finite-difference-verified implementations, with R-side FD checks of
  the link chain and the assembled gradient/Hessian.
