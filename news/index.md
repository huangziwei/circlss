# Changelog

## circlss 0.0.1

Initial release.

- [`vmlss()`](https://circstat.github.io/circlss/reference/vmlss.md):
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
