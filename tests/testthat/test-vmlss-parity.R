# Differential test against pycircstat2 (the plan's promotion gate): refit
# the frozen battery and compare against the frozen Python results in
# fixtures/ (written by dev/parity/compare.py --freeze). Hermetic: no
# Python needed at test time.
#
# Cross-engine conventions normalized here (see dev/parity/compare.py):
# - hea's REML criterion is -2*log REML; mgcv's gcv.ubre is -log REML --
#   compare R against Python/2.
# - thin-plate basis column signs are arbitrary (eigen decomposition), so
#   coefficients compare sign-aligned: min(|a-b|, |a+b|). The fitted
#   function is invariant; parametric/cc coefficients match without it.
#
# Tolerance classes pinned at ~10x the noise observed at v0.0.1:
# TIGHT (deterministic bases: parametric, cc) and EIGEN (thin plate).

skip_if_not_installed("jsonlite")

TIGHT <- list(coef = 1e-7, sp_log = 1e-7, edf_total = 1e-7,
              loglik = 1e-8, reml = 1e-7, mu = 1e-7, kappa_rel = 1e-9)
EIGEN <- list(coef = 1e-4, sp_log = 3e-4, edf_total = 3e-4,
              loglik = 3e-4, reml = 1e-6, mu = 3e-5, kappa_rel = 3e-5)

cases <- list(
  lin = list(formula = list(y ~ x, ~ x), knots = NULL, var = "x",
             tol = TIGHT),
  smooth = list(formula = list(y ~ s(x, k = 10), ~ s(x, k = 10)),
                knots = NULL, var = "x", tol = EIGEN),
  cyclic = list(formula = list(y ~ s(phi, bs = "cc", k = 10),
                               ~ s(phi, bs = "cc", k = 10)),
                knots = list(phi = c(-pi, pi)), var = "phi", tol = TIGHT),
  small = list(formula = list(y ~ s(x, k = 8), ~ s(x, k = 8)),
               knots = NULL, var = "x", tol = EIGEN)
)

wrap_abs <- function(d) abs(atan2(sin(d), cos(d)))

for (nm in names(cases)) {
  test_that(paste0("parity vs pycircstat2: ", nm), {
    cs <- cases[[nm]]
    csv <- test_path("fixtures", paste0(nm, ".csv"))
    js <- test_path("fixtures", paste0(nm, "_py.json"))
    skip_if(!file.exists(csv) || !file.exists(js),
            "parity fixtures not generated")
    dat <- read.csv(csv)
    py <- jsonlite::read_json(js, simplifyVector = TRUE)

    b <- mgcv::gam(cs$formula, family = vmlss(), data = dat,
                   method = "REML", knots = cs$knots)
    expect_true(is.null(b$outer.info$conv) ||
                  identical(b$outer.info$conv, "full convergence"))

    co <- unname(coef(b))
    expect_equal(length(co), length(py$coef))
    expect_lt(max(pmin(abs(co - py$coef), abs(co + py$coef))), cs$tol$coef)

    expect_equal(length(b$sp), length(py$sp))
    if (length(py$sp) > 0)
      expect_lt(max(abs(log(unname(b$sp) / py$sp))), cs$tol$sp_log)
    expect_lt(abs(sum(b$edf) - py$edf_total), cs$tol$edf_total)
    expect_lt(abs(as.numeric(logLik(b)) - py$loglik), cs$tol$loglik)
    expect_lt(abs(as.numeric(b$gcv.ubre) - py$reml / 2), cs$tol$reml)

    nd <- stats::setNames(data.frame(py$grid), cs$var)
    pr <- predict(b, newdata = nd, type = "response")
    expect_lt(max(wrap_abs(pr[, 1] - py$mu_grid)), cs$tol$mu)
    expect_lt(max(abs(pr[, 2] - py$kappa_grid) / abs(py$kappa_grid)),
              cs$tol$kappa_rel)
  })
}
