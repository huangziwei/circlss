# Differential test for wclss against pycircstat2, on the frozen battery
# (conventions as in test-vmlss-parity.R: REML/2, sign-aligned
# coefficients).
#
# The parity instrument runs both engines at tightened convergence
# (conv.tol = 1e-11) so stopping points don't masquerade as disagreement;
# a single tolerance class (2e-6, ~10x the worst observed cross-engine
# floor) covers every quantity. Coefficients compare sign-aligned (tp
# basis signs are eigen-arbitrary). See dev/parity/compare.py.

skip_if_not_installed("jsonlite")

TOL <- list(coef = 2e-6, sp_log = 2e-6, edf_total = 2e-6,
            loglik = 2e-6, reml = 2e-6, mu = 2e-6, rho = 2e-6)
ctl <- mgcv::gam.control(epsilon = 1e-10, newton = list(conv.tol = 1e-11))

cases <- list(
  wc_lin = list(formula = list(y ~ x, ~ x), knots = NULL, var = "x"),
  wc_smooth = list(formula = list(y ~ s(x, k = 10), ~ s(x, k = 10)),
                   knots = NULL, var = "x"),
  wc_cyclic = list(formula = list(y ~ s(phi, bs = "cc", k = 10),
                                  ~ s(phi, bs = "cc", k = 10)),
                   knots = list(phi = c(-pi, pi)), var = "phi"),
  wc_small = list(formula = list(y ~ s(x, k = 8), ~ s(x, k = 8)),
                  knots = NULL, var = "x")
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

    b <- suppressWarnings(mgcv::gam(cs$formula, family = wclss(),
                   data = dat, method = "REML", knots = cs$knots,
                   control = ctl))
    expect_true(is.null(b$outer.info$conv) ||
                  identical(b$outer.info$conv, "full convergence"))

    co <- unname(coef(b))
    expect_equal(length(co), length(py$coef))
    expect_lt(max(pmin(abs(co - py$coef), abs(co + py$coef))), TOL$coef)

    expect_equal(length(b$sp), length(py$sp))
    if (length(py$sp) > 0)
      expect_lt(max(abs(log(unname(b$sp) / py$sp))), TOL$sp_log)
    expect_lt(abs(sum(b$edf) - py$edf_total), TOL$edf_total)
    expect_lt(abs(as.numeric(logLik(b)) - py$loglik), TOL$loglik)
    expect_lt(abs(as.numeric(b$gcv.ubre) - py$reml / 2), TOL$reml)

    nd <- stats::setNames(data.frame(py$grid), cs$var)
    pr <- predict(b, newdata = nd, type = "response")
    expect_lt(max(wrap_abs(pr[, 1] - py$mu_grid)), TOL$mu)
    expect_lt(max(abs(pr[, 2] - py$rho_grid)), TOL$rho)
  })
}
