# End-to-end smoke for wclss: REML and EFS converge, truth is recovered,
# and the mgcv hook surface functions.

conv_ok <- function(b) {
  is.null(b$outer.info$conv) || identical(b$outer.info$conv, "full convergence")
}

rwc <- function(n, mu, rho) {
  # exact: wrapped Cauchy = wrapped Cauchy(mu, gamma = -log(rho))
  th <- mu + (-log(rho)) * rcauchy(n)
  atan2(sin(th), cos(th))
}

test_that("wclss REML fit converges, recovers truth, hooks work", {
  set.seed(22)
  n <- 400
  x <- runif(n)
  mu_true <- 2 * atan(1.4 * sin(2 * pi * x))
  rho_true <- plogis(0.5 + 0.8 * cos(2 * pi * x))
  dat <- data.frame(y = rwc(n, mu_true, rho_true), x = x)

  b <- mgcv::gam(list(y ~ s(x, k = 10), ~ s(x, k = 10)),
                 family = wclss(), data = dat, method = "REML")
  expect_true(conv_ok(b))
  expect_true(is.finite(logLik(b)))
  expect_true(is.finite(AIC(b)))

  fv <- fitted(b)
  expect_equal(ncol(fv), 2)
  expect_true(all(abs(fv[, 1]) < pi))
  expect_true(all(fv[, 2] > 0 & fv[, 2] < 1))
  ang_err <- abs(atan2(sin(fv[, 1] - mu_true), cos(fv[, 1] - mu_true)))
  expect_lt(mean(ang_err), 0.2)
  expect_lt(median(abs(fv[, 2] - rho_true)), 0.1)

  sm <- summary(b)
  expect_equal(length(sm$chi.sq), 2L)
  expect_gt(sm$dev.expl, 0)
  expect_lte(sm$dev.expl, 1)
  nd <- data.frame(x = seq(0.05, 0.95, length.out = 11))
  pr <- predict(b, newdata = nd, type = "response")
  expect_equal(dim(pr), c(11L, 2L))
  expect_true(all(is.finite(residuals(b))))
  expect_true(all(is.finite(residuals(b, type = "pearson"))))
  expect_true(all(abs(residuals(b, type = "response")) <= pi))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(b, pages = 1))
  expect_no_error(invisible(capture.output(mgcv::gam.check(b))))

  sim <- b$family$rd(fv, wt = rep(1, n), scale = 1)
  expect_true(all(abs(sim) <= pi))

  b_efs <- mgcv::gam(list(y ~ s(x, k = 10), ~ s(x, k = 10)),
                     family = wclss(), data = dat, optimizer = "efs")
  expect_true(conv_ok(b_efs))
})

test_that("wclss parametric-only fit works", {
  set.seed(23)
  n <- 400
  x <- runif(n)
  y <- rwc(n, 2 * atan(-0.3 + 1.1 * x), plogis(0.2 + 0.9 * x))
  b <- mgcv::gam(list(y ~ x, ~ x), family = wclss(),
                 data = data.frame(y = y, x = x), method = "REML")
  expect_true(conv_ok(b))
  expect_equal(length(coef(b)), 4L)
  expect_true(all(is.finite(coef(b))))
})
