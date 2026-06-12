# End-to-end smoke for pnlss: REML and EFS converge, the mgcv hook surface
# functions, and -- the distinguishing capability -- the fitted mean
# direction can wind around the circle (winding number 1), which the
# tan-half vmlss family cannot represent.

conv_ok <- function(b) {
  is.null(b$outer.info$conv) || identical(b$outer.info$conv, "full convergence")
}

test_that("pnlss REML fit converges and the hook surface works", {
  set.seed(11)
  n <- 300
  x <- runif(n)
  m1 <- 1.8 * sin(2 * pi * x) + 0.5
  m2 <- 1.2 * cos(2 * pi * x) - 0.3
  y <- atan2(m2 + rnorm(n), m1 + rnorm(n))
  dat <- data.frame(y = y, x = x)

  b <- mgcv::gam(list(y ~ s(x, k = 10), ~ s(x, k = 10)),
                 family = pnlss(), data = dat, method = "REML")
  expect_true(conv_ok(b))
  expect_true(is.finite(logLik(b)))
  expect_true(is.finite(AIC(b)))

  # fitted values: the two Cartesian mean components
  fv <- fitted(b)
  expect_equal(ncol(fv), 2)
  dirhat <- atan2(fv[, 2], fv[, 1])
  dirtrue <- atan2(m2, m1)
  ang_err <- abs(atan2(sin(dirhat - dirtrue), cos(dirhat - dirtrue)))
  expect_lt(mean(ang_err), 0.2)

  sm <- summary(b)
  expect_equal(length(sm$chi.sq), 2L)
  expect_gt(sm$dev.expl, 0)
  expect_lte(sm$dev.expl, 1)
  nd <- data.frame(x = seq(0.05, 0.95, length.out = 11))
  pr <- predict(b, newdata = nd, type = "response")
  expect_equal(dim(pr), c(11L, 2L))
  expect_true(all(is.finite(residuals(b))))
  expect_true(all(abs(residuals(b, type = "response")) <= pi))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(b, pages = 1))
  expect_no_error(invisible(capture.output(mgcv::gam.check(b))))

  sim <- b$family$rd(fv, wt = rep(1, n), scale = 1)
  expect_true(all(abs(sim) <= pi))

  b_efs <- mgcv::gam(list(y ~ s(x, k = 10), ~ s(x, k = 10)),
                     family = pnlss(), data = dat, optimizer = "efs")
  expect_true(conv_ok(b_efs))
})

test_that("pnlss can fit a winding mean direction (theta ~ phi)", {
  set.seed(12)
  n <- 300
  phi <- runif(n, -pi, pi)
  y <- atan2(2 * sin(phi) + rnorm(n), 2 * cos(phi) + rnorm(n))
  dat <- data.frame(y = y, phi = phi)

  b <- mgcv::gam(list(y ~ s(phi, bs = "cc", k = 10),
                      ~ s(phi, bs = "cc", k = 10)),
                 family = pnlss(), data = dat, method = "REML",
                 knots = list(phi = c(-pi, pi)))
  expect_true(conv_ok(b))

  phig <- seq(-pi, pi, length.out = 73)
  pr <- predict(b, newdata = data.frame(phi = phig), type = "response")
  dirhat <- atan2(pr[, 2], pr[, 1])
  # tracks the identity map around the whole circle...
  err <- abs(atan2(sin(dirhat - phig), cos(dirhat - phig)))
  expect_lt(max(err), 0.35)
  expect_lt(mean(err), 0.12)
  # ...i.e. winding number 1 (sum of wrapped increments = 2*pi)
  dd <- diff(dirhat)
  winding <- sum(atan2(sin(dd), cos(dd))) / (2 * pi)
  expect_equal(winding, 1, tolerance = 0.05)
})

test_that("pnlss parametric-only fit works", {
  set.seed(13)
  n <- 400
  x <- runif(n)
  y <- atan2(-0.4 + 1.5 * x + rnorm(n), 0.8 + 1.2 * x + rnorm(n))
  b <- mgcv::gam(list(y ~ x, ~ x), family = pnlss(),
                 data = data.frame(y = y, x = x), method = "REML")
  expect_true(conv_ok(b))
  expect_equal(length(coef(b)), 4L)
  expect_true(all(is.finite(coef(b))))
})
