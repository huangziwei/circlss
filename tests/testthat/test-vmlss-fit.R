# End-to-end smoke: a REML fit on simulated data converges, recovers the
# true curves to within statistical noise, and the whole mgcv hook surface
# (summary, predict, plot, gam.check, AIC, residuals, rd) functions.

# general-family fits record outer convergence in outer.info$conv; fully
# parametric fits have no outer loop (outer.info is NULL)
conv_ok <- function(b) {
  is.null(b$outer.info$conv) || identical(b$outer.info$conv, "full convergence")
}

rvm <- function(n, mu, kappa) {
  # Best & Fisher (1979) rejection sampler, scalar-recycled
  mu <- rep_len(mu, n)
  kappa <- rep_len(kappa, n)
  out <- numeric(n)
  for (i in seq_len(n)) {
    k <- kappa[i]
    if (k < 1e-9) {
      out[i] <- runif(1, -pi, pi)
      next
    }
    a <- 1 + sqrt(1 + 4 * k * k)
    b <- (a - sqrt(2 * a)) / (2 * k)
    r <- (1 + b * b) / (2 * b)
    repeat {
      z <- cos(pi * runif(1))
      f <- (1 + r * z) / (r + z)
      cc <- k * (r - f)
      u2 <- runif(1)
      if (cc * (2 - cc) - u2 > 0 || log(cc / u2) + 1 - cc >= 0) {
        out[i] <- sign(runif(1) - 0.5) * acos(max(min(f, 1), -1)) + mu[i]
        break
      }
    }
  }
  atan2(sin(out), cos(out))
}

test_that("vmlss REML fit converges and recovers the truth", {
  set.seed(20260612)
  n <- 400
  x <- runif(n)
  mu_true <- 2 * atan(1.5 * sin(2 * pi * x))
  kappa_true <- exp(1.2 + 0.8 * cos(2 * pi * x))
  y <- rvm(n, mu_true, kappa_true)
  dat <- data.frame(y = y, x = x)

  b <- mgcv::gam(list(y ~ s(x, k = 10), ~ s(x, k = 10)),
                 family = vmlss(), data = dat, method = "REML")

  expect_true(conv_ok(b))
  expect_true(is.finite(logLik(b)))
  expect_true(is.finite(AIC(b)))

  # fitted values: two response-scale columns, mu in (-pi, pi), kappa > 0
  fv <- fitted(b)
  expect_equal(ncol(fv), 2)
  expect_true(all(abs(fv[, 1]) < pi))
  expect_true(all(fv[, 2] > 0))

  # curve recovery: mean angular error of mu-hat small; kappa-hat within a
  # factor band (kappa is the harder parameter at this n)
  ang_err <- abs(atan2(sin(fv[, 1] - mu_true), cos(fv[, 1] - mu_true)))
  expect_lt(mean(ang_err), 0.15)
  expect_lt(median(abs(log(fv[, 2] / kappa_true))), 0.5)

  # hook surface
  sm <- summary(b)
  expect_equal(length(sm$chi.sq), 2L)
  expect_gt(sm$dev.expl, 0)
  expect_lte(sm$dev.expl, 1)
  pr_link <- predict(b)
  expect_equal(ncol(pr_link), 2)
  nd <- data.frame(x = seq(0.05, 0.95, length.out = 11))
  pr_resp <- predict(b, newdata = nd, type = "response")
  expect_equal(dim(pr_resp), c(11L, 2L))
  expect_true(all(abs(pr_resp[, 1]) < pi))
  expect_true(all(pr_resp[, 2] > 0))

  rsd <- residuals(b)            # deviance
  expect_true(all(is.finite(rsd)))
  expect_true(all(is.finite(residuals(b, type = "pearson"))))
  expect_true(all(abs(residuals(b, type = "response")) <= pi))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_silent(plot(b, pages = 1))
  expect_no_error(invisible(capture.output(mgcv::gam.check(b))))

  # rd: simulates plausible angles
  sim <- b$family$rd(fv, wt = rep(1, n), scale = 1)
  expect_true(all(abs(sim) <= pi))
})

test_that("parametric-only and efs fits work", {
  set.seed(3)
  n <- 250
  x <- runif(n)
  y <- rvm(n, 2 * atan(-0.3 + 1.1 * x), exp(0.9 + 0.6 * x))
  dat <- data.frame(y = y, x = x)

  b_par <- mgcv::gam(list(y ~ x, ~ x), family = vmlss(), data = dat,
                     method = "REML")
  expect_true(conv_ok(b_par))
  expect_equal(length(coef(b_par)), 4L)
  expect_true(all(is.finite(coef(b_par))))

  b_efs <- mgcv::gam(list(y ~ s(x, k = 8), ~ s(x, k = 8)),
                     family = vmlss(), data = dat, optimizer = "efs")
  expect_true(conv_ok(b_efs))
  expect_identical(b_efs$outer.info$conv, "full convergence")
})
