# FD checks of the tan-half link derivative chain. Each analytic derivative
# is checked against a central difference of the previous one, so every
# level gets one FD step (accurate to ~1e-9), not a noisy high-order stencil.

fd <- function(f, x, h = 1e-5) (f(x + h) - f(x - h)) / (2 * h)

test_that("tanhalf link derivatives are mutually consistent", {
  lk <- circlss:::tanhalf.link()
  mu <- seq(-2.8, 2.8, length.out = 41)

  # linkinv inverts linkfun
  expect_equal(lk$linkinv(lk$linkfun(mu)), mu, tolerance = 1e-12)

  # mu.eta(eta) = d linkinv / d eta
  eta <- seq(-8, 8, length.out = 41)
  expect_equal(lk$mu.eta(eta), fd(lk$linkinv, eta), tolerance = 1e-8)

  # g'(mu) implied by mu.eta: g' = 1 / mu.eta(g(mu))
  g1 <- function(m) 1 / lk$mu.eta(lk$linkfun(m))
  expect_equal(g1(mu), fd(lk$linkfun, mu), tolerance = 1e-7)

  # d2link = d g'/d mu, d3link = d d2link/d mu, d4link = d d3link/d mu
  expect_equal(lk$d2link(mu), fd(g1, mu), tolerance = 1e-6)
  expect_equal(lk$d3link(mu), fd(lk$d2link, mu), tolerance = 1e-6)
  expect_equal(lk$d4link(mu), fd(lk$d3link, mu), tolerance = 1e-5)
})

test_that("A1 derivative chain is FD-consistent, including the series branch", {
  a1 <- circlss:::A1
  d1 <- function(k) circlss:::A1prime(k)
  d2 <- function(k) circlss:::A1prime2(k)
  d3 <- function(k) circlss:::A1prime3(k)
  # d1 down to tiny kappa (step kept below kappa so FD stays in-domain)
  kappa <- c(1e-6, 0.005, 0.009, 0.011, 0.05, 0.5, 2, 10, 100, 700)
  expect_equal(d1(kappa), fd(a1, kappa, 1e-7), tolerance = 1e-5)
  # d2/d3 across the series/recurrence switch at kappa = 0.01
  kappa <- c(0.005, 0.009, 0.011, 0.05, 0.5, 2, 10, 100, 700)
  expect_equal(d2(kappa), fd(d1, kappa, 1e-5), tolerance = 1e-6)
  expect_equal(d3(kappa), fd(d2, kappa, 1e-4), tolerance = 1e-5)
})
