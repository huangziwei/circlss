## Bessel-ratio helpers for the von Mises family, transcribed from the
## FD-verified reference implementations in pycircstat2/utils.py (A1,
## A1prime, A1prime2, A1prime3). Exponentially scaled Bessel functions keep
## every expression finite at large kappa, where I0/I1 themselves overflow.

A1 <- function(kappa) {
  besselI(kappa, 1, expon.scaled = TRUE) / besselI(kappa, 0, expon.scaled = TRUE)
}

A1prime <- function(kappa, a1 = A1(kappa)) {
  ## A1'(k) = 1 - A1(k)/k - A1(k)^2; A1(k)/k -> 1/2 as k -> 0 is the
  ## removable singularity, filled in to avoid 0/0.
  r <- a1 / kappa
  r[kappa == 0] <- 0.5
  1 - r - a1 * a1
}

A1prime2 <- function(kappa, a1 = A1(kappa), d1 = A1prime(kappa, a1)) {
  ## A1''(k) = -A1'/k + A1/k^2 - 2 A1 A1'. The two 1/k terms cancel
  ## catastrophically as k -> 0 (true value ~ -3k/8), so below k < 0.01 a
  ## Maclaurin series is used; both branches agree to ~1e-12 at the switch.
  rec <- -d1 / kappa + a1 / kappa^2 - 2 * a1 * d1
  k2 <- kappa * kappa
  series <- kappa * (-3 / 8 + k2 * (5 / 24 - k2 * (77 / 1024)))
  out <- ifelse(kappa < 0.01, series, rec)
  out
}

A1prime3 <- function(kappa, a1 = A1(kappa), d1 = A1prime(kappa, a1),
                     d2 = A1prime2(kappa, a1, d1)) {
  ## A1'''(k) = -A1''/k + 2 A1'/k^2 - 2 A1/k^3 - 2 A1'^2 - 2 A1 A1'';
  ## same removable cancellation as A1prime2, series below k < 0.01.
  rec <- -d2 / kappa + 2 * d1 / kappa^2 - 2 * a1 / kappa^3 -
    2 * d1 * d1 - 2 * a1 * d2
  k2 <- kappa * kappa
  series <- -3 / 8 + k2 * (5 / 8 - k2 * (385 / 1024))
  out <- ifelse(kappa < 0.01, series, rec)
  out
}
