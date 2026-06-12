# Circular-circular regression showcase: theta ~ s(phi, bs="cc") with
# vmlss(), drawn on a torus (predictor angle around the ring, response
# angle around the tube). Writes man/figures/README-torus.png; the same
# code lives in vignettes/articles/circular-circular-regression.Rmd.

library(mgcv)
library(circlss)

## ---- simulate circular-circular data (same DGP family as the parity
## battery's cyclic case) ----
rvm <- function(n, mu, kappa) {
  mu <- rep_len(mu, n); kappa <- rep_len(kappa, n); out <- numeric(n)
  for (i in seq_len(n)) {
    k <- kappa[i]
    a <- 1 + sqrt(1 + 4 * k * k); b <- (a - sqrt(2 * a)) / (2 * k)
    r <- (1 + b * b) / (2 * b)
    repeat {
      z <- cos(pi * runif(1)); f <- (1 + r * z) / (r + z)
      cc <- k * (r - f); u2 <- runif(1)
      if (cc * (2 - cc) - u2 > 0 || log(cc / u2) + 1 - cc >= 0) {
        out[i] <- sign(runif(1) - 0.5) * acos(max(min(f, 1), -1)) + mu[i]
        break
      }
    }
  }
  atan2(sin(out), cos(out))
}

set.seed(20260612)
n <- 200
phi <- runif(n, -pi, pi)
mu_true <- 2 * atan(1.6 * sin(phi))
kappa_true <- exp(1.4 + 0.6 * cos(phi))
theta <- rvm(n, mu_true, kappa_true)
dat <- data.frame(theta = theta, phi = phi)

## ---- fit: distributional von Mises, cyclic smooths in both parameters ----
b <- gam(list(theta ~ s(phi, bs = "cc", k = 10),
                    ~ s(phi, bs = "cc", k = 10)),
         family = vmlss(), data = dat, method = "REML",
         knots = list(phi = c(-pi, pi)))

phig <- seq(-pi, pi, length.out = 400)
pred <- predict(b, newdata = data.frame(phi = phig), type = "response")
mu_hat <- pred[, 1]

## ---- torus drawing helpers ----
torus_xyz <- function(phi, theta, R = 1.9, r = 0.95) {
  list(x = (R + r * cos(theta)) * cos(phi),
       y = (R + r * cos(theta)) * sin(phi),
       z = r * sin(theta))
}

# view-space depth from the persp transformation matrix (trans3d gives
# only the projected x/y; the z/w component orders front vs back)
depth3d <- function(x, y, z, pm) {
  p <- cbind(x, y, z, 1) %*% pm
  p[, 3] / p[, 4]
}

draw_torus <- function(dat, mu_hat, phig, R = 1.9, r = 0.95,
                       theta_view = 45, phi_view = 38) {
  op <- par(mar = c(0.2, 0.2, 0.2, 0.2))
  on.exit(par(op))
  lim <- R + r + 0.1
  pm <- persp(x = c(-lim, lim), y = c(-lim, lim),
              z = matrix(c(-r, -r, r, r), 2, 2),
              zlim = c(-r - 0.5, r + 0.5),
              theta = theta_view, phi = phi_view, d = 4,
              scale = FALSE, expand = 1,
              col = NA, border = NA, box = FALSE, axes = FALSE)

  ## wireframe: tube circles at fixed phi, ring circles at fixed theta
  thd <- seq(-pi, pi, length.out = 80)
  for (p in seq(-pi, pi, length.out = 49)[-1]) {
    w <- torus_xyz(rep(p, length(thd)), thd, R, r)
    lines(trans3d(w$x, w$y, w$z, pm), col = "gray88", lwd = 0.6)
  }
  phd <- seq(-pi, pi, length.out = 160)
  for (t in seq(-pi, pi, length.out = 25)[-1]) {
    w <- torus_xyz(phd, rep(t, length(phd)), R, r)
    lines(trans3d(w$x, w$y, w$z, pm), col = "gray88", lwd = 0.6)
  }

  ## outer equator (theta = 0): the zero-response reference line
  w <- torus_xyz(phd, rep(0, length(phd)), R, r)
  lines(trans3d(w$x, w$y, w$z, pm), col = "gray55", lty = 3, lwd = 0.9)

  ## data on the torus surface, depth-faded (near = opaque, far = faint)
  w <- torus_xyz(dat$phi, dat$theta, R, r)
  dp <- depth3d(w$x, w$y, w$z, pm)
  a <- 0.15 + 0.7 * (dp - min(dp)) / diff(range(dp))
  pt <- trans3d(w$x, w$y, w$z, pm)
  ord <- order(dp)  # paint far points first
  points(pt$x[ord], pt$y[ord], pch = 19,
         cex = 0.4 + 0.25 * a[ord], col = adjustcolor("#2c5aa0", 1)[1] |>
           (\(cc) sapply(a[ord], function(ai) adjustcolor(cc, ai)))())

  ## fitted mean-direction curve mu-hat(phi), winding over the tube:
  ## far half faint, near half full -- painted in depth order
  w <- torus_xyz(phig, mu_hat, R, r)
  dp <- depth3d(w$x, w$y, w$z, pm)
  cv <- trans3d(w$x, w$y, w$z, pm)
  a <- 0.25 + 0.75 * (dp - min(dp)) / diff(range(dp))
  seg <- data.frame(x0 = head(cv$x, -1), y0 = head(cv$y, -1),
                    x1 = tail(cv$x, -1), y1 = tail(cv$y, -1),
                    d = head(dp, -1), a = head(a, -1))
  seg <- seg[order(seg$d), ]
  segments(seg$x0, seg$y0, seg$x1, seg$y1,
           col = sapply(seg$a, function(ai) adjustcolor("#c0392b", ai)),
           lwd = 2.2 + 1.6 * seg$a, lend = 1)
  invisible(pm)
}

## ---- render ----
dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
png("man/figures/README-torus.png", width = 1500, height = 1180, res = 220)
draw_torus(dat, mu_hat, phig)
invisible(dev.off())
cat("wrote man/figures/README-torus.png\n")
cat("fit: conv =", b$outer.info$conv, " edf =", round(sum(b$edf), 2), "\n")
