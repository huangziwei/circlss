# README headline figure: circular-circular regression with pnlss() --
# the fitted mean direction winds once around the torus, with a 95%
# delta-method ribbon. Writes man/figures/README-torus.png; the same code
# (and the vmlss oscillation companion) lives in
# vignettes/articles/circular-circular-regression.Rmd.

library(mgcv)
library(circlss)

## ---- simulate: winding C-C data (projected normal draws) ----
set.seed(20260612)
n <- 300
phi <- runif(n, -pi, pi)
dir_true <- phi + 0.6 * sin(phi)            # winds once per cycle
gamma_true <- exp(0.6 + 0.5 * cos(phi))     # concentration scale
theta <- atan2(gamma_true * sin(dir_true) + rnorm(n),
               gamma_true * cos(dir_true) + rnorm(n))
dat <- data.frame(theta = theta, phi = phi)

## ---- fit: projected normal, cyclic smooths in both mean components ----
b <- gam(list(theta ~ s(phi, bs = "cc", k = 10),
                    ~ s(phi, bs = "cc", k = 10)),
         family = pnlss(), data = dat, method = "REML",
         knots = list(phi = c(-pi, pi)))

phig <- seq(-pi, pi, length.out = 400)
pr <- predict(b, newdata = data.frame(phi = phig), type = "response")
dir_hat <- atan2(pr[, 2], pr[, 1])

## 95% band by the delta method: joint Vp through atan2(mu2, mu1)
Xp <- predict(b, newdata = data.frame(phi = phig), type = "lpmatrix")
lpi <- attr(Xp, "lpi")
X1 <- Xp[, lpi[[1]], drop = FALSE]
X2 <- Xp[, lpi[[2]], drop = FALSE]
V <- b$Vp
v11 <- rowSums((X1 %*% V[lpi[[1]], lpi[[1]]]) * X1)
v22 <- rowSums((X2 %*% V[lpi[[2]], lpi[[2]]]) * X2)
v12 <- rowSums((X1 %*% V[lpi[[1]], lpi[[2]]]) * X2)
m1 <- pr[, 1]; m2 <- pr[, 2]; r2 <- m1^2 + m2^2
se_dir <- sqrt(pmax(m2^2 * v11 + m1^2 * v22 - 2 * m1 * m2 * v12, 0)) / r2
dir_lo <- dir_hat - 1.96 * se_dir
dir_hi <- dir_hat + 1.96 * se_dir

## ---- torus drawing helpers (as in the article) ----
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

draw_torus <- function(dat, mu_hat, phig, lo = NULL, hi = NULL,
                       R = 1.9, r = 0.95, theta_view = 45, phi_view = 38) {
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
         cex = 0.4 + 0.25 * a[ord],
         col = sapply(a[ord], function(ai) adjustcolor("#2c5aa0", ai)))

  ## 95% ribbon between lo(phi) and hi(phi): translucent quads painted in
  ## depth order, subdivided across the band so it follows the tube
  if (!is.null(lo) && !is.null(hi)) {
    K <- 4
    i0 <- seq_len(length(phig) - 1)
    quads <- list()
    for (k in seq_len(K)) {
      t0 <- lo + (hi - lo) * (k - 1) / K
      t1 <- lo + (hi - lo) * k / K
      for (i in i0) {
        th <- c(t0[i], t0[i + 1], t1[i + 1], t1[i])
        ph <- c(phig[i], phig[i + 1], phig[i + 1], phig[i])
        w <- torus_xyz(ph, th, R, r)
        quads[[length(quads) + 1]] <-
          list(p = trans3d(w$x, w$y, w$z, pm),
               d = mean(depth3d(w$x, w$y, w$z, pm)))
      }
    }
    dq <- vapply(quads, `[[`, numeric(1), "d")
    aq <- 0.05 + 0.13 * (dq - min(dq)) / diff(range(dq))
    for (j in order(dq)) {
      polygon(quads[[j]]$p$x, quads[[j]]$p$y, border = NA,
              col = adjustcolor("#c0392b", aq[j]))
    }
  }

  ## fitted mean-direction curve, painted in depth order, far half faint
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
draw_torus(dat, dir_hat, phig, lo = dir_lo, hi = dir_hi)
invisible(dev.off())
cat("wrote man/figures/README-torus.png\n")
cat("fit: conv =", b$outer.info$conv, " edf =", round(sum(b$edf), 2),
    " se_dir range:", paste(round(range(se_dir), 3), collapse = ".."), "\n")
