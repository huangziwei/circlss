## von Mises location-scale general family for mgcv, structured on the
## gaulss() template (mgcv 1.9-4). The log-likelihood derivative algebra is
## transcribed from pycircstat2's vonmises_gen (dlogpdf/d2logpdf/d3logpdf/
## d4logpdf), which is FD-verified there; differential tests against
## pycircstat2 gate this port (see tests/testthat/test-vmlss-parity.R).

vmlss <- function(link = list("tanhalf", "log")) {
  if (length(link) != 2) stop("vmlss requires 2 links specified as character strings")
  okLinks <- list("tanhalf", "log")
  if (!(link[[1]] %in% okLinks[[1]]))
    stop(link[[1]], " link not available for the location parameter of vmlss")
  if (!(link[[2]] %in% okLinks[[2]]))
    stop(link[[2]], " link not available for the concentration parameter of vmlss")

  stats <- list()
  stats[[1]] <- tanhalf.link()
  stats[[2]] <- stats::make.link("log")
  fam <- structure(list(link = "log", canonical = "none",
                        linkfun = stats[[2]]$linkfun,
                        mu.eta = stats[[2]]$mu.eta),
                   class = "family")
  fam <- mgcv::fix.family.link(fam)
  stats[[2]]$d2link <- fam$d2link
  stats[[2]]$d3link <- fam$d3link
  stats[[2]]$d4link <- fam$d4link

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    mu <- object$fitted[, 1]
    kappa <- object$fitted[, 2]
    d <- object$y - mu
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## Var(sin(y - mu)) = (1 - A2(kappa))/2 = A1(kappa)/kappa
      v <- A1(kappa) / kappa
      v[kappa == 0] <- 0.5
      return(sind / sqrt(v))
    }
    ## deviance: per-obs saturated log lik (y = mu) minus attained, doubled:
    ## 2*kappa*(1 - cos(y - mu)), signed by the sine of the residual angle
    sign(sind) * sqrt(pmax(2 * kappa * (1 - cos(d)), 0))
  }

  postproc <- expression({
    ## Deviance and null deviance against a COMMON saturated reference, the
    ## per-observation saturated log-likelihood at the fitted kappa-hat:
    ##   deviance      = 2*(l_sat - l_hat)  = sum of squared dev. residuals
    ##   null.deviance = 2*(l_sat - l_null), null = intercept-only vM MLE
    ## so summary.gam's "deviance explained" is the fraction of the
    ## saturated-vs-null log-likelihood gap closed by the model, in [0, 1].
    ## Evaluated inside mgcv's frames: all work happens in local() so no
    ## host variable is clobbered.
    .vmlss.dev <- local({
      muh <- object$fitted[, 1]
      kah <- object$fitted[, 2]
      sy <- mean(sin(object$y))
      cy <- mean(cos(object$y))
      mu0 <- atan2(sy, cy)
      Rbar <- min(sqrt(sy * sy + cy * cy), 1 - 1e-08)
      kappa0 <- if (Rbar < 0.53) 2 * Rbar + Rbar^3 + 5 * Rbar^5 / 6 else
        if (Rbar < 0.85) -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar) else
          1 / (Rbar^3 - 4 * Rbar^2 + 3 * Rbar)
      dev <- sum(2 * kah * (1 - cos(object$y - muh)))
      lsat <- sum(-log(2 * pi * besselI(kah, 0, expon.scaled = TRUE)))
      lnull <- sum(kappa0 * (cos(object$y - mu0) - 1) -
                     log(2 * pi * besselI(kappa0, 0, expon.scaled = TRUE)))
      c(dev, 2 * (lsat - lnull))
    })
    object$deviance <- .vmlss.dev[1]
    object$null.deviance <- .vmlss.dev[2]
    rm(.vmlss.dev)
  })

  ll <- function(y, X, coef, wt, family, offset = NULL, deriv = 0,
                 d1b = 0, d2b = 0, Hp = NULL, rank = 0, fh = NULL,
                 D = NULL, eta = NULL, ncv = FALSE, sandwich = FALSE) {
    ## deriv: 0 - eval; 1 - grad and Hess; 2 - diagonal of first deriv of
    ## Hess; 3 - first deriv of Hess; 4 - everything (gaulss convention).
    if (!is.null(offset)) offset[[3]] <- 0
    jj <- attr(X, "lpi")
    if (is.null(eta)) {
      eta <- X[, jj[[1]], drop = FALSE] %*% coef[jj[[1]]]
      if (!is.null(offset[[1]])) eta <- eta + offset[[1]]
      eta1 <- X[, jj[[2]], drop = FALSE] %*% coef[jj[[2]]]
      if (!is.null(offset[[2]])) eta1 <- eta1 + offset[[2]]
    } else {
      eta1 <- eta[, 2]
      eta <- eta[, 1]
    }
    eta <- drop(eta); eta1 <- drop(eta1)
    mu <- family$linfo[[1]]$linkinv(eta)
    kappa <- family$linfo[[2]]$linkinv(eta1)
    n <- length(y)
    d <- y - mu
    cosd <- cos(d)
    sind <- sin(d)
    ## scaled-Bessel log density: kappa*(cos d - 1) - log(2*pi*i0e(kappa)),
    ## finite at every kappa (the naive kappa*cos d - log(2*pi*I0) pair
    ## overflows from kappa ~ 713)
    l0 <- kappa * (cosd - 1) - log(2 * pi * besselI(kappa, 0, expon.scaled = TRUE))
    l <- sum(l0)

    if (deriv) {
      a1 <- A1(kappa)
      ## l1: d l / d(mu, kappa); l2 columns ordered (mm, mk, kk)
      l1 <- cbind(kappa * sind, cosd - a1)
      l2 <- cbind(-kappa * cosd, sind, -A1prime(kappa, a1))
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),
                   family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(kappa))
    }
    l3 <- l4 <- g3 <- g4 <- 0
    if (deriv > 1) {
      ## l3 columns ordered (mmm, mmk, mkk, kkk); the mkk triple vanishes
      l3 <- cbind(-kappa * sind, -cosd, 0, -A1prime2(kappa))
      g3 <- cbind(family$linfo[[1]]$d3link(mu),
                  family$linfo[[2]]$d3link(kappa))
    }
    if (deriv > 3) {
      ## l4 columns ordered (mmmm, mmmk, mmkk, mkkk, kkkk); mixed quadruples
      ## with two or more kappa derivatives vanish
      l4 <- cbind(kappa * cosd, -sind, 0, 0, -A1prime3(kappa))
      g4 <- cbind(family$linfo[[1]]$d4link(mu),
                  family$linfo[[2]]$d4link(kappa))
    }
    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4
      de <- mgcv::gamlss.etamu(l1, l2, l3, l4, ig1, g2, g3, g4,
                               i2, i3, i4, deriv - 1)
      ret <- mgcv::gamlss.gH(X, jj, de$l1, de$l2, i2, l3 = de$l3, i3 = i3,
                             l4 = de$l4, i4 = i4, d1b = d1b, d2b = d2b,
                             deriv = deriv - 1, fh = fh, D = D,
                             sandwich = sandwich)
      if (ncv) {
        ret$l1 <- de$l1
        ret$l2 <- de$l2
        ret$l3 <- de$l3
      }
    } else ret <- list()
    ret$l <- l
    ret$l0 <- l0
    ret
  }

  sandwich <- function(y, X, coef, wt, family, offset = NULL) {
    ll(y, X, coef, wt, family, offset = NULL, deriv = 1, sandwich = TRUE)$lbb
  }

  initialize <- expression({
    ## Start from the global (intercept-only) von Mises MLE: mean direction
    ## mu0 and Fisher's A1-inverse approximation for kappa0, pushed through
    ## each link as a constant target and fitted to each linear predictor's
    ## design block by penalized least squares. Evaluated inside mgcv's
    ## frames (gam.fit5, initial.spg): everything except the contract
    ## variables n, use.unscaled and start stays inside local(), so host
    ## locals (S, off, rank, ...) are never clobbered.
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      start <- local({
        jj <- attr(x, "lpi")
        sy <- mean(sin(y))
        cy <- mean(cos(y))
        mu0 <- atan2(sy, cy)
        Rbar <- min(sqrt(sy * sy + cy * cy), 1 - 1e-08)
        kappa0 <- if (Rbar < 0.53) 2 * Rbar + Rbar^3 + 5 * Rbar^5 / 6 else
          if (Rbar < 0.85) -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar) else
            1 / (Rbar^3 - 4 * Rbar^2 + 3 * Rbar)
        kappa0 <- min(max(kappa0, 0.01), 500)
        st <- rep(0, ncol(x))
        yt1 <- rep(family$linfo[[1]]$linkfun(mu0), nobs)
        if (!is.null(offset) && length(offset) >= 1 && !is.null(offset[[1]]))
          yt1 <- yt1 - offset[[1]]
        x1 <- x[, jj[[1]], drop = FALSE]
        e1 <- E[, jj[[1]], drop = FALSE]
        stj <- qr.coef(qr(rbind(x1, e1)), c(yt1, rep(0, nrow(e1))))
        stj[!is.finite(stj)] <- 0
        st[jj[[1]]] <- stj
        yt2 <- rep(family$linfo[[2]]$linkfun(kappa0), nobs)
        if (!is.null(offset) && length(offset) >= 2 && !is.null(offset[[2]]))
          yt2 <- yt2 - offset[[2]]
        x1 <- x[, jj[[2]], drop = FALSE]
        e1 <- E[, jj[[2]], drop = FALSE]
        stj <- qr.coef(qr(rbind(x1, e1)), c(yt2, rep(0, nrow(e1))))
        stj[!is.finite(stj)] <- 0
        st[jj[[2]]] <- stj
        st
      })
    }
  })

  rd <- function(mu, wt, scale) {
    ## von Mises deviates by Best & Fisher (1979) rejection sampling;
    ## kappa ~ 0 falls back to the circular uniform. Angles return wrapped
    ## to (-pi, pi], the tan-half branch.
    mu0 <- mu[, 1]
    kappa <- mu[, 2]
    n <- length(mu0)
    th <- numeric(n)
    small <- kappa < 1e-09
    th[small] <- stats::runif(sum(small), -pi, pi)
    todo <- which(!small)
    if (length(todo)) {
      k <- kappa[todo]
      a <- 1 + sqrt(1 + 4 * k * k)
      b <- (a - sqrt(2 * a)) / (2 * k)
      r <- (1 + b * b) / (2 * b)
      out <- numeric(length(todo))
      left <- seq_along(todo)
      while (length(left)) {
        m <- length(left)
        z <- cos(pi * stats::runif(m))
        f <- (1 + r[left] * z) / (r[left] + z)
        cc <- k[left] * (r[left] - f)
        u2 <- stats::runif(m)
        ok <- (cc * (2 - cc) - u2 > 0) | (log(cc / u2) + 1 - cc >= 0)
        if (any(ok)) {
          u3 <- stats::runif(sum(ok))
          out[left[ok]] <- sign(u3 - 0.5) * acos(pmax(pmin(f[ok], 1), -1))
        }
        left <- left[!ok]
      }
      th[todo] <- out
    }
    th <- th + mu0
    atan2(sin(th), cos(th))
  }

  structure(list(family = "vmlss", ll = ll, link = paste(link), nlp = 2,
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 2, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
