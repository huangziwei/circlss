## Projected normal general family for mgcv, structured on the same
## template as vmlss (and mgcv's gaulss). The response y (radians) is the
## angle of a bivariate normal with mean (mu1, mu2) and identity
## covariance; BOTH linear predictors carry identity links, so the mean
## direction atan2(mu2, mu1) can sweep the whole circle and wind -- no
## tan-half branch cut. The derivative algebra is transcribed from
## pycircstat2's projectednormal_gen (dlogpdf/d2logpdf/d3logpdf/d4logpdf),
## FD-verified there; differential tests against pycircstat2 gate this
## port (tests/testthat/test-pnlss-parity.R).
##
## Rotated coordinates: t = mu1*cos(y) + mu2*sin(y) (radial),
##                      s = mu2*cos(y) - mu1*sin(y) (tangent), and with
## R(t) = dnorm(t)/pnorm(t) (inverse Mills ratio, computed stably via
## log-scale pnorm):
##   l  = -s^2/2 - log(2*pi)/2 + log(Phi(t)) + log(t + R)
##   G' = 1/(t + R),  P = R*G',  G'' = P - G'^2,
##   G''' = -t*P - 3*P*G' + 2*G'^3,
##   G'''' = P*(t^2 - 1) + 4*t*P*G' + 12*P*G'^2 - 3*P^2 - 6*G'^4,
## and the (mu1, mu2) derivatives are G^(k) times powers of cos(y)/sin(y),
## plus the tangent contributions at first/second order.

pn_mills_inv <- function(t) {
  ## R(t) = dnorm(t)/pnorm(t), stable at both extremes
  ## (t -> -Inf: R -> |t|; t -> +Inf: R -> 0)
  exp(-0.5 * t * t - 0.5 * log(2 * pi) - pnorm(t, log.p = TRUE))
}

pnlss <- function(link = list("identity", "identity")) {
  if (length(link) != 2) stop("pnlss requires 2 links specified as character strings")
  if (!(link[[1]] %in% "identity") || !(link[[2]] %in% "identity"))
    stop("pnlss supports only identity links for mu1 and mu2")

  stats <- list()
  for (i in 1:2) {
    stats[[i]] <- stats::make.link("identity")
    fam <- structure(list(link = "identity", canonical = "none",
                          linkfun = stats[[i]]$linkfun,
                          mu.eta = stats[[i]]$mu.eta),
                     class = "family")
    fam <- mgcv::fix.family.link(fam)
    stats[[i]]$d2link <- fam$d2link
    stats[[i]]$d3link <- fam$d3link
    stats[[i]]$d4link <- fam$d4link
  }

  residuals <- function(object, type = c("deviance", "pearson", "response")) {
    type <- match.arg(type)
    mu1 <- object$fitted[, 1]
    mu2 <- object$fitted[, 2]
    dirhat <- atan2(mu2, mu1)
    d <- object$y - dirhat
    if (type == "response") return(atan2(sin(d), cos(d)))
    ## deviance (pearson aliases it, as in pycircstat2): twice the log-lik
    ## gap to the fitted-direction mode, signed by the residual angle.
    ## At the mode t = gamma = ||mu||, s = 0.
    gam_ <- sqrt(mu1 * mu1 + mu2 * mu2)
    lmode <- -0.5 * log(2 * pi) + pnorm(gam_, log.p = TRUE) +
      log(gam_ + pn_mills_inv(gam_))
    t <- mu1 * cos(object$y) + mu2 * sin(object$y)
    s <- mu2 * cos(object$y) - mu1 * sin(object$y)
    l <- -0.5 * s * s - 0.5 * log(2 * pi) + pnorm(t, log.p = TRUE) +
      log(t + pn_mills_inv(t))
    sign(sin(d)) * sqrt(pmax(2 * (lmode - l), 0))
  }

  postproc <- expression({
    ## Deviance pair against the common saturated reference (per-obs
    ## log-lik at the fitted-direction mode), as in vmlss: "deviance
    ## explained" is the fraction of the saturated-vs-null gap closed.
    ## Null = moment fit (mean direction + A1inv concentration scale, the
    ## same start pycircstat2 uses). local() so mgcv's frames are safe.
    .pnlss.dev <- local({
      mu1 <- object$fitted[, 1]
      mu2 <- object$fitted[, 2]
      y <- object$y
      mills <- function(t) exp(-0.5 * t * t - 0.5 * log(2 * pi) -
                                 pnorm(t, log.p = TRUE))
      lpdf <- function(y, m1, m2) {
        t <- m1 * cos(y) + m2 * sin(y)
        s <- m2 * cos(y) - m1 * sin(y)
        -0.5 * s * s - 0.5 * log(2 * pi) + pnorm(t, log.p = TRUE) +
          log(t + mills(t))
      }
      gam_ <- sqrt(mu1 * mu1 + mu2 * mu2)
      lsat <- sum(-0.5 * log(2 * pi) + pnorm(gam_, log.p = TRUE) +
                    log(gam_ + mills(gam_)))
      lhat <- sum(lpdf(y, mu1, mu2))
      sy <- mean(sin(y)); cy <- mean(cos(y))
      dir0 <- atan2(sy, cy)
      Rbar <- min(sqrt(sy * sy + cy * cy), 1 - 1e-08)
      g0 <- if (Rbar < 0.53) 2 * Rbar + Rbar^3 + 5 * Rbar^5 / 6 else
        if (Rbar < 0.85) -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar) else
          1 / (Rbar^3 - 4 * Rbar^2 + 3 * Rbar)
      g0 <- min(max(g0, 1e-3), 500)
      lnull <- sum(lpdf(y, g0 * cos(dir0), g0 * sin(dir0)))
      c(2 * (lsat - lhat), 2 * (lsat - lnull))
    })
    object$deviance <- .pnlss.dev[1]
    object$null.deviance <- .pnlss.dev[2]
    rm(.pnlss.dev)
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
    mu1 <- family$linfo[[1]]$linkinv(eta)
    mu2 <- family$linfo[[2]]$linkinv(eta1)
    cy <- cos(y)
    sy <- sin(y)
    t <- mu1 * cy + mu2 * sy
    s <- mu2 * cy - mu1 * sy
    Rm <- pn_mills_inv(t)
    l0 <- -0.5 * s * s - 0.5 * log(2 * pi) + pnorm(t, log.p = TRUE) +
      log(t + Rm)
    l <- sum(l0)

    if (deriv) {
      g1 <- 1 / (t + Rm)
      P <- Rm * g1
      g2 <- P - g1 * g1
      ## l1: d l / d(mu1, mu2); l2 columns ordered (11, 12, 22)
      l1 <- cbind(g1 * cy + s * sy, g1 * sy - s * cy)
      l2 <- cbind(g2 * cy * cy - sy * sy,
                  (g2 + 1) * cy * sy,
                  g2 * sy * sy - cy * cy)
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),
                   family$linfo[[2]]$mu.eta(eta1))
      g2l <- cbind(family$linfo[[1]]$d2link(mu1),
                   family$linfo[[2]]$d2link(mu2))
    }
    l3 <- l4 <- g3l <- g4l <- 0
    if (deriv > 1) {
      ## only the radial direction survives at third order:
      ## l3 columns (111, 112, 122, 222) = G''' times cos/sin powers
      G3 <- -t * P - 3 * P * g1 + 2 * g1^3
      l3 <- cbind(G3 * cy^3, G3 * cy * cy * sy, G3 * cy * sy * sy,
                  G3 * sy^3)
      g3l <- cbind(family$linfo[[1]]$d3link(mu1),
                   family$linfo[[2]]$d3link(mu2))
    }
    if (deriv > 3) {
      ## l4 columns (1111, 1112, 1122, 1222, 2222)
      G4 <- P * (t * t - 1) + 4 * t * P * g1 + 12 * P * g1 * g1 -
        3 * P * P - 6 * g1^4
      l4 <- cbind(G4 * cy^4, G4 * cy^3 * sy, G4 * cy * cy * sy * sy,
                  G4 * cy * sy^3, G4 * sy^4)
      g4l <- cbind(family$linfo[[1]]$d4link(mu1),
                   family$linfo[[2]]$d4link(mu2))
    }
    if (deriv) {
      i2 <- family$tri$i2
      i3 <- family$tri$i3
      i4 <- family$tri$i4
      de <- mgcv::gamlss.etamu(l1, l2, l3, l4, ig1, g2l, g3l, g4l,
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
    ## Start from the moment fit pycircstat2 uses: mean direction d0 and a
    ## von Mises-scale concentration gamma0 = A1inv(Rbar), pushed through
    ## the identity links as constant targets per linear predictor
    ## (mu1 = gamma0*cos d0, mu2 = gamma0*sin d0), fitted by penalized LS.
    ## Everything except the contract variables stays inside local().
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      start <- local({
        jj <- attr(x, "lpi")
        sy <- mean(sin(y))
        cy <- mean(cos(y))
        dir0 <- atan2(sy, cy)
        Rbar <- min(sqrt(sy * sy + cy * cy), 1 - 1e-08)
        gamma0 <- if (Rbar < 0.53) 2 * Rbar + Rbar^3 + 5 * Rbar^5 / 6 else
          if (Rbar < 0.85) -0.4 + 1.39 * Rbar + 0.43 / (1 - Rbar) else
            1 / (Rbar^3 - 4 * Rbar^2 + 3 * Rbar)
        gamma0 <- min(max(gamma0, 1e-3), 500)
        targets <- c(gamma0 * cos(dir0), gamma0 * sin(dir0))
        st <- rep(0, ncol(x))
        for (j in 1:2) {
          ytj <- rep(targets[j], nobs)
          if (!is.null(offset) && length(offset) >= j &&
              !is.null(offset[[j]])) ytj <- ytj - offset[[j]]
          x1 <- x[, jj[[j]], drop = FALSE]
          e1 <- E[, jj[[j]], drop = FALSE]
          stj <- qr.coef(qr(rbind(x1, e1)), c(ytj, rep(0, nrow(e1))))
          stj[!is.finite(stj)] <- 0
          st[jj[[j]]] <- stj
        }
        st
      })
    }
  })

  rd <- function(mu, wt, scale) {
    ## exact: the angle of a bivariate normal draw around (mu1, mu2)
    n <- nrow(mu)
    atan2(mu[, 2] + stats::rnorm(n), mu[, 1] + stats::rnorm(n))
  }

  structure(list(family = "pnlss", ll = ll, link = paste(link), nlp = 2,
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 2, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
