## Wrapped Cauchy general family for mgcv, on the vmlss/gaulss template:
## y ~ WC(mu, rho), location mu on the Fisher-Lee tan-half link, mean
## resultant length rho in (0, 1) on a logit link, one linear predictor
## (and smooths) each. The derivative algebra is transcribed from
## pycircstat2's wrapcauchy_gen (dlogpdf/d2logpdf/d3logpdf/d4logpdf),
## FD-verified there; differential tests against pycircstat2 gate this
## port (tests/testthat/test-wclss-parity.R).
##
## Numerical backbone: the density denominator is computed as
##   D = (1 - rho)^2 + 4*rho*sin(d/2)^2   (== 1 + rho^2 - 2*rho*cos d)
## which does not cancel at the peak as rho -> 1, and the normalizer uses
## log1p so the log-density stays exact near rho = 1:
##   l = log1p(-rho) + log1p(rho) - log(2*pi) - log(D).
## Derivatives combine D's partials (Dm = -2*rho*s, Dr = 2*rho - 2*c,
## Dmm = 2*rho*c, Dmr = -2*s, Drr = 2, Dmmm = 2*rho*s, Dmmr = 2*c,
## Dmmmm = -2*rho*c, Dmmmr = 2*s) through the log-derivative partition
## formulas, plus the rho-only chain from log(1 - rho^2).

wclss <- function(link = list("tanhalf", "logit")) {
  if (length(link) != 2) stop("wclss requires 2 links specified as character strings")
  if (!(link[[1]] %in% "tanhalf"))
    stop(link[[1]], " link not available for the location parameter of wclss")
  if (!(link[[2]] %in% "logit"))
    stop(link[[2]], " link not available for the concentration parameter of wclss")

  stats <- list()
  stats[[1]] <- tanhalf.link()
  stats[[2]] <- stats::make.link("logit")
  fam <- structure(list(link = "logit", canonical = "none",
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
    rho <- object$fitted[, 2]
    d <- object$y - mu
    if (type == "response") return(atan2(sin(d), cos(d)))
    sind <- sin(d)
    if (type == "pearson") {
      ## WC trig moments are rho^p, so Var(sin(y - mu)) = (1 - rho^2)/2
      return(sind / sqrt((1 - rho * rho) / 2))
    }
    ## deviance: 2*(l_mode - l) = 2*log(D / (1 - rho)^2), signed
    D <- (1 - rho)^2 + 4 * rho * sin(0.5 * d)^2
    sign(sind) * sqrt(pmax(2 * (log(D) - 2 * log1p(-rho)), 0))
  }

  postproc <- expression({
    ## Deviance pair against the common saturated reference (per-obs
    ## log-lik at the fitted mode), as in vmlss; null = moment fit
    ## (mean direction, rho = Rbar -- the WC moment estimator). local()
    ## so mgcv's frames are safe.
    .wclss.dev <- local({
      mu <- object$fitted[, 1]
      rho <- object$fitted[, 2]
      y <- object$y
      D <- (1 - rho)^2 + 4 * rho * sin(0.5 * (y - mu))^2
      dev <- sum(2 * (log(D) - 2 * log1p(-rho)))
      lsat <- sum(log1p(rho) - log1p(-rho) - log(2 * pi))
      sy <- mean(sin(y)); cy <- mean(cos(y))
      mu0 <- atan2(sy, cy)
      rho0 <- min(max(sqrt(sy * sy + cy * cy), 1e-3), 1 - 1e-3)
      D0 <- (1 - rho0)^2 + 4 * rho0 * sin(0.5 * (y - mu0))^2
      lnull <- sum(log1p(-rho0) + log1p(rho0) - log(2 * pi) - log(D0))
      c(dev, 2 * (lsat - lnull))
    })
    object$deviance <- .wclss.dev[1]
    object$null.deviance <- .wclss.dev[2]
    rm(.wclss.dev)
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
    rho <- family$linfo[[2]]$linkinv(eta1)
    d <- y - mu
    s <- sin(d)
    cd <- cos(d)
    Dq <- (1 - rho)^2 + 4 * rho * sin(0.5 * d)^2
    l0 <- log1p(-rho) + log1p(rho) - log(2 * pi) - log(Dq)
    l <- sum(l0)

    if (deriv) {
      w <- 1 / Dq
      w2 <- w * w
      om <- 1 - rho * rho            # 1 - rho^2
      Dm <- -2 * rho * s
      Dr <- 2 * rho - 2 * cd
      ## l1: d l / d(mu, rho); l2 columns ordered (mm, mr, rr)
      l1 <- cbind(2 * rho * s * w,
                  -2 * rho / om - Dr * w)
      l2 <- cbind(-2 * rho * cd * w + Dm * Dm * w2,
                  2 * s * w + Dm * Dr * w2,
                  -2 * (1 + rho * rho) / om^2 - 2 * w + Dr * Dr * w2)
      ig1 <- cbind(family$linfo[[1]]$mu.eta(eta),
                   family$linfo[[2]]$mu.eta(eta1))
      g2 <- cbind(family$linfo[[1]]$d2link(mu),
                  family$linfo[[2]]$d2link(rho))
    }
    l3 <- l4 <- g3 <- g4 <- 0
    if (deriv > 1) {
      ## l3 columns (mmm, mmr, mrr, rrr)
      w3 <- w2 * w
      Dmm <- 2 * rho * cd
      Dmr <- -2 * s
      Drr <- 2
      Dmmm <- 2 * rho * s
      Dmmr <- 2 * cd
      l3 <- cbind(
        -Dmmm * w + 3 * Dmm * Dm * w2 - 2 * Dm^3 * w3,
        -Dmmr * w + (Dmm * Dr + 2 * Dmr * Dm) * w2 - 2 * Dm * Dm * Dr * w3,
        (2 * Dmr * Dr + Drr * Dm) * w2 - 2 * Dm * Dr * Dr * w3,
        -4 * rho * (3 + rho * rho) / om^3 + 3 * Drr * Dr * w2 -
          2 * Dr^3 * w3
      )
      g3 <- cbind(family$linfo[[1]]$d3link(mu),
                  family$linfo[[2]]$d3link(rho))
    }
    if (deriv > 3) {
      ## l4 columns (mmmm, mmmr, mmrr, mrrr, rrrr)
      w4 <- w3 * w
      Dmmmm <- -2 * rho * cd
      Dmmmr <- 2 * s
      l4 <- cbind(
        -Dmmmm * w + (4 * Dmmm * Dm + 3 * Dmm * Dmm) * w2 -
          12 * Dmm * Dm * Dm * w3 + 6 * Dm^4 * w4,
        -Dmmmr * w + (Dmmm * Dr + 3 * Dmmr * Dm + 3 * Dmm * Dmr) * w2 -
          6 * (Dmm * Dm * Dr + Dmr * Dm * Dm) * w3 + 6 * Dm^3 * Dr * w4,
        (2 * Dmmr * Dr + Dmm * Drr + 2 * Dmr * Dmr) * w2 -
          2 * (Dmm * Dr * Dr + Drr * Dm * Dm + 4 * Dmr * Dm * Dr) * w3 +
          6 * Dm * Dm * Dr * Dr * w4,
        3 * Dmr * Drr * w2 - 6 * (Drr * Dm * Dr + Dmr * Dr * Dr) * w3 +
          6 * Dm * Dr^3 * w4,
        -12 * (1 + 6 * rho * rho + rho^4) / om^4 + 3 * Drr * Drr * w2 -
          12 * Drr * Dr * Dr * w3 + 6 * Dr^4 * w4
      )
      g4 <- cbind(family$linfo[[1]]$d4link(mu),
                  family$linfo[[2]]$d4link(rho))
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
    ## Start from the moment fit: mean direction mu0 through the tan-half
    ## link, and rho0 = Rbar (the WC moment estimator, E cos(y - mu) =
    ## rho) through the logit link; constant targets per linear predictor
    ## fitted by penalized LS. Everything except the contract variables
    ## stays inside local().
    n <- rep(1, nobs)
    use.unscaled <- if (!is.null(attr(E, "use.unscaled"))) TRUE else FALSE
    if (is.null(start)) {
      start <- local({
        jj <- attr(x, "lpi")
        sy <- mean(sin(y))
        cy <- mean(cos(y))
        mu0 <- atan2(sy, cy)
        rho0 <- min(max(sqrt(sy * sy + cy * cy), 0.01), 0.95)
        targets <- c(family$linfo[[1]]$linkfun(mu0),
                     family$linfo[[2]]$linkfun(rho0))
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
    ## exact: the wrapped Cauchy is the wrapping of a Cauchy with scale
    ## gamma = -log(rho)
    n <- nrow(mu)
    gam_ <- -log(pmin(pmax(mu[, 2], 1e-12), 1 - 1e-12))
    th <- mu[, 1] + gam_ * stats::rcauchy(n)
    atan2(sin(th), cos(th))
  }

  structure(list(family = "wclss", ll = ll, link = paste(link), nlp = 2,
                 sandwich = sandwich, tri = mgcv::trind.generator(2),
                 initialize = initialize, postproc = postproc,
                 residuals = residuals, linfo = stats, rd = rd,
                 d2link = 1, d3link = 1, d4link = 1, ls = 1,
                 available.derivs = 2, discrete.ok = FALSE),
            class = c("general.family", "extended.family", "family"))
}
