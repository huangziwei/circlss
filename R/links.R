## Fisher-Lee tan-half link for a circular location parameter:
##   eta = g(mu) = tan(mu/2),  mu = g^{-1}(eta) = 2*atan(eta) in (-pi, pi).
## Transcribed from pycircstat2's TanHalfLink. Derivative conventions are
## mgcv's: mu.eta is d mu / d eta as a function of eta; d2link, d3link and
## d4link are the 2nd-4th derivatives of the link g with respect to mu, as
## functions of mu. In t = tan(mu/2) (so dt/dmu = (1+t^2)/2):
##   g'  = (1+t^2)/2,            g'' = t(1+t^2)/2,
##   g''' = (1+t^2)(1+3t^2)/4,   g'''' = t(1+t^2)(2+3t^2)/2.
## The antipode mu = pi (mod 2pi) is the one singularity: unrepresentable
## under this link, the standard Fisher-Lee caveat.
tanhalf.link <- function() {
  list(
    link = "tanhalf",
    linkfun = function(mu) tan(0.5 * mu),
    linkinv = function(eta) 2 * atan(eta),
    mu.eta = function(eta) 2 / (1 + eta * eta),
    d2link = function(mu) {
      t <- tan(0.5 * mu)
      0.5 * t * (1 + t * t)
    },
    d3link = function(mu) {
      t2 <- tan(0.5 * mu)^2
      0.25 * (1 + t2) * (1 + 3 * t2)
    },
    d4link = function(mu) {
      t <- tan(0.5 * mu)
      t2 <- t * t
      0.5 * t * (1 + t2) * (2 + 3 * t2)
    },
    valideta = function(eta) TRUE
  )
}
