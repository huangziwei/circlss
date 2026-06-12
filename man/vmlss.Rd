\name{vmlss}
\alias{vmlss}
\title{von Mises location-scale family for mgcv GAMs}
\usage{
vmlss(link = list("tanhalf", "log"))
}
\arguments{
\item{link}{Two-element list of link names, for the mean direction and the
concentration. Currently only the defaults are available:
\code{"tanhalf"} for the location and \code{"log"} for the concentration.}
}
\value{
An object of class \code{c("general.family", "extended.family", "family")}
for use with \code{\link[mgcv]{gam}}.
}
\description{
A general family implementing distributional regression for a circular
(angular) response \eqn{y} in radians, with
\deqn{y \sim \mathrm{vM}(\mu, \kappa),}
where both the mean direction \eqn{\mu} and the (log) concentration
\eqn{\kappa} get their own linear predictor, each of which may contain
smooth terms. Used with \code{mgcv::gam} and a list of two formulas: the
first specifies the response and the model for \eqn{\mu}, the second the
model for \eqn{\log\kappa}.
}
\details{
The mean direction uses the Fisher-Lee tan-half link,
\eqn{\mu = 2\arctan(\eta) \in (-\pi, \pi)}: the linear predictor's
intercept sets the reference direction and the antipode of that direction
is unrepresentable, the standard caveat of this link. The concentration
uses a log link.

Log-likelihood derivatives up to fourth order are implemented, so the
family supports full Newton smoothing-parameter estimation with
\code{method = "REML"} (and the other general-family criteria, including
\code{optimizer = "efs"}).

The derivative implementations are transcriptions of the finite-difference
verified implementations in the Python package \code{pycircstat2}, and the
fitted results are differentially tested against it.

The response should be supplied in radians; any branch (for example
\eqn{[0, 2\pi)} or \eqn{(-\pi, \pi]}) is acceptable since the density is
periodic. Fitted values are a two-column matrix: the mean direction in
\eqn{(-\pi, \pi)} and the concentration.
}
\examples{
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
mu <- 2 * atan(1.5 * sin(2 * pi * x))
kappa <- exp(1 + 0.8 * cos(2 * pi * x))
## von Mises deviates via wrapped rejection-free approximation for the
## example only: use circular::rvonmises or the family's rd in practice
y <- mu + rnorm(n) / sqrt(kappa)  # high-kappa approximation, example only
b <- gam(list(y ~ s(x), ~ s(x)), family = vmlss(), method = "REML")
summary(b)
plot(b, pages = 1)
}
\references{
Fisher, N. I. and Lee, A. J. (1992) Regression models for an angular
response. Biometrics 48, 665-677.

Wood, S. N., Pya, N. and Saefken, B. (2016) Smoothing parameter and model
selection for general smooth models. Journal of the American Statistical
Association 111, 1548-1575.
}
\seealso{
\code{\link[mgcv]{gaulss}}, \code{\link[mgcv]{gam}}
}
