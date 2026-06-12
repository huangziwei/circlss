# Wrapped Cauchy location-scale family for mgcv GAMs

A general family implementing distributional regression for a circular
(angular) response \\y\\ in radians under the wrapped Cauchy law
\$\$f(y) = \frac{1 - \rho^2}{2\pi\\(1 + \rho^2 - 2\rho\cos(y -
\mu))},\$\$ with both the mean direction \\\mu\\ and the mean resultant
length \\\rho\\ getting their own linear predictor, each of which may
contain smooth terms. Used with
[`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) and a list of two
formulas: the first names the response and models \\\mu\\, the second
models \\\mathrm{logit}(\rho)\\.

## Usage

``` r
wclss(link = list("tanhalf", "logit"))
```

## Arguments

- link:

  Two-element list of link names, for the mean direction and the mean
  resultant length. Currently only the defaults are available:
  `"tanhalf"` for the location and `"logit"` for the concentration
  parameter \\\rho \in (0, 1)\\.

## Value

An object of class `c("general.family", "extended.family", "family")`
for use with [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html).

## Details

The wrapped Cauchy is the heavy-tailed counterpart of the von Mises
([`vmlss`](https://huangziwei.github.io/circlss/reference/vmlss.md)):
sharply peaked with fat circular tails, so it is the more robust choice
when the data contain angular outliers. Its trigonometric moments are
simply \\\rho^p\\, which gives clean residual conventions: Pearson
residuals standardize by \\\mathrm{Var}\\\sin(y-\mu)\\ = (1-\rho^2)/2\\.

The mean direction uses the Fisher-Lee tan-half link (\\\mu \in (-\pi,
\pi)\\, antipode unrepresentable, winding number zero – see
[`pnlss`](https://huangziwei.github.io/circlss/reference/pnlss.md) when
the mean direction must wind). Log-likelihood derivatives are
implemented to fourth order, so full Newton `method = "REML"` works (as
does `optimizer = "efs"`). Internally the density denominator is
computed in the cancellation-free form \\(1-\rho)^2 +
4\rho\sin^2((y-\mu)/2)\\ so the log-likelihood stays exact as \\\rho \to
1\\.

The derivative implementations are transcriptions of the
finite-difference-verified implementations in the Python package
`pycircstat2`, and fitted results are differentially tested against it.

## References

Fisher, N. I. and Lee, A. J. (1992) Regression models for an angular
response. Biometrics 48, 665-677.

Wood, S. N., Pya, N. and Saefken, B. (2016) Smoothing parameter and
model selection for general smooth models. Journal of the American
Statistical Association 111, 1548-1575.

## See also

[`vmlss`](https://huangziwei.github.io/circlss/reference/vmlss.md),
[`pnlss`](https://huangziwei.github.io/circlss/reference/pnlss.md),
[`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Examples

``` r
library(mgcv)
set.seed(1)
n <- 300
x <- runif(n)
mu <- 2 * atan(1.4 * sin(2 * pi * x))
rho <- plogis(0.5 + 0.8 * cos(2 * pi * x))
y <- atan2(sin(mu - log(rho) * rcauchy(n)), cos(mu - log(rho) * rcauchy(n)))
b <- gam(list(y ~ s(x), ~ s(x)), family = wclss(), method = "REML")
summary(b)
#> 
#> Family: wclss 
#> Link function: tanhalf logit 
#> 
#> Formula:
#> y ~ s(x)
#> <environment: 0x55bbe9babe40>
#> ~s(x)
#> <environment: 0x55bbe9babe40>
#> 
#> Parametric coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)    0.05089    0.04567   1.114    0.265    
#> (Intercept).1  0.49627    0.10869   4.566 4.98e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Approximate significance of smooth terms:
#>          edf Ref.df Chi.sq p-value    
#> s(x)   6.764  7.789 590.40 < 2e-16 ***
#> s.1(x) 2.837  3.495  20.92 0.00025 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Deviance explained = 35.4%
#> -REML = 417.49  Scale est. = 1         n = 300
```
