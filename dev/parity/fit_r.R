# Fit the parity battery with circlss (the port under test).
# Mirrors fit_python.py case for case: same CSV in, same JSON shape out,
# Newton-REML both sides. Run from dev/parity/ (or pass the dir as arg 1).

suppressMessages({
  library(mgcv)
  library(circlss)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
here <- if (length(args) >= 1) args[1] else "."
out_dir <- file.path(here, "out")
dir.create(out_dir, showWarnings = FALSE)

cases <- list(
  lin = list(
    formula = list(y ~ x, ~ x),
    knots = NULL, var = "x", grid = seq(0, 1, length.out = 101)
  ),
  smooth = list(
    formula = list(y ~ s(x, k = 10), ~ s(x, k = 10)),
    knots = NULL, var = "x", grid = seq(0, 1, length.out = 101)
  ),
  cyclic = list(
    formula = list(y ~ s(phi, bs = "cc", k = 10), ~ s(phi, bs = "cc", k = 10)),
    knots = list(phi = c(-pi, pi)), var = "phi",
    grid = seq(-pi, pi, length.out = 101)
  ),
  small = list(
    formula = list(y ~ s(x, k = 8), ~ s(x, k = 8)),
    knots = NULL, var = "x", grid = seq(0, 1, length.out = 101)
  )
)

for (nm in names(cases)) {
  cs <- cases[[nm]]
  dat <- read.csv(file.path(here, "data", paste0(nm, ".csv")))
  b <- gam(cs$formula, family = vmlss(), data = dat, method = "REML",
           knots = cs$knots)
  conv <- is.null(b$outer.info$conv) ||
    identical(b$outer.info$conv, "full convergence")
  Xp <- predict(b, type = "lpmatrix")
  lpi <- attr(Xp, "lpi")
  nd <- stats::setNames(data.frame(cs$grid), cs$var)
  pr <- predict(b, newdata = nd, type = "response")
  edf_by <- lapply(b$smooth, function(s) sum(b$edf[s$first.para:s$last.para]))
  names(edf_by) <- vapply(b$smooth, function(s) s$label, character(1))

  out <- list(
    case = nm,
    coef_names = names(coef(b)),
    coef = unname(coef(b)),
    lpi = lapply(lpi, function(i) as.integer(i) - 1L),  # 0-based like Python
    sp = unname(b$sp),
    edf_total = sum(b$edf),
    edf_by_smooth = edf_by,
    loglik = as.numeric(logLik(b)),
    reml = as.numeric(b$gcv.ubre),
    converged = conv,
    grid = cs$grid,
    mu_grid = unname(pr[, 1]),
    kappa_grid = unname(pr[, 2])
  )
  path <- file.path(out_dir, paste0(nm, "_r.json"))
  write_json(out, path, digits = NA, auto_unbox = TRUE)
  cat(sprintf("%s: converged=%s loglik=%.6f edf=%.4f sp=%s\n",
              nm, conv, out$loglik, out$edf_total,
              paste(signif(out$sp, 5), collapse = ",")))
}
