## partialling out the included exogenous controls, which reduces the model to
## the one without them at effective sample size n = N - q, q = rank(W). Every
## average E_n[.] therefore divides by n, not by the number of rows N. Residuals
## are kept as N-vectors, since a robust or clustered variance is not invariant
## to a rotation onto an orthonormal basis of the complement.


## QR decomposition of the controls, or NULL if there are none
control_qr <- function(w) {
  if (is.null(w) || ncol(w) == 0L) return(NULL)
  qr(w)
}


## residual from projecting on the controls, M_W a
partial_out <- function(a, w_qr) {
  if (is.null(w_qr)) return(a)
  qr.resid(w_qr, a)
}


## Residualize a dictionary on the controls, drop columns the projection
## annihilates, and scale the rest to unit empirical second moment,
## sum(z_j^2) / n == 1. This is Assumption 2(i) of Ma, Navjeevan and Salahub;
## the penalty is assigned after the scaling. A column in the span of the
## controls contributes nothing to the fit and cannot be scaled.
prepare_dictionary <- function(z, w_qr, n) {
  raw_scale <- sqrt(colSums(z^2) / n)
  z <- partial_out(z, w_qr)
  scale <- sqrt(colSums(z^2) / n)

  ## test relative to the column's own size, so it does not depend on the units
  keep <- is.finite(scale) & raw_scale > 0 & scale > 1e-10 * raw_scale
  if (!any(keep)) {
    stop("every dictionary column lies in the span of the controls")
  }

  z <- z[, keep, drop = FALSE]
  scale <- scale[keep]
  list(z = sweep(z, 2L, scale, "/"), scale = scale, kept = which(keep))
}
