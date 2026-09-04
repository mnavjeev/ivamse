## effective dimension of a lasso fit, defined through the fitted values and the
## residual, which are unique when the coefficient vector is not

#' How Many Instruments a Lasso Fit Uses
#'
#' @description Counts the instruments a lasso first stage is really using. This
#' is the complexity measure that enters the selection criterion, and it plays
#' the part that the number of regressors plays in the degrees of freedom of an
#' ordinary regression.
#'
#' When no instrument can be written as a combination of the others, it is the
#' number of nonzero first-stage coefficients. When some can, the coefficients
#' themselves are not unique and two solvers may disagree about which
#' instruments were used; the count returned here does not depend on which
#' solution was computed.
#'
#' \code{\link{ivamse}} reports it for every candidate in the \code{d} column,
#' so it is seldom called directly.
#'
#' @param z instrument set, with columns already scaled so that
#'   \code{colSums(z^2) / n} equals one. A fitted model has done this scaling
#'   internally.
#' @param residual the lasso residual: the endogenous regressor minus the fitted
#'   first stage.
#' @param lambda the penalty that produced the fit, on the scale described in
#'   \code{\link{bcch_lambda}}. Supplying \pkg{glmnet}'s penalty instead gives a
#'   wrong answer without an error, since the two differ by a factor
#'   \code{n / nrow(z)} whenever there are controls.
#' @param n effective sample size.
#' @param pi_hat optional first-stage coefficient vector. When supplied, every
#'   instrument with a nonzero coefficient is counted, which guards against a
#'   solver leaving a coefficient just below the threshold.
#' @param dictionary_rank optional rank of \code{z}. When it equals
#'   \code{ncol(z)} no instrument is a combination of the others, and the count
#'   is simply the number of nonzero coefficients.
#' @param tolerance relative tolerance used to decide which instruments the
#'   penalty is holding at its boundary; \code{0.001} by default, the value used
#'   in the replication code for Ma, Navjeevan and Salahub. Because instruments
#'   with a nonzero coefficient are counted regardless, a smaller value can only
#'   lower the count.
#'
#' @details Formally the count is the rank of the columns in the
#' \emph{equicorrelation set}
#' \deqn{E(x) = \{ j : |E_n[z_{ij} u_i]| = \lambda \},}{%
#'   E(x) = { j : |mean(z_ij u_i)| = lambda },}
#' that is, the instruments whose sample covariance with the lasso residual is
#' as large as the penalty permits: those the lasso used, together with any it
#' came as close to using as the penalty allows. Every lasso solution is
#' supported inside this set, which is why the count does not depend on the
#' solution computed.
#'
#' @return A single non-negative integer, at most the number of instruments. A
#'   value of zero means the lasso shrank the whole first stage away, in which
#'   case the candidate defines no IV estimate.
#'
#' @seealso \code{\link{feasible_criterion}}, which uses it, and
#'   \code{\link{ivamse}}.
#'
#' @examples
#' ## a fit that keeps a handful of the twenty available instruments
#' set.seed(1)
#' N <- 300
#' z <- matrix(rnorm(N * 20), N, 20)
#' z <- sweep(z, 2, sqrt(colMeans(z^2)), "/")
#' x <- z[, 1] + 0.6 * z[, 2] + rnorm(N)
#' fit <- glmnet::glmnet(z, x, alpha = 1, lambda = 0.1,
#'                       standardize = FALSE, intercept = FALSE)
#' pi_hat <- as.vector(fit$beta)
#' effective_dimension(z, x - predict(fit, newx = z), 0.1, N, pi_hat)
#' sum(pi_hat != 0)
#'
#' ## a fitted model reports it for every candidate
#' y <- 2 * x + rnorm(N)
#' m <- ivamse_fit(y, x, z)
#' m$candidates[, c("kappa", "d")]
#'
#' @export
effective_dimension <- function(z, residual, lambda, n, pi_hat = NULL,
                                dictionary_rank = NULL,
                                tolerance = 0.001) {
  ## 0.001 is the value used in the replication code for Ma, Navjeevan and
  ## Salahub. The active set is unioned in below, so a column with a nonzero
  ## coefficient is never missed and the tolerance can only add columns; the
  ## shortcut and the general branch can therefore disagree by a column or two.
  active <- integer(0)
  if (!is.null(pi_hat)) {
    if (length(pi_hat) != ncol(z) || any(!is.finite(pi_hat))) {
      stop("the coefficient vector does not match the dictionary")
    }
    active <- which(pi_hat != 0)

    ## with full column rank the solution is unique and the equicorrelation set
    ## is generically the support
    if (!is.null(dictionary_rank) && dictionary_rank == ncol(z)) {
      return(length(active))
    }
  }

  score <- abs(crossprod(z, residual)) / n
  ## union with the active set: at small penalties floating-point error in the
  ## solver residual can push an active column below the equality threshold
  equicorrelation <- union(active, which(score >= lambda * (1 - tolerance)))

  if (!length(equicorrelation)) {
    0L
  } else {
    qr(z[, equicorrelation, drop = FALSE])$rank
  }
}
