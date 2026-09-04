## fitting one dictionary's candidate first stages. The problem is
##
##   minimize  (1/(2n)) * ||x - Z pi||^2  +  lambda * ||pi||_1
##
## glmnet's Gaussian objective normalizes the loss by the number of rows N, not
## by n = N - q, so the two agree at glmnet_lambda = lambda * n / N

GLMNET_THRESHOLD <- 1e-13
GLMNET_MAXIT <- 500000L


## match requested penalties back to the columns glmnet returned; a duplicate
## index means glmnet dropped a path point and two penalties map to one column
match_lambda <- function(target, returned, tolerance = 1e-10) {
  index <- vapply(target, function(v) which.min(abs(returned - v)), integer(1))
  ok <- abs(returned[index] - target) <= tolerance * pmax(1, abs(target))
  if (any(!ok) || anyDuplicated(index)) {
    stop("glmnet did not return the requested penalty grid")
  }
  index
}


## fit every candidate penalty for one dictionary and return the quantities the
## criterion is built from. The outcome y does not enter; only (x, Z).
fit_dictionary <- function(z, x, lambda, n, dictionary_rank) {
  N <- nrow(z)
  glmnet_lambda <- lambda * n / N
  path <- sort(glmnet_lambda, decreasing = TRUE)

  fit <- glmnet::glmnet(
    z, x,
    alpha = 1, lambda = path,
    standardize = FALSE, intercept = FALSE,
    thresh = GLMNET_THRESHOLD, maxit = GLMNET_MAXIT
  )

  pi_hat <- as.matrix(fit$beta)[, match_lambda(glmnet_lambda, fit$lambda), drop = FALSE]
  Pi_hat <- z %*% pi_hat

  J <- length(lambda)
  out <- list(Pi_hat = Pi_hat, pi_hat = pi_hat, h = numeric(J),
              fitted_moment = numeric(J), d = numeric(J),
              nonzero = integer(J), pi_l1 = numeric(J))

  for (j in seq_len(J)) {
    fitted <- Pi_hat[, j]
    coefficients <- pi_hat[, j]

    out$h[j] <- sum(fitted * x) / n
    out$fitted_moment[j] <- sum(fitted^2) / n
    out$nonzero[j] <- sum(coefficients != 0)
    out$pi_l1[j] <- sum(abs(coefficients))
    ## the penalty passed here is on the objective's scale, not glmnet's:
    ## glmnet's |z'r| / N <= glmnet_lambda is exactly |z'r| / n <= lambda
    out$d[j] <- effective_dimension(
      z, x - fitted, lambda[j], n,
      pi_hat = coefficients, dictionary_rank = dictionary_rank
    )
  }

  out
}


## tolerance below which a fitted second moment is numerical dust rather than a
## fit, scaled to the size of x
EMPTY_FIT_TOLERANCE <- 100 * .Machine$double.eps

## empty-fit check: the fitted values, the second moment and the effective
## dimension must all be non-null, since the solver can leave numerical dust at
## the largest penalties. A merely nonzero h is not part of the test.
is_empty_fit <- function(Pi_hat, fitted_moment, d, n, x) {
  tolerance <- EMPTY_FIT_TOLERANCE * max(1, sum(x^2) / n)

  fit_has_a_nonzero <- all(is.finite(Pi_hat)) && any(Pi_hat != 0)
  moment_above_dust <- is.finite(fitted_moment) && fitted_moment > tolerance
  dimension_positive <- is.finite(d) && d > 0

  !(fit_has_a_nonzero && moment_above_dust && dimension_positive)
}
