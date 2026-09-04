## Penalty levels and candidate grids. Two rules are available. "path" is a
## geometric grid running down from the smallest penalty that shrinks every
## coefficient to zero, as in glmnet; it needs no preliminary estimate of
## sigma_v. "bcch" is fixed multiples of the homoscedastic plug-in penalty of
## Belloni, Chen, Chernozhukov and Hansen (2012), reproducing the candidate list
## of Ma, Navjeevan and Salahub. Grids are built separately for each dictionary.


#' Plug-In Lasso Penalty of Belloni, Chen, Chernozhukov and Hansen
#'
#' @description Returns the penalty level that the plug-in rule of Belloni,
#' Chen, Chernozhukov and Hansen (2012) prescribes for a lasso first stage. It
#' is computed from the number of instruments, the sample size and the scale of
#' the first-stage error, rather than by cross-validation, and is set just high
#' enough that the lasso discards instruments whose apparent correlation with
#' the endogenous regressor is only sampling noise.
#'
#' Most users do not need to call this. \code{\link{ivamse}} uses it internally
#' when \code{penalty = "bcch"} and reports every candidate's penalty as a
#' multiple of it in the \code{kappa} column. Call it directly only to build a
#' penalty grid by hand.
#'
#' @param p number of columns in the instrument set.
#' @param n effective sample size. When controls are partialled out this is the
#'   number of observations less the rank of the controls, not
#'   \code{nrow(data)}; a fitted model reports it as \code{fit$n}.
#' @param sigma_v standard deviation of the first-stage error. The penalty is
#'   proportional to it, so the default of \code{1} returns the penalty per unit
#'   of that standard deviation and will not match the penalties a fitted model
#'   used. A fitted model reports the value it used as \code{fit$sigma_v}.
#' @param c_lambda inflation constant, greater than one. Belloni et al. and Ma
#'   et al. both use \code{1.1}, the default.
#'
#' @details The penalty is
#' \deqn{\lambda = c_\lambda \sigma_v \Phi^{-1}(1 - \gamma / (2p)) / \sqrt{n},
#'       \qquad \gamma = 0.1 / \log(\max(p, n)),}{%
#'   lambda = c_lambda * sigma_v * qnorm(1 - gamma / (2p)) / sqrt(n),
#'   gamma = 0.1 / log(max(p, n)),}
#' with \eqn{\Phi}{Phi} the standard normal distribution function. It grows with
#' the number of instruments, since more instruments give noise more chances to
#' look like signal, and shrinks with the sample size.
#'
#' @return A single number: a penalty on the scale of the objective
#'   \eqn{(1/2n)\|x - Z\pi\|^2 + \lambda\|\pi\|_1}{(1/2n)||x - Z pi||^2 + lambda
#'   ||pi||_1}, which is the scale of the \code{lambda} argument of
#'   \code{\link{ivamse_fit}}. It is not \pkg{glmnet}'s scale, which normalizes
#'   by the number of rows rather than by the effective sample size;
#'   \code{\link{ivamse_fit}} converts internally.
#'
#' @references
#' Belloni, A., Chen, D., Chernozhukov, V., and Hansen, C. (2012). Sparse models
#' and methods for optimal instruments with an application to eminent domain.
#' \emph{Econometrica} \bold{80}, 2369--2429.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{implied_alpha}}.
#'
#' @examples
#' ## the plug-in penalty for 50 instruments and 500 observations
#' bcch_lambda(p = 50, n = 500, sigma_v = 1)
#'
#' ## a grid of multiples of it, which is what penalty = "bcch" fits
#' 2^seq(-4, 2, by = 0.5) * bcch_lambda(p = 50, n = 500, sigma_v = 1)
#'
#' ## recovering the level a fitted model used: kappa = 1 by construction
#' set.seed(1)
#' N <- 300
#' z <- matrix(rnorm(N * 20), N, 20)
#' v <- rnorm(N)
#' x <- z[, 1] + v
#' y <- 2 * x + 0.7 * v + rnorm(N)
#' m <- ivamse_fit(y, x, z, penalty = "bcch", kappa = c(0.5, 1, 2))
#' m$candidates[, c("lambda", "kappa")]
#' bcch_lambda(p = 20, n = m$n, sigma_v = m$sigma_v)
#'
#' @export
bcch_lambda <- function(p, n, sigma_v = 1, c_lambda = 1.1) {
  stopifnot(p >= 1, n >= 1, sigma_v > 0, c_lambda > 1)
  gamma <- 0.1 / log(max(p, n))
  c_lambda * sigma_v * stats::qnorm(1 - gamma / (2 * p)) / sqrt(n)
}


#' Whether a Penalty Satisfies the Condition the Theory Requires
#'
#' @description Reports whether a penalty is large enough for the guarantee
#' behind the selection criterion to apply. The theory of Ma, Navjeevan and
#' Salahub requires the penalty to exceed a threshold that grows with the number
#' of instruments; \code{implied_alpha} returns a number that is positive when
#' the penalty clears that threshold and negative when it does not. Its
#' magnitude has no interpretation, only its sign.
#'
#' \code{\link{ivamse}} reports this for every candidate in the
#' \code{implied_alpha} column, so it rarely needs to be called directly.
#'
#' @param lambda penalty level, or a vector of penalty levels, on the scale
#'   described in \code{\link{bcch_lambda}}.
#' @param p number of columns in the instrument set.
#' @param n effective sample size; see \code{\link{bcch_lambda}}.
#' @param sigma_v standard deviation of the first-stage error. The verdict
#'   depends on it, so a value that does not match the one used in fitting can
#'   flip the sign. A fitted model reports it as \code{fit$sigma_v}.
#' @param c_lambda inflation constant; \code{1.1} by default, as in
#'   \code{\link{bcch_lambda}}.
#'
#' @details The condition is
#' \deqn{\lambda \ge c_\lambda \sqrt{2} \sigma_v
#'       \sqrt{(\log(2p) + \alpha \log\log(p \vee n)) / n}}{%
#'   lambda >= c_lambda * sqrt(2) * sigma_v *
#'   sqrt((log(2p) + alpha * log(log(max(p, n)))) / n)}
#' for some \eqn{\alpha > 0}{alpha > 0}, and the function returns the largest
#' \eqn{\alpha}{alpha} the penalty supports.
#'
#' A negative value is not an error and does not invalidate the fit. Lightly
#' penalized candidates routinely fail the condition; they are still fitted,
#' still scored, and may still be selected. The value is reported so that a
#' reader can see which candidates lie outside the range the theory covers, and
#' \code{summary} says so when the selected one does.
#'
#' @return A numeric vector the same length as \code{lambda}. Positive entries
#'   satisfy the condition; negative entries do not.
#'
#' @seealso \code{\link{bcch_lambda}}, \code{\link{ivamse}}.
#'
#' @examples
#' ## across a grid of multiples of the plug-in penalty, the lightly penalized
#' ## candidates fail the condition and the heavily penalized ones satisfy it
#' lambda <- 2^seq(-4, 2, by = 0.5) * bcch_lambda(50, 500)
#' data.frame(kappa = 2^seq(-4, 2, by = 0.5),
#'            alpha = round(implied_alpha(lambda, p = 50, n = 500), 2))
#'
#' @export
implied_alpha <- function(lambda, p, n, sigma_v = 1, c_lambda = 1.1) {
  stopifnot(p >= 1, n >= 1, sigma_v > 0, c_lambda > 1)
  ## log(log(.)) is the rate the assumption is stated at, and is positive only
  ## once max(p, n) exceeds e
  if (max(p, n) <= exp(1)) return(rep(NA_real_, length(lambda)))
  (n * lambda^2 / (2 * c_lambda^2 * sigma_v^2) - log(2 * p)) /
    log(log(max(p, n)))
}


## smallest penalty at which every coefficient is zero: the subgradient
## condition at pi = 0 is |z_j'x| / n <= lambda for every column. Columns are
## already scaled to unit second moment, so no further standardization enters.
lambda_max <- function(z, x, n) {
  max(abs(crossprod(z, x))) / n
}


## geometric grid of nlambda penalties descending from lambda_max, with glmnet's
## default ratio. Returned increasing in lambda, the candidate order used
## throughout. The grid includes lambda_max itself, whose fit is null, so one
## candidate per dictionary is always empty and ineligible.
lambda_path <- function(z, x, n, nlambda = 25, lambda.min.ratio = NULL) {
  p <- ncol(z)
  if (is.null(lambda.min.ratio)) {
    lambda.min.ratio <- if (n < p) 1e-2 else 1e-4
  }
  stopifnot(nlambda >= 1, lambda.min.ratio > 0, lambda.min.ratio < 1)

  top <- lambda_max(z, x, n)
  if (!is.finite(top) || top <= 0) {
    stop("the dictionary has no variation left after partialling out the controls")
  }
  if (nlambda == 1) return(top)
  sort(top * lambda.min.ratio^seq(0, 1, length.out = nlambda))
}


## penalty grid for one dictionary under either rule, with the multiplier kappa
## expressing each penalty as a fraction of the plug-in level. kappa is reported
## whichever rule built the grid.
dictionary_lambda <- function(z, x, n, penalty, lambda, kappa, nlambda,
                              lambda.min.ratio, sigma_v, c_lambda) {
  p <- ncol(z)

  grid <- if (!is.null(lambda)) {
    sort(as.numeric(lambda))
  } else if (penalty == "bcch") {
    sort(kappa) * bcch_lambda(p, n, sigma_v, c_lambda)
  } else {
    lambda_path(z, x, n, nlambda, lambda.min.ratio)
  }

  if (any(!is.finite(grid)) || any(grid <= 0)) {
    stop("penalties must be finite and strictly positive")
  }
  if (anyDuplicated(grid)) {
    stop("the penalty grid for a dictionary contains duplicate values")
  }

  plugin <- bcch_lambda(p, n, sigma_v, c_lambda)
  list(lambda = grid, kappa = grid / plugin, plugin = plugin,
       implied_alpha = implied_alpha(grid, p, n, sigma_v, c_lambda))
}


## degrees-of-freedom corrected residual scale from least squares of x on a
## dictionary, used as sigma_v for the "bcch" grid. Uses (x, Z) only.
residual_scale <- function(z, x, n) {
  fit <- qr(z)
  residual <- qr.resid(fit, x)
  df <- n - fit$rank
  if (df <= 0) {
    ## dictionary saturates the sample: fall back to the total variation, as in
    ## the original BCCH programs
    return(sqrt(sum(x^2) / n))
  }
  sqrt(sum(residual^2) / df)
}
