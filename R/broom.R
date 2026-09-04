## tidiers for the broom package

#' \pkg{broom} Methods for an \code{"ivamse"} Fit
#'
#' @description Put a fitted model into data frames, for tables and plots.
#' \code{tidy} gives one row per coefficient, \code{glance} one row summarizing
#' the fit and the selection, and \code{augment} returns the data with fitted
#' values and residuals attached.
#'
#' These are methods for generics in the \pkg{broom} package rather than
#' functions to be called by name. Call the generic on the fitted model, as in
#' \code{broom::tidy(m)}. \pkg{broom} must be installed.
#'
#' @param x an object of class \code{"ivamse"}.
#' @param conf.int logical. For \code{tidy}, if \code{TRUE} confidence limits
#'   are added as \code{conf.low} and \code{conf.high}.
#' @param conf.level the confidence level for those limits.
#' @param data for \code{augment}, the data frame to attach the fitted values
#'   to; the model frame by default.
#' @param newdata for \code{augment}, an optional data frame of new
#'   observations, for which fitted values are computed.
#' @param ... for \code{tidy}, further arguments passed to
#'   \code{\link{vcov.ivamse}}, so that \code{cluster = NULL} gives unclustered
#'   standard errors. Not used by \code{glance} or \code{augment}.
#'
#' @return \code{tidy} returns a data frame with columns \code{term},
#'   \code{estimate}, \code{std.error}, \code{statistic} and \code{p.value}.
#'
#'   \code{glance} returns a one-row data frame with the following columns.
#'   \item{nobs}{number of observations.}
#'   \item{n.effective}{observations less the rank of the controls.}
#'   \item{n.controls}{rank of the controls, counting the intercept.}
#'   \item{n.candidates}{number of candidates scored.}
#'   \item{rule}{the selection rule used.}
#'   \item{dictionary}{name of the instrument set selected.}
#'   \item{lambda, kappa}{the selected penalty, and that penalty as a multiple
#'     of the plug-in penalty.}
#'   \item{dimension}{roughly how many instruments the selected fit uses.}
#'   \item{criterion}{the score of the selected candidate.}
#'   \item{rho.hat}{estimated correlation between the structural and
#'     first-stage errors.}
#'   \item{sigma}{residual standard error.}
#'
#'   \code{augment} returns \code{data} with \code{.fitted} and \code{.resid}
#'   added, or \code{newdata} with \code{.fitted} added.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{selected}},
#'   \code{\link{vcov.ivamse}}.
#'
#' @examples
#' set.seed(1)
#' N <- 300
#' z <- matrix(rnorm(N * 10), N, 10)
#' v <- rnorm(N)
#' x <- z[, 1] + v
#' y <- 2 * x + 0.7 * v + rnorm(N)
#' m <- ivamse_fit(y, x, z)
#'
#' broom::tidy(m, conf.int = TRUE)
#' broom::glance(m)
#'
#' @name ivamse-tidiers
#' @usage \method{tidy}{ivamse}(x, conf.int = FALSE, conf.level = 0.95, ...)
#' @rdname ivamse-tidiers
#' @exportS3Method broom::tidy
tidy.ivamse <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  table <- summary(x, ...)$coefficients
  out <- data.frame(
    term = rownames(table),
    estimate = table[, 1L],
    std.error = table[, 2L],
    statistic = table[, 3L],
    p.value = table[, 4L],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  if (conf.int) {
    ci <- stats::confint(x, level = conf.level, ...)
    out$conf.low <- ci[, 1L]
    out$conf.high <- ci[, 2L]
  }
  out
}

#' @usage \method{glance}{ivamse}(x, ...)
#' @rdname ivamse-tidiers
#' @exportS3Method broom::glance
glance.ivamse <- function(x, ...) {
  data.frame(
    nobs = x$nobs,
    n.effective = x$n,
    n.controls = x$q,
    n.candidates = nrow(x$candidates),
    rule = x$select,
    dictionary = x$selected$dictionary,
    lambda = x$selected$lambda,
    kappa = x$selected$kappa,
    dimension = x$selected$d,
    criterion = x$selected$criterion,
    rho.hat = x$rho_hat,
    sigma = sqrt(sum(x$residuals^2) / x$df.residual),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @usage \method{augment}{ivamse}(x, data = x$model, newdata = NULL, ...)
#' @rdname ivamse-tidiers
#' @exportS3Method broom::augment
augment.ivamse <- function(x, data = x$model, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    newdata$.fitted <- stats::predict(x, newdata = newdata)
    return(newdata)
  }
  if (is.null(data)) {
    stop("no model frame stored; refit with model = TRUE or supply 'data'")
  }
  data$.fitted <- stats::fitted(x)
  data$.resid <- stats::residuals(x)
  data
}
