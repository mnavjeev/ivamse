## standard errors: given the selected candidate the estimator is just-identified
## IV, with the usual variance. The selected first stage is treated as fixed, so
## the standard errors reflect neither its estimation error nor the selection
## step. The homoscedastic variance is invariant to how the controls are
## partialled out; a robust or clustered one is not.


#' Covariance Matrix for an \code{"ivamse"} Fit
#'
#' @description Returns the estimated covariance matrix of the coefficients.
#' The square roots of its diagonal are the standard errors reported by
#' \code{summary}. If the model was fitted with a \code{cluster} argument, the
#' clustered covariance is returned; otherwise the ordinary one, which assumes
#' the errors have constant variance.
#'
#' @param object an object of class \code{"ivamse"}.
#' @param cluster the clustering to use: a one-sided formula, a vector of
#'   cluster identifiers, or a data frame of them. It defaults to whatever was
#'   supplied when the model was fitted, so passing \code{cluster = NULL} forces
#'   the ordinary covariance.
#' @param ... further arguments passed to \code{\link[sandwich]{vcovCL}} when
#'   clustering, such as \code{type}. Without clustering there is nothing to
#'   pass them to and supplying any is a warning, since a \code{type} argument
#'   silently ignored would leave a non-robust standard error looking robust.
#'
#' @details For heteroscedasticity-consistent covariances without clustering,
#' use \code{\link{vcovHC.ivamse}}. Any other estimator in the \pkg{sandwich}
#' package may also be applied directly to the fitted object.
#'
#' Standard errors from any of these treat the selected first stage as fixed.
#' They account for neither the estimation error in that first stage nor the
#' fact that the candidate was chosen using the same data, and they are not
#' valid under weak identification.
#'
#' @return A square covariance matrix with one row and column per coefficient.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{vcovHC.ivamse}},
#'   \code{\link{ivamse-sandwich}}, \code{\link[sandwich]{vcovCL}}.
#'
#' @examples
#' set.seed(1)
#' N <- 400
#' d <- data.frame(firm = factor(rep(1:40, each = 10)))
#' Z <- matrix(rnorm(N * 12), N, 12)
#' colnames(Z) <- paste0("z", 1:12)
#' d <- cbind(d, Z)
#' v <- rnorm(N)
#' d$x <- d$z1 + 0.6 * d$z2 + v
#' d$y <- 2 * d$x + 0.7 * v + rnorm(N)
#'
#' m <- ivamse(y ~ 1 | x | z1 + z2 + z3 + z4, data = d, cluster = ~ firm)
#'
#' ## standard errors are the square roots of the diagonal
#' sqrt(diag(vcov(m)))
#'
#' ## the same model without clustering
#' sqrt(diag(vcov(m, cluster = NULL)))
#'
#' @export
vcov.ivamse <- function(object, cluster = object$cluster, ...) {
  if (is.null(cluster)) {
    ## the homoscedastic branch takes no further arguments; silently dropping
    ## something like type = "HC1" would report a non-robust standard error as
    ## though it were robust
    if (...length()) {
      warning("arguments ", paste(sQuote(...names()), collapse = ", "),
              " are ignored without clustering; use sandwich::vcovHC()",
              call. = FALSE)
    }
    ## this residual variance carries the usual degrees-of-freedom correction,
    ## so it is not the sigma_eps2 the criterion uses, which divides by n
    sigma2 <- sum(object$residuals^2) / object$df.residual
    v <- sigma2 * object$cov.unscaled
    dimnames(v) <- list(names(object$coefficients), names(object$coefficients))
    return(v)
  }
  if (!requireNamespace("sandwich", quietly = TRUE)) {
    stop("clustered standard errors require the 'sandwich' package")
  }
  sandwich::vcovCL(object, cluster = cluster, ...)
}


#' Sandwich Building Blocks for an \code{"ivamse"} Fit
#'
#' @description Methods that let the \pkg{sandwich} package compute robust and
#' clustered covariance matrices for a fitted model. \code{estfun} returns each
#' observation's contribution to the estimating equations, and \code{bread} the
#' matrix that scales them.
#'
#' These are methods for generics in other packages rather than functions to be
#' called by name. Call the generic on the fitted model, as in
#' \code{sandwich::vcovCL(m, cluster = ~ firm)}. Most users need neither
#' directly: \code{\link{vcov.ivamse}} and \code{\link{vcovHC.ivamse}} cover the
#' usual cases.
#'
#' @param x,object an object of class \code{"ivamse"}.
#' @param ... currently not used.
#'
#' @details Both are built from the projected regressors, the fitted first stage
#' together with the controls, and the structural residuals. This is what makes
#' the \pkg{sandwich} estimators give the same answers here as they would for
#' the equivalent two-stage least-squares fit.
#'
#' @return \code{estfun} returns a matrix with one row per observation and one
#'   column per coefficient. \code{bread} returns a square matrix with one row
#'   and column per coefficient. \code{model.matrix} for the internal
#'   \code{"ivamse_projected"} class returns the projected regressors.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{vcov.ivamse}},
#'   \code{\link{vcovHC.ivamse}}, \code{\link[sandwich]{vcovCL}}.
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
#' dim(sandwich::estfun(m))
#' sqrt(diag(sandwich::vcovCL(m, cluster = rep(1:30, each = 10))))
#'
#' @name ivamse-sandwich
#' @usage \method{estfun}{ivamse}(x, ...)
#' @rdname ivamse-sandwich
#' @exportS3Method sandwich::estfun
estfun.ivamse <- function(x, ...) {
  ef <- x$residuals * x$projected
  colnames(ef) <- names(x$coefficients)
  ef
}

#' @usage \method{bread}{ivamse}(x, ...)
#' @rdname ivamse-sandwich
#' @exportS3Method sandwich::bread
bread.ivamse <- function(x, ...) {
  b <- x$cov.unscaled * x$nobs
  dimnames(b) <- list(names(x$coefficients), names(x$coefficients))
  b
}


## resolve the vcov. argument, which may be a matrix, a function, or NULL
resolve_vcov <- function(object, vcov., ...) {
  if (is.null(vcov.)) return(stats::vcov(object, ...))
  if (is.function(vcov.)) return(vcov.(object, ...))
  vcov.
}


#' @rdname ivamse
#' @param parm parameters for which confidence intervals are to be computed; a
#'   vector of numbers or names. The default is all parameters.
#' @param level confidence level; the default is \code{0.95}.
#' @param vcov. optionally either a coefficient covariance matrix or a function
#'   to compute one from a fitted \code{ivamse} object. If \code{NULL} (the
#'   default) \code{\link{vcov.ivamse}} is used. If \code{vcov.} is a function,
#'   the \code{...} argument can be used to pass on further arguments to it.
#' @param df optional residual degrees of freedom for the reference
#'   distribution. \code{NULL} and \code{Inf} both use the normal.
#' @param object,x an object of class \code{"ivamse"}.
#' @export
confint.ivamse <- function(object, parm, level = 0.95, vcov. = NULL,
                           df = NULL, ...) {
  estimate <- stats::coef(object)
  se <- sqrt(diag(resolve_vcov(object, vcov., ...)))
  if (is.null(df)) df <- Inf

  a <- (1 - level) / 2
  quantile <- if (is.finite(df)) stats::qt(1 - a, df) else stats::qnorm(1 - a)

  ci <- cbind(estimate - quantile * se, estimate + quantile * se)
  colnames(ci) <- paste(format(100 * c(a, 1 - a), trim = TRUE), "%")
  rownames(ci) <- names(estimate)

  if (missing(parm)) return(ci)
  if (is.numeric(parm)) parm <- names(estimate)[parm]
  ci[parm, , drop = FALSE]
}


#' Robust Covariance Matrix for an \code{"ivamse"} Fit
#'
#' @description Heteroscedasticity-consistent ("robust") covariance matrix,
#' the method behind \code{\link[sandwich]{vcovHC}} for fitted models of class
#' \code{"ivamse"}. The square roots of its diagonal are robust standard errors.
#'
#' Call it as \code{sandwich::vcovHC(m)} rather than by name.
#'
#' @param x an object of class \code{"ivamse"}.
#' @param ... further arguments passed to \code{\link[sandwich]{vcovHC}},
#'   notably \code{type}, which selects the small-sample adjustment and takes
#'   \code{"HC3"} by default. Users coming from Stata's \code{robust} option
#'   will want \code{type = "HC1"}.
#'
#' @details Any clustering supplied when the model was fitted is \emph{not} used
#' here; this is the unclustered robust covariance. For clustered standard
#' errors use \code{\link{vcov.ivamse}} or
#' \code{\link[sandwich]{vcovCL}}.
#'
#' As with every covariance this package reports, the selected first stage is
#' treated as fixed, so the standard errors reflect neither its estimation error
#' nor the selection step.
#'
#' @return A square covariance matrix with one row and column per coefficient.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{vcov.ivamse}},
#'   \code{\link{ivamse-sandwich}}, \code{\link[sandwich]{vcovHC}}.
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
#' sqrt(diag(sandwich::vcovHC(m, type = "HC1")))
#'
#' @usage \method{vcovHC}{ivamse}(x, ...)
#' @exportS3Method sandwich::vcovHC
vcovHC.ivamse <- function(x, ...) {
  class(x) <- c("ivamse_projected", "ivamse")
  sandwich::vcovHC.default(x, ...)
}

#' @rdname ivamse-sandwich
#' @export
model.matrix.ivamse_projected <- function(object, ...) object$projected


#' Leverage Values for an \code{"ivamse"} Fit
#'
#' @description Leverage of each observation, measuring how far its regressors
#' lie from the average. \code{\link[sandwich]{vcovHC}} needs these for its
#' \code{"HC2"} and \code{"HC3"} small-sample adjustments.
#'
#' @param model an object of class \code{"ivamse"}.
#' @param ... currently not used.
#'
#' @return A numeric vector with one value per observation. The values are
#'   between zero and one and sum to the number of coefficients.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{vcovHC.ivamse}}.
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
#' summary(hatvalues(m))
#' sum(hatvalues(m))
#'
#' @export
hatvalues.ivamse <- function(model, ...) {
  ## leverage from the projected regressors, needed by HC2 and HC3
  as.vector(rowSums((model$projected %*% model$cov.unscaled) * model$projected))
}
