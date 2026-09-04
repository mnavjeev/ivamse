## The formula interface, using the two shapes ivreg accepts:
##
##   y ~ x + w | z + w    regressors | instruments
##   y ~ w | x | z        exogenous | endogenous | excluded instruments
##
## The three-part form is preferable here, since a dictionary of technical
## instruments can have hundreds of columns and the two-part form requires
## writing all of them twice.


## text of one part of a formula, so that parts can be reassembled; deparse()
## splits long expressions across strings, so collapse the pieces
part_text <- function(f, side, k) {
  paste(deparse(stats::formula(f, lhs = if (side == "lhs") 1 else 0,
                               rhs = if (side == "lhs") 0 else k)[[2]]),
        collapse = " ")
}


## assemble the model formula, the dictionary formulas and any grouping
## formulas into one multi-part formula, so that a single model frame covers
## every variable and 'subset' and 'na.action' apply to all of them alike
combine_parts <- function(formula, dictionaries, extra = list()) {
  f <- Formula::as.Formula(formula)
  ## for a Formula, length() returns c(number of LHS parts, number of RHS parts)
  n_parts <- length(f)[2L]
  if (!n_parts %in% c(2L, 3L)) {
    stop("'formula' must have two or three parts on the right-hand side, ",
         "as in y ~ x + w | z + w or y ~ w | x | z")
  }
  if (n_parts == 2L && length(dictionaries)) {
    stop("'dictionaries' requires a three-part formula, as in y ~ w | x | z")
  }

  parts <- vapply(seq_len(n_parts), function(k) part_text(f, "rhs", k), character(1))
  for (d in dictionaries) {
    parts <- c(parts, paste(deparse(d[[length(d)]]), collapse = " "))
  }
  ## append cluster and grouping variables as final parts, so that 'subset' and
  ## 'na.action' apply to them too
  for (e in extra) {
    if (inherits(e, "formula")) {
      parts <- c(parts, paste(deparse(e[[length(e)]]), collapse = " "))
    }
  }

  list(
    formula = Formula::as.Formula(stats::as.formula(
      paste(part_text(f, "lhs", 1L), "~", paste(parts, collapse = " | ")),
      env = environment(formula)
    )),
    n_parts = n_parts
  )
}


## one part of the model frame as a matrix, without its intercept column
part_matrix <- function(f, k, mf, contrasts) {
  m <- stats::model.matrix(stats::formula(f, lhs = 0, rhs = k), mf,
                           contrasts.arg = contrasts)
  m[, colnames(m) != "(Intercept)", drop = FALSE]
}


## split the regressor parts of a model frame into the endogenous column and the
## included exogenous controls. Used when fitting and when predicting, so that
## both build the same matrices in the same column order.
split_regressors <- function(f, n_parts, mf, contrasts) {
  if (n_parts == 3L) {
    return(list(endogenous = part_matrix(f, 2L, mf, contrasts),
                w = part_matrix(f, 1L, mf, contrasts)))
  }
  regressors <- part_matrix(f, 1L, mf, contrasts)
  instruments <- part_matrix(f, 2L, mf, contrasts)
  is_endogenous <- !colnames(regressors) %in% colnames(instruments)
  list(
    endogenous = regressors[, is_endogenous, drop = FALSE],
    w = regressors[, !is_endogenous, drop = FALSE],
    excluded = instruments[, !colnames(instruments) %in% colnames(regressors),
                           drop = FALSE]
  )
}


#' Instrumental-Variable Regression with a Lasso First Stage
#'
#' @description Fits an instrumental-variables regression in which the first
#' stage is estimated by the lasso. Where 2SLS uses one fixed set of
#' instruments, \code{ivamse} is given several candidate sets and several
#' penalty levels, fits the first stage once for each combination, and keeps the
#' one that minimizes an estimate of the mean squared error of the coefficient
#' on the endogenous regressor. The estimator is otherwise an ordinary
#' just-identified IV regression, and the usual extractor and inference methods
#' apply to the fitted object.
#'
#' @param formula a model formula with two or three parts on the right-hand
#'   side. In the two-part form \code{y ~ x + w | z + w} the first part gives
#'   the regressors and the second the instruments, with the exogenous
#'   regressors appearing in both, as in \code{\link[ivreg]{ivreg}}. In the
#'   three-part form \code{y ~ w | x | z} the parts give the exogenous
#'   regressors, the endogenous regressor, and the excluded instruments. The
#'   three-part form is usually more convenient here, since an instrument set
#'   may have hundreds of columns. Exactly one endogenous regressor is allowed.
#' @param data an optional data frame containing the variables in the model. By
#'   default the variables are taken from the environment of the
#'   \code{formula}.
#' @param subset an optional vector specifying a subset of observations to be
#'   used in fitting the model.
#' @param na.action a function that indicates what should happen when the data
#'   contain \code{NA}s. The default is set by the \code{na.action} option.
#' @param contrasts an optional list. See the \code{contrasts.arg} of
#'   \code{\link[stats:model.matrix]{model.matrix.default}}.
#' @param model logical. If \code{TRUE}, the default, the model frame is
#'   returned as a component of the fitted object.
#' @param dictionaries an optional named list of one-sided formulas, one per
#'   candidate set of instruments, evaluated in \code{data}. When supplied these
#'   replace the excluded instruments given in \code{formula}, whose third part
#'   is then ignored, and the three-part form is required. The order of the list
#'   is the order in which candidates are scored, and ties in the criterion are
#'   broken toward the earlier one, so the simpler sets should be listed first.
#'   To supply matrices instead of formulas, use \code{\link{ivamse_fit}}.
#' @param cluster an optional one-sided formula, vector of cluster identifiers,
#'   or name of a column of \code{data}, giving the clustering used for the
#'   reported standard errors. It does not affect which candidate is selected.
#' @param cv_group an optional one-sided formula, vector, or column name. Rows
#'   sharing a value are kept in the same cross-validation fold. It has an
#'   effect only when \code{select = "cv"}.
#' @param ... further arguments passed to \code{\link{ivamse_fit}}, which
#'   documents them: \code{lambda}, \code{nlambda}, \code{lambda.min.ratio},
#'   \code{penalty}, \code{kappa}, \code{c_lambda}, \code{sigma_v},
#'   \code{select}, \code{pilot}, \code{screen}, and the cross-validation
#'   settings \code{nfolds}, \code{foldid} and \code{cv_s}. The fold rule is
#'   named \code{cv_s} rather than \code{s} because \code{s} would partially
#'   match \code{subset}.
#'
#' @details
#' A candidate is one instrument set paired with one penalty level. Every
#' candidate is a consistent first stage, and first-order asymptotics therefore
#' cannot separate them: they all leave the estimator of the structural
#' coefficient with the same limiting distribution. The criterion separates them
#' by the next terms in the expansion, which trade the fit of the first stage
#' against a many-instrument bias of the kind studied by Donald and Newey
#' (2001). Because the weight on that bias grows with the correlation between
#' the structural and first-stage errors, a strongly endogenous regressor is met
#' with a heavier penalty and a simpler first stage than a nearly exogenous one.
#' \code{\link{feasible_criterion}} gives the score itself, and
#' \code{\link{ivamse_fit}} documents the alternatives \code{select = "cv"} and
#' \code{select = "bcch"}, which score the first stage alone.
#'
#' Included exogenous controls are removed from the outcome, from the endogenous
#' regressor and from every instrument by least squares before anything else
#' happens. They are never penalized, so they are always kept and their
#' coefficients are returned alongside the endogenous one. The effective sample
#' size \code{n} is the number of observations less the rank of the controls,
#' and it, rather than \code{nrow(data)}, is what the penalties are scaled by.
#'
#' \code{vcov} uses the clustering supplied at fitting, and the homoscedastic
#' form otherwise. Because \code{bread} and \code{estfun} methods are provided,
#' the variance estimators of the \pkg{sandwich} package apply directly, and
#' \code{summary}, \code{confint} and \code{predict} take a \code{vcov.}
#' argument that is either a covariance matrix or a function returning one.
#' Standard errors treat the selected first stage as fixed. They therefore
#' account for neither the estimation error in that first stage nor the
#' selection step, and they are not valid under weak identification.
#'
#' The theory is derived for one endogenous regressor and homoscedastic Gaussian
#' errors. Clustering and heteroscedasticity affect the reported standard errors
#' only; the criterion is always evaluated in its homoscedastic form.
#'
#' @return \code{ivamse} returns an object of class \code{"ivamse"} with the
#'   following components.
#'   \item{coefficients}{the coefficient on the endogenous regressor, followed
#'     by the coefficients on the controls.}
#'   \item{residuals}{vector of structural residuals.}
#'   \item{fitted.values}{vector of fitted values for the response.}
#'   \item{candidates}{data frame with one row per candidate; see Details of
#'     \code{\link{ivamse_fit}} for its columns.}
#'   \item{selected}{list describing the candidate used: \code{index} is the one
#'     the rule chose, \code{estimation_index} the one actually used, and the
#'     remaining fields describe the latter. The two differ only when
#'     \code{fallback} is \code{TRUE}.}
#'   \item{instrument}{the fitted first stage, of length \code{nobs}, used as
#'     the single instrument.}
#'   \item{first_stage, scale}{the first-stage coefficients on the normalized
#'     instrument set, and the column scales used to normalize it.}
#'   \item{sigma_eps2, sigma_epsv}{the variance of the structural error and its
#'     covariance with the first-stage error, both estimated from the pilot.}
#'   \item{rho_hat}{the implied correlation between the two errors.}
#'   \item{sigma_v}{the first-stage error scale; see \code{\link{ivamse_fit}}.}
#'   \item{first_stage_F, instrument_F}{least-squares first-stage F statistic on
#'     the reference instrument set, and the F statistic on the fitted
#'     instrument. The latter is computed from an instrument built using the
#'     endogenous regressor and is biased upward.}
#'   \item{pilot}{list naming the pilot candidate and its instrument set.}
#'   \item{n, N, q}{the effective sample size, the number of observations, and
#'     the rank of the controls.}
#'   \item{df.residual}{residual degrees of freedom.}
#'   \item{cov.unscaled}{unscaled covariance matrix for the coefficients.}
#'   \item{foldid}{the cross-validation folds used, or \code{NULL}.}
#'   \item{call, formula, terms, levels, contrasts, na.action, model}{the usual
#'     model-fitting components.}
#'
#' @references
#' Belloni, A., Chen, D., Chernozhukov, V., and Hansen, C. (2012). Sparse models
#' and methods for optimal instruments with an application to eminent domain.
#' \emph{Econometrica} \bold{80}, 2369--2429.
#'
#' Donald, S. G., and Newey, W. K. (2001). Choosing the number of instruments.
#' \emph{Econometrica} \bold{69}, 1161--1191.
#'
#' Ma, Y., Navjeevan, M., and Salahub, B. Choosing the dictionary and penalty
#' for IV-LASSO. Working paper.
#'
#' @seealso \code{\link{ivamse_fit}} for the arguments controlling the candidate
#'   list and the selection rule, \code{\link{selected}} for what was chosen,
#'   \code{\link{plot.ivamse}} for the criterion along the penalty grid,
#'   \code{\link{vcov.ivamse}} and \code{\link{ivamse-sandwich}} for standard
#'   errors, \code{\link{ivamse-tidiers}} for \pkg{broom} methods,
#'   \code{\link{feasible_criterion}}, \code{\link{effective_dimension}},
#'   \code{\link{bcch_lambda}} and \code{\link{implied_alpha}} for the
#'   quantities involved, and \code{vignette("ivamse")}.
#'
#' @examples
#' ## a first stage in which only the first few instruments matter
#' set.seed(1)
#' N <- 400
#' d <- data.frame(w = rnorm(N))
#' Z <- matrix(rnorm(N * 15), N, 15)
#' colnames(Z) <- paste0("z", 1:15)
#' d <- cbind(d, Z)
#' v <- rnorm(N)
#' d$x <- d$z1 + 0.5 * d$z2 + d$w + v
#' d$y <- 2 * d$x + d$w + 0.7 * v + rnorm(N)
#'
#' ## exogenous | endogenous | excluded instruments; the true coefficient is 2
#' m <- ivamse(y ~ w | x | z1 + z2 + z3 + z4 + z5, data = d)
#' coef(m)["x"]
#' summary(m)
#'
#' ## offer three instrument sets and let the criterion choose among them
#' m2 <- ivamse(y ~ w | x | z1,
#'   dictionaries = list(
#'     small  = ~ z1 + z2,
#'     medium = ~ z1 + z2 + z3 + z4 + z5,
#'     wide   = ~ (z1 + z2 + z3) * (z4 + z5)
#'   ),
#'   data = d)
#' selected(m2)
#'
#' @export
ivamse <- function(formula, data, subset, na.action, contrasts = NULL,
                   model = TRUE, dictionaries = NULL, cluster = NULL,
                   cv_group = NULL, ...) {
  cl <- match.call()

  ## s = "lambda.1se" partially matches 'subset', which would subset the data
  ## instead of setting the cross-validation rule
  if (!missing(subset)) {
    guess <- tryCatch(subset, error = function(e) NULL)
    if (is.character(guess) && length(guess) == 1L &&
        guess %in% c("lambda.min", "lambda.1se")) {
      stop("'subset' received \"", guess,
           "\"; the cross-validation rule is set by 'cv_s'")
    }
  }
  dots <- list(...)
  clash <- intersect(names(dots), c("x_name", "intercept", "cluster", "w", "z"))
  if (length(clash)) {
    stop("argument", if (length(clash) > 1L) "s " else " ",
         paste(sQuote(clash), collapse = ", "),
         " ", if (length(clash) > 1L) "are" else "is",
         " determined by the formula; use ivamse_fit()")
  }

  if (!is.null(dictionaries)) {
    if (!is.list(dictionaries) || is.null(names(dictionaries)) ||
        any(!nzchar(names(dictionaries)))) {
      stop("'dictionaries' must be a named list")
    }
    if (!all(vapply(dictionaries, inherits, logical(1), "formula"))) {
      stop("'dictionaries' must be one-sided formulas; to supply matrices ",
           "directly, use ivamse_fit()")
    }
  }

  ## a bare column name becomes a one-sided formula, so that it reaches the
  ## model frame like any other variable
  cluster <- as_group_formula(cluster)
  cv_group <- as_group_formula(cv_group)

  combined <- combine_parts(formula, dictionaries, list(cluster, cv_group))
  f <- combined$formula
  n_parts <- combined$n_parts

  ## ---- Model frame ---------------------------------------------------------
  mf_call <- match.call(expand.dots = FALSE)
  keep <- match(c("data", "subset", "na.action"), names(mf_call), 0L)
  mf_call <- mf_call[c(1L, keep)]
  mf_call$formula <- f
  mf_call$drop.unused.levels <- TRUE
  mf_call[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf_call, parent.frame())

  y <- stats::model.response(mf, "numeric")
  intercept <- attr(stats::terms(stats::formula(f, lhs = 0, rhs = 1)),
                    "intercept") == 1

  ## ---- Split the regressors into endogenous and exogenous ------------------
  split <- split_regressors(f, n_parts, mf, contrasts)
  endogenous <- split$endogenous
  w <- split$w

  if (ncol(endogenous) != 1L) {
    stop("exactly one endogenous regressor is required, but ",
         ncol(endogenous), " were found",
         if (n_parts == 2L) {
           "; in the two-part form these are the regressors not among the instruments"
         } else "")
  }
  if (ncol(w) == 0L) w <- NULL

  ## ---- Dictionaries --------------------------------------------------------
  z <- if (is.null(dictionaries)) {
    list(instruments = if (n_parts == 3L) part_matrix(f, 3L, mf, contrasts)
                       else split$excluded)
  } else {
    stats::setNames(
      lapply(seq_along(dictionaries), function(i) {
        part_matrix(f, n_parts + i, mf, contrasts)
      }),
      names(dictionaries)
    )
  }
  if (any(vapply(z, ncol, integer(1)) == 0L)) {
    stop("a dictionary has no columns")
  }

  ## ---- Fit -----------------------------------------------------------------
  fit <- ivamse_fit(
    y = y, x = endogenous[, 1L], z = z, w = w, intercept = intercept,
    cluster = resolve_variable(cluster, mf),
    cv_group = resolve_variable(cv_group, mf),
    x_name = colnames(endogenous), ...
  )

  fit$call <- cl
  fit$formula <- formula
  fit$parts <- list(f = f, n_parts = n_parts)
  fit$cluster_label <- variable_label(cluster)
  fit$terms <- stats::terms(mf)
  fit$levels <- stats::.getXlevels(stats::terms(mf), mf)
  fit$contrasts <- contrasts
  fit$na.action <- attr(mf, "na.action")
  if (model) fit$model <- mf
  fit
}


## a grouping variable may be a one-sided formula, a column name, or a vector.
## The first two are looked up in the model frame, so they are subset and
## filtered for missingness along with everything else.
resolve_variable <- function(value, mf) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "formula")) {
    frame <- stats::model.frame(value, data = mf, na.action = stats::na.pass)
    return(if (ncol(frame) == 1L) frame[[1L]] else frame)
  }
  value
}


## a bare column name given as a length-one character is turned into a
## one-sided formula; anything else is passed through
as_group_formula <- function(value) {
  if (is.character(value) && length(value) == 1L) {
    return(stats::as.formula(paste("~", value)))
  }
  value
}


## a short name for a grouping variable, for printing
variable_label <- function(value) {
  if (is.null(value)) return(NULL)
  if (inherits(value, "formula")) {
    return(paste(all.vars(value), collapse = " + "))
  }
  if (is.character(value) && length(value) == 1L) return(value)
  NULL
}
