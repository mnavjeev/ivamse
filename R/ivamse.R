## The estimator. ivamse_fit() is the engine; ivamse() in formula.R parses a
## model formula and calls it. The work is split into four steps:
## prepare_model(), prepare_dictionaries(), build_candidates() and
## choose_index(). ivamse_fit() runs them in order and assembles the result.


## dictionaries arrive either as one matrix or as a named list of matrices
as_dictionary_list <- function(z) {
  if (is.list(z) && !is.data.frame(z)) {
    if (!length(z)) stop("at least one dictionary is required")
    if (is.null(names(z)) || any(!nzchar(names(z)))) {
      stop("every dictionary in the list must be named")
    }
    if (anyDuplicated(names(z))) stop("dictionary names must be unique")
    return(lapply(z, as.matrix))
  }
  list(instruments = as.matrix(z))
}


## an argument that may be given once for all dictionaries or as a named list
per_dictionary <- function(value, name) {
  if (is.null(value)) return(NULL)
  if (is.list(value)) {
    if (!name %in% names(value)) {
      stop("no entry named '", name, "' in the per-dictionary argument")
    }
    return(value[[name]])
  }
  value
}


## partial the controls out of the outcome and the endogenous regressor.
## Aliased control columns are dropped rather than only warned about, since
## keeping them leaves a singular cross-product later.
prepare_model <- function(y, x, w, intercept) {
  N <- length(y)
  if (intercept) w <- cbind("(Intercept)" = rep(1, N), w)

  if (!is.null(w) && ncol(w) > 0L) {
    decomposition <- qr(w)
    if (decomposition$rank < ncol(w)) {
      keep <- sort(decomposition$pivot[seq_len(decomposition$rank)])
      warning("some controls are collinear, dropping ",
              paste(colnames(w)[-keep], collapse = ", "), call. = FALSE)
      w <- w[, keep, drop = FALSE]
    }
  } else {
    w <- NULL
  }

  w_qr <- control_qr(w)
  q <- if (is.null(w_qr)) 0L else w_qr$rank
  n <- N - q
  if (n <= 1L) stop("the controls leave no residual variation")

  list(w = w, w_qr = w_qr, q = q, n = n, N = N,
       y_res = partial_out(y, w_qr), x_res = partial_out(x, w_qr))
}


## residualize and normalize every dictionary, and refuse any that contains the
## endogenous regressor itself
prepare_dictionaries <- function(z, model) {
  prepared <- lapply(z, prepare_dictionary, w_qr = model$w_qr, n = model$n)

  ## a dictionary column that is numerically the endogenous regressor makes the
  ## first stage trivially perfect and the IV estimate meaningless. This happens
  ## most easily by writing a dot on the instrument side of a formula.
  scale_x <- sqrt(sum(model$x_res^2) / model$n)
  if (scale_x > 0) {
    for (name in names(prepared)) {
      score <- abs(crossprod(prepared[[name]]$z, model$x_res) / model$n) / scale_x
      offending <- which(score > 1 - 1e-8)
      if (length(offending)) {
        label <- colnames(prepared[[name]]$z)[offending[1L]]
        stop("column ", if (is.null(label)) offending[1L] else paste0("'", label, "'"),
             " of dictionary '", name, "' is collinear with the endogenous ",
             "regressor; check for a '.' on the instrument side of the formula")
      }
    }
  }

  dropped <- vapply(seq_along(z), function(i) {
    ncol(z[[i]]) - ncol(prepared[[i]]$z)
  }, integer(1))
  if (any(dropped > 0L)) {
    message("dropped ", sum(dropped), " dictionary column",
            if (sum(dropped) == 1L) "" else "s",
            " lying in the span of the controls")
  }

  prepared
}


## build the penalty grid for every dictionary, fit every candidate, and collect
## the quantities the criterion needs into one table in candidate order
build_candidates <- function(prepared, model, penalty, lambda, kappa, nlambda,
                             lambda.min.ratio, sigma_v, c_lambda) {
  rows <- list()
  fits <- list()

  for (name in names(prepared)) {
    dictionary <- prepared[[name]]$z
    grid <- dictionary_lambda(
      z = dictionary, x = model$x_res, n = model$n, penalty = penalty,
      lambda = per_dictionary(lambda, name),
      kappa = per_dictionary(kappa, name),
      nlambda = per_dictionary(nlambda, name),
      lambda.min.ratio = per_dictionary(lambda.min.ratio, name),
      sigma_v = sigma_v, c_lambda = c_lambda
    )
    fit <- fit_dictionary(dictionary, model$x_res, grid$lambda, model$n,
                          qr(dictionary)$rank)

    rows[[name]] <- data.frame(
      dictionary = name,
      p = ncol(dictionary),
      ## position within this dictionary, so the selected candidate's
      ## first-stage coefficients can be recovered from 'fits' later
      column = seq_along(grid$lambda),
      lambda = grid$lambda,
      kappa = grid$kappa,
      implied_alpha = grid$implied_alpha,
      h = fit$h,
      fitted_moment = fit$fitted_moment,
      d = fit$d,
      nonzero = fit$nonzero,
      stringsAsFactors = FALSE
    )
    fits[[name]] <- fit
  }

  candidates <- do.call(rbind, rows)
  rownames(candidates) <- NULL
  candidates$index <- seq_len(nrow(candidates))

  ## 'rows' and 'fits' were filled in the same order, so column j of Pi_hat is
  ## row j of the table; everything downstream relies on that alignment
  list(candidates = candidates, fits = fits,
       Pi_hat = do.call(cbind, lapply(fits, function(f) f$Pi_hat)))
}


## apply the selection rule and return the chosen candidate index
choose_index <- function(select, candidates, eligible, prepared, model,
                         reference, c_lambda, cv_s) {
  switch(select,
    amse = select_candidate(candidates$criterion, eligible),
    cv = select_by_cv(candidates$cvm, candidates$cvsd,
                      candidates$lambda, eligible)[[cv_s]],
    bcch = {
      at <- which(candidates$dictionary == reference)
      level <- bcch_plugin_level(prepared[[reference]]$z, model$x_res,
                                 model$n, c_lambda)
      ## nearest grid point on the log scale; use penalty = "bcch" with
      ## kappa = 1 for the plug-in penalty itself
      at[which.min(abs(log(candidates$lambda[at]) - log(level)))]
    }
  )
}


## conventional first-stage F: least squares of x on the reference dictionary,
## after the controls are partialled out. It does not use the selected first
## stage, so it is not a post-selection quantity.
first_stage_f <- function(z, x_res, n) {
  p <- ncol(z)
  if (p >= n) return(NA_real_)
  residual <- qr.resid(qr(z), x_res)
  explained <- sum(x_res^2) - sum(residual^2)
  denominator <- sum(residual^2) / (n - p)
  if (denominator <= 0) return(NA_real_)
  (explained / p) / denominator
}


#' Fitting an IV-Lasso Model from Design Matrices
#'
#' @description \code{ivamse_fit} is the work-horse function behind
#' \code{\link{ivamse}}, taking design matrices rather than a model formula. It
#' fits the first stage by the lasso once for every candidate, scores them, and
#' estimates the structural coefficient from the best. \code{\link{ivamse}} is
#' the same estimator with a formula interface and is what most users want.
#'
#' This page documents the arguments that build the candidate list and choose
#' among them, which \code{\link{ivamse}} passes through unchanged.
#'
#' @param y outcome, a numeric vector of length \code{N}.
#' @param x the single endogenous regressor, a numeric vector of length \code{N}.
#' @param z excluded instruments: either one matrix, or a named list of
#'   matrices, one per candidate dictionary. Candidates are ordered by
#'   dictionary and then by increasing penalty, and ties in the criterion are
#'   broken toward the earlier candidate, so the order of this list matters.
#' @param w included exogenous controls, a matrix or \code{NULL}. These enter
#'   both equations and are partialled out; they are never penalized. The same
#'   controls apply to every dictionary.
#' @param intercept logical. If \code{TRUE}, the default, an intercept column
#'   is added to \code{w}.
#' @param lambda penalty grid. \code{NULL} (the default) builds one
#'   automatically. May be a numeric vector used for every dictionary, or a
#'   named list with one vector per dictionary. Supplying it makes
#'   \code{penalty}, \code{nlambda}, \code{lambda.min.ratio} and \code{kappa}
#'   inert.
#' @param nlambda number of penalties per dictionary, used only when
#'   \code{penalty = "path"}. May be a single number or a named list with one
#'   per dictionary. The default is 25, rather than the 100 of
#'   \code{\link[glmnet]{glmnet}}, since the criterion ranks a fixed list of
#'   candidates; the grids of Ma et al. have 13 to 29 points per dictionary.
#' @param lambda.min.ratio smallest penalty as a fraction of the largest, used
#'   only when \code{penalty = "path"}. Defaults to \code{1e-2} if the
#'   dictionary is wider than the effective sample size and \code{1e-4}
#'   otherwise, as in \pkg{glmnet}.
#' @param penalty how to build the grid. \code{"path"} (the default) runs
#'   geometrically down from the smallest penalty giving a null fit;
#'   \code{"bcch"} uses fixed multiples \code{kappa} of the plug-in penalty.
#' @param kappa multipliers of the plug-in penalty, used only when
#'   \code{penalty = "bcch"}. May be a numeric vector used for every dictionary
#'   or a named list with one vector per dictionary.
#' @param c_lambda inflation constant in the plug-in penalty; \code{1.1} in
#'   Ma et al. and in Belloni et al.
#' @param sigma_v first-stage error scale used by \code{penalty = "bcch"}.
#'   \code{NULL} estimates it by least squares of \code{x} on the reference
#'   dictionary. It affects the penalty grid only, never the reported error
#'   correlation.
#' @param select which rule picks the candidate. \code{"amse"} (the default)
#'   minimizes the criterion; \code{"cv"} and \code{"bcch"} are
#'   first-stage-only comparators, see Details.
#' @param pilot which candidate supplies the two error moments, and which
#'   dictionary is the reference. \code{"auto"} takes the strongest first stage
#'   within the first dictionary. May also be the name of a dictionary or an
#'   integer candidate index.
#' @param screen logical. If \code{TRUE}, selection is restricted to candidates
#'   whose IV denominator satisfies \code{abs(h) >= 1/sqrt(n)}. The default is
#'   \code{FALSE}, in which case every candidate with a non-degenerate first
#'   stage is eligible. The pilot candidate is screened either way.
#' @param cluster clustering for the reported standard errors: a vector of
#'   cluster identifiers of length \code{N}, or \code{NULL}. Selection is
#'   unaffected.
#' @param nfolds,foldid,cv_group,cv_s cross-validation settings, used only when
#'   \code{select = "cv"}. \code{cv_group} keeps rows sharing a value in the
#'   same fold. \code{cv_s} is \code{"lambda.min"} or \code{"lambda.1se"}, with
#'   the meanings these have in \code{\link[glmnet]{cv.glmnet}}; it is named
#'   \code{cv_s} rather than \code{s} because \code{s} would partially match the
#'   \code{subset} argument of \code{\link{ivamse}}.
#' @param x_name name to give the endogenous coefficient.
#'
#' @details
#' Both comparators score the first stage alone. Cross-validation ranks
#' candidates by how well they predict \code{x}; \code{"bcch"} takes the
#' candidate closest to the plug-in penalty of Belloni et al., with the
#' first-stage scale refined iteratively. Both ignore the outcome
#' equation and are therefore unaffected by the endogeneity of the regressor.
#' Note that
#' \code{select = "bcch"} refines \code{sigma_v} on the LASSO residual, as the
#' original programs of Belloni et al. do, whereas the \code{penalty = "bcch"}
#' grid uses the least-squares \code{sigma_v}; the two levels therefore differ.
#'
#' The reference dictionary, set by \code{pilot}, supplies the pilot candidate,
#' the default \code{sigma_v} behind the \code{"bcch"} grid, and the range
#' searched by \code{select = "bcch"}. It defaults to the first dictionary and
#' should be one fixed before the outcome is examined.
#'
#' Clustering and heteroscedasticity affect the reported standard errors only.
#' The criterion is derived under homoscedastic Gaussian errors and is always
#' evaluated in that form.
#'
#' @return An object of class \code{"ivamse"}, whose components are listed
#'   under \code{\link{ivamse}}. Its \code{candidates} component is a data frame
#'   with one row per candidate and the following columns.
#'   \item{dictionary}{name of the instrument set the candidate belongs to.}
#'   \item{p}{number of columns in that set, after any column lying in the span
#'     of the controls has been dropped.}
#'   \item{column}{position of the candidate within its own instrument set.}
#'   \item{lambda}{the penalty level.}
#'   \item{kappa}{that penalty as a multiple of the plug-in penalty, so
#'     \code{1} is the level the plug-in rule would choose and smaller values
#'     keep more instruments.}
#'   \item{implied_alpha}{positive when the penalty satisfies the condition the
#'     theory requires; see \code{\link{implied_alpha}}.}
#'   \item{h}{the IV denominator, \code{mean(Pi_hat * x)}. A value near zero
#'     means a weak first stage.}
#'   \item{fitted_moment}{second moment of the fitted first stage.}
#'   \item{d}{roughly how many instruments the fit uses; see
#'     \code{\link{effective_dimension}}.}
#'   \item{nonzero}{number of nonzero first-stage coefficients.}
#'   \item{index}{the candidate's position in the table.}
#'   \item{empty}{\code{TRUE} if the lasso shrank the whole first stage away.}
#'   \item{eligible}{\code{TRUE} if the candidate may be selected.}
#'   \item{criterion}{the score, smaller being better; see
#'     \code{\link{feasible_criterion}}.}
#'   \item{cvm, cvsd}{cross-validated prediction error of the first stage and
#'     its standard error, computed only when \code{select = "cv"}.}
#'   \item{used}{\code{TRUE} for the one candidate the reported estimate came
#'     from.}
#'
#' @seealso \code{\link{ivamse}} for the formula interface.
#'
#' @examples
#' set.seed(1)
#' N <- 400
#' z <- matrix(rnorm(N * 20), N, 20)
#' colnames(z) <- paste0("z", 1:20)
#' w <- cbind(control = rnorm(N))
#' v <- rnorm(N)
#' x <- as.vector(z[, 1] + 0.5 * z[, 2] + w + v)
#' y <- as.vector(2 * x + w + 0.7 * v + rnorm(N))
#'
#' ## the default: a grid of 25 penalties, chosen by the criterion
#' m <- ivamse_fit(y, x, z, w = w)
#' coef(m)
#' m$selected$d
#'
#' ## two instrument sets, each with its own grid length
#' m2 <- ivamse_fit(y, x, list(small = z[, 1:5], wide = z), w = w,
#'                  nlambda = list(small = 10, wide = 30))
#' table(m2$candidates$dictionary)
#' m2$selected$dictionary
#'
#' ## multiples of the plug-in penalty instead of an automatic grid
#' m3 <- ivamse_fit(y, x, z, w = w, penalty = "bcch",
#'                  kappa = 2^seq(-4, 2, by = 0.5))
#' m3$candidates[, c("kappa", "d", "criterion")]
#'
#' ## the two first-stage-only comparators
#' sapply(c("amse", "cv", "bcch"), function(rule)
#'   ivamse_fit(y, x, z, w = w, select = rule)$selected$d)
#'
#' @export
ivamse_fit <- function(y, x, z, w = NULL, intercept = TRUE,
                       lambda = NULL, nlambda = 25, lambda.min.ratio = NULL,
                       penalty = c("path", "bcch"),
                       kappa = 2^seq(-4, 2, by = 0.5),
                       c_lambda = 1.1, sigma_v = NULL,
                       select = c("amse", "cv", "bcch"),
                       pilot = "auto", screen = FALSE, cluster = NULL,
                       nfolds = 10, foldid = NULL, cv_group = NULL,
                       cv_s = c("lambda.min", "lambda.1se"),
                       x_name = "x") {
  ## record which grid arguments the caller supplied, before match.arg() and
  ## the defaults make everything look supplied
  supplied <- list(kappa = !missing(kappa), nlambda = !missing(nlambda),
                   lambda.min.ratio = !missing(lambda.min.ratio),
                   cv_s = !missing(cv_s))

  penalty <- match.arg(penalty)
  select <- match.arg(select)
  cv_s <- match.arg(cv_s)

  ## ---- Validate ------------------------------------------------------------
  if (is.matrix(x) && ncol(x) > 1L) {
    stop("only one endogenous regressor is supported, but 'x' has ", ncol(x),
         " columns")
  }
  if (!is.numeric(y) || !is.numeric(x)) {
    stop("'y' and 'x' must be numeric")
  }

  y <- as.vector(y)
  x <- as.vector(x)
  z <- as_dictionary_list(z)
  N <- length(y)

  if (length(x) != N) stop("'x' and 'y' must have the same length")
  if (any(vapply(z, nrow, integer(1)) != N)) {
    stop("every dictionary must have as many rows as 'y'")
  }
  if (!is.null(w) && nrow(as.matrix(w)) != N) {
    stop("'w' must have as many rows as 'y'")
  }
  if (!all(is.finite(y)) || !all(is.finite(x))) {
    stop("'y' and 'x' must be finite")
  }
  if (!is.logical(screen) || length(screen) != 1L || is.na(screen)) {
    stop("'screen' must be TRUE or FALSE")
  }
  if (any(vapply(z, function(m) !all(is.finite(m)), logical(1)))) {
    stop("dictionaries must be finite; remove missing values first")
  }
  if (!is.null(w) && !all(is.finite(as.matrix(w)))) {
    stop("'w' must be finite; remove missing values first")
  }
  if (!is.null(cluster) && NROW(cluster) != N) {
    stop("'cluster' has ", NROW(cluster), " observations but the data has ", N)
  }
  warn_inert_arguments(penalty, lambda, select, supplied)

  ## ---- Prepare -------------------------------------------------------------
  model <- prepare_model(y, x, if (is.null(w)) NULL else as.matrix(w), intercept)
  prepared <- prepare_dictionaries(z, model)

  reference <- if (is.character(pilot) && !identical(pilot, "auto")) {
    if (!pilot %in% names(z)) {
      stop("'pilot' does not name a dictionary; available: ",
           paste(names(z), collapse = ", "))
    }
    pilot
  } else {
    names(z)[1L]
  }
  ## an integer pilot names a candidate but cannot move the reference
  ## dictionary, since sigma_v and the penalty grid are built from the reference
  ## before any candidate exists

  ## sigma_v is estimated from (x, Z) alone, never from the outcome
  sigma_v_hat <- residual_scale(prepared[[reference]]$z, model$x_res, model$n)
  if (is.null(sigma_v)) sigma_v <- sigma_v_hat

  ## ---- Fit every candidate -------------------------------------------------
  used_foldid <- NULL
  built <- build_candidates(prepared, model, penalty, lambda, kappa, nlambda,
                            lambda.min.ratio, sigma_v, c_lambda)
  candidates <- built$candidates
  Pi_hat <- built$Pi_hat
  J <- nrow(candidates)
  n <- model$n

  if (J > 100) {
    message("scoring ", J, " candidates; consider a smaller 'nlambda'")
  }

  ## ---- Which candidates define an estimator at all -------------------------
  candidates$empty <- vapply(seq_len(J), function(j) {
    is_empty_fit(Pi_hat[, j], candidates$fitted_moment[j], candidates$d[j],
                 n, model$x_res)
  }, logical(1))

  ## an empty fit is treated as having no denominator even when numerical dust
  ## leaves h arithmetically nonzero. Selection uses this zeroed copy; the table
  ## keeps the raw value.
  h <- candidates$h
  h[candidates$empty] <- 0

  eligible <- eligible_candidates(h, n, screen)
  candidates$eligible <- seq_len(J) %in% eligible
  if (!length(eligible)) {
    stop("no candidate has a usable first stage; every lasso fit is empty",
         if (screen) " or fails the denominator screen" else "")
  }

  beta <- candidate_beta(Pi_hat, model$y_res, h, n)

  ## ---- Pilot and the two error moments -------------------------------------
  pilot_choice <- if (is.numeric(pilot)) {
    if (length(pilot) != 1L || pilot < 1L || pilot > J) {
      stop("'pilot' is candidate ", pilot, " but there are only ", J,
           " candidates")
    }
    list(index = as.integer(pilot), screened = NA)
  } else {
    choose_pilot(h, candidates$dictionary, reference, n)
  }
  pilot_index <- pilot_choice$index
  if (is.na(pilot_index) || !is.finite(beta[pilot_index])) {
    stop("the pilot candidate has no usable first stage; widen the penalty ",
         "grid or set 'pilot'")
  }

  moments <- nuisance_moments(model$y_res, model$x_res, beta[pilot_index], n)

  ## ---- Score and select ----------------------------------------------------
  candidates$criterion <- feasible_criterion(
    candidates$fitted_moment, candidates$d, h,
    moments$sigma_eps2, moments$sigma_epsv, n
  )
  ## scores are kept for every candidate; eligibility is applied separately
  candidates$cvm <- NA_real_
  candidates$cvsd <- NA_real_

  if (select == "cv") {
    if (is.null(foldid)) foldid <- assign_folds(N, nfolds, cv_group)
    nfolds <- length(unique(foldid))
    used_foldid <- foldid
    for (name in names(prepared)) {
      at <- candidates$dictionary == name
      ## training fits use the full-sample penalty. The conversion uses the
      ## full-sample n and N, since the controls were partialled out before the
      ## folds were formed.
      scores <- cv_scores(prepared[[name]]$z, model$x_res,
                          candidates$lambda[at] * n / N, foldid, nfolds)
      candidates$cvm[at] <- scores$cvm
      candidates$cvsd[at] <- scores$cvsd
    }
  }

  chosen <- choose_index(select, candidates, eligible, prepared, model,
                         reference, c_lambda, cv_s)

  resolution <- resolve_estimate(chosen, beta, h, pilot_index)
  if (resolution$unresolved) {
    stop("no usable candidate and no usable pilot")
  }
  used <- resolution$index

  ## ---- Estimate ------------------------------------------------------------
  estimate <- iv_estimate(y, x, model$w, model$w_qr, Pi_hat[, used], n, x_name)
  candidates$used <- seq_len(J) == used

  used_dictionary <- candidates$dictionary[used]
  first_stage <- built$fits[[used_dictionary]]$pi_hat[, candidates$column[used]]
  names(first_stage) <- colnames(prepared[[used_dictionary]]$z)

  out <- list(
    coefficients = estimate$coefficients,
    residuals = estimate$residuals,
    fitted.values = estimate$fitted.values,
    regressors = estimate$regressors,
    projected = estimate$projected,
    cov.unscaled = estimate$cov.unscaled,
    candidates = candidates,
    beta = beta,
    selected = list(
      index = chosen,
      estimation_index = used,
      fallback = resolution$fallback,
      reason = resolution$reason,
      dictionary = used_dictionary,
      lambda = candidates$lambda[used],
      kappa = candidates$kappa[used],
      d = candidates$d[used],
      h = candidates$h[used],
      fitted_moment = candidates$fitted_moment[used],
      implied_alpha = candidates$implied_alpha[used],
      criterion = candidates$criterion[used]
    ),
    instrument = as.vector(Pi_hat[, used]),
    first_stage = first_stage,
    scale = prepared[[used_dictionary]]$scale,
    sigma_eps2 = moments$sigma_eps2,
    sigma_epsv = moments$sigma_epsv,
    sigma_v = sigma_v,
    ## the reported correlation always uses the least-squares scale, never a
    ## sigma_v the user supplied to shape the penalty grid
    rho_hat = moments$sigma_epsv / sqrt(moments$sigma_eps2 * sigma_v_hat^2),
    ## joint F on the reference dictionary, and F on the fitted instrument used
    first_stage_F = first_stage_f(prepared[[reference]]$z, model$x_res, n),
    instrument_F = first_stage_f(cbind(Pi_hat[, used]), model$x_res, n),
    foldid = if (select == "cv") used_foldid else NULL,
    pilot = list(index = pilot_index, dictionary = reference,
                 beta = beta[pilot_index], screened = pilot_choice$screened),
    select = select,
    screen = screen,
    cluster = cluster,
    endogenous = x_name,
    controls = colnames(model$w),
    n = n, N = N, q = model$q, nobs = N,
    df.residual = N - sum(!is.na(estimate$coefficients))
  )
  class(out) <- "ivamse"
  out
}


## warn when an argument cannot have any effect given the others, rather than
## letting it be ignored silently
warn_inert_arguments <- function(penalty, lambda, select, supplied) {
  inert <- character(0)
  if (!is.null(lambda)) {
    if (supplied$kappa) inert <- c(inert, "kappa")
    if (supplied$nlambda) inert <- c(inert, "nlambda")
    if (supplied$lambda.min.ratio) inert <- c(inert, "lambda.min.ratio")
    if (length(inert)) {
      warning("'lambda' was supplied, so ", paste(sQuote(inert), collapse = ", "),
              " ", if (length(inert) == 1L) "is" else "are", " ignored",
              call. = FALSE)
      return(invisible(NULL))
    }
  } else if (penalty == "path" && supplied$kappa) {
    warning("'kappa' is ignored when penalty = \"path\"", call. = FALSE)
  } else if (penalty == "bcch") {
    if (supplied$nlambda) inert <- c(inert, "nlambda")
    if (supplied$lambda.min.ratio) inert <- c(inert, "lambda.min.ratio")
    if (length(inert)) {
      warning(paste(sQuote(inert), collapse = " and "), " ",
              if (length(inert) == 1L) "is" else "are",
              " ignored when penalty = \"bcch\"", call. = FALSE)
    }
  }
  if (select != "cv" && supplied$cv_s) {
    warning("'cv_s' is ignored unless select = \"cv\"", call. = FALSE)
  }
  invisible(NULL)
}
