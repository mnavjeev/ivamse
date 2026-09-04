## methods for fitted "ivamse" objects

#' @rdname ivamse
#' @param digits minimal number of significant digits for printing.
#' @export
print.ivamse <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")
  cat("First stage: dictionary \"", x$selected$dictionary,
      "\", lambda = ", format(x$selected$lambda, digits = digits),
      " (kappa = ", format(x$selected$kappa, digits = 2),
      "), effective dimension ", x$selected$d, "\n\n", sep = "")
  cat("Coefficients:\n")
  print.default(format(x$coefficients, digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\n")
  invisible(x)
}


## best candidate within each dictionary, ordered by criterion
dictionary_ranking <- function(candidates) {
  usable <- candidates[candidates$eligible & is.finite(candidates$criterion),
                       , drop = FALSE]
  if (!nrow(usable)) return(NULL)
  best <- do.call(rbind, lapply(
    split(usable, factor(usable$dictionary, levels = unique(usable$dictionary))),
    function(g) g[which.min(g$criterion), , drop = FALSE]
  ))
  best[order(best$criterion), c("dictionary", "p", "kappa", "d", "criterion",
                                "implied_alpha"), drop = FALSE]
}


#' @rdname ivamse
#' @export
summary.ivamse <- function(object, vcov. = NULL, df = NULL, ...) {
  estimate <- object$coefficients
  se <- sqrt(diag(resolve_vcov(object, vcov., ...)))
  if (is.null(df)) df <- Inf

  statistic <- estimate / se
  p <- if (is.finite(df)) {
    2 * stats::pt(abs(statistic), df, lower.tail = FALSE)
  } else {
    2 * stats::pnorm(abs(statistic), lower.tail = FALSE)
  }

  table <- cbind(estimate, se, statistic, p)
  dimnames(table) <- list(
    names(estimate),
    c("Estimate", "Std. Error",
      if (is.finite(df)) "t value" else "z value",
      if (is.finite(df)) "Pr(>|t|)" else "Pr(>|z|)")
  )

  ## relative gap between the best and the second-best score
  usable <- object$candidates$eligible & is.finite(object$candidates$criterion)
  scores <- sort(object$candidates$criterion[usable])
  margin <- if (length(scores) > 1L) {
    (scores[2L] - scores[1L]) / max(abs(scores[1L]), .Machine$double.eps)
  } else {
    NA_real_
  }

  structure(
    list(call = object$call, coefficients = table, df = df,
         selected = object$selected, select = object$select,
         screen = object$screen,
         candidates = nrow(object$candidates),
         eligible = sum(object$candidates$eligible),
         outside_assumption = sum(object$candidates$implied_alpha <= 0,
                                  na.rm = TRUE),
         ranking = dictionary_ranking(object$candidates),
         margin = margin,
         sigma = sqrt(sum(object$residuals^2) / object$df.residual),
         df.residual = object$df.residual,
         sigma_eps2 = object$sigma_eps2, sigma_epsv = object$sigma_epsv,
         rho_hat = object$rho_hat, pilot = object$pilot,
         first_stage_F = object$first_stage_F,
         instrument_F = object$instrument_F,
         cluster_label = object$cluster_label,
         n_clusters = if (is.null(object$cluster)) NA_integer_ else
           length(unique(if (is.data.frame(object$cluster))
             do.call(paste, object$cluster) else object$cluster)),
         n = object$n, N = object$N, q = object$q),
    class = "summary.ivamse"
  )
}


#' @rdname ivamse
#' @param signif.stars show "significance stars" in summary output?
#' @export
print.summary.ivamse <- function(x, digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
  cat("\nCall:\n", paste(deparse(x$call), collapse = "\n"), "\n\n", sep = "")

  cat("Coefficients:\n")
  stats::printCoefmat(x$coefficients, digits = digits,
                      signif.stars = signif.stars, na.print = "NA", ...)

  cat("\nResidual standard error: ", format(x$sigma, digits = digits),
      " on ", x$df.residual, " degrees of freedom", sep = "")
  if (!is.na(x$n_clusters)) {
    cat("\nStandard errors clustered",
        if (!is.null(x$cluster_label)) paste0(" by ", x$cluster_label) else "",
        ": ", x$n_clusters, " clusters", sep = "")
  }
  cat("\n")

  rule <- switch(x$select,
    amse = "minimum AMSE criterion",
    cv   = "cross-validated first-stage fit",
    bcch = "plug-in penalty"
  )
  cat("\nFirst-stage selection by ", rule, ":\n", sep = "")
  cat("  dictionary:          ", x$selected$dictionary, "\n", sep = "")
  cat("  penalty:             ", format(x$selected$lambda, digits = digits),
      " (kappa = ", format(x$selected$kappa, digits = 3), ")\n", sep = "")
  cat("  effective dimension: ", x$selected$d, "\n", sep = "")
  cat("  IV denominator h:    ", format(x$selected$h, digits = digits), "\n",
      sep = "")
  cat("  first-stage F:       ",
      if (is.na(x$first_stage_F)) "NA" else format(x$first_stage_F, digits = 4),
      " (reference dictionary), ",
      if (is.na(x$instrument_F)) "NA" else format(x$instrument_F, digits = 4),
      " (fitted instrument)\n", sep = "")
  cat("  chosen from ", x$eligible, " eligible of ", x$candidates,
      " candidates\n", sep = "")

  if (!is.null(x$ranking) && nrow(x$ranking) > 1L) {
    cat("\nBest candidate in each dictionary:\n")
    display <- x$ranking
    display$kappa <- signif(display$kappa, 3)
    display$criterion <- signif(display$criterion, 6)
    display$implied_alpha <- signif(display$implied_alpha, 3)
    print(display, row.names = FALSE, digits = digits)
  }

  cat("\nError moments from pilot candidate ", x$pilot$index, " (dictionary \"",
      x$pilot$dictionary, "\"):\n", sep = "")
  cat("  sigma_eps^2: ", format(x$sigma_eps2, digits = digits),
      ",  sigma_epsv: ", format(x$sigma_epsv, digits = digits),
      ",  implied correlation: ", format(x$rho_hat, digits = 3), "\n", sep = "")

  cat("\n", x$N, " observations, ", x$q, " control column",
      if (x$q == 1L) "" else "s", ", effective sample size ", x$n, "\n", sep = "")

  if (x$selected$fallback) {
    cat("Fallback to the pilot candidate (", x$selected$reason, ")\n", sep = "")
  }
  if (!is.na(x$margin) && x$margin < 1e-3) {
    cat("Runner-up within ", format(100 * x$margin, digits = 2),
        "% of the selected criterion\n", sep = "")
  }
  if (!is.na(x$selected$implied_alpha) && x$selected$implied_alpha <= 0) {
    cat("Assumption 2(iii) not satisfied at the selected penalty (implied alpha ",
        format(x$selected$implied_alpha, digits = 3), "); ",
        x$outside_assumption, " of ", x$candidates, " candidates affected\n",
        sep = "")
  }
  if (x$selected$d > sqrt(x$n)) {
    cat("Effective dimension ", x$selected$d, " exceeds sqrt(n) = ",
        round(sqrt(x$n), 1), "\n", sep = "")
  }
  if (abs(x$selected$h) < 2 / sqrt(x$n)) {
    cat("Weak first stage: |h| = ", format(abs(x$selected$h), digits = digits),
        " < 2/sqrt(n); standard errors not robust to weak identification\n",
        sep = "")
  }
  cat("\n")
  invisible(x)
}


#' Which Candidate Was Selected
#'
#' @description Reports which candidate the selection rule chose, how the best
#' candidate in each instrument set compared, and which individual instruments
#' the chosen first stage uses. It answers the question a fitted model does not
#' show directly: what did the procedure actually pick, and was it a close call?
#'
#' @param object an object of class \code{"ivamse"}.
#' @param x an object of class \code{"selected.ivamse"}.
#' @param digits minimal number of significant digits for printing.
#' @param ... currently not used.
#'
#' @details The printed table has one row per instrument set, giving the best
#' candidate within each, so that the margin between them can be judged. Its
#' columns are
#' \describe{
#'   \item{\code{kappa}}{the penalty as a multiple of the plug-in penalty, so
#'     \code{1} is the level that rule would choose and smaller values keep more
#'     instruments.}
#'   \item{\code{d}}{roughly how many instruments that fit uses; see
#'     \code{\link{effective_dimension}}.}
#'   \item{\code{h}}{the IV denominator; a value near zero means a weak first
#'     stage.}
#'   \item{\code{criterion}}{the score, smaller being better. Only differences
#'     within one fit are meaningful.}
#'   \item{\code{implied_alpha}}{positive when the penalty satisfies the
#'     condition the theory requires; see \code{\link{implied_alpha}}.}
#' }
#'
#' The reported first-stage coefficients are on the instrument set after the
#' controls have been partialled out and each column rescaled to unit second
#' moment. They say which instruments the fit leans on and how heavily relative
#' to one another; they are not coefficients on the original variables and
#' should not be read as such.
#'
#' @return An object of class \code{"selected.ivamse"}: a list with the chosen
#'   candidate (\code{candidate}), the rule used (\code{rule}), whether
#'   estimation fell back to the pilot and why (\code{fallback},
#'   \code{reason}), the per-instrument-set comparison
#'   (\code{dictionaries}), and the nonzero first-stage coefficients
#'   (\code{instruments}), largest in absolute value first.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{plot.ivamse}},
#'   \code{\link{effective_dimension}}, \code{\link{implied_alpha}}.
#'
#' @examples
#' set.seed(1)
#' N <- 400
#' d <- data.frame(w = rnorm(N))
#' Z <- matrix(rnorm(N * 12), N, 12)
#' colnames(Z) <- paste0("z", 1:12)
#' d <- cbind(d, Z)
#' v <- rnorm(N)
#' d$x <- d$z1 + 0.6 * d$z2 + d$w + v
#' d$y <- 2 * d$x + d$w + 0.7 * v + rnorm(N)
#'
#' m <- ivamse(y ~ w | x | z1,
#'   dictionaries = list(small = ~ z1 + z2,
#'                       medium = ~ z1 + z2 + z3 + z4 + z5,
#'                       wide = ~ (z1 + z2 + z3) * (z4 + z5)),
#'   data = d)
#'
#' selected(m)
#'
#' ## the pieces are available individually
#' selected(m)$candidate$dictionary
#' names(selected(m)$instruments)
#'
#' @export
selected <- function(object, ...) UseMethod("selected")

#' @rdname selected
#' @export
selected.ivamse <- function(object, ...) {
  nonzero <- object$first_stage[object$first_stage != 0]
  structure(
    list(
      candidate = data.frame(
        dictionary = object$selected$dictionary,
        lambda = object$selected$lambda,
        kappa = object$selected$kappa,
        d = object$selected$d,
        h = object$selected$h,
        criterion = object$selected$criterion,
        implied_alpha = object$selected$implied_alpha,
        stringsAsFactors = FALSE
      ),
      rule = object$select,
      fallback = object$selected$fallback,
      reason = object$selected$reason,
      dictionaries = dictionary_ranking(object$candidates),
      instruments = nonzero[order(abs(nonzero), decreasing = TRUE)]
    ),
    class = "selected.ivamse"
  )
}

#' @rdname selected
#' @export
print.selected.ivamse <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  ...) {
  cat("Selected by ", x$rule, ":\n", sep = "")
  print(x$candidate, row.names = FALSE, digits = digits)
  if (x$fallback) {
    cat("\nFallback to the pilot (", x$reason, ")\n", sep = "")
  }
  if (!is.null(x$dictionaries) && nrow(x$dictionaries) > 1L) {
    cat("\nBest candidate in each dictionary:\n")
    print(x$dictionaries, row.names = FALSE, digits = digits)
  }
  cat("\nFirst stage uses ", length(x$instruments), " instrument",
      if (length(x$instruments) == 1L) "" else "s",
      " on the normalized dictionary:\n", sep = "")
  print(signif(x$instruments, digits))
  invisible(x)
}


#' @rdname ivamse
#' @param complete logical. If \code{TRUE}, the default, the returned
#'   coefficient vector includes elements for aliased regressors.
#' @export
coef.ivamse <- function(object, complete = TRUE, ...) {
  if (complete) object$coefficients else object$coefficients[!is.na(object$coefficients)]
}

#' @rdname ivamse
#' @export
fitted.ivamse <- function(object, ...) object$fitted.values

#' @rdname ivamse
#' @export
residuals.ivamse <- function(object, ...) object$residuals

#' @rdname ivamse
#' @export
nobs.ivamse <- function(object, ...) object$nobs

#' @rdname ivamse
#' @export
sigma.ivamse <- function(object, ...) {
  sqrt(sum(object$residuals^2) / object$df.residual)
}

#' @rdname ivamse
#' @export
formula.ivamse <- function(x, ...) x$formula

#' @rdname ivamse
#' @export
terms.ivamse <- function(x, ...) x$terms


#' @rdname ivamse
#' @param component character indicating \code{"regressors"},
#'   \code{"projected"} (the projected regressors used by the variance
#'   estimator), or \code{"instrument"}.
#' @export
model.matrix.ivamse <- function(object,
                                component = c("regressors", "projected",
                                              "instrument"),
                                ...) {
  component <- match.arg(component)
  switch(component,
    regressors = object$regressors,
    projected = object$projected,
    instrument = cbind(instrument = object$instrument)
  )
}


#' @rdname ivamse
#' @param newdata optionally, a data frame of new observations.
#' @param type character. Only \code{"response"} is supported.
#' @param se.fit logical. If \code{TRUE}, standard errors of the fitted
#'   values are returned alongside them.
#' @export
predict.ivamse <- function(object, newdata = NULL, type = "response",
                           se.fit = FALSE, vcov. = NULL, ...) {
  type <- match.arg(type, "response")
  X <- if (is.null(newdata)) {
    object$regressors
  } else {
    if (is.null(object$parts)) {
      stop("prediction from new data requires a model fitted with ivamse()")
    }
    mf <- stats::model.frame(stats::terms(object), newdata,
                             na.action = stats::na.pass, xlev = object$levels)
    split <- split_regressors(object$parts$f, object$parts$n_parts, mf,
                              object$contrasts)
    m <- cbind(split$endogenous, split$w)
    if ("(Intercept)" %in% names(object$coefficients)) {
      m <- cbind(split$endogenous, "(Intercept)" = 1, split$w)
    }
    m[, names(object$coefficients), drop = FALSE]
  }

  beta <- ifelse(is.na(object$coefficients), 0, object$coefficients)
  fit <- as.vector(X %*% beta)
  if (!se.fit) return(fit)

  v <- resolve_vcov(object, vcov., ...)
  list(fit = fit, se.fit = sqrt(rowSums((X %*% v) * X)))
}


## one panel of plot.ivamse: 'value' against log(lambda), one line per
## dictionary, with the candidate used marked. Arguments the caller supplied in
## '...' override the defaults, so they must be merged rather than passed
## alongside them, which would give plot() two values for the same formal.
criterion_panel <- function(candidates, value, keep, dictionaries, colours,
                            defaults, dots) {
  graphics::plot.new()
  do.call(graphics::plot,
          c(list(log(candidates$lambda[keep]), value[keep], type = "n"),
            utils::modifyList(defaults, dots)))
  for (i in seq_along(dictionaries)) {
    at <- candidates$dictionary == dictionaries[i] & keep
    graphics::lines(log(candidates$lambda[at]), value[at],
                    col = colours[i], lwd = 2)
  }
  at <- which(candidates$used)
  graphics::abline(v = log(candidates$lambda[at]), lty = 3, col = "grey40")
  graphics::points(log(candidates$lambda[at]), value[at], pch = 19, cex = 1.2)
}


#' Plot the Criterion along the Penalty Grid
#'
#' @description Draws two panels side by side. The left shows the selection
#' criterion against the penalty, with one line per instrument set and the
#' chosen candidate marked; the right shows how many instruments the fit uses at
#' each penalty. Together they show what the selection traded off.
#'
#' @param x an object of class \code{"ivamse"}.
#' @param main a character vector of length two giving the titles of the two
#'   panels.
#' @param ... further graphical arguments passed to
#'   \code{\link[graphics]{plot}} for both panels, such as \code{xlab},
#'   \code{ylab} or \code{cex.axis}. These override the method's own defaults.
#'   Line colours are chosen by the method, one per instrument set, and are not
#'   taken from \code{col}.
#'
#' @details The horizontal axis is the natural logarithm of the penalty, so the
#' left of each panel is a lightly penalized fit using many instruments and the
#' right a heavily penalized one using few. The dotted vertical line and the
#' filled point mark the candidate used for the reported estimate.
#'
#' Reading the two panels together: the criterion is usually flat over a range
#' of penalties, and a minimum reached at the very edge of the grid suggests
#' extending it with \code{nlambda} or \code{lambda.min.ratio}. Candidates whose
#' first stage is empty are omitted from the left panel, since their score is
#' infinite.
#'
#' @return \code{x}, invisibly. Called for its side effect.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{selected}}.
#'
#' @examples
#' set.seed(1)
#' N <- 400
#' d <- data.frame(w = rnorm(N))
#' Z <- matrix(rnorm(N * 12), N, 12)
#' colnames(Z) <- paste0("z", 1:12)
#' d <- cbind(d, Z)
#' v <- rnorm(N)
#' d$x <- d$z1 + 0.6 * d$z2 + d$w + v
#' d$y <- 2 * d$x + d$w + 0.7 * v + rnorm(N)
#'
#' m <- ivamse(y ~ w | x | z1,
#'   dictionaries = list(small = ~ z1 + z2,
#'                       wide = ~ (z1 + z2 + z3) * (z4 + z5)),
#'   data = d)
#' plot(m)
#'
#' @export
plot.ivamse <- function(x, main = c("Criterion", "Complexity"), ...) {
  candidates <- x$candidates
  dictionaries <- unique(candidates$dictionary)
  colours <- grDevices::hcl.colors(max(length(dictionaries), 2L), "Dark 3")
  main <- rep_len(main, 2L)

  old <- graphics::par(mfrow = c(1L, 2L), mar = c(4.2, 4.2, 2.5, 1))
  on.exit(graphics::par(old))

  dots <- list(...)
  finite <- is.finite(candidates$criterion)
  criterion_panel(candidates, candidates$criterion, finite, dictionaries, colours,
                  list(xlab = expression(log(lambda)), ylab = "AMSE criterion",
                       main = main[1L]), dots)
  if (length(dictionaries) > 1L) {
    graphics::legend("topright", legend = dictionaries, col = colours,
                     lwd = 2, bty = "n", cex = 0.85)
  }
  criterion_panel(candidates, candidates$d, rep(TRUE, nrow(candidates)),
                  dictionaries, colours,
                  list(xlab = expression(log(lambda)),
                       ylab = "effective dimension", main = main[2L]), dots)
  invisible(x)
}


#' @rdname ivamse
#' @export
deviance.ivamse <- function(object, ...) sum(object$residuals^2)


#' @rdname ivamse
#' @param formula. changes to the model formula, as in
#'   \code{\link[stats]{update}}. Multi-part formulas are updated part by part,
#'   so \code{. ~ . | . | . + z9} adds an instrument.
#' @param evaluate logical. If \code{TRUE}, the default, the updated call is
#'   evaluated; otherwise the call itself is returned.
#' @export
update.ivamse <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call
  if (is.null(call)) stop("need an object with call component")

  ## the model formula has two or three parts, so it must be updated through
  ## Formula's method rather than the default one, which would collapse them
  if (!missing(formula.)) {
    call$formula <- stats::formula(
      stats::update(Formula::as.Formula(stats::formula(object)), formula.)
    )
  }

  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) call <- as.call(c(as.list(call), extras[!existing]))
  }

  if (evaluate) eval(call, parent.frame()) else call
}
