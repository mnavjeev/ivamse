## The feasible AMSE criterion of Ma, Navjeevan and Salahub. Squared error
## decomposes as
##
##   n * (beta_hat_c - beta)^2 = C_0 + Q_c + remainder,
##
## with C_0 common to every candidate, so ranking candidates by approximate mean
## squared error means ranking E[Q_c]. That expectation is approximated by
##
##   S_c = sigma_eps2 * A_c / h_c^2
##         + sigma_epsv^2 * E[d_c^2 + d_c] / (n h_c^2),
##
## where A_c is the error in approximating the optimal instrument by the best
## rescaling of the fit and d_c is the effective dimension. The first term is
## first-stage approximation error, the second a many-instrument bias.
##
## A_c = E[Pi_hat_c^2] - h_c^2 / H, so replacing A_c by the observed fitted
## second moment adds sigma_eps2 / H, which is common to all candidates and
## cannot change the ranking. Replacing E[d_c^2 + d_c] by its realized value is
## justified by ratio consistency of the effective dimension, giving the
## feasible criterion below.
##
## The level of S_hat_c carries that common offset, so only differences and
## rankings across candidates are meaningful.


#' The Selection Criterion
#'
#' @description Computes the score that \code{\link{ivamse}} minimizes when
#' choosing among candidates. Each candidate is one instrument set paired with
#' one lasso penalty; the score estimates the mean squared error of the
#' coefficient on the endogenous regressor that the candidate would produce, so
#' smaller is better.
#'
#' \code{\link{ivamse}} computes this for every candidate and reports it in the
#' \code{criterion} column, so it is seldom called directly. It is exported so
#' that the score can be recomputed by hand, for instance with error moments
#' fixed at known values.
#'
#' @param fitted_moment second moment of the fitted first stage,
#'   \code{mean(Pi_hat^2)}. A numeric vector with one entry per candidate.
#' @param d effective dimension of each fit, roughly how many instruments it
#'   uses; see \code{\link{effective_dimension}}.
#' @param h IV denominator of each candidate, \code{mean(Pi_hat * x)}. A value
#'   near zero means a weak first stage and gives a large score.
#' @param sigma_eps2 estimated variance of the structural error. A single
#'   number, the same for every candidate; a fitted model reports it as
#'   \code{fit$sigma_eps2}.
#' @param sigma_epsv estimated covariance of the structural and first-stage
#'   errors, not its square. A single number; a fitted model reports it as
#'   \code{fit$sigma_epsv}.
#' @param n effective sample size, the number of observations less the rank of
#'   the controls.
#'
#' @details The score is
#' \deqn{\hat S_c = \frac{\hat\sigma_\varepsilon^2 E_n[\hat\Pi_{c,i}^2] +
#'   n^{-1}\hat\sigma_{\varepsilon v}^2 (d_c^2 + d_c)}{\hat h_c^2}}{%
#'   S_c = (sigma_eps2 * fitted_moment + sigma_epsv^2 * (d^2 + d) / n) / h^2}
#'
#' The first term in the numerator falls as the first stage fits better. The
#' second is a many-instrument bias, of the kind studied by Donald and Newey
#' (2001), which grows with the complexity of the fit. Its weight is the squared
#' covariance of the two errors, so the more endogenous the regressor, the more
#' heavily complexity is punished and the simpler the selected first stage. This
#' is the term that cross-validation, which scores the first stage alone,
#' leaves out.
#'
#' The score carries an offset that is the same for every candidate. Differences
#' and rankings within one fit are therefore meaningful, but the level is not an
#' estimate of any coefficient's mean squared error, and scores from different
#' fits cannot be compared.
#'
#' @return A numeric vector of scores, one per candidate, and \code{Inf} for any
#'   candidate whose first stage is degenerate. The chosen candidate is the one
#'   with the smallest value, so \code{which.min} recovers it.
#'
#' @references
#' Donald, S. G., and Newey, W. K. (2001). Choosing the number of instruments.
#' \emph{Econometrica} \bold{69}, 1161--1191.
#'
#' Ma, Y., Navjeevan, M., and Salahub, B. Choosing the dictionary and penalty
#' for IV-LASSO. Working paper.
#'
#' @seealso \code{\link{ivamse}}, \code{\link{effective_dimension}}.
#'
#' @examples
#' ## recomputing the scores of a fitted model from its own candidate table
#' set.seed(1)
#' N <- 300
#' z <- matrix(rnorm(N * 20), N, 20)
#' v <- rnorm(N)
#' x <- z[, 1] + 0.6 * z[, 2] + v
#' y <- 2 * x + 0.7 * v + rnorm(N)
#' m <- ivamse_fit(y, x, z)
#'
#' cand <- m$candidates
#' by_hand <- feasible_criterion(cand$fitted_moment, cand$d, cand$h,
#'                               m$sigma_eps2, m$sigma_epsv, m$n)
#' all.equal(by_hand, cand$criterion)
#' which.min(by_hand)
#'
#' ## a more endogenous regressor is met with a simpler first stage
#' strong <- feasible_criterion(cand$fitted_moment, cand$d, cand$h,
#'                              m$sigma_eps2, 0.9, m$n)
#' c(selected = cand$d[which.min(by_hand)],
#'   more_endogenous = cand$d[which.min(strong)])
#'
#' @export
feasible_criterion <- function(fitted_moment, d, h, sigma_eps2, sigma_epsv, n) {
  value <- (sigma_eps2 * fitted_moment + sigma_epsv^2 * (d^2 + d) / n) / h^2
  ## a candidate whose first stage is identically zero defines no IV estimator
  value[!is.finite(h) | h == 0 | !is.finite(value)] <- Inf
  value
}


## finite-sample screen on the IV denominator, at its root-n fluctuation scale
denominator_threshold <- function(n) 1 / sqrt(n)


## candidates that define an IV estimator, optionally after the screen
eligible_candidates <- function(h, n, screen) {
  usable <- is.finite(h) & h != 0
  if (screen) usable <- usable & abs(h) >= denominator_threshold(n)
  which(usable)
}


## minimizer, with ties broken toward the smallest candidate index. Candidates
## are ordered by dictionary in the order the user listed them and, within a
## dictionary, by increasing penalty, so a tie resolves toward the earlier
## dictionary and the lighter penalty.
select_candidate <- function(criterion, eligible) {
  eligible <- eligible[is.finite(criterion[eligible])]
  if (!length(eligible)) return(NA_integer_)
  best <- min(criterion[eligible])
  tolerance <- 1e-12 * max(1, abs(best))
  min(eligible[abs(criterion[eligible] - best) <= tolerance])
}


## Pilot candidate, whose structural estimate supplies the two error moments:
## the strongest observed first stage within a reference dictionary fixed in
## advance. Restricting to one dictionary keeps a large interacted dictionary
## from contaminating the preliminary estimate and makes the pilot a function of
## (x, Z) alone.
##
## the pilot is always screened, whatever 'screen' says, which governs selection
## only; if nothing survives the unscreened set is used, and that is recorded
choose_pilot <- function(h, dictionary, reference, n) {
  usable <- which(dictionary == reference & is.finite(h) & h != 0)
  if (!length(usable)) return(list(index = NA_integer_, screened = NA))

  screened <- usable[abs(h[usable]) >= denominator_threshold(n)]
  candidates <- if (length(screened)) screened else usable

  strongest <- max(abs(h[candidates]))
  list(
    index = min(candidates[abs(h[candidates]) >=
                             strongest - 1e-12 * max(1, strongest)]),
    screened = length(screened) > 0L
  )
}


## the two error moments, from the pilot's structural residual. The second uses
## x itself rather than a first-stage residual, so no estimate of v is needed.
nuisance_moments <- function(y, x, beta_pilot, n) {
  residual <- y - beta_pilot * x
  list(sigma_eps2 = sum(residual^2) / n,
       sigma_epsv = sum(residual * x) / n)
}


## which candidate is used for estimation: the selected one, or the pilot if the
## selected one gives no usable estimate. The selected index is kept either way;
## if the pilot is also unusable the fit is unresolved.
resolve_estimate <- function(selected, beta, h, pilot) {
  usable <- !is.na(selected) && is.finite(h[selected]) && h[selected] != 0 &&
    is.finite(beta[selected])
  if (usable) {
    return(list(index = selected, fallback = FALSE, unresolved = FALSE,
                reason = "none"))
  }

  reason <- if (is.na(selected)) {
    "selection_unresolved"
  } else if (!is.finite(h[selected])) {
    "nonfinite_denominator"
  } else if (h[selected] == 0) {
    "zero_denominator"
  } else {
    "nonfinite_estimate"
  }

  if (is.na(pilot) || !is.finite(beta[pilot])) {
    list(index = NA_integer_, fallback = FALSE, unresolved = TRUE, reason = reason)
  } else {
    list(index = pilot, fallback = TRUE, unresolved = FALSE, reason = reason)
  }
}
