#' @description Instrumental-variables regression in which the first stage is
#' fitted by the lasso, and both the set of instruments and the amount of
#' penalization are chosen from the data. The main fitting function is
#' \code{\link{ivamse}}, which takes a model formula in the manner of
#' \code{\link[ivreg]{ivreg}}; \code{\link{ivamse_fit}} is the same estimator
#' taking design matrices. See \code{vignette("ivamse")} for a worked example.
#'
#' @details
#' Suppose an endogenous regressor has many available instruments, or a few
#' instruments that could be entered as polynomials and interactions. Using all
#' of them makes the first stage fit well but biases the second stage; using few
#' of them wastes information. The lasso offers a middle course, but it leaves
#' two choices open, and standard asymptotic theory settles neither: any
#' reasonable first stage gives the same limiting distribution. This package
#' scores the available choices by an estimate of the mean squared error each
#' would produce, and keeps the best one.
#'
#' @section Terminology:
#' The following terms are used throughout the help pages.
#'
#' \describe{
#'   \item{instrument set (a \emph{dictionary})}{one set of columns offered to
#'     the first stage: the excluded instruments as they stand, or
#'     transformations of them such as polynomials and interactions. A richer
#'     set gives the first stage more shapes to fit, at the cost of more columns
#'     to fit them with. The \code{dictionaries} argument takes several.}
#'   \item{candidate}{one instrument set paired with one lasso penalty level.
#'     \code{\link{ivamse}} fits the first stage once per candidate, scores them
#'     all, and keeps one.}
#'   \item{the criterion}{the score being minimized: an estimate of the mean
#'     squared error of the coefficient on the endogenous regressor that each
#'     candidate would produce. Smaller is better. Every score carries the same
#'     unknown constant, so scores can be compared within one fit but are not
#'     estimates of any one coefficient's mean squared error, and scores from
#'     different fits are not comparable. See \code{\link{feasible_criterion}}.}
#'   \item{effective dimension (\code{d})}{roughly, how many instruments the
#'     lasso fit is really using. When no instrument can be written as a
#'     combination of the others it is the number of nonzero first-stage
#'     coefficients. See \code{\link{effective_dimension}}.}
#'   \item{the IV denominator (\code{h})}{the sample average of the fitted first
#'     stage times the endogenous regressor, which is the denominator of the IV
#'     estimate. A value near zero means a weak first stage.}
#'   \item{pilot candidate}{one preliminary first stage, fitted before the
#'     candidates are compared, used to estimate the two error moments the
#'     criterion needs. Every candidate's score depends on it.}
#'   \item{reference instrument set}{the set the pilot is chosen within. It is
#'     the first one listed unless \code{pilot} names another.}
#'   \item{plug-in penalty}{a penalty computed directly from the number of
#'     instruments, the sample size and the first-stage error scale, rather than
#'     by cross-validation, and set just high enough that the lasso discards
#'     instruments whose apparent correlation with the endogenous regressor is
#'     only sampling noise. Due to Belloni, Chen, Chernozhukov and Hansen (2012);
#'     see \code{\link{bcch_lambda}}.}
#'   \item{\code{kappa}}{a candidate's penalty as a multiple of the plug-in
#'     penalty, so \code{kappa = 1} is the level that rule would choose and
#'     smaller values keep more instruments.}
#'   \item{effective sample size (\code{n})}{the number of observations less the
#'     number of linearly independent controls, counting the intercept. This,
#'     not \code{nrow(data)}, is what the penalties are scaled by.}
#' }
#'
#' @section Notation:
#' Object names in the source follow Ma, Navjeevan and Salahub.
#'
#' \tabular{lll}{
#'   \strong{Paper} \tab \strong{Code} \tab \strong{Meaning} \cr
#'   \eqn{y}{y} \tab \code{y} \tab outcome \cr
#'   \eqn{x}{x} \tab \code{x} \tab the single endogenous regressor \cr
#'   \eqn{W}{W} \tab \code{w} \tab included exogenous controls, with the intercept \cr
#'   \eqn{Z_c}{Z_c} \tab \code{z} \tab one instrument set \cr
#'   \eqn{N}{N} \tab \code{N} \tab number of observations \cr
#'   \eqn{q}{q} \tab \code{q} \tab rank of the controls \cr
#'   \eqn{n = N - q}{n = N - q} \tab \code{n} \tab effective sample size \cr
#'   \eqn{p_c}{p_c} \tab \code{p} \tab number of columns in an instrument set \cr
#'   \eqn{\lambda_c}{lambda_c} \tab \code{lambda} \tab penalty level \cr
#'   \eqn{\kappa_c}{kappa_c} \tab \code{kappa} \tab penalty as a multiple of the plug-in level \cr
#'   \eqn{\hat\pi_c}{pi_hat} \tab \code{pi_hat} \tab first-stage lasso coefficients \cr
#'   \eqn{\hat\Pi_c}{Pi_hat} \tab \code{Pi_hat} \tab fitted instrument \eqn{Z_c \hat\pi_c}{Z_c pi_hat} \cr
#'   \eqn{\hat h_c}{h} \tab \code{h} \tab \eqn{E_n[\hat\Pi_{c,i} x_i]}{mean(Pi_hat * x)}, the IV denominator \cr
#'   \eqn{E_n[\hat\Pi_{c,i}^2]}{mean(Pi_hat^2)} \tab \code{fitted_moment} \tab second moment of the fitted instrument \cr
#'   \eqn{d_c}{d} \tab \code{d} \tab effective dimension of the lasso fit \cr
#'   \eqn{\hat\sigma_\varepsilon^2}{sigma_eps2} \tab \code{sigma_eps2} \tab structural error variance \cr
#'   \eqn{\hat\sigma_{\varepsilon v}}{sigma_epsv} \tab \code{sigma_epsv} \tab covariance of the two errors \cr
#'   \eqn{\hat S_c}{S_c} \tab \code{criterion} \tab the criterion \cr
#'   \eqn{\hat c}{c_hat} \tab \code{selected} \tab index of the chosen candidate
#' }
#'
#' @section Scope:
#' The model may contain only one endogenous regressor. The criterion is derived
#' under homoscedastic Gaussian errors; simulations with heavier-tailed errors
#' leave the ranking of candidates essentially unchanged. Clustering and
#' heteroscedasticity affect the reported standard errors only, never which
#' candidate is selected, and the standard errors treat the selected first stage
#' as fixed.
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
#' @seealso \code{\link{ivamse}}, \code{\link{ivamse_fit}}.
#'
#' @examples
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
#' m <- ivamse(y ~ w | x | z1 + z2 + z3 + z4 + z5, data = d)
#' coef(m)["x"]
#'
#' @importFrom stats coef confint fitted formula hatvalues model.matrix nobs
#' @importFrom stats predict residuals sigma terms vcov
#'
#' @aliases ivamse-package
#' @keywords internal
"_PACKAGE"
