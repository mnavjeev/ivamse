## forming the structural estimate from a chosen first stage. The fitted
## instrument g = M_W Z pi_hat is the single excluded instrument in a
## just-identified IV regression of y on (x, W) instrumented by (g, W); by
## Frisch-Waugh-Lovell its coefficient on x is E_n[g_i y_i] / E_n[g_i x_i],
## written as a regression so that the standard variance methods apply.
##
## The instrument is the lasso fit itself, not a refit on the selected support.


## structural coefficient implied by each candidate's fitted instrument
candidate_beta <- function(Pi_hat, y, h, n) {
  beta <- (crossprod(Pi_hat, y)[, 1L] / n) / h
  beta[!is.finite(h) | h == 0] <- NA_real_
  beta
}


## coefficients, residuals and the matrices inference needs, given the selected
## fitted instrument g, which is already orthogonal to the controls
iv_estimate <- function(y, x, w, w_qr, g, n, x_name = "x") {
  h <- sum(g * x) / n
  fitted_moment <- sum(g^2) / n
  beta <- (sum(g * y) / n) / h

  gamma <- if (is.null(w_qr)) numeric(0) else qr.coef(w_qr, y - beta * x)
  coefficients <- c(beta, gamma)
  names(coefficients) <- c(x_name, colnames(w))

  regressors <- cbind(x, w)
  colnames(regressors) <- names(coefficients)

  ## projected regressors P_G X with instruments G = (g, W). The controls are
  ## already in the span of G and g is orthogonal to them, so only the
  ## endogenous column changes.
  x_projected <- (x - partial_out(x, w_qr)) + g * h / fitted_moment
  projected <- cbind(x_projected, w)
  colnames(projected) <- names(coefficients)

  fitted <- as.vector(regressors %*% ifelse(is.na(coefficients), 0, coefficients))

  list(
    coefficients = coefficients,
    residuals = y - fitted,
    fitted.values = fitted,
    regressors = regressors,
    projected = projected,
    ## (X_hat' X_hat)^{-1}; its leading element is fitted_moment / (n h^2)
    cov.unscaled = chol2inv(chol(crossprod(projected))),
    h = h,
    fitted_moment = fitted_moment
  )
}
