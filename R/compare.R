## first-stage-only selection rules, for comparison with the AMSE criterion:
## cross-validation of the first stage, and the plug-in penalty of Belloni,
## Chen, Chernozhukov and Hansen. Both ignore the outcome equation.


## fold assignment, keeping whole groups together when 'group' is supplied
assign_folds <- function(N, nfolds, group) {
  if (is.null(group)) return(sample(rep_len(seq_len(nfolds), N)))
  levels <- unique(group)
  by_group <- sample(rep_len(seq_len(nfolds), length(levels)))
  by_group[match(group, levels)]
}


## K-fold cross-validated prediction error of the first stage, per candidate.
## Training fits use the full-sample penalties. 'cvm' is the mean over folds and
## 'cvsd' its standard error, as in cv.glmnet.
## controls are partialled out and columns normalized on the full sample, so a
## held-out fold enters the projection and the column scales.
cv_scores <- function(z, x, glmnet_lambda, foldid, nfolds) {
  path <- sort(glmnet_lambda, decreasing = TRUE)
  fold_error <- matrix(NA_real_, nfolds, length(glmnet_lambda))

  for (f in seq_len(nfolds)) {
    test <- foldid == f
    fit <- glmnet::glmnet(
      z[!test, , drop = FALSE], x[!test],
      alpha = 1, lambda = path, standardize = FALSE, intercept = FALSE,
      thresh = GLMNET_THRESHOLD, maxit = GLMNET_MAXIT
    )
    columns <- match_lambda(glmnet_lambda, fit$lambda)
    prediction <- z[test, , drop = FALSE] %*% as.matrix(fit$beta)[, columns, drop = FALSE]
    fold_error[f, ] <- colMeans((prediction - x[test])^2)
  }

  list(cvm = colMeans(fold_error),
       cvsd = apply(fold_error, 2L, stats::sd) / sqrt(nfolds))
}


## the cross-validation choice, in glmnet's terms: "lambda.min" is the largest
## penalty attaining the minimum error and "lambda.1se" the largest penalty
## within one standard error of it. Ties use the same relative tolerance as the
## AMSE rule.
select_by_cv <- function(cvm, cvsd, lambda, eligible) {
  eligible <- eligible[is.finite(cvm[eligible])]
  if (!length(eligible)) {
    return(list(lambda.min = NA_integer_, lambda.1se = NA_integer_))
  }

  smallest <- min(cvm[eligible])
  at_min <- eligible[cvm[eligible] <= smallest + 1e-12 * max(1, abs(smallest))]
  best <- at_min[which.max(lambda[at_min])]

  list(
    lambda.min = best,
    lambda.1se = {
      within <- eligible[cvm[eligible] <= cvm[best] + cvsd[best]]
      if (!length(within)) best else within[which.max(lambda[within])]
    }
  )
}


## plug-in penalty level, with sigma_v refined on the LASSO residual as in the
## original programs of Belloni, Chen, Chernozhukov and Hansen
bcch_plugin_level <- function(z, x, n, c_lambda = 1.1, iterations = 8L,
                              tolerance = 1e-6) {
  N <- nrow(z)
  ## penalty per unit of sigma_v; the returned level is this times the converged
  ## scale
  unit_lambda <- bcch_lambda(ncol(z), n, sigma_v = 1, c_lambda = c_lambda)
  sigma <- sqrt(sum(x^2) / n)

  for (iteration in seq_len(iterations)) {
    fit <- glmnet::glmnet(
      z, x, alpha = 1, lambda = unit_lambda * sigma * n / N,
      standardize = FALSE, intercept = FALSE,
      thresh = GLMNET_THRESHOLD, maxit = GLMNET_MAXIT
    )
    updated <- sqrt(sum((x - z %*% as.vector(fit$beta))^2) / n)
    if (!is.finite(updated) || updated <= 0) break
    if (abs(updated - sigma) <= tolerance * sigma) break
    sigma <- updated
  }

  unit_lambda * sigma
}
