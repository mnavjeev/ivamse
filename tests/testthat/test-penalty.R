test_that("bcch_lambda reproduces the published plug-in formula", {
  p <- 50; n <- 500; c_lambda <- 1.1
  gamma <- 0.1 / log(max(p, n))
  expect_equal(
    bcch_lambda(p, n, sigma_v = 2, c_lambda = c_lambda),
    c_lambda * 2 * qnorm(1 - gamma / (2 * p)) / sqrt(n)
  )
})

test_that("the penalty scales with sigma_v and shrinks with n", {
  expect_equal(bcch_lambda(50, 500, 3), 3 * bcch_lambda(50, 500, 1))
  expect_lt(bcch_lambda(50, 5000), bcch_lambda(50, 500))
})

test_that("implied_alpha is positive exactly when the penalty dominates", {
  p <- 50; n <- 500
  plugin <- bcch_lambda(p, n)
  expect_gt(implied_alpha(4 * plugin, p, n), 0)
  expect_lt(implied_alpha(plugin / 64, p, n), 0)
})

test_that("lambda_max is the boundary of the null fit", {
  design <- make_design(N = 300, p = 15)
  n <- 300
  z <- sweep(design$z, 2, sqrt(colSums(design$z^2) / n), "/")
  top <- ivamse:::lambda_max(z, design$x, n)

  fit <- ivamse:::fit_dictionary(z, design$x, c(top * 0.9, top), n, qr(z)$rank)
  # At lambda_max the fit is numerically zero; just below it is not.
  expect_true(ivamse:::is_empty_fit(fit$Pi_hat[, 2], fit$fitted_moment[2],
                                    fit$d[2], n, design$x))
  expect_false(ivamse:::is_empty_fit(fit$Pi_hat[, 1], fit$fitted_moment[1],
                                     fit$d[1], n, design$x))
})

test_that("the automatic grid is increasing and spans the requested range", {
  design <- make_design(N = 300, p = 15)
  n <- 300
  z <- sweep(design$z, 2, sqrt(colSums(design$z^2) / n), "/")
  grid <- ivamse:::lambda_path(z, design$x, n, nlambda = 10,
                               lambda.min.ratio = 1e-3)

  expect_length(grid, 10L)
  expect_false(is.unsorted(grid, strictly = TRUE))
  expect_equal(grid[10] / grid[1], 1000, tolerance = 1e-8)
  expect_equal(grid[10], ivamse:::lambda_max(z, design$x, n))
})
