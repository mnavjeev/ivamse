test_that("effective dimension counts the support when the fit is unique", {
  design <- make_design(N = 300, p = 15)
  n <- 300
  z <- sweep(design$z, 2, sqrt(colSums(design$z^2) / n), "/")
  lambda <- 0.1
  fit <- glmnet::glmnet(z, design$x, alpha = 1, lambda = lambda,
                        standardize = FALSE, intercept = FALSE, thresh = 1e-13)
  pi_hat <- as.vector(fit$beta)
  residual <- design$x - z %*% pi_hat

  expect_equal(
    effective_dimension(z, residual, lambda, n, pi_hat, dictionary_rank = ncol(z)),
    sum(pi_hat != 0)
  )
})

test_that("effective dimension is a rank when the dictionary is collinear", {
  design <- make_design(N = 300, p = 8)
  n <- 300
  z <- cbind(design$z, design$z[, 1:3])          # exact duplicate columns
  z <- sweep(z, 2, sqrt(colSums(z^2) / n), "/")
  lambda <- 0.05
  fit <- glmnet::glmnet(z, design$x, alpha = 1, lambda = lambda,
                        standardize = FALSE, intercept = FALSE, thresh = 1e-13)
  pi_hat <- as.vector(fit$beta)
  residual <- design$x - z %*% pi_hat

  d <- effective_dimension(z, residual, lambda, n, pi_hat,
                           dictionary_rank = qr(z)$rank)
  expect_lte(d, qr(z)$rank)
  expect_lte(d, 8L)          # cannot exceed the rank of the duplicated block
  expect_gt(d, 0L)
})

test_that("effective dimension is zero for a null fit", {
  z <- matrix(rnorm(100 * 5), 100, 5)
  z <- sweep(z, 2, sqrt(colSums(z^2) / 100), "/")
  expect_equal(effective_dimension(z, rep(0, 100), 1, 100, rep(0, 5),
                                   dictionary_rank = 4L), 0L)
})

test_that("a mismatched coefficient vector is rejected", {
  z <- matrix(rnorm(50 * 4), 50, 4)
  expect_error(effective_dimension(z, rnorm(50), 0.1, 50, pi_hat = rep(0, 3)),
               "does not match")
})
