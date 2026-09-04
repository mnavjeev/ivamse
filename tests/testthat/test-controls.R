test_that("controls reduce the effective sample size by their rank", {
  design <- make_design()
  fit <- ivamse_fit(design$y, design$x, design$z, w = design$w)

  expect_equal(fit$q, 3L)                 # intercept plus two controls
  expect_equal(fit$n, fit$N - fit$q)
  expect_equal(fit$N, length(design$y))
})

test_that("partialling out matches rotating onto the orthogonal complement", {
  # The corollary on controls says the model with controls reduces to the model
  # without them, at sample size n = N - q, after rotating onto an orthonormal
  # basis of the complement of the controls.  Rotating gives n-vectors, which is
  # tidier algebraically; the package instead keeps N-vectors so that robust and
  # clustered variances stay meaningful.  The two must agree exactly.
  design <- make_design()
  w <- cbind(1, design$w)
  N <- length(design$y)
  q <- qr(w)$rank
  n <- N - q

  partialled <- ivamse_fit(design$y, design$x, design$z, w = design$w)
  lambda <- partialled$candidates$lambda

  U <- qr.Q(qr(w), complete = TRUE)[, (q + 1L):N, drop = FALSE]
  rotated <- ivamse_fit(
    as.vector(crossprod(U, design$y)), as.vector(crossprod(U, design$x)),
    crossprod(U, design$z), w = NULL, intercept = FALSE, lambda = lambda
  )

  expect_equal(rotated$n, n)
  expect_equal(partialled$n, n)
  expect_equal(partialled$candidates$h, rotated$candidates$h, tolerance = 1e-10)
  expect_equal(partialled$candidates$fitted_moment,
               rotated$candidates$fitted_moment, tolerance = 1e-10)
  expect_equal(partialled$candidates$d, rotated$candidates$d)
  expect_equal(partialled$candidates$criterion, rotated$candidates$criterion,
               tolerance = 1e-8)
  expect_equal(unname(coef(partialled)[1]), unname(coef(rotated)[1]),
               tolerance = 1e-10)
})

test_that("dictionary columns inside the span of the controls are dropped", {
  design <- make_design()
  # Append a column that is an exact linear combination of the controls.
  z <- cbind(design$z, bad = design$w %*% c(2, -3))
  fit <- ivamse_fit(design$y, design$x, z, w = design$w)

  expect_equal(unique(fit$candidates$p), ncol(design$z))
  expect_false("bad" %in% names(fit$first_stage))
})

test_that("a dictionary entirely absorbed by the controls is an error", {
  design <- make_design()
  z <- cbind(a = design$w[, 1], b = design$w[, 2])
  expect_error(ivamse_fit(design$y, design$x, z, w = design$w),
               "span of the controls")
})

test_that("dictionary columns are normalized to unit second moment", {
  design <- make_design()
  w_qr <- ivamse:::control_qr(cbind(1, design$w))
  n <- length(design$y) - w_qr$rank
  prepared <- ivamse:::prepare_dictionary(design$z * 1000, w_qr, n)

  expect_equal(colSums(prepared$z^2) / n, rep(1, ncol(prepared$z)),
               ignore_attr = TRUE)
})
