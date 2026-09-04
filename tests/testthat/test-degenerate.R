test_that("a dictionary unrelated to x still produces an estimate or a clear error", {
  set.seed(9)
  N <- 200
  z <- matrix(rnorm(N * 10), N, 10)
  x <- rnorm(N)                       # no first stage at all
  y <- 2 * x + rnorm(N)

  # Either it estimates from a weak fit, or it says plainly that it cannot.
  result <- tryCatch(ivamse_fit(y, x, z), error = function(e) e)
  if (inherits(result, "error")) {
    expect_match(conditionMessage(result), "usable first stage|pilot")
  } else {
    expect_true(is.finite(coef(result)[1]))
  }
})

test_that("a wide dictionary with p greater than n is handled", {
  set.seed(10)
  N <- 80
  z <- matrix(rnorm(N * 200), N, 200)
  v <- rnorm(N)
  x <- z[, 1] + 0.8 * z[, 2] + v
  y <- 2 * x + 0.7 * v + rnorm(N)

  fit <- ivamse_fit(y, x, z)
  expect_true(is.finite(coef(fit)[1]))
  expect_lte(max(fit$candidates$d), fit$n)
})

test_that("a rank-deficient dictionary is scored by rank, not column count", {
  design <- make_design(N = 300, p = 8)
  z <- cbind(design$z, design$z[, 1:4])       # duplicated columns
  fit <- ivamse_fit(design$y, design$x, z, w = design$w)

  expect_true(is.finite(coef(fit)[1]))
  expect_lte(max(fit$candidates$d), 8L)
})

test_that("empty fits are marked and excluded from selection", {
  design <- make_design()
  fit <- ivamse_fit(design$y, design$x, design$z, w = design$w)

  # The largest penalty on the grid gives a null fit by construction.
  expect_true(any(fit$candidates$empty))
  expect_true(all(fit$candidates$criterion[fit$candidates$empty] == Inf))
  expect_false(fit$candidates$used[which(fit$candidates$empty)[1]])
})

test_that("the screen restricts the eligible set", {
  design <- make_design()
  loose <- ivamse_fit(design$y, design$x, design$z, w = design$w, screen = FALSE)
  tight <- ivamse_fit(design$y, design$x, design$z, w = design$w, screen = TRUE)

  expect_gte(sum(loose$candidates$eligible), sum(tight$candidates$eligible))
  expect_true(all(abs(tight$candidates$h[tight$candidates$eligible]) >=
                    1 / sqrt(tight$n)))
})

test_that("mismatched inputs give informative errors", {
  design <- make_design()
  expect_error(ivamse_fit(design$y, design$x[-1], design$z), "same length")
  expect_error(ivamse_fit(design$y, design$x, design$z[-1, ]), "as many rows")
  expect_error(ivamse_fit(design$y, design$x, list(design$z)), "must be named")
  expect_error(ivamse_fit(design$y, cbind(design$x, design$x), design$z),
               "one endogenous regressor")
})

test_that("a missing value is caught rather than propagated", {
  design <- make_design()
  y <- design$y; y[1] <- NA
  expect_error(ivamse_fit(y, design$x, design$z), "finite")
})

test_that("an instrument that is the endogenous regressor is refused", {
  design <- make_design()
  z <- cbind(design$z, x = design$x)
  expect_error(ivamse_fit(design$y, design$x, z, w = design$w),
               "collinear with the endogenous")

  # This is what `.` on the instrument side of a formula does.
  d <- design$data
  expect_error(ivamse(y ~ w1 + w2 | x | ., data = d),
               "collinear with the endogenous")
})

test_that("non-finite and malformed inputs are refused up front", {
  design <- make_design()
  z <- design$z; z[1, 1] <- NA
  expect_error(ivamse_fit(design$y, design$x, z), "must be finite")

  w <- design$w; w[1, 1] <- Inf
  expect_error(ivamse_fit(design$y, design$x, design$z, w = w), "must be finite")

  expect_error(ivamse_fit(design$y, design$x, design$z, screen = 2),
               "must be TRUE or FALSE")
  expect_error(ivamse_fit(design$y, design$x, design$z, cluster = 1:3),
               "observations but the data has")
  expect_error(ivamse_fit(design$y, design$x, design$z, pilot = 9999),
               "only .* candidates")
})

test_that("grouping variables may be named as columns", {
  design <- make_design()
  by_name <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data,
                    cluster = "firm")
  by_formula <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data,
                       cluster = ~ firm)
  expect_equal(vcov(by_name), vcov(by_formula))
})

test_that("s is caught before it silently becomes subset", {
  design <- make_design()
  expect_error(
    ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data,
           select = "cv", s = "lambda.1se"),
    "'cv_s'"
  )
})
