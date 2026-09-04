test_that("the estimator recovers the structural coefficient", {
  design <- make_design(N = 800, seed = 4)
  fit <- ivamse_fit(design$y, design$x, design$z, w = design$w)
  expect_equal(unname(coef(fit)[1]), 2, tolerance = 0.15)
})

test_that("the KKT denominator identity holds along the whole grid", {
  design <- make_design()
  w_qr <- ivamse:::control_qr(cbind(1, design$w))
  n <- length(design$y) - w_qr$rank
  prepared <- ivamse:::prepare_dictionary(design$z, w_qr, n)
  x_res <- ivamse:::partial_out(design$x, w_qr)
  lambda <- ivamse:::lambda_path(prepared$z, x_res, n, nlambda = 20)
  fit <- ivamse:::fit_dictionary(prepared$z, x_res, lambda, n, qr(prepared$z)$rank)

  # h = E_n[Pi_hat^2] + lambda * ||pi||_1, so in particular h is never negative.
  expect_equal(fit$h, fit$fitted_moment + lambda * fit$pi_l1, tolerance = 1e-7)
  expect_true(all(fit$h >= 0))
})

test_that("the two formula shapes describe the same model", {
  design <- make_design()
  d <- design$data
  three <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = d)
  two <- ivamse(y ~ x + w1 + w2 | z1 + z2 + z3 + z4 + z5 + w1 + w2, data = d)

  expect_equal(coef(three), coef(two))
  expect_equal(three$candidates$criterion, two$candidates$criterion)
})

test_that("the formula and matrix interfaces agree", {
  design <- make_design()
  by_formula <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5,
                       data = design$data)
  by_matrix <- ivamse_fit(design$y, design$x, design$z[, 1:5], w = design$w,
                          x_name = "x")

  expect_equal(unname(coef(by_formula)), unname(coef(by_matrix)))
})

test_that("several dictionaries are scored together and one is chosen", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1, data = design$data,
                dictionaries = list(small = ~ z1 + z2 + z3,
                                    large = ~ z1 + z2 + z3 + z4 + z5 + z6 + z7))

  expect_equal(unique(fit$candidates$dictionary), c("small", "large"))
  expect_equal(nrow(fit$candidates), 50L)
  expect_true(fit$selected$dictionary %in% c("small", "large"))
  expect_equal(sum(fit$candidates$used), 1L)
})

test_that("more than one endogenous regressor is refused", {
  design <- make_design()
  d <- design$data
  d$x2 <- rnorm(nrow(d))
  expect_error(
    ivamse(y ~ w1 | x + x2 | z1 + z2 + z3, data = d),
    "one endogenous regressor"
  )
})

test_that("the BCCH grid puts kappa = 1 at the plug-in penalty", {
  design <- make_design()
  fit <- ivamse_fit(design$y, design$x, design$z, w = design$w,
                    penalty = "bcch", kappa = c(0.5, 1, 2))

  expect_equal(fit$candidates$kappa, c(0.5, 1, 2))
  expect_equal(
    fit$candidates$lambda[2],
    bcch_lambda(ncol(design$z), fit$n, sigma_v = fit$sigma_v)
  )
})

test_that("the comparators run and pick a candidate", {
  design <- make_design()
  for (rule in c("amse", "cv", "bcch")) {
    fit <- ivamse_fit(design$y, design$x, design$z, w = design$w, select = rule)
    expect_equal(fit$select, rule)
    expect_true(is.finite(coef(fit)[1]))
  }
})

test_that("selection is invariant to rescaling the outcome", {
  design <- make_design()
  base <- ivamse_fit(design$y, design$x, design$z, w = design$w)
  scaled <- ivamse_fit(10 * design$y, design$x, design$z, w = design$w)

  expect_equal(base$selected$index, scaled$selected$index)
  expect_equal(unname(coef(scaled)[1]), 10 * unname(coef(base)[1]))
})

test_that("penalty grids can differ across dictionaries", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1, data = design$data,
                dictionaries = list(small = ~ z1 + z2, large = ~ z1 + z2 + z3),
                nlambda = list(small = 10, large = 4))
  expect_equal(unname(table(fit$candidates$dictionary)[c("small", "large")]),
               c(10L, 4L), ignore_attr = TRUE)

  # the asymmetric grid of the empirical application
  bcch <- ivamse(y ~ w1 + w2 | x | z1, data = design$data,
                 dictionaries = list(small = ~ z1 + z2, large = ~ z1 + z2 + z3),
                 penalty = "bcch",
                 kappa = list(small = 2^seq(-6, 2, by = 0.5),
                              large = 2^seq(-4, 2, by = 0.5)))
  expect_equal(nrow(bcch$candidates), 17L + 13L)
})

test_that("inert arguments are reported rather than ignored", {
  design <- make_design()
  expect_warning(ivamse_fit(design$y, design$x, design$z, kappa = c(1, 2)),
                 "ignored when penalty")
  expect_warning(ivamse_fit(design$y, design$x, design$z, penalty = "bcch",
                            nlambda = 5), "ignored when penalty")
  expect_warning(ivamse_fit(design$y, design$x, design$z,
                            lambda = c(0.05, 0.1), nlambda = 5), "ignored")
})

test_that("cross-validation records the folds it used", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4, data = design$data,
                select = "cv", nfolds = 5)
  expect_length(fit$foldid, nobs(fit))
  expect_equal(length(unique(fit$foldid)), 5L)

  # supplying them back reproduces the selection exactly
  again <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4, data = design$data,
                  select = "cv", foldid = fit$foldid)
  expect_equal(again$selected$index, fit$selected$index)
})
