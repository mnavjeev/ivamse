test_that("the criterion matches the formula in the paper", {
  fitted_moment <- c(0.5, 0.8, 1.2)
  d <- c(3, 10, 25)
  h <- c(0.4, 0.6, 0.9)
  sigma_eps2 <- 1.3
  sigma_epsv <- -0.6
  n <- 100

  expect_equal(
    feasible_criterion(fitted_moment, d, h, sigma_eps2, sigma_epsv, n),
    (sigma_eps2 * fitted_moment + sigma_epsv^2 * (d^2 + d) / n) / h^2
  )
})

test_that("a vanishing denominator gives an infinite criterion", {
  expect_equal(feasible_criterion(1, 2, 0, 1, 1, 100), Inf)
  expect_equal(feasible_criterion(1, 2, NA_real_, 1, 1, 100), Inf)
})

test_that("more endogeneity pushes the criterion toward simpler fits", {
  # Two candidates: a rich fit and a parsimonious one with slightly worse fit.
  fitted_moment <- c(1.00, 0.95)
  d <- c(40, 5)
  h <- c(1.0, 0.95)

  weak <- feasible_criterion(fitted_moment, d, h, 1, 0.01, 500)
  strong <- feasible_criterion(fitted_moment, d, h, 1, 0.90, 500)

  expect_equal(which.min(weak), 1L)   # near exogeneity: fit wins
  expect_equal(which.min(strong), 2L) # strong endogeneity: complexity is costly
})

test_that("ties are broken toward the smallest candidate index", {
  criterion <- c(2, 1, 1, 3)
  expect_equal(ivamse:::select_candidate(criterion, 1:4), 2L)
  expect_true(is.na(ivamse:::select_candidate(c(Inf, Inf), 1:2)))
})

test_that("the fallback taxonomy distinguishes its cases", {
  beta <- c(1, 2, NA)
  h <- c(1, 2, 0.5)

  expect_equal(ivamse:::resolve_estimate(1L, beta, h, 2L)$reason, "none")
  expect_false(ivamse:::resolve_estimate(1L, beta, h, 2L)$fallback)

  bad <- ivamse:::resolve_estimate(3L, beta, h, 2L)
  expect_equal(bad$reason, "nonfinite_estimate")
  expect_true(bad$fallback)
  expect_equal(bad$index, 2L)

  unresolved <- ivamse:::resolve_estimate(3L, beta, h, 3L)
  expect_true(unresolved$unresolved)
  expect_true(is.na(unresolved$index))

  zero <- ivamse:::resolve_estimate(1L, beta, c(0, 2, 0.5), 2L)
  expect_equal(zero$reason, "zero_denominator")
})
