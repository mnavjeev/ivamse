skip_if_not_installed("ivreg")
skip_if_not_installed("sandwich")

equivalent_ivreg <- function(fit, design) {
  d <- design$data
  d$ghat <- fit$instrument
  ivreg::ivreg(y ~ w1 + w2 | x | ghat, data = d)
}

test_that("the point estimate is the IV estimate with the fitted instrument", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data)
  reference <- equivalent_ivreg(fit, design)

  expect_equal(unname(coef(fit)),
               unname(coef(reference)[c("x", "(Intercept)", "w1", "w2")]))
})

test_that("homoskedastic and robust variances match ivreg", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data)
  reference <- equivalent_ivreg(fit, design)
  order <- c("x", "(Intercept)", "w1", "w2")

  expect_equal(unname(vcov(fit)),
               unname(vcov(reference)[order, order]), tolerance = 1e-10)

  for (type in c("HC0", "HC1", "HC2", "HC3")) {
    expect_equal(
      unname(sandwich::vcovHC(fit, type = type)),
      unname(sandwich::vcovHC(reference, type = type)[order, order]),
      tolerance = 1e-10
    )
  }
})

test_that("clustered variances match ivreg", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data,
                cluster = ~ firm)
  reference <- equivalent_ivreg(fit, design)
  order <- c("x", "(Intercept)", "w1", "w2")

  expect_equal(
    unname(sandwich::vcovCL(fit, cluster = design$data$firm)),
    unname(sandwich::vcovCL(reference, cluster = design$data$firm)[order, order]),
    tolerance = 1e-10
  )

  # vcov() uses the clustering given at fitting; cluster = NULL turns it off.
  expect_equal(vcov(fit), sandwich::vcovCL(fit, cluster = design$data$firm))
  expect_false(isTRUE(all.equal(vcov(fit), vcov(fit, cluster = NULL))))
})

test_that("the homoskedastic variance equals the formula in the paper", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data)

  sigma2 <- sum(residuals(fit)^2) / fit$df.residual
  by_hand <- sigma2 * fit$selected$fitted_moment /
    (fit$n * fit$selected$h^2)
  expect_equal(unname(vcov(fit)[1, 1]), by_hand, tolerance = 1e-10)
})

test_that("confint and summary route through vcov", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data)

  ci <- confint(fit)
  se <- sqrt(diag(vcov(fit)))
  expect_equal(ci[, 1], coef(fit) - qnorm(0.975) * se)

  s <- summary(fit)
  expect_s3_class(s, "summary.ivamse")
  expect_equal(s$coefficients[, "Std. Error"], se, ignore_attr = TRUE)

  # a supplied function is used, with ... passed on
  robust <- summary(fit, vcov. = sandwich::vcovHC, type = "HC1")
  expect_equal(robust$coefficients[, "Std. Error"],
               sqrt(diag(sandwich::vcovHC(fit, type = "HC1"))),
               ignore_attr = TRUE)
})

test_that("variance arguments are not silently discarded", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data)

  # a `type` argument means nothing without clustering; ignoring it silently
  # would report a non-robust standard error as though it were robust
  expect_warning(vcov(fit, type = "HC1"), "ignored without clustering")

  clustered <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data,
                      cluster = ~ firm)
  expect_false(isTRUE(all.equal(vcov(clustered),
                                vcov(clustered, cluster = NULL))))
})

test_that("tidy passes its arguments through to vcov", {
  skip_if_not_installed("broom")
  design <- make_design()
  clustered <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data,
                      cluster = ~ firm)

  expect_equal(broom::tidy(clustered, cluster = NULL)$std.error,
               sqrt(diag(vcov(clustered, cluster = NULL))), ignore_attr = TRUE)
  expect_false(isTRUE(all.equal(broom::tidy(clustered)$std.error,
                                broom::tidy(clustered, cluster = NULL)$std.error)))
})
