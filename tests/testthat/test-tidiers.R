skip_if_not_installed("broom")

test_that("tidy, glance and augment follow the broom contract", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data)

  tidied <- broom::tidy(fit, conf.int = TRUE)
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value",
                    "conf.low", "conf.high") %in% names(tidied)))
  expect_equal(nrow(tidied), length(coef(fit)))
  expect_equal(tidied$estimate, unname(coef(fit)))

  glanced <- broom::glance(fit)
  expect_equal(nrow(glanced), 1L)
  expect_equal(glanced$nobs, nobs(fit))
  expect_equal(glanced$dictionary, fit$selected$dictionary)

  augmented <- broom::augment(fit)
  expect_true(all(c(".fitted", ".resid") %in% names(augmented)))
  expect_equal(augmented$.resid, unname(residuals(fit)))
})
