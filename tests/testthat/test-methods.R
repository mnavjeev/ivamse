test_that("the standard extractors behave", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data)

  expect_equal(nobs(fit), nrow(design$data))
  expect_length(residuals(fit), nobs(fit))
  expect_length(fitted(fit), nobs(fit))
  expect_equal(fitted(fit) + residuals(fit), design$data$y, ignore_attr = TRUE)
  expect_named(coef(fit), c("x", "(Intercept)", "w1", "w2"))
  expect_s3_class(formula(fit), "formula")
})

test_that("summary and print return invisibly and print nothing unasked", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data)

  expect_silent(s <- summary(fit))
  expect_silent(ci <- confint(fit))     # must not print, unlike hdm
  expect_true(is.matrix(ci))
  expect_output(print(fit), "First stage")
  expect_output(print(summary(fit)), "First-stage selection")
})

test_that("predict reproduces the fitted values and accepts new data", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data)

  expect_equal(predict(fit), fitted(fit))
  expect_equal(predict(fit, newdata = design$data), fitted(fit),
               ignore_attr = TRUE)

  with_se <- predict(fit, newdata = head(design$data, 5), se.fit = TRUE)
  expect_length(with_se$fit, 5L)
  expect_true(all(with_se$se.fit > 0))
})

test_that("selected() reports the chosen candidate and its instruments", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3 + z4 + z5, data = design$data)
  s <- selected(fit)

  expect_equal(s$candidate$dictionary, fit$selected$dictionary)
  expect_true(all(s$instruments != 0))
  expect_equal(length(s$instruments), sum(fit$first_stage != 0))
  expect_false(is.unsorted(rev(abs(s$instruments))))
})

test_that("model.matrix returns each requested component", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data)

  expect_equal(colnames(model.matrix(fit)), names(coef(fit)))
  expect_equal(dim(model.matrix(fit, "projected")), dim(model.matrix(fit)))
  expect_equal(ncol(model.matrix(fit, "instrument")), 1L)
})

test_that("plot runs", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1, data = design$data,
                dictionaries = list(small = ~ z1 + z2, large = ~ z1 + z2 + z3 + z4))
  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  expect_invisible(plot(fit))
  grDevices::dev.off()
  expect_true(file.exists(path))
})

test_that("update refits, including part-wise formula changes", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data)

  expect_equal(coef(update(fit, . ~ .)), coef(fit))

  wider <- update(fit, . ~ . | . | . + z4 + z5)
  expect_equal(unique(wider$candidates$p), unique(fit$candidates$p) + 2L)

  smaller <- update(fit, data = design$data[1:200, ])
  expect_equal(nobs(smaller), 200L)
})

test_that("deviance is the residual sum of squares", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data)
  expect_equal(deviance(fit), sum(residuals(fit)^2))
})

test_that("plot accepts graphical arguments", {
  design <- make_design()
  fit <- ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = design$data)
  path <- tempfile(fileext = ".png")
  grDevices::png(path)
  on.exit({grDevices::dev.off(); unlink(path)}, add = TRUE)

  # these all collided with the panel helper's own formals before
  expect_silent(plot(fit))
  expect_silent(plot(fit, main = c("a", "b")))
  expect_silent(plot(fit, col = "red"))
  expect_silent(plot(fit, ylab = "criterion"))
})
