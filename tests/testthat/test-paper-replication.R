# Numerical agreement with the replication code for the paper.  The paper's data
# is not redistributable, so this compares the two implementations on a
# simulated design.  Point IVAMSE_PAPER_REPO at a checkout of the paper
# repository to run it.
test_that("candidate quantities match the paper's replication code", {
  repo <- Sys.getenv("IVAMSE_PAPER_REPO", "")
  skip_if(repo == "", "IVAMSE_PAPER_REPO is not set")
  reference <- file.path(repo, "simulations", "bcch_candidates.R")
  skip_if_not(file.exists(reference), "bcch_candidates.R not found")

  environment <- new.env()
  suppressWarnings(sys.source(reference, envir = environment))

  set.seed(101)
  N <- 500
  p <- 40
  z <- matrix(rnorm(N * p), N, p)
  z <- sweep(z, 2, colMeans(z), "-")
  z <- sweep(z, 2, sqrt(colMeans(z^2)), "/")
  v <- rnorm(N)
  x <- as.vector(z[, 1] + 0.7 * z[, 2] + 0.4 * z[, 3] + v)
  y <- 2 * x + 0.7 * v + rnorm(N)

  kappa <- 2^seq(-4, 2, by = 0.5)
  sigma_v <- 1

  # The reference code works without controls, so n equals N there.
  theirs <- environment$fit_candidate_list(
    Z_full = z, dictionaries = list(only = seq_len(p)),
    kappa_grids = list(only = kappa), x = x, n = N, scale = sigma_v
  )

  ours <- ivamse_fit(y, x, z, intercept = FALSE, penalty = "bcch",
                     kappa = kappa, sigma_v = sigma_v)

  expect_equal(ours$candidates$lambda, theirs$lambda, tolerance = 1e-12)
  expect_equal(ours$candidates$h, theirs$h, tolerance = 1e-10)
  expect_equal(ours$candidates$fitted_moment, theirs$fitted_moment,
               tolerance = 1e-10)
  expect_equal(ours$candidates$d, theirs$dimension)

  # And the criterion itself, computed the way the drivers compute it.
  pilot <- which.max(abs(theirs$h))
  beta <- colSums(theirs$fitted * y) / N / theirs$h
  residual <- y - beta[pilot] * x
  sigma_eps2 <- mean(residual^2)
  sigma_epsv <- mean(residual * x)
  their_criterion <- (sigma_eps2 * theirs$fitted_moment +
    sigma_epsv^2 * (theirs$dimension^2 + theirs$dimension) / N) / theirs$h^2

  expect_equal(ours$candidates$criterion, their_criterion, tolerance = 1e-10)
  expect_equal(ours$selected$index, which.min(their_criterion))
})
