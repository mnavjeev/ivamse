# A small IV design with a sparse first stage, controls, and clustering.
make_design <- function(N = 400, p = 20, rho = 0.7, clusters = 40, seed = 1) {
  set.seed(seed)
  z <- matrix(rnorm(N * p), N, p)
  colnames(z) <- paste0("z", seq_len(p))
  w <- cbind(w1 = rnorm(N), w2 = rnorm(N))
  v <- rnorm(N)
  e <- rho * v + sqrt(1 - rho^2) * rnorm(N)
  x <- as.vector(z[, 1] + 0.7 * z[, 2] + 0.4 * z[, 3] + w %*% c(1, -1) + v)
  y <- as.vector(2 * x + w %*% c(1, -0.5) + e)
  list(
    y = y, x = x, z = z, w = w,
    firm = factor(rep_len(seq_len(clusters), N)),
    data = data.frame(y = y, x = x, w1 = w[, 1], w2 = w[, 2],
                      firm = factor(rep_len(seq_len(clusters), N)), z)
  )
}
