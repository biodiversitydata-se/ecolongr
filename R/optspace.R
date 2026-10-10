# OptSpace matrix completion (Keshavan, Montanari and Oh 2010) as in vegan::optspace(), but
# on the observed cells only. The steps, constants and the line search follow vegan's code
# so that the same fixed number of iterations gives the same matrix.

optspace_g <- function(x, m0, r) {
  z <- rowSums(x^2) / (2 * m0 * r)
  y <- exp((z - 1)^2) - 1
  y[z < 1] <- 0
  sum(y)
}

optspace_gp <- function(x, m0, r) {
  z <- rowSums(x^2) / (2 * m0 * r)
  z <- 2 * exp((z - 1)^2) / (z - 1)
  z[z < 0] <- 0
  x * z / (m0 * r)
}

optspace_triplets <- function(i, j, v, n, p, r = 3, niter = 5, tol = 1e-5, threads = 1L) {
  nnz <- length(v)
  eps <- nnz / sqrt(p * n)
  rho <- eps * n
  rescale <- sqrt(nnz * r / sum(v^2))
  v <- v * rescale

  observed <- Matrix::sparseMatrix(i = i + 1L, j = j + 1L, x = v, dims = c(n, p))
  sv <- if (as.numeric(n) * p <= 2e6 || r >= min(n, p) - 1L) {
    s <- svd(as.matrix(observed), nu = r, nv = r)
    list(u = s$u, v = s$v)
  } else {
    irlba::irlba(observed, nv = r, nu = r, tol = 1e-12, maxit = 5000)
  }
  reverse <- rev(seq_len(r))
  X <- sv$u[, reverse, drop = FALSE] * sqrt(n)
  Y <- sv$v[, reverse, drop = FALSE] * sqrt(p)

  get_s <- function(X, Y) {
    c_vec <- as.vector(crossprod(X, cell_multiply(i, j, v, Y, n, threads)))
    a_mat <- cell_gram(X, Y, i, j, threads)
    matrix(solve(a_mat, c_vec), nrow = r)
  }
  fit_error <- function(X, Y, S) cell_fitted(X %*% S, Y, i, j, threads) - v
  objective <- function(X, Y, S, m0) {
    0.5 * sum(fit_error(X, Y, S)^2) + rho * optspace_g(Y, m0, r) + rho * optspace_g(X, m0, r)
  }

  m0 <- 10000
  S <- get_s(X, Y)
  for (iter in seq_len(niter)) {
    XS <- X %*% S
    YS <- Y %*% t(S)
    resid <- cell_fitted(XS, Y, i, j, threads) - v
    RY <- cell_multiply(i, j, resid, YS, n, threads)
    RX <- cell_multiply(j, i, resid, XS, p, threads)
    Qx <- -crossprod(X, RY) / n
    Qy <- -crossprod(Y, RX) / p
    W <- RY + X %*% Qx + rho * optspace_gp(X, m0, r)
    Z <- RX + Y %*% Qy + rho * optspace_gp(Y, m0, r)

    norm2 <- sum(W^2) + sum(Z^2)
    f0 <- objective(X, Y, S, m0)
    t <- -0.1
    for (step in seq_len(21)) {
      if (objective(X + t * W, Y + t * Z, S, m0) - f0 <= 0.5 * t * norm2) break
      t <- t / 2
    }
    X <- X + t * W
    Y <- Y + t * Z
    S <- get_s(X, Y)
    if (sqrt(sum(fit_error(X, Y, S)^2)) / sqrt(nnz) < tol) break
  }

  S <- S / rescale
  keep <- order(diag(S), decreasing = TRUE)
  list(X = X[, keep, drop = FALSE], S = S[keep, keep, drop = FALSE], Y = Y[, keep, drop = FALSE])
}
