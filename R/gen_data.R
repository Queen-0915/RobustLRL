#' @export
gen_data <- function(p = 20,
                     q = 20,
                     N = 1e3,
                     m = 10,
                     pc = .3,
                     r = 3,
                     sigma2 = 1,
                     rho = 0.1,
                     noise_type = "Cauchy",
                     ishomo = TRUE,
                     hetercase = 1) {
  if (!(noise_type %in% c("Cauchy", "Normal", "T2")))
    stop("The error type must be one of Cauchy, Normal and T2.")
  U <- matrix(rnorm(p * r), p, r)
  U.qr <- qr(U)
  U <- qr.Q(U.qr)
  V <- matrix(rnorm(q * r), q, r)
  V.qr <- qr(V)
  V <- qr.Q(V.qr)
  B <- U %*% t(V)
  betaT <- as.vector(B)
  if (ishomo) {

    if (noise_type == "Cauchy") {
      noise <- rcauchy(N)
    }
    if (noise_type == "Normal") {
      noise <- rnorm(N)
    }
    if (noise_type == "T2") {
      noise <- rt(N, 2)
    }
    # X <- matrix(rnorm(p*q*N), N, p*q)
    Sigma <- sigma2 * toeplitz(rho^seq(0, p*q - 1, by = 1))
    X <- MASS::mvrnorm(N, rep(0, p*q), Sigma)
    y <- X%*%as.vector(B) + noise
    # The following code is wrong.
    # X <- array(rnorm(p*q*N), dim = c(p, q, N))
    # y <- apply(X, 3, function(x) {
    #   sum(diag(t(x)%*%B))
    #   }) + noise
    graph <- igraph::sample_gnp(m, pc)
    while(!igraph::is_connected(graph)) {
      graph <- igraph::sample_gnp(m, pc)
    }


    return(list(
      U = U,
      V = V,
      B = B,
      X = X,
      y = y,
      graph = graph
    ))
  } else {
    graph <- igraph::sample_gnp(m, pc)
    while(!igraph::is_connected(graph)) {
      graph <- igraph::sample_gnp(m, pc)
    }

    # Generate heterogeneous data
    n <- N/m
    if (hetercase == 1) {
      X <- matrix(rep(NA, N * p*q), ncol = p*q)
      y <- matrix(rep(NA, N), ncol = 1)
      if (noise_type == "Cauchy")
        noise <- rcauchy(N)
      if (noise_type == "Normal")
        noise <- rnorm(N)
      if (noise_type == "T2")
        noise <- rt(N, 2)
      sigma2_set <- c(1, 3)
      rho_set <- c(0.1, 0.3)
      for (j in 1:m) {
        idx <- (1 + (j - 1) * n):(j * n)
        sigma2 <- sample(sigma2_set, 1)
        rho <- sample(rho_set, 1)
        Sigma <- sigma2 * toeplitz(rho^seq(0, p*q - 1, by = 1))
        X[idx,] <- MASS::mvrnorm(n, rep(0, p*q), Sigma)
        y[idx] <- X[idx,] %*% betaT + noise[idx]
      }
    }
    if (hetercase == 2) {
      # Generate heterogeneous data
      Sigma <- sigma2 * toeplitz(rho^seq(0, p*q - 1, by = 1))
      X <- MASS::mvrnorm(N, rep(0, p*q), Sigma)
      idx <- ceiling(runif(m, 0, 3))
      noise <- (rep(idx == 1, each = n))*rcauchy(N) + (rep(idx == 2, each = n))*rnorm(N) + (rep(idx == 3, each = n))*rt(N, 2)
      y <- X %*% betaT + noise
    }

    return(list(
      U = U,
      V = V,
      B = B,
      X = X,
      y = y,
      graph = graph
    ))
  }

}
