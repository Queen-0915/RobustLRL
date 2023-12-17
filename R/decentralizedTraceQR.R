decentralizedTraceQR <- function(X,
                                 y,
                                 adjacency_matrix,
                                 B_init,
                                 BT,
                                 tau = 1 / 2,
                                 r = 3,
                                 T_outer = 10,
                                 T_inner = 20,
                                 c0 = 0.1,
                                 tau_penalty_factor = 1 / 6,
                                 nlambda = 100,
                                 lambda_factor = 1e-4,
                                 lambda_max = 1,
                                 quiet = TRUE) {
  m <- ncol(adjacency_matrix)
  p <- nrow(X)
  q <- ncol(X)
  d <- p + q
  N <- dim(X)[3]
  n <- N / m
  svd.out <- svd(matrix(X[, , 1:n], p * q, n), nu = 0, nv = 0)$d
  tau_penalty <- svd.out[1]^2 / n * tau_penalty_factor

  yt <- matrix(rep(NA, N), N, 1)

  # cache omega and rho
  rho <- rep(NA, m)
  omega <- rep(NA, m)
  for (j in 1:m) {
    idx <- calN_j_cpp(n, j - 1) + 1
    s_vec <- svd(matrix(X[, , idx], p * q, n), nu = 0, nv = 0)$d
    rho[j] <- s_vec[1]^2/n
    omega[j] <-
      1 / (2 * tau_penalty * sum(adjacency_matrix[j, ] != 0) + rho[j])

  }


  # === Main routine === #
  B_outer <- B_init
  bic_array <- rep(NA, nlambda - 1)
  for (v in 1:T_outer) {
    hv <- (r * d * log(N) / N) ^ {
      1 / 2
    } + (r * d / n) ^ {
      1 / 2
    } * (r * d / n * log(N)) ^ (v / 2)
    if (!quiet) {
      cat("Outer iteration: v =", v, "\n")
    }

    # Construct the pseudo-responses
    for (j in 1:m) {
      idx <- calN_j_cpp(n, j - 1) + 1
      XtB <- apply(X[, , idx], 3, function(x) {
        sum(diag(t(x) %*% B_outer[, , j]))
      })
      E <-
        y[idx] - XtB
      fhat <- as.numeric(kernel(matrix(E, n, 1), matrix(hv, 1, 1)))
      if (!quiet) {
        cat("fhat =", fhat, "\n")
      }
      yt[idx] <- XtB - 1 / fhat * ((E <= 0) - tau)
    }

    # Tuning
    lambda_max <-
      norm(matrix(-matrix(X, p * q, N) %*% yt / m / n, p, q), "2")
    lambda_array <-
      exp(seq(log(1), log(lambda_factor), length.out = nlambda))
    lambda_array <- lambda_array[2:nlambda]
    lambda_array <- lambda_array * lambda_max
    if (!quiet) {
      cat("lambda_max =", lambda_max, "\n")
    }
    bic_array <- rep(0, nlambda - 1)

    for (ilambda in 1:(nlambda - 1)) {
      lambda <- lambda_array[ilambda]
      out.inner <- decentralizedTraceQR_inner(
        X,
        yt,
        adjacency_matrix = adjacency_matrix,
        rho = rho,
        omega = omega,
        B_outer = B_outer,
        tau_penalty = tau_penalty,
        lambda = lambda,
        T_inner = T_inner
      )

      # bic
      r_hat <- 0
      for (j in 1:m) {
        idx <- calN_j_cpp(n, j - 1) + 1
        bic_array[ilambda] <-
          bic_array[ilambda] + quantile_lossCPP(y[idx] - apply(X[, , idx], 3, function(x) {
            sum(diag(t(x) %*% out.inner$B_inner[, , j]))
          }), tau = tau)
        r_hat <- r_hat + pracma::Rank(out.inner$B_inner[, , j])
      }
      r_hat <- r_hat / m
      # if(!quiet) {
      #   cat("r_hat =", r_hat, "\n")
      # }
      bic_array[ilambda] <-
        bic_array[ilambda] / N + log(N) * log(log(N)) / N * r_hat
    }
    ilambda <- which.min(bic_array)
    lambda <- lambda_array[ilambda]
    if (!quiet) {
      cat("lambda =", lambda, "\n")
    }

    B_outer <- decentralizedTraceQR_inner(
      X,
      yt,
      adjacency_matrix = adjacency_matrix,
      rho = rho,
      omega = omega,
      B_outer = B_outer,
      tau_penalty = tau_penalty,
      lambda = lambda,
      T_inner = T_inner
    )$B_inner
  }

  return(B_outer)

}
