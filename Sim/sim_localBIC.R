# ============================================================== #
# Simulation: LocalBIC V.S. Decentralized BIC
# ============================================================== #

simulation_name <- "localBIC"


# ============================================================== #
# LOAD LIBRARY
# ============================================================== #

library(MASS)
library(pracma)
# for renyi graph
library(igraph)
library(tictoc)
library(glmnet)
library(ggplot2)
library(ggpubr)
# for parallel computing
library(doRNG)
library(foreach)
library(doFuture)
library(parallel)
library(RobustLRL)
library(xtable)

# ============================================================== #
# PREPROCESSING
# ============================================================== #

Platform <- Sys.info()['sysname']
fig_dir <- "Output/figs"
LOG_dir <- "Output/LOG"
# createdir(fig_dir)
# createdir(LOG_dir)

if (Platform == "Linux") {
  Nreps <- 100
  T_outer <- 10
  T_inner <- 80
  noise_type_arr <- c("Cauchy", "Normal", "T2")
  registerDoFuture()
  # use multiple workers to accelerate the time of replication, change
  # the number 123 to smaller number, e.g., 4 or 8 to accommodate your machine.
  # plan(multisession, workers = 100)    ## on MS Windows
  plan(multisession, workers = 100)     ## on Linux, Solaris, and macOS
}
if (Platform == "Darwin") {
  Nreps <- 8
  T_outer <- 10
  T_inner <- 80
  noise_type_arr <- c("Cauchy")
  registerDoFuture()
  plan(multisession, workers = 8)    ## on MS Windows
}


# ============================================================== #
# PARAMETERS
# ============================================================== #

m <- 10 # the number of machines
# n <- 2e2 # local sample size
# N <- m*n # sample size
# p <- 10 # row dimension
# q <- 10 # column dimension
n_p_arr <- list(c(100, 10), c(200, 10), c(200, 20))
r <- 3 # rank
pc <- .3 # the connection probability
tau <- 1 / 2 # quatile level
rho <- .1
sigma2 <- 1
ishomo <- T
hetercase <- 1
c0 <- 0.045
tau_penalty_factor <- 1 / 6
nlambda = 100
lambda_factor <- 1e-3
# lambda_factor = 1e-4
lambda_max <- .1
quiet = T
MAXIT <- 2e3
eps <- 1e-3
result <- vector(mode = "list", length = length(noise_type_arr))


# ============================================================== #
# MAIN ROUNTINE
# ============================================================== #
n_methods <- 1
tic()
for (inoise_type_arr in 1:length(noise_type_arr)) { # inoise_type_arr <- 1
  noise_type <- noise_type_arr[inoise_type_arr]
  cat("The noise distribution is ", noise_type, "\n")

  set.seed(2022) # fix the seed
  r_out <- foreach(
    iNreps = 1:Nreps,
    .init = replicate(length(n_p_arr), list()),
    .combine = "comb"
  ) %dorng% {

    output_list <-
      vector(mode = "list", length = length(n_p_arr))
    for (in_p_arr in 1:length(n_p_arr)) { # in_p_arr <- 1
      c0 <- 0.045
      if(in_p_arr == 3) {
        c0 <- c0/2
      }
      # Load the parameter
      n <- n_p_arr[[in_p_arr]][1]
      p <- q <- n_p_arr[[in_p_arr]][[2]]
      N <- m * n
      cat("n = ", n, "p = ", p, "q = ", q, "\n")

      # Generate data
      # fix seed
      # RNGkind("L'Ecuyer-CMRG")
      # .Random.seed <- attr(r0, "rng")[[98]]
      data <- gen_data(
        p = p,
        q = q,
        r = r,
        N = N,
        m = m,
        pc = pc,
        noise_type = noise_type,
        rho = rho,
        sigma2 = sigma2,
        ishomo = ishomo,
        hetercase = hetercase
      )
      X <- data$X
      y <- data$y
      B <- data$B
      betaT <- matrix(as.numeric(B), p * q, 1)
      graph <- data$graph
      adjacency_matrix <- as.matrix(as_adjacency_matrix(graph))

      # Estimate


      # Initial
      # MQR
      # Local estimate by BIC
      B_init_local <- matrix(rnorm(p * q), p, q)
      B_init <- matrix(rep(NA, p * q * m), p * q, m)
      for (im in 1:m) {
        idx <- ((im - 1) * n + 1):(im * n)
        out_local <- bic.quantile_trace_regression(
          X[idx,],
          y[idx],
          tau = tau,
          nlambda = nlambda,
          B_init = B_init_local,
          B = B,
          eps = eps,
          MAXIT = MAXIT,
          quiet = T
        )
        B_init[, im] <- c(out_local$B)
      }
     # error_init <- computeError(B_init, betaT)

      # # Trace Reg
      # # Local estimate by BIC
      # B_init_local <- matrix(rnorm(p * q), p, q)
      # B_init_LS <- matrix(rep(NA, p * q * m), p * q, m)
      # for (im in 1:m) {
      #   idx <- ((im - 1) * n + 1):(im * n)
      #   out_local <- bic.trace_regression(
      #     X[idx,],
      #     y[idx],
      #     nlambda = nlambda,
      #     B_init = B_init_local,
      #     B = B,
      #     eps = eps,
      #     MAXIT = MAXIT,
      #     quiet = F
      #   )
      #   B_init_LS[, im] <- c(out_local$B)
      # }
      # error_init_LS <- computeError(B_init_LS, betaT)


      out_deMQR <- decentralizedTraceQR_local_cpp(
        X,
        y,
        d1 = p,
        d2 = q,
        adjacency_matrix,
        B_init,
        betaT = betaT,
        s  = r,
        tau = tau,
        T_outer = T_outer,
        T_inner = T_inner,
        c0 = c0,
        tau_penalty_factor = tau_penalty_factor,
        nlambda = nlambda,
        lambda_factor = lambda_factor,
        lambda_max = lambda_max,
        quiet = T
      )
      error_deMQR <-
        out_deMQR$history$errors_inner[length(out_deMQR$history$errors_inner)]
      rank_deMQR <-
        mean(apply(
          out_deMQR$B,
          2,
          FUN = function(x)
            compute_rank(matrix(x, p, q), cutoff = 1e-1)
        ))

      # out_deLR <- deTraceReg(
      #   X,
      #   y,
      #   d1 = p,
      #   d2 = q,
      #   adjacency_matrix = adjacency_matrix,
      #   B_init = B_init_LS,
      #   betaT = betaT,
      #   T_inner = T_inner * T_outer,
      #   tau_penalty_factor = tau_penalty_factor,
      #   nlambda = nlambda,
      #   lambda_factor = lambda_factor,
      #   lambda_max = lambda_max,
      #   quiet = quiet
      # )
      # error_deLR <-
      #   out_deLR$history$errors_inner[length(out_deLR$history$errors_inner)]
      # rank_deLR <-
      #   mean(apply(
      #     out_deLR$B,
      #     2,
      #     FUN = function(x)
      #       compute_rank(matrix(x, p, q), cutoff = 1e-1)
      #   ))


      output_list[[in_p_arr]] <- c(error_deMQR, rank_deMQR)
    }
    output_list
  }
  result[[inoise_type_arr]] <- r_out

  # output table
  # output_table <- vector(mode = "list", length = 1)
  # r0 <- result[[inoise_type_arr]]
  # output_table[[1]] <-
  #   t(sapply(r0, function(x)
  #     rowMeans(do.call(cbind, x), na.rm = T)))
  # xtable(do.call(rbind, output_table), digits = 4)

}
toc()


# ============================================================== #
# OUPUT TABLE
# ============================================================== #
output_table <- vector(mode = "list", length = length(result))
for (inoise_type_arr in seq_along(noise_type_arr)) {
  noise_type <- noise_type_arr[inoise_type_arr]
  cat("The noise distribution is ", noise_type, "\n")
  r0 <- result[[inoise_type_arr]]
  output_table[[inoise_type_arr]] <-
    t(sapply(r0, function(x)
      rowMeans(do.call(cbind, x), na.rm = T)))
}
xtable(do.call(rbind, output_table), digits = 4)

# ============================================================== #
# SAVE DATA
# ============================================================== #

save.image(file = paste0(
  LOG_dir,
  "/sim_",
  simulation_name,
  "_",
  format(Sys.time(), "%Y%m%d%H%M%S"),
  ".RData"
))
