
# Simulation: number_of_nodes
simulation_name <- "number_of_nodes"


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
createdir(fig_dir)
createdir(LOG_dir)

if (Platform == "Linux") {
  Nreps <- 100
  T_outer <- 10
  T_inner <- 80
  noise_case_arr <- c("Normal", "T2", "Cauchy")
  # noise_case_arr <- c("Cauchy")
  m_arr <- c(5, 10, 20)
  # m_arr <- c(15)
  registerDoFuture()
  # use multiple workers to accelerate the time of replication, change
  # the number 123 to smaller number, e.g., 4 or 8 to accommodate your machine.
  # plan(multisession, workers = 50)    ## on MS Windows
  # plan(multicore, workers = 123)     ## on Linux, Solaris, and macOS
  plan(multisession, workers = 100)     ## on Linux, Solaris, and macOS
}
if (Platform == "Darwin") {
  Nreps <- 8
  T_outer <- 5
  T_inner <- 10
  # noise_case_arr <- c("T2")
  noise_case_arr <- c("Cauchy")
  m_arr <- c(15)
  registerDoFuture()
  plan(multisession, workers = 8)    ## on MS Windows
}


# ============================================================== #
# PARAMETERS
# ============================================================== #

# m <- 10 # the number of machines
# n <- 2e2 # local sample size
# N <- m*n # sample size
# N <- 4200
# N <- 3000
# N <- 4600
N <- 6000
p <- 10 # row dimension
q <- 10 # column dimension
r <- 3 # rank
# pc <- .3 # the connection probability
pc <- 1 # the connection probability
tau = 1 / 2 # quatile level
rho <- .1
sigma2 <- 1
ishomo <- T
# c0 <- 0.04
# c0 <- 0.013
c0 <- 0.045
tau_penalty_factor <- 1 / 6
# tau_penalty_factor <- 0.05
# tau_penalty_factor <- 0.1
# tau_penalty_factor <- 0.05/4
# tau_penalty_factor <- 0.05/6
nlambda = 100
lambda_factor <- 1e-3
# lambda_factor = 1e-4
lambda_max <- .1
quiet = T
MAXIT <- 2e3
eps <- 1e-3
result <- vector(mode = "list", length = length(noise_case_arr))


# ============================================================== #
# MAIN ROUNTINE
# ============================================================== #
n_methods <- 4
t_start <- tic()
for (inoise_case_arr in 1:length(noise_case_arr)) {
  noise_case <- noise_case_arr[inoise_case_arr]
  cat("The noise_case type is ", noise_case, "\n")

  set.seed(2022) # fix the seed
  r_out <- foreach(
    iNreps = 1:Nreps,
    .init = replicate(length(m_arr), list()),
    .combine = "comb"
  ) %dorng% {

    output_list <-
      vector(mode = "list", length = length(m_arr))
    for (im_arr in 1:length(m_arr)) {
      # Load the parameter
      m <- m_arr[[im_arr]]
      n <- N/m
      # c0 <- log(2000)/200*0.045/log(N)*n
      cat("n = ", n, "p = ", p, "q = ", q, "\n")

      # RNGkind("L'Ecuyer-CMRG")
      # .Random.seed <- attr(r, "rng")[[30]]

      # Generate data
      data <- gen_data(
        p = p,
        q = q,
        r = r,
        N = N,
        m = m,
        pc = pc,
        noise_type = noise_case,
        rho = rho,
        sigma2 = sigma2,
        ishomo = ishomo
      )
      X <- data$X
      y <- data$y
      B <- data$B
      betaT <- matrix(as.numeric(B), p * q, 1)
      graph <- data$graph
      adjacency_matrix <- as.matrix(as_adjacency_matrix(graph))

      # Estimate

      # Pooled estimate by BIC
      # B_init_pooled <- matrix(rnorm(p * q), p, q)
      # out_pooled <- bic.quantile_trace_regression(
      #   X,
      #   y,
      #   tau = tau,
      #   nlambda = nlambda,
      #   B_init = B_init_pooled,
      #   B = B,
      #   eps = eps,
      #   MAXIT = MAXIT,
      #   quiet = T
      # )
      # B_pooled <- out_pooled$B
      # error_pooled <- computeError(matrix(rep(B_pooled, m), p*q, m), betaT)
      # rank_pooled <- compute_rank(matrix(B_pooled, p, q), cutoff = 1e-1)


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
        )
        B_init[, im] <- c(out_local$B)
      }
      error_local <- computeError(B_init, betaT)
      rank_local <- mean(apply(
        B_init,
        2,
        FUN = function(x)
          compute_rank(matrix(x, p, q), cutoff = 1e-1)
      ))

      # Avg estimate
      B_avg <- rowMeans(B_init)
      error_avg <- computeError(B_avg, betaT)
      rank_avg <- compute_rank(B_avg, cutoff = 1e-1)

      # deSubGD
      out_deSubGD <-  deSubGD(X,
                              y,
                              d1 = p,
                              d2 = q,
                              adjacency_matrix,
                              # B_init = B_init,
                              B_init = matrix(rnorm(p * q * m), p * q, m),
                              betaT,
                              T = T_inner*T_outer,
                              lambda_max = 20,
                              eps = eps,
                              quiet = quiet)
      B_deSubGD <- out_deSubGD$B
      error_deSubGD <- computeError(B_deSubGD, betaT)
      rank_deSubGD <- mean(apply(
        B_deSubGD,
        2,
        FUN = function(x)
          compute_rank(matrix(x, p, q), cutoff = 1e-1)
      ))

      # deSCQR
      out_deMQR <- decentralizedTraceQR_cpp(
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
        quiet = quiet
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

      output_list[[im_arr]] <- c(error_local, error_avg, error_deSubGD,
                                   error_deMQR,
                                   rank_local, rank_avg, rank_deSubGD,
                                   rank_deMQR)
    }
    output_list
  }
  result[[inoise_case_arr]] <- r_out

}
t_end <- toc(t_start)

# ============================================================== #
# OUPUT TABLE
# ============================================================== #
output_table_error <- vector(mode = "list", length = length(result))
output_table_rank <- vector(mode = "list", length = length(result))
output_table <- vector(mode = "list", length = length(result))
odds <- seq(1, 2*n_methods - 1, by = 2)
even <- seq(2, 2*n_methods, by = 2)
for(inoise_case_arr in seq_along(noise_case_arr)) {
  noise_type <- noise_case_arr[inoise_case_arr]
  cat("The noise distribution is ", noise_type, "\n")
  r <- result[[inoise_case_arr]]
  output_table_error[[inoise_case_arr]] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))[,1:n_methods]
  output_table_rank[[inoise_case_arr]] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))[,(1 + n_methods):(2*n_methods)]

  output_table[[inoise_case_arr]] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))
  output_table[[inoise_case_arr]][,odds] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))[,1:n_methods]
  output_table[[inoise_case_arr]][,even] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))[,(1 + n_methods):(2*n_methods)]
}
xtable(do.call(rbind, output_table_error), digits = 4)

xtable(do.call(rbind, output_table_rank), digits = 4)

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
