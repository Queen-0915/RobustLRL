# ============================================================== #
# Simulation: Heterogeneity
# ============================================================== #

simulation_name <- "heterogeneity"

#' @param Nreps the independent replications
#' @param T_outer the outer iteration
#' @param T_inner the inner iteration
#' @param n the local sample size
#' @param m the number of machines
#' @param N the whole sample size
#' @param p the row dimension
#' @param q the column dimension
#' @param r the ture rank of generated matrix 
#' @param pc the connection probability of the network
#' @param tau the quatile level
#' @param hetercase  hetercase = 1: data generation follows the setting in Section 4.4; hetercase = 2: data generation follows the setting in Section 4.1.
#' @param noise_type_arr different type of noise: Cauchy, Normal, Student's t(2)
#' @param X the input p*q matrix 
#' @param Y the response vector
#' @param B the coefficient matrix
#' @param tau_penalty_factor the penalty parameter in the augmented Lagrangian 
#' @param nlambda the length of tuning lambda
#' @function decentralizedTraceQR_cpp  Our deSMQR method


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
  n_p_arr <- list(c(100, 10), c(200, 10), c(200, 20))
  heter_case_arr <- c(1, 2)
  registerDoFuture()
  # use multiple workers to accelerate the time of replication, change
  # the number 123 to smaller number, e.g., 4 or 8 to accommodate your machine.
  plan(multisession, workers = 50)   ## on Linux, Solaris, and macOS
}
if (Platform == "Darwin") {
  Nreps <- 8
  T_outer <- 5
  T_inner <- 20
  n_p_arr <- list(c(100, 10))
  heter_case_arr <- c(1)
  registerDoFuture()
  plan(multisession, workers = 8)    ## on MS Windows
}


# ============================================================== #
# PARAMETERS
# ============================================================== #

m <- 10 # the number of machines
r <- 3 # rank
pc <- .3 # the connection probability
tau = 1 / 2 # quatile level
rho <- .1
sigma2 <- 1
ishomo <- F
hetercase <- 1
c0 <- 0.045
nlambda = 100
lambda_factor <- 1e-3
lambda_max <- .1
quiet = T
MAXIT <- 2e3
eps <- 1e-3
result <- vector(mode = "list", length = length(heter_case_arr))


# ============================================================== #
# MAIN ROUNTINE
# ============================================================== #
n_methods <- 5
t_start <- tic()
for (iheter_case_arr in 1:length(heter_case_arr)) {
  hetercase <- heter_case_arr[iheter_case_arr]
  if (hetercase == 2) {
    tau_penalty_factor <- 0.05/6
  }
  if (hetercase == 1) {
    tau_penalty_factor <- 0.05/2
  }
  cat("The heterogenous type is ", hetercase, "\n")

  set.seed(2022) # fix the seed
  r_out <- foreach(
    iNreps = 1:Nreps,
    .init = replicate(length(n_p_arr), list()),
    .combine = "comb"
  ) %dorng% {

    output_list <-
      vector(mode = "list", length = length(n_p_arr))
    for (in_p_arr in 1:length(n_p_arr)) {
      # Load the parameter
      n <- n_p_arr[[in_p_arr]][1]
      p <- q <- n_p_arr[[in_p_arr]][[2]]
      N <- m * n
      cat("n = ", n, "p = ", p, "q = ", q, "\n")

      # Generate data
      data <- gen_data(
        p = p,
        q = q,
        r = r,
        N = N,
        m = m,
        pc = pc,
        noise_type = "Cauchy",
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
      adjacency_matrix <- as.matrix(as_adjacency_matrix(graph))  ##adjacency matrix of the network

      # Estimate
      # Pooled estimate by BIC
      B_init_pooled <- matrix(rnorm(p * q), p, q)
      out_pooled <- bic.quantile_trace_regression(
        X,
        y,
        tau = tau,
        nlambda = nlambda,
        B_init = B_init_pooled,
        B = B,
        eps = eps,
        MAXIT = MAXIT,
        quiet = T
      )
      B_pooled <- out_pooled$B
      error_pooled <- computeError(matrix(rep(B_pooled, m), p*q, m), betaT)
      rank_pooled <- compute_rank(matrix(B_pooled, p, q), cutoff = 1e-1)


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

      # Our deSMQR method
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

      output_list[[in_p_arr]] <- c(error_pooled, error_local, error_avg, error_deSubGD,
                                   error_deMQR,
                                   rank_pooled, rank_local, rank_avg, rank_deSubGD,
                                   rank_deMQR)
    }
    output_list
  }
  result[[iheter_case_arr]] <- r_out

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
for(iheter_case_arr in seq_along(heter_case_arr)) {
  noise_type <- heter_case_arr[iheter_case_arr]
  cat("The noise distribution is ", noise_type, "\n")
  r <- result[[iheter_case_arr]]
  output_table_error[[iheter_case_arr]] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))[,1:n_methods]
  output_table_rank[[iheter_case_arr]] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))[,(1 + n_methods):(2*n_methods)]

  output_table[[iheter_case_arr]] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))
  output_table[[iheter_case_arr]][,odds] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))[,1:n_methods]
  output_table[[iheter_case_arr]][,even] <- t(sapply(r, function(x) rowMeans(do.call(cbind, x), na.rm = T)))[,(1 + n_methods):(2*n_methods)]
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
