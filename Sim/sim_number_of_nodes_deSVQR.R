
# ============================================================== #
# Simulation: Number of nodes with deSVQR
# ============================================================== #

simulation_name <- "number_of_nodes_deSVQR"

#' @param Nreps the independent replications
#' @param T_outer the outer iteration
#' @param T_inner the inner iteration
#' @param n the local sample size
#' @param m the number of machines
#‘ @param N the whole sample size
#‘ @param p the row dimension
#‘ @param q the column dimension
#‘ @param r the ture rank of generated matrix 
#‘ @param pc the connection probability of the network
#‘ @param tau the quatile level
#‘ @param hetercase  hetercase = 1: data generation follows the setting in Section 4.4; hetercase = 2: data generation follows the setting in Section 4.1.
#‘ @param noise_type_arr different type of noise: Cauchy, Normal, Student's t(2)
#' @param X the input p*q matrix 
#' @param Y the response vector
#' @param B the coefficient matrix
#' @param tau_penalty_factor the penalty parameter in the augmented Lagrangian 
#' @param nlambda the length of tuning lambda
#' @function  decentralizedCQR_cpp  Decentralized surrogate vector quantile regression (deSVQR)



# ============================================================== #
# LOAD LIBRARY
# ============================================================== #

# library(FHDQR)
# library(FHDCQR)
library(cqrReg)

library(MASS)
library(pracma)
# for renyi graph
library(igraph)
library(tictoc)
library(glmnet)
library(ggplot2)
library(ggpubr)
library(decentralizedCQR)
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
  m_arr <- c(5, 10, 20)
  registerDoFuture()
  # use multiple workers to accelerate the time of replication, change
  # the number 123 to smaller number, e.g., 4 or 8 to accommodate your machine.
  plan(multisession, workers = 100)     ## on Linux, Solaris, and macOS
}
if (Platform == "Darwin") {
  Nreps <- 8
  T_outer <- 5
  T_inner <- 10
  noise_case_arr <- c("Cauchy")
  m_arr <- c(15)
  registerDoFuture()
  plan(multisession, workers = 8)    ## on MS Windows
}


# ============================================================== #
# PARAMETERS
# ============================================================== #


N <- 6000
p <- 10 # row dimension
q <- 10 # column dimension
r <- 3 # rank
pc <- 1 # the connection probability
tau = 1 / 2 # quatile level
rho <- .1
sigma2 <- 1
ishomo <- T
c0 <- 0.045
tau_penalty_factor <- 1 / 6
nlambda = 100
lambda_factor <- 1e-3
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
      cat("n = ", n, "p = ", p, "q = ", q, "\n")

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
      X <- data$X; #dim(data$X)
      y <- data$y
      betaT <- c(data$B)
      suppT <- which(abs(betaT)>0) ##support set
      graph <- data$graph
      adjacency_matrix <- as.matrix(as_adjacency_matrix(graph))


      # Estimate
      # Local estimate by BIC
      A_init_CQR <- matrix(rep(0, m * K), nrow = K)
      B_init_CQR <- matrix(rep(0, m * p), nrow = p)
      for (j in 1:m) { # j <- 1
        # get the index set for node j
        idx <- calN_j(n, j)
        B_init_CQR[, j] <- bic.hqreg(X[idx,], y[idx], tau_K = tau_K, nlambda = nlambda)
        A_init_CQR[, j] <-
          matrix(quantile(y[idx] - X[idx,] %*% B_init_CQR[, j], tau_K))
      }

       # deSCQR
      if (hetercase == 2) {
        tau_penalty_factor <- 0.05/2 
      }
      if (hetercase == 1) {
        tau_penalty_factor <- 0.1  
      }
      out_beta_deSCQR <- decentralizedCQR_cpp(
        X = X,
        y = y,
        adjacency_matrix = adjacency_matrix,
        A_init = A_init_CQR,
        B_init = B_init_CQR,
        betaT = betaT,
        T_outer = T_outer,
        T_inner = T_inner,
        s = 57,
        K = K,
        c0 = c0,
        tau_penalty_factor = tau_penalty_factor,
        lambda_factor = lambda_factor,
        nlambda = nlambda,
        lambda_max = lambda_max, # tunning parameter is very important; how to determine the upper bound
        quiet = quiet
      )
      error_deSCQR <- out_beta_deSCQR$history$errors_outer[length(out_beta_deSCQR$history$errors_outer)]
      rank_deSCQR <-
        mean(apply(
          out_beta_deSCQR$B,
          2,
          FUN = function(x)
            compute_rank(matrix(x, sqrt(p), sqrt(p)), cutoff = 1e-1)
        ))

      output_list[[in_p_arr]] <- c( error_deSCQR,
                                    rank_deSCQR)
    }
    output_list
  }
  result[[inoise_case_arr]] <- r_out

}
t_end <- toc(t_start)

# ============================================================== #
# OUPUT TABLE
# ============================================================== #
n_methods <- 1
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
