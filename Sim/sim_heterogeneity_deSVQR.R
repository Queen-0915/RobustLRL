# ============================================================== #
# Simulation: Heterogeneity with deSVQR
# ============================================================== #

simulation_name <- "heterogeneity_deSVQR"


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
#' @function  decentralizedCQR_cpp  The Decentralized surrogate vector quantile regression (deSVQR) method


# ============================================================== #
# LOAD LIBRARY
# ============================================================== #

library(rqPen)
# library(FHDQR)
# library(FHDCQR)
library(cqrReg)
library(MASS)
library(pracma)
library(igraph) # for graph
library(tictoc)
library(hqreg)
library(glmnet)
library(ggplot2)
library(ggpubr)
# for parallel computing
library(doRNG)
library(foreach)
library(doFuture)
library(parallel)
library(decentralizedCQR)
library(xtable)
library(RobustLRL)


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
  n_p_arr <- list(c(100, 10), c(200, 10), c(200, 20))  ## 
  heter_case_arr <- c(1, 2)
  registerDoFuture()
  plan(multisession, workers = 50)   ## on Linux, Solaris, and macOS
}
if (Platform == "Darwin") {
  Nreps <- 8
  T_outer <- 10
  T_inner <- 50
  heter_case_arr <- c(1)
  registerDoFuture()
  plan(multisession, workers = 8)    ## on MS Windows
}


# ============================================================== #
# PARAMETERS
# ============================================================== #

m <- 10 # the number of machines
r <- 3
pc <- 0.3 # the connection probability
K <- 1  # the number of quantile levels 
tau_K <- seq(1, K) / (K + 1)
rho <- .1
sigma2 <- 1
ishomo <- F
c0 <- 0.045

nlambda <- 100L
lambda_factor <- 1e-3


lambda_max <- 0.1 #1
quiet <- TRUE


MAXIT <- 2e3
eps <- 1e-3
result <- vector(mode = "list", length = length(heter_case_arr))



# ============================================================== #
# MAIN ROUNTINE
# ============================================================== #
t_start <- tic()
for (iheter_case_arr in 1:length(heter_case_arr)) { # iheter_case_arr <- 1
  hetercase <- heter_case_arr[iheter_case_arr]
  cat("The heterogenous type is ", hetercase, "\n")

  set.seed(2022) # fix the seed
  r <- foreach(
    iNreps = 1:Nreps,
    .init = replicate(3, list()),
    .combine = "comb"
  ) %dorng% {
    output_list <-
      vector(mode = "list", length = length(n_p_arr))
    for (in_p_arr in 1:length(n_p_arr)) { #in_p_arr <- 1
      # Load the parameter
      n <- n_p_arr[[in_p_arr]][[1]]
      p <-  n_p_arr[[in_p_arr]][[2]]^2
      N <- m * n
      cat("n = ", n, "p = ", p, "\n")

      # Generate data
      data <- RobustLRL::gen_data(
        p = sqrt(p),
        q = sqrt(p),
        r = 3,
        N = N,
        m = m,
        pc = pc,
        noise_type = "Cauchy",
        rho = rho,
        sigma2 = sigma2,
        ishomo = F,
        hetercase = hetercase
      )

      X <- data$X;
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
    
       # The Decentralized surrogate vector quantile regression (deSVQR) method
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
  result[[iheter_case_arr]] <- r

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

xtable(do.call(rbind, output_table_f1), digits = 4)

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
  "2.RData"
))



