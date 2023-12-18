# ============================================================== #
# Simulation: Iteration
# ============================================================== #

simulation_name <- "iteration"

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

# ============================================================== #
# PREPROCESSING
# ============================================================== #

Platform <- Sys.info()['sysname']
fig_dir <- "Output/figs"
LOG_dir <- "Output/LOG"
createdir(fig_dir)
createdir(LOG_dir)
if (Platform == "Linux") {
  Nreps <- 100  #the independent replications
  T_outer <- 50 #Outer iteration
  noise_type_arr <- c("Cauchy", "Normal", "T2")
  registerDoFuture()
  # use multiple workers to accelerate the time of replication, change
  # the number 123 to smaller number, e.g., 4 or 8 to accommodate your machine.
  plan(multisession, workers = 50)     ## on Linux, Solaris, and macOS
}
if (Platform == "Darwin") {
  Nreps <- 10
  T_outer <- 10
  noise_type_arr <- c("Cauchy")
  registerDoFuture()
  plan(multisession, workers = 10)    ## on MS Windows
}


# ============================================================== #
# PARAMETERS
# ============================================================== #

m <- 10 # the number of machines
n <- 2e2 # local sample size
N <- m * n # sample size
p <- 10 # row dimension
q <- 10 # column dimension
r <- 3 # rank
pc <- .3 # the connection probability
tau <- 1 / 2 # quatile level
rho <- .1
sigma2 <- 1
ishomo <- T
hetercase <- 1
c0 <- 0.045
tau_penalty_factor <-  1 / 6

T_inner_arr <- c(40, 80, 120) # array of the inner iteration
nlambda <- 100
lambda_factor <- 1e-3
lambda_max <- .1
quiet = F
MAXIT <- 2e3
eps <- 1e-3
result <- vector(mode = "list", length = length(noise_type_arr))


# ============================================================== #
# MAIN ROUNTINE
# ============================================================== #


t_start <- tictoc::tic()
for (inoise_type_arr in 1:length(noise_type_arr)) {
  noise_type <- noise_type_arr[inoise_type_arr]
  cat("The noise distribution is ", noise_type)

  set.seed(2022) # fix the seed
  r_out <- foreach(
    iNreps = 1:Nreps,
    .init = list(list(), list(), list()),
    .combine = "comb"
  ) %dorng% {
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


    # Initial
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
    error_init <- computeError(B_init, betaT)

    output_list <-
      vector(mode = "list", length = length(T_inner_arr))
    for (iT_inner_arr in 1:length(T_inner_arr)) {
      # Estimate
      T_inner <- T_inner_arr[iT_inner_arr]

      tic() ##calculate computation time
      ### Our deSMQR method
      out_deRRR <- decentralizedTraceQR_cpp(
        X,
        y,
        d1 = p,
        d2 = q,
        adjacency_matrix,
        B_init,
        betaT,
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
      toc()

      output_list[[iT_inner_arr]] <-
        c(error_init, out_deRRR$history$errors_outer)
    }
    output_list
  }
  result[[inoise_type_arr]] <- r_out

}
t_end <- tictoc::toc(t_start)


# ============================================================== #
# PLOT-- logY
# ============================================================== #

for (inoise_type_arr in 1:length(noise_type_arr)) {
  noise_type <- noise_type_arr[inoise_type_arr]
  cat("The noise distribution is ", noise_type, "\n")

  r <- result[[inoise_type_arr]]


  pdf(paste0(fig_dir, "/fig_iterations_log_", noise_type, ".pdf"))

  op <- par(no.readonly = TRUE)
  fontsize <- 12
  type <- "l"
  lwd <- 3
  par(
    font = fontsize,
    font.axis = fontsize,
    font.lab = fontsize,
    font.main = fontsize,
    font.sub = fontsize,
    cex = 1.5,
    mar = c(2, 2, 1, 1)
  )
  plot(
    log(rowMeans(do.call(cbind, r[[1]]), na.rm = T)),
    type = type,
    ylim = c(
      min(
        log(rowMeans(do.call(cbind, r[[1]]))), ### Using log-Y
        log(rowMeans(do.call(cbind, r[[2]]))),
        log(rowMeans(do.call(cbind, r[[3]]))),
        na.rm = T
      ),
      max(
        log(rowMeans(do.call(cbind, r[[1]]))),
        log(rowMeans(do.call(cbind, r[[2]]))),
        log(rowMeans(do.call(cbind, r[[3]]))),
        na.rm = T
      )
    ),
    lwd = lwd,
    xlab = "",
    ylab = ""
  )
  lines(
    log(rowMeans(do.call(cbind, r[[2]]), na.rm = T)),
    type = type,
    pch = 2,
    col = 2,
    lty = 2,
    lwd = lwd
  )
  lines(
    log(rowMeans(do.call(cbind, r[[3]]), na.rm = T)),
    type = type,
    pch = 3,
    col = 3,
    lty = 3,
    lwd = lwd
  )

  par(op)

  dev.off()
}






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
