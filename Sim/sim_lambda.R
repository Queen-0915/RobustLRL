# ============================================================== #
# Simulation: Sensitivity Study of Tuning Parameter under Distributed BIC
# ============================================================== #

simulation_name <- "lambda"
#rm(list = ls())

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
#createdir(fig_dir)
#createdir(LOG_dir)
if (Platform == "Linux") {
  # Nreps <- 100
  Nreps <- 100
  T_outer <- 1
  noise_type_arr <- c("Cauchy", "Normal", "T2")
  #noise_type_arr <- "Cauchy"
  # noise_type_arr <- c("Normal")
  registerDoFuture()
  # use multiple workers to accelerate the time of replication, change
  # the number 123 to smaller number, e.g., 4 or 8 to accommodate your machine.
  # plan(multisession, workers = 100)    ## on MS Windows
  plan(multisession, workers = 100)     ## on Linux, Solaris, and macOS
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

m <- 20 # the number of machines
n <- 2e2 # local sample size
N <- m * n # sample size
p <- 10 # row dimension
q <- 10 # column dimension
r <- 3 # rank
# r <- 3 # rank
pc <- .3 # the connection probability
tau <- 1 / 2 # quatile level
rho <- .1
# sigma2 <- 1
sigma2 <- 1
ishomo <- T
#hetercase <- 1
c0 <- 0.045
# c0 <- 0.04
# c0 <- 0.045
# c0 <- 0.013
# c0 <- 0.034
# c0*r^2*(p + q)*log(N)/n
# tau_penalty_factor <-  1
# tau_penalty_factor <-  1 / 12
tau_penalty_factor <-  1 / 6
# tau_penalty_factor <-  1 / 3
T_inner_arr <- 80 #c(40, 80, 120)
# T_inner_arr <- c(30, 60, 90)
# nlambda <- 100
# nlambda <- 100
# # lambda_factor = 1e-4
# lambda_factor <- 1e-3
# lambda_max <- .1
quiet = F
MAXIT <- 2e3
eps <- 1e-3
result <- vector(mode = "list", length = length(noise_type_arr))


# ============================================================== #
# MAIN ROUNTINE
# ============================================================== #
# 定义参数
nlambda <- 100
lambda_min <- 1e-4
lambda_max <- 1e-1

# 使用seq函数生成对数网格的指数幂次值
exponent <- seq(log10(lambda_min), log10(lambda_max), length.out = nlambda)

# 计算对应的lambda值
lambda_array <- 10^exponent

t_start <- tictoc::tic()
for (inoise_type_arr in 1:length(noise_type_arr)) { # inoise_type_arr <- 1
  noise_type <- noise_type_arr[inoise_type_arr]
  cat("The noise distribution is ", noise_type)

  set.seed(2022) # fix the seed
  r_out <- foreach(
    iNreps = 1:Nreps #,
    # .init = replicate(length(lambda_array), list()),
    # .combine = "comb"
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
      ishomo = ishomo
    )
    X <- data$X
    y <- data$y
    B <- data$B
    betaT <- matrix(as.numeric(B), p * q, 1)
    graph <- data$graph
    adjacency_matrix <- as.matrix(as_adjacency_matrix(graph))


    # c("Cauchy", "Normal", "T2")
    if(noise_type == 'Cauchy'){ f0 <- dcauchy(0, location = 0, scale = 1)}
    if(noise_type == 'Normal'){ f0 <- dnorm(0)}
    if(noise_type == "T2"){ f0 <- dt(0, 2)}


    # Local estimate by BIC
    B_init_local <- matrix(0, p, q)
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
    #error_init <- computeError(B_init, betaT)


    output_list <- list() #vector(mode = "list", length = length(lambda_array))#matrix(NA, nrow = nlambda,  ncol = 3)
      #vector(mode = "list", length = length(lambda_array))
    for (iT_lambda in 1:length(lambda_array)) {
      # Estimate
      lambda_i <- lambda_array[iT_lambda]

     # tic()
      out_lambda <- decentralizedTraceQR_lambda_cpp(
        X,
        y,
        d1 = p,
        d2 = q,
        adjacency_matrix,
        B_init,
        betaT,
        s  = r,
        tau = tau,
        T_outer = 1,
        T_inner = 80,
        c0 = c0,
        tau_penalty_factor = tau_penalty_factor,
        f0 = f0,
        lambda = lambda_i,
        quiet = quiet
      )
     # toc()
      error_deMQR <-
        out_lambda$history$errors_inner[length(out_lambda$history$errors_inner)]
      rank_deMQR <-
        mean(apply(
          out_lambda$B,
          2,
          FUN = function(x)
            compute_rank(matrix(x, p, q), cutoff = 1e-1)
        ))
      BIC <- out_lambda$BIC

      output_list[[iT_lambda]] <- c(lambda_i, error_deMQR, rank_deMQR, BIC)
    }
    output_list
  }

  result[[inoise_type_arr]] <- r_out

}
t_end <- toc(t_start)




# ============================================================== #
# PLOT
# ============================================================== #

plot_type <- c('Error', 'Rank', 'BIC')


out_error <- matrix(NA, nrow = length(lambda_array), ncol = length(noise_type_arr))
for (inoise_type_arr in 1:length(noise_type_arr)) { # inoise_type_arr <- 1
  noise_type <- noise_type_arr[inoise_type_arr]
  cat("The noise type is ", noise_type, "\n")

  r <- result[[inoise_type_arr]]
  ll <- matrix(0, nrow = length(lambda_array), ncol = 1)
  for(i in 1:length(lambda_array)){ # i <- 1
    for (j in 1:100){ # j <- 1
      ll[i,] <- ll[i,] + do.call(rbind, r[[j]])[ i, 2]
    }
  }
  ll_lambda <- ll/length(lambda_array)
  out_error[,inoise_type_arr] <- ll_lambda
}


out_rank <- matrix(NA, nrow = length(lambda_array), ncol = length(noise_type_arr))
for (inoise_type_arr in 1:length(noise_type_arr)) { # inoise_type_arr <- 1
  noise_type <- noise_type_arr[inoise_type_arr]
  cat("The noise type is ", noise_type, "\n")

  r <- result[[inoise_type_arr]]
  ll <- matrix(0, nrow = length(lambda_array), ncol = 1)
  for(i in 1:length(lambda_array)){ # i <- 1
    for (j in 1:100){ # j <- 1
      ll[i,] <- ll[i,] + do.call(rbind, r[[j]])[ i, 3]
    }
  }
  ll_lambda <- ll/length(lambda_array)
  out_rank[,inoise_type_arr] <- ll_lambda
}


out_bic <- matrix(NA, nrow = length(lambda_array), ncol = length(noise_type_arr))
for (inoise_type_arr in 1:length(noise_type_arr)) { # inoise_type_arr <- 1
  noise_type <- noise_type_arr[inoise_type_arr]
  cat("The noise type is ", noise_type, "\n")

  r <- result[[inoise_type_arr]]
  ll <- matrix(0, nrow = length(lambda_array), ncol = 1)
  for(i in 1:length(lambda_array)){ # i <- 1
    for (j in 1:100){ # j <- 1
      ll[i,] <- ll[i,] + do.call(rbind, r[[j]])[ i, 4]
    }
  }
  ll_lambda <- ll/length(lambda_array)
  out_bic[,inoise_type_arr] <- ll_lambda
}


out_list_plot <- list(out_error, out_rank, out_bic )

for(iplot_type in 1:length(plot_type)){ #iplot_type <- 3
  pdf(paste0(fig_dir, "/fig_lambda_", plot_type[iplot_type], ".pdf"))
  cat("The noise distribution is ", plot_type[iplot_type], "\n")

  outcome_plot <- data.frame(out_list_plot[[iplot_type]])
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
    outcome_plot[,1],
    type = type,
    ylim = c(
      min(
        outcome_plot[,1],
        outcome_plot[,2],
        outcome_plot[,3],
        na.rm = T
      ),
      max(
        outcome_plot[,1],
        outcome_plot[,2],
        outcome_plot[,3],
        na.rm = T
      )
    ),
    x = lambda_array, #1:length(lambda_array),
    lwd = lwd,
    xlab = "",
    ylab = ""
  )
  lines(
    outcome_plot[,2],
    type = type,
    pch = 2,
    col = 2,
    lty = 2,
    lwd = lwd
  )
  lines(
    outcome_plot[,3],
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
  "/sim_lambda",
  #simulation_name,
  "_",
  format(Sys.time(), "%Y%m%d%H%M%S"),
  ".RData"
))

#
# error <- do.call(rbind, result[[1]][[1]])[,2]
#
# plot(x = 1:length(lambda_array), y = error)

