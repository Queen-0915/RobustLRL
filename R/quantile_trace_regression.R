#' @title cv Quantile Trace Regression
#' @param gamma the step size
#' @export
bic.quantile_trace_regression <- function(X,
                                          y,
                                          tau,
                                          nlambda,
                                          B_init,
                                          B,
                                          gamma = 0.1,
                                          MAXIT = 1e4L,
                                          eps = 1e-2,
                                          quiet = T) {
  p <- nrow(B_init)
  q <- nrow(B_init)
  n <- length(y)
  lambda_vec <- exp(seq(log(1e-6), log(1), length.out = nlambda))*sqrt((p + q)*log(n)/n)
  bic_vec <- rep(NA, nlambda)
  for(i in 1:nlambda) {
    lambda <- lambda_vec[i]
    if(!quiet) {
    cat("lambda = ", lambda, "\n")}
    out_local <- mQR(
      X = X,
      y = y,
      tau = tau,
      lambda = lambda,
      B_init = B_init,
      B = B,
      eps = eps,
      gamma = gamma,
      MAXIT = MAXIT,
      quiet = quiet
    )
    if(!quiet) {
    cat(norm(out_local$B - B, "F"), "\t")}
    bic_vec[i] <- bic(X, y, out_local$B, sum(svd(out_local$B)$d>1e-4))
  }
  lambda <- lambda_vec[which.min(bic_vec)]
  out_local <- quantile_trace_regression(
    X = X,
    y = y,
    tau = tau,
    lambda = lambda,
    B_init = B_init,
    B = B,
    eps = eps,
    gamma = gamma,
    MAXIT = MAXIT,
    quiet = quiet
  )
  return(list(
    "B" = out_local$B,
    "real_error" = out_local$real_error,
    "lambda" = lambda,
    "lambda_vec" = lambda_vec,
    "bic_vec" = bic_vec
  ))
}

#' @title Quantile Trace Regression
#' @param gamma the step size
#' @export
quantile_trace_regression <- function(X,
                                      y,
                                      tau,
                                      lambda,
                                      B_init,
                                      B,
                                      gamma = 0.1,
                                      MAXIT = 1e4L,
                                      eps = 1e-2,
                                      quiet = T) {
  stopifnot(MAXIT >= 1)
  p <- nrow(B_init)
  q <- nrow(B_init)
  n <- length(y)
  B_prev <- B_init
  real_error = rep(NA, MAXIT)
  for (t in 1:MAXIT) {
    B_cur <- SVT(B_prev, gamma*lambda)
    gradient_vec <- t(X)%*%((y - X%*%c(B_cur))<= 0 - tau)/n
    # B_cur <- B_cur - gamma*1/t*matrix(gradient_vec, p, q)
    B_cur <- B_cur - gamma*matrix(gradient_vec, p, q)
    # B_cur <- (t(B_cur) + B_cur)/2
    if(!quiet) {
      cat(paste0("t = ", t, "relerror = ", norm(B_cur - B_prev, "F")/norm(B_prev, "F"), "\n"))
    }

    real_error[t] <- norm(B_cur - B, "F")
    if(norm(B_cur - B_prev, "F")/norm(B_prev, "F") < eps) {
      return(list(
        "B" = B_cur,
        "real_error" = na.omit(real_error)
      ))
    } else {
      B_prev <- B_cur
    }
  }
  return(list(
    "B" = B_cur,
    "real_error" = na.omit(real_error)
  ))
}



SVT <-  function(B, lambda) {
  svdProxJ <- svd(B)
  svdu <- svdProxJ$u
  svdd <- svdProxJ$d
  svdv <- svdProxJ$v

  nzD = sum(pmax(svdd - lambda, 0) > 1e-10)

  if (nzD > 1) {
    svtx <- svdu[, seq(nzD)] %*% (pmax(svdd - lambda,
                                       0)[seq(nzD)] * t(svdv[, seq(nzD)]))
  } else {
    svtx <- pmax(svdd - lambda, 0)[1] * svdu[, 1] %*%
      t(svdv[, 1])
  }

  return(svtx)

}


# Backtrcking line search
backtracking <- function(Y,
                         X,
                         parameters,
                         grad,
                         grad.g,
                         lambda,
                         step_size,
                         constant) {
  # list2env(Grads, envir = environment()) list2env(grad, envir
  # = environment()) list2env(grad.g, envir = environment())

  while (Obj.func(Y,
                  X,
                  beta = parameters - step_size * grad.g,
                  step_size = step_size) * step_size > Obj.func(Y, X, parameters,
                                                                step_size = step_size) * step_size - step_size * t(grad) %*%
         grad.g + (step_size / 2) * tcrossprod(t(grad.g))) {
    step_size <- constant * step_size
    grad.g <- (parameters - S(lambda, (parameters - step_size *
                                         grad))) / step_size
  }
  return(step_size)
}
