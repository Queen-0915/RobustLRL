#' @title cv Trace Regression
#' @param gamma the step size
#' @export
bic.trace_regression <- function(X,
                                 y,
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
    out_local <- TraceReg(
      X = X,
      y = y,
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
    bic_vec[i] <- bic_traceReg(X, y, out_local$B, sum(svd(out_local$B)$d>1e-4))
  }
  lambda <- lambda_vec[which.min(bic_vec)]
  out_local <- TraceReg(
    X = X,
    y = y,
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
