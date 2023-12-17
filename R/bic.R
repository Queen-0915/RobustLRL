#' @export
bic <- function(X, y, B, r) {
  N <- length(y)
  bic_value <- quantile_lossCPP(y - X%*%c(B), tau) / N  + log(N)*log(log(N))/N*r
  return(bic_value)
}


#' @export
bic_traceReg <- function(X, y, B, r) {
  N <- length(y)
  bic_value <- norm(matrix(y - X%*%c(B)), "F")^2 / N  + log(N)*log(log(N))/N*r
  return(bic_value)
}
