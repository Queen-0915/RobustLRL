#' @export
createdir <- function(path, recursive = T) {
  ifelse(!dir.exists(path), dir.create(path, recursive = recursive), FALSE)
}

# for `foreach`
#' @export
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#' @export
computeError <- function(x, y, type = "F") {
  stopifnot(is.numeric(x) && is.numeric(y))
  if((length(x) != length(y)) && is.matrix(x)) {
    m <- ncol(x)
    sqrt(norm(x - matrix(rep(y, ncol(x)), ncol = ncol(x)), type = type)^2/m)
  } else {
    stopifnot(length(x)==length(y))
    return(drop(norm(matrix(x - y), type =  type)))
  }
}

#' @export
compute_rank <- function(B, cutoff = 1e-1) {
  return(sum(svd(B)$d>cutoff))
}
