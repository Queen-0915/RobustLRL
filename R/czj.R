n <- 1e5
Y <- sapply(1:n, function(x) rbinom(1, 1, 1/x))
X <- rnorm(n)
S <- sum((1:n)*Y*X)
eps <- 1e-6
abs(S)/n/(log(log(n)))^eps
