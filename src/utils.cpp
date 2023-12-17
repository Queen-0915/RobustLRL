#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "decentralizedTraceQR.h"

//' RcppArmadillo implementation of the R pmax()
 //' function
 //'
 //' @title RcppArmadillo implemenation of the r pmax function
 //' @param x the (aramdillo) vector
 //' @param bound the (double) bound
 //'
 //' @export
 // [[Rcpp::export]]
 arma::vec pmax_arma(arma::vec x, double bound)
 {
   int length_x = x.n_elem;
   for (int i = 0; i < length_x; i++)
   {
     if (x(i) < bound)
       x(i) = bound;
   }
   return x;
 }

 //' Kernel fitting
 //'
 //' Conduct the Nadaraya-Watson kernel regression to estimate the density of the
 //' errs.
 //'
 //' @param E A n-by-K residual matrix, where n is sample size and K is the
 //' quantile levels
 //' @param h A n-by-K bandwidth
 //'
 //' @return
 //' A K-by-1 density estimate vector.
 //' @examples
 //' set.seed(1)
 //' library(pracma)
 //' p <- 600
 //' n <- 100
 //' case <- "cauchy"
 //' betaT <- rep(0, p)
 //' betaT[1] <- 3
 //' betaT[2] <- 1.5
 //' betaT[5] <- 2
 //' s <- length(which(abs(betaT)>0))
 //' rho <- 0.5
 //' data <- gen_data(n, p, betaT, case, rho)
 //' x <- data$x
 //' y <- data$y
 //' data_val <- gen_data(n, p, betaT, case, rho)
 //' x_val <- data_val$x
 //' y_val <- data_val$y
 //' K <- 19
 //' tau <- 0.05 * (1:K)
 //' m <- qraenet(x = x, y = y, tau = 0.5, standardize = FALSE,
 //' intercept = FALSE, lambda2 = 0, sigma = 1, method = "padmm")
 //' nlambda <- length(m$lambda)
 //' validation_err <- rep(NA, nlambda)
 //' for(ivalidation_err in 1:nlambda) {
 //'   validation_err[ivalidation_err] <-
 //'     sum((y_val - x_val%*%m$beta[,ivalidation_err])^2)
 //' }
 //' beta_init <-  matrix(as.numeric(m$beta[,which.min(validation_err)]))
 //' alpha_init <- quantile(y - x%*%beta_init, tau)
 //' E <- repmat(y - x%*%beta_init, 1, K) - repmat(alpha_init, n, 1);
 //' h <- rep(sqrt(s*log(n)/n), K)
 //' fh <- kernel(E, h)
 //' @export
 // [[Rcpp::export]]
 arma::mat kernel(arma::mat E, arma::vec h)
 {
   arma::mat f(arma::size(E));
   E = E / arma::repmat(h.t(), E.n_rows, 1);
   // E.each_row() /= h.t();
   f = (-315 * arma::pow(E, 6) + 735 * arma::pow(E, 4) -
     525 * arma::pow(E, 2) + 105) /
       64 % (arma::abs(E) < 1);
   return arma::trans(arma::mean(f, 0)) / h;
 }

 //' @title Mode of a given vector
 //' @description
 //' Find the mode of a given numerical vector.
 //' @name mod
 //' @param v A vector of real numbers
 //' @return the mode of the vector
 //' @examples
 //' v <- c(rep(1,10),rep(2,9))
 //' mod(v)
 //' @export
 // [[Rcpp::export]]
 double mod(arma::vec v)
 {
   std::map<int, size_t> frequencyCount;
   using pair_type = decltype(frequencyCount)::value_type;

   for (auto i : v)
     frequencyCount[i]++;

   auto pr = std::max_element(
     std::begin(frequencyCount), std::end(frequencyCount),
     [](const pair_type &p1, const pair_type &p2)
     {
       return p1.second < p2.second;
     });
   return pr->first;
 }

 //' @export
 // [[Rcpp::export]]
 double quantile_lossCPP(arma::vec x, double tau)
 {
   // return arma::accu(x % (tau - (x < 0))); // why this is wrong?
   // return arma::accu(x % ((x < 0) - tau)); // shall be this?
   return arma::accu((arma::abs(x) + (2 * tau - 1) * x) * 0.5);
 }

 //' @export
 // [[Rcpp::export]]
 double cqr_loss_cpp(arma::mat X, arma::vec y, arma::vec beta, arma::vec alpha, arma::vec &tau)
 {
   double s = 0;
   for (int i = 0; i < tau.n_elem; i++)
   {
     s += quantile_lossCPP(y - X * beta - alpha(i), tau(i));
   }
   return s;
 }

 //' @export
 // [[Rcpp::export]]
 arma::uvec calN_j_cpp(int n, int j)
 {
   return arma::regspace<arma::uvec>(n * (j), 1, n * (j + 1) - 1);
 }


 // arma::uvec calN_j_cpp2(int n, int n1, int m)
 // {
 //   return arma::regspace<arma::uvec>(n * m, 1, n * m + n1  - 1);
 // }

 //' @export
 // [[Rcpp::export]]
 arma::vec soft_thresholding_cpp(arma::vec x, double t)
 {
   return pmax_arma(x - t, 0) - pmax_arma(-x - t, 0);
 }

 //' @export
 // [[Rcpp::export]]
 arma::mat SVT_cpp(arma::mat &B, double lambda)
 {
   arma::mat U, V;

   arma::vec s;

   arma::svd(U, s, V, B);

   s = pmax_arma(s - lambda, 0);

   int nzD = arma::accu(s > 1e-10);

   if (nzD > 1) {
     return U.cols(0, nzD - 1) * arma::diagmat(s.head(nzD)) *
       arma::trans(V.cols(0, nzD - 1));
   } else {
     return U.col(0) * s.head(1) *
       arma::trans(V.col(0)); //matrix
   }
 }


 //' @export
 // [[Rcpp::export]]
 arma::vec SVT_vec_cpp(arma::vec &beta, int d1, int d2, double lambda)
 {
   arma::mat B = arma::reshape(beta, d1, d2);
   arma::mat U, V;

   arma::vec s;

   arma::svd(U, s, V, B);

   s = pmax_arma(s - lambda, 0);

   int nzD = arma::accu(s > 1e-10);

   if (nzD > 1) {
     return arma::vectorise(U.cols(0, nzD - 1) * arma::diagmat(s.head(nzD)) *
                            arma::trans(V.cols(0, nzD - 1)));
   } else {
     return arma::vectorise(U.col(0) * arma::diagmat(s.head(1)) *
                            arma::trans(V.col(0))); //vector
   }
   // return arma::vectorise(U.cols(0, nzD - 1) * arma::diagmat(s.head(nzD)) *
   //                        arma::trans(V.cols(0, nzD - 1)));
 }
