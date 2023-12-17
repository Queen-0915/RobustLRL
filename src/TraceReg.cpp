#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include "decentralizedTraceQR.h"


//' @export
// [[Rcpp::export]]
Rcpp::List TraceReg(arma::mat &X,
               arma::vec &y,
               double lambda,
               arma::mat &B_init,
               arma::mat &B,
               double gamma = 0.1,
               int MAXIT = 1e4L,
               double eps = 1e-2,
               bool quiet = true) {
  int p = B_init.n_rows, q = B_init.n_cols, n = y.n_elem;

  arma::mat B_prev = B_init, B_cur = B_init;
  arma::vec real_error = arma::zeros(MAXIT);
  arma::vec gradient_vec(p*q);

  int t;
  for(t = 0; t < MAXIT; t++) {
    Rcpp::checkUserInterrupt();

    B_cur = SVT_cpp(B_prev, gamma*lambda);
    gradient_vec = -X.t()*(y - X*arma::vectorise(B_cur))/n;
    B_cur = B_cur - gamma*arma::reshape(gradient_vec, p, q);
    if(!quiet) {
      Rcpp::Rcout<<"t = "<<t<<", relerror = " <<
        arma::norm(B_cur - B_prev, "fro")/arma::norm(B_prev, "fro") << "\n";
    }

    real_error[t] = arma::norm(B_cur - B, "fro");

    // Check convergence
    if(arma::norm(B_cur - B_prev, "fro") < eps*arma::norm(B_prev, "fro")) {
      return Rcpp::List::create(
        Rcpp::Named("B") = B_cur,
        Rcpp::Named("real_error") = real_error.subvec(0, t)
      );
    } else {
      B_prev = B_cur;
    }

  }
  return Rcpp::List::create(
    Rcpp::Named("B") = B_cur,
    Rcpp::Named("real_error") = real_error
  );


}
