// #include <RcppArmadillo.h>
// // [[Rcpp::depends(RcppArmadillo)]]
//
// #include <cmath>
// #include "decentralizedCQR.h"
//
// //' @export
// // [[Rcpp::export]]
// Rcpp::List SubGD_single(arma::mat &X,
//                         arma::vec &y,
//                         arma::mat &b_init,
//                         arma::vec &betaT,
//                         arma::mat &tau_K_mat, // quantile level
//                         int T = 20,
//                         double a_lambda = 3.5, // for stepsize
//                         double b_lambda = 0.51, // for stepsize
//                         double lambda = 1e-2, // penalty parameter
//                         double eps = 1e-4,
//                         bool quiet = true) {
//
//   int p = X.n_cols,
//     N = X.n_rows;
//   int K = tau_K_mat.n_cols;
//   arma::vec g_mat(p); // subgradient matrix
//   arma::vec b(p), b_old(p);
//   arma::vec errors(T, arma::fill::zeros);
//   double eta;
//   arma::umat E(N, K, arma::fill::zeros);
//
//
//   b = b_init;
//   int t, j;
//   for(t = 0; t < T; t++) { // Iteration
//     // Rcpp::Rcout << "t = " << t << "\n";
//     b_old = b;
//     // Computation
//     // cqr loss
//     E.each_col() = ((y - X*b_old)<=0);
//     g_mat = 1.0/N * arma::mean(arma::trans(X)*(E - tau_K_mat), 1) +
//       lambda*arma::sign(b_old); // optimize repmat
//
//     // Estimation
//     eta = std::pow(a_lambda/(t + 1), b_lambda);
//     b = b_old - eta*g_mat;
//
//     errors(t) = std::sqrt(std::pow(arma::norm(b - betaT, "fro"), 2));
//     if((t>0) && (std::abs(errors(t) - errors(t - 1))<eps)) {
//       break;
//     }
//   }
//
//   if ((t + 1) <= (T - 1)) {
//     errors.shed_rows(t + 1, T - 1);
//   }
//   return Rcpp::List::create(
//     Rcpp::Named("b") = b,
//     Rcpp::Named("history") = Rcpp::List::create(Rcpp::Named("errors") = errors));
//
// }
//
//
// //' @export
// // [[Rcpp::export]]
// Rcpp::List SubGD(arma::mat &X,
//                  arma::vec &y,
//                  int d1,
//                  int d2,
//                  arma::vec &b_init,
//                  arma::vec &betaT,
//                  double tau = 1/2,
//                  int T = 20,
//                  double a_lambda = 3.5, // for stepsize
//                  double b_lambda = 0.51, // for stepsize
//                  int nlambda = 100L,
//                  double lambda_factor = 1e-4,
//                  double lambda_max = 1,
//                  double eps = 1e-4,
//                  bool quiet = true) {
//
//
//
//   int N = X.n_rows;
//   arma::vec b;
//
//   double lambda;
//   arma::vec bic_array(nlambda - 1);
//   arma::vec lambda_array(nlambda);
//   Rcpp::List L;
//
//   lambda_max = std::min(
//     arma::norm(arma::reshape(X.t() * y / N, d1, d2)),
//     lambda_max);
//   lambda_array = arma::exp(arma::linspace(std::log(1), std::log(lambda_factor), nlambda));
//   lambda_array.shed_row(0);
//   lambda_array *= lambda_max;
//
//   if (!quiet)
//   {
//     Rcpp::Rcout << "lambda_max = " << lambda_max << std::endl;
//   }
//   bic_array = arma::zeros(arma::size(lambda_array));
//   arma::vec shat_array = arma::zeros(arma::size(lambda_array));
//
//   // BIC
//   int shat;
//   for (int ilambda = 0; ilambda < lambda_array.n_elem; ilambda++) {
//     Rcpp::checkUserInterrupt(); // checking interruption
//     lambda = lambda_array(ilambda);
//     L = SubGD_single(X,
//                      y,
//                      b_init,
//                      betaT,
//                      tau,
//                      T,
//                      a_lambda,
//                      b_lambda,
//                      lambda,
//                      eps,
//                      quiet);
//     b = Rcpp::as<arma::vec>(L["b"]);
//     // bic
//     bic_array(ilambda) +=  std::pow(arma::norm(y - X*b, "fro"), 2);
//     shat = int(arma::as_scalar(arma::accu(b != 0)));
//     bic_array(ilambda) = bic_array(ilambda) / N + std::log(N) * std::log(std::log(N)) / N * std::max(shat, 1);
//   }
//
//
//   arma::uword ilambda = arma::index_min(bic_array);
//   lambda = lambda_array(ilambda);
//
//
//   if (!quiet)
//   {
//     Rcpp::Rcout << "lambda = " << lambda << std::endl;
//   }
//
//   L = SubGD_single(X,
//                    y,
//                    b_init,
//                    betaT,
//                    tau_K_mat,
//                    T,
//                    a_lambda,
//                    b_lambda,
//                    lambda,
//                    eps,
//                    quiet);
//   return L;
// }
