#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include "decentralizedTraceQR.h"

//' @export
// [[Rcpp::export]]
Rcpp::List deSubGD_single(arma::mat &X,
                          arma::vec &y,
                          int d1,
                          int d2,
                          arma::mat &adjacency_matrix,
                          arma::mat &B_init,
                          arma::vec &betaT,
                          double tau, // quantile level
                          int T = 20,
                          double a = 3.5, // for stepsize
                          double b = 0.51, // for stepsize
                          double lambda = 1e-2, // penalty parameter
                          double eps = 1e-2,
                          bool quiet = true) {

  int m = adjacency_matrix.n_cols,
    p = X.n_cols,
    N = X.n_rows,
    n = N / m;
  arma::mat g_mat(p, m); // subgradient matrix
  arma::mat B(size(B_init)), B_old(size(B_init));
  arma::mat Phi(size(B_init));

  arma::vec errors(T, arma::fill::zeros);
  arma::uvec idx(n);
  double eta;
  arma::uvec E(n, 1, arma::fill::zeros);

  arma::mat C(m, m, arma::fill::zeros);
  arma::uvec n_arr(m, arma::fill::zeros);
  for(int j = 0; j < m; j++) {
    arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));
    n_arr(j) = neighbors_j.n_elem;
  }
  for(int j = 0; j < m; j++) {
    arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));
    for(int l = 0; l<neighbors_j.n_elem; l++) {
      C(neighbors_j(l),j) = 1.0/std::max(n_arr(j), n_arr(neighbors_j(l)));
    }
  }
  for(int j = 0; j<m;j++) {
    arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));
    arma::vec tmp_C = C.col(j);
    C(j,j) = 1 - arma::sum(tmp_C(neighbors_j));
  }
  // if(!quiet) {
  //   Rcpp::Rcout << C << std::endl;
  // }

  B = B_init;
  int t, j, k;
  for(t = 0; t < T; t++) { // Iteration
    // Rcpp::Rcout << "t = " << t << "\n";
    // Computation
    for(j = 0; j < m; j++) {
      idx = calN_j_cpp(n, j);
      // quantile loss
      // g_mat.col(j) = 1.0/N * arma::trans(X.rows(idx)) *((y(idx) - X.rows(idx)*B.col(j))<=0 - tau) +
      //   lambda*arma::sign(B.col(j));

      // cqr loss
      E = ((y(idx) - X.rows(idx)*B.col(j))<=0);
      // E.each_col() = ((y(idx) - X.rows(idx)*B.col(j))<=0);
      // g_mat.col(j) = 1.0/N * arma::mean(arma::trans(X.rows(idx))*(E - tau), 1) +
      //   lambda/1.0/m*arma::sign(B.col(j)); // optimize repmat
      arma::mat U_tmp, V_tmp;
      arma::vec s_tmp;
      arma::svd(U_tmp, s_tmp, V_tmp, arma::reshape(B.col(j), d1, d2));
      g_mat.col(j) = 1.0/N * arma::mean(arma::trans(X.rows(idx))*(E - tau), 1) +
        lambda/1.0/m*arma::vectorise(U_tmp*V_tmp.t()); // optimize repmat
    }

    // Local Estimation
    eta = a*std::pow(1.0/(t + 1), b);
    // for(j = 0; j < m; j++) {
    //   Phi.col(j) = B.col(j) - eta*g_mat.col(j);
    // }
    Phi = B - eta*g_mat;

    // Information Exchange and Combination
    // for(j = 0; j < m; j++) {
    //   arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));
    //   B.col(j) = sum(B.cols(neighbors_j), 1)/neighbors_j.n_elem;
    // }
    B_old = B;
    B = Phi*C;

    errors(t) = std::sqrt(std::pow(arma::norm(B - B_old, "fro"), 2) / m);
    if((t>0) && ((std::abs(errors(t) - errors(t - 1)))/errors(t - 1)<eps)) {
      break;
    }
  }

  if ((t + 1) <= (T - 1)) {
    errors.shed_rows(t + 1, T - 1);
  }
  return Rcpp::List::create(
    Rcpp::Named("B") = B,
    Rcpp::Named("history") = Rcpp::List::create(Rcpp::Named("errors") = errors));

}


//' @export
// [[Rcpp::export]]
Rcpp::List deSubGD(arma::mat &X,
                   arma::vec &y,
                   int d1, int d2,
                   arma::mat &adjacency_matrix,
                   arma::mat &B_init,
                   arma::vec &betaT,
                   int T = 20,
                   double a = 3.5, // for stepsize
                   double b = 0.51, // for stepsize
                   double tau = 0.5, // quantile level
                   int nlambda = 100L,
                   double lambda_factor = 1e-4,
                   double lambda_max = 1,
                   double eps = 1e-2,
                   bool quiet = true) {



  int m = adjacency_matrix.n_cols,
    N = X.n_rows,
    n = N / m;
  arma::mat B;
  arma::uvec idx(n);

  double lambda;
  arma::vec bic_array(nlambda - 1);
  arma::vec lambda_array(nlambda);
  Rcpp::List L;

  lambda_max = std::min(
    arma::max(arma::abs(X.t() * y / N)),
    lambda_max);
  lambda_array = arma::exp(arma::linspace(std::log(1),
                                          std::log(lambda_factor), nlambda));
  lambda_array.shed_row(0);
  lambda_array *= lambda_max;

  if (!quiet)
  {
    Rcpp::Rcout << "lambda_max = " << lambda_max << std::endl;
  }
  bic_array = arma::zeros(arma::size(lambda_array));
  arma::vec shat_array = arma::zeros(arma::size(lambda_array));

  // BIC
  int shat;
  for (int ilambda = 0; ilambda < lambda_array.n_elem; ilambda++) {
    Rcpp::checkUserInterrupt(); // checking interruption
    lambda = lambda_array(ilambda);
    L = deSubGD_single(X,
                       y,
                       d1, d2,
                       adjacency_matrix,
                       B_init,
                       betaT,
                       tau,
                       T,
                       a,
                       b,
                       lambda,
                       eps,
                       quiet);
    B = Rcpp::as<arma::mat>(L["B"]);

    // bic
    shat = 0;
    for (int j = 0; j < m; j++)
    {
      idx = calN_j_cpp(n, j);
      // bic_array(ilambda) +=  std::pow(arma::norm(y(idx) - X.rows(idx)*B.col(j), "fro"), 2);
      // bic_array(ilambda) += cqr_loss_cpp(X.rows(idx), y(idx), B.col(j), A.col(j), tau_K);
      bic_array(ilambda) += quantile_lossCPP(y(idx) - X.rows(idx) * B.col(j), tau);
      shat += int(arma::rank(arma::reshape(B.col(j), d1, d2)));
    }
    shat = 1.0 * shat / m;
    bic_array(ilambda) = bic_array(ilambda) / N +
      std::log(N) * std::log(std::log(N)) / N * std::max(shat, 1);
    shat_array(ilambda) = shat;
  }


  arma::uword ilambda = arma::index_min(bic_array);
  lambda = lambda_array(ilambda);


  if (!quiet)
  {
    Rcpp::Rcout << "lambda = " << lambda << std::endl;
  }

  L = deSubGD_single(X,
                     y,
                     d1, d2,
                     adjacency_matrix,
                     B_init,
                     betaT,
                     tau,
                     T,
                     a,
                     b,
                     lambda,
                     eps,
                     quiet);
  return L;
}
