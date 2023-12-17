#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include "decentralizedTraceQR.h"


//' @export
// [[Rcpp::export]]
Rcpp::List deTraceReg(arma::mat &X, arma::vec &y,
                                    int d1, int d2,
                                    arma::mat &adjacency_matrix,
                                    arma::mat &B_init,
                                    arma::vec &betaT,
                                    int T_inner = 20,
                                    double tau_penalty_factor = 1 / 6,
                                    int nlambda = 100L,
                                    double lambda_factor = 1e-4,
                                    double lambda_max = 1,
                                    bool quiet = true)
{

  // Parameters
  int m = adjacency_matrix.n_cols, p = d1*d2, N = X.n_rows;
  int n = N / m;
  arma::vec s_vec;
  s_vec = arma::eig_sym(arma::trans(X.rows(0, n - 1))*X.rows(0, n - 1));
  double tau_penalty_beta = s_vec.back() / n * tau_penalty_factor;

  // Results
  int K = 1;
  arma::vec errors_inner(T_inner, arma::fill::zeros);
  arma::mat B_inner(p, m, arma::fill::zeros),
  B_inner_old(p, m, arma::fill::zeros),
  P_beta(p, m, arma::fill::zeros);
  arma::uvec idx(n);

  // cache omega and rho
  arma::vec rho(m), omega(m);
  for (int j = 0; j < m; j++)
  {
    idx = calN_j_cpp(n, j);
    s_vec = arma::eig_sym(arma::trans(X.rows(idx))*X.rows(idx)/n);
    rho(j) = s_vec.back();
    omega(j) = 1 / (2 * tau_penalty_beta * arma::accu(adjacency_matrix.row(j)!=0) + rho(j));
  }

  // === Main routine === //
  // double lambda_max;
  arma::vec bic_array(nlambda - 1);
  arma::vec tmp(p);


    // ===== Tune parameter ===== //
    arma::vec lambda_array(nlambda);
    lambda_max = std::min(
      arma::norm(arma::reshape(X.t() * y / N, d1, d2)),
      lambda_max
    );
    // lambda_max = .5;
    lambda_array = arma::exp(arma::linspace(std::log(1),
                                            std::log(lambda_factor), nlambda));
    lambda_array.shed_row(0);
    lambda_array *= lambda_max;
    if(!quiet) {
      Rcpp::Rcout << "lambda_max = " << lambda_max << std::endl;
    }
    bic_array = arma::zeros(nlambda - 1, 1);
    arma::vec shat_array = arma::zeros(nlambda - 1, 1);

    double lambda;

    for (int ilambda = 0; ilambda < lambda_array.n_elem; ilambda++)
    {
      lambda = lambda_array(ilambda);

      B_inner = B_init;
      P_beta = arma::zeros(p, m);
      for (int t = 0; t < T_inner; t++)
      {
        B_inner_old = B_inner;
        for (int j = 0; j < m; j++)
        {
          idx = calN_j_cpp(n, j);
          arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));
          P_beta.col(j) = P_beta.col(j) + tau_penalty_beta *
            arma::sum(arma::repmat(B_inner_old.col(j), 1, neighbors_j.n_elem) -
            B_inner_old.cols(neighbors_j), 1);

          tmp = omega(j) * (rho(j) * B_inner_old.col(j) - 1.0 / m / n *
            arma::trans(X.rows(idx)) * (X.rows(idx) * B_inner_old.col(j) - y(idx)) -
            P_beta.col(j) + tau_penalty_beta *
            arma::sum(arma::repmat(B_inner_old.col(j), 1, neighbors_j.n_elem) +
            B_inner_old.cols(neighbors_j), 1));
          // Update B
          // B_inner.col(j) = soft_thresholding_cpp(tmp, lambda * omega(j));
          B_inner.col(j) = SVT_vec_cpp(tmp, d1, d2, lambda * omega(j));
        }
      }

      // bic
      int shat = 0;
      for (int j = 0; j < m; j++)
      {
        idx = calN_j_cpp(n, j);
        bic_array(ilambda) += std::pow(arma::norm(y(idx) -
          X.rows(idx)*B_inner.col(j), "fro"), 2);
        shat += int(arma::rank(arma::reshape(B_inner.col(j), d1, d2)));
        // shat = std::max(shat, int(arma::as_scalar(arma::accu(B_inner.col(j) != 0))));
        // shat += int(arma::as_scalar(arma::accu(B_inner.col(j) != 0)));
      }
      shat = 1.0*shat/m;

      // bic_array(ilambda) = std::log(bic_array(ilambda) / K) + std::log(N) * std::log(std::log(N)) / N * std::max(shat, 1); // not perform well under Cauchy noise; see gu2022
      bic_array(ilambda) = bic_array(ilambda) / N + std::log(N) * std::log(std::log(N)) / N * std::max(shat, 1);
      // bic_array(ilambda) = std::log(bic_array(ilambda) / K) + std::log(N) * std::log(N) / N * shat;
      shat_array(ilambda) = shat;

    }

    arma::uword ilambda = arma::index_min(bic_array);
    lambda = lambda_array(ilambda);

    // arma::uvec imode = arma::find(shat_array == mod(shat_array));
    // int ilambda = imode(imode.n_elem - 1); // get the smallest possible penalty parameter
    // lambda = lambda_array(ilambda);

    // lambda = 0.01488635;

    if (!quiet)
    {
      Rcpp::Rcout << "lambda = " << lambda << std::endl;
    }

    B_inner = B_init;
    P_beta = arma::zeros(p, m);
    for (int t = 0; t < T_inner; t++)
    {
      B_inner_old = B_inner;
      for (int j = 0; j < m; j++)
      {
        idx = calN_j_cpp(n, j);
        arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));
        P_beta.col(j) = P_beta.col(j) + tau_penalty_beta *
          arma::sum(arma::repmat(B_inner_old.col(j), 1, neighbors_j.n_elem) - B_inner_old.cols(neighbors_j), 1);

        tmp = omega(j) * (rho(j) * B_inner_old.col(j) - 1.0 / m / n *
          arma::trans(X.rows(idx)) * (X.rows(idx) * B_inner_old.col(j) - y(idx)) -
          P_beta.col(j) + tau_penalty_beta * arma::sum(arma::repmat(B_inner_old.col(j),
                                                       1, neighbors_j.n_elem) +
                                                         B_inner_old.cols(neighbors_j), 1));
        // Update B
        B_inner.col(j) = SVT_vec_cpp(tmp, d1, d2, lambda * omega(j));
      }
      errors_inner(t) = std::sqrt(std::pow(arma::norm(B_inner - arma::repmat(betaT, 1, m), "fro"), 2) / m);
      if (!quiet)
      {
        Rcpp::Rcout << errors_inner(t) << "\t";
      }
    }



  return Rcpp::List::create(
    Rcpp::Named("B") = B_inner,
    Rcpp::Named("history") = Rcpp::List::create(Rcpp::Named("errors_inner") = arma::vectorise(errors_inner)
    ));
}
