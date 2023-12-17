#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include "decentralizedTraceQR.h"

//' @export
// [[Rcpp::export]]
Rcpp::List decentralizedTraceQR_inner(arma::cube &X,
                                      arma::vec &yt,
                                      arma::mat &adjacency_matrix,
                                      arma::vec &rho,
                                      arma::vec &omega,
                                      arma::cube &B_outer,
                                      double tau_penalty,
                                      double lambda = 0.015,
                                      int T_inner = 20)
{
  int p = B_outer.n_rows,
    q = B_outer.n_cols,
    m = B_outer.n_slices;
  int N = yt.n_elem;
  int n = N / m;
  arma::cube B_inner(p, q, m),
  B_inner_old(p, q, m), P(p, q, m, arma::fill::zeros);
  arma::uvec idx(n);
  arma::mat tmp(p, q);
  arma::mat grad(p, q);

  B_inner = B_outer;
  for (int t = 0; t < T_inner; t++)
  {
    B_inner_old = B_inner;
    for (int j = 0; j < m; j++)
    {
      grad.fill(0); // reset the gradient
      idx = calN_j_cpp(n, j);
      arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));
      P.slice(j) -= tau_penalty * arma::sum(B_inner_old.each_slice(neighbors_j) -
                                                B_inner_old.slice(j),
                                            2);
      for (int i : idx)
      {
        grad += X.slice(i) * (arma::trace(arma::trans(X.slice(i)) *
          B_inner_old.slice(j)) - yt(i)) / m / n;
      }
      tmp = rho(j) * B_inner_old.slice(j) - grad - P.slice(j);
      tmp += tau_penalty * arma::sum(B_inner_old.each_slice(neighbors_j) +
        B_inner_old.slice(j), 2);
      tmp *= omega(j);
      B_inner.slice(j) = SVT_cpp(tmp, lambda * omega(j));
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("B_inner") = B_inner);
}
