#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include "decentralizedTraceQR.h"

//' @export
// [[Rcpp::export]]
Rcpp::List decentralizedTraceQR_lambda_cpp(arma::mat &X, arma::vec &y,
                                            int d1, int d2,
                                            arma::mat &adjacency_matrix,
                                            arma::mat &B_init,
                                            arma::vec &betaT,
                                            int s,
                                            double tau = 0.5,
                                            int T_outer = 1,
                                            int T_inner = 20, double c0 = 0.1,
                                            double tau_penalty_factor = 1 / 6,
                                            double f0 = 1,
                                            double lambda = 0.05,
                                            bool quiet = true)
 {
   // int nlambda = 100L,
   // double lambda_factor = 1e-4,
   // Parameters
   int m = adjacency_matrix.n_cols, p = d1 * d2, N = X.n_rows;
   int n = N / m;
   arma::vec s_vec;
   // if(!quiet) {
   //   Rcpp::Rcout<<arma::trans(s_vec)<<std::endl;
   // }
   s_vec = arma::eig_sym(arma::trans(X.rows(0, n - 1)) * X.rows(0, n - 1));
   double tau_penalty_beta = s_vec.back() / n * tau_penalty_factor;

   // Results
   int K = 1;
   arma::mat errors_inner(T_inner, T_outer, arma::fill::zeros),
     B_out(p, m, arma::fill::zeros),
     B_inner(p, m, arma::fill::zeros),
     B_inner_old(p, m, arma::fill::zeros),
     P_beta(p, m, arma::fill::zeros);
   arma::vec yt(N, arma::fill::zeros);
   double fhat;
   double bicj = 0;
   arma::mat E(n, K);
   arma::uvec idx(n);
   // arma::mat shat_mat(nlambda - 1, T_outer);
   // arma::mat bic_mat(nlambda - 1, T_outer);
   // arma::vec lambda_max_array(T_outer);
   arma::mat Anew(K, m, arma::fill::zeros);

   // cache omega and rho
   arma::vec rho(m), omega(m);
   for (int j = 0; j < m; j++)
   {
     idx = calN_j_cpp(n, j);
     s_vec = arma::eig_sym(arma::trans(X.rows(idx)) * X.rows(idx) / n);
     rho(j) = s_vec.back();
     omega(j) = 1 / (2 * tau_penalty_beta * arma::accu(adjacency_matrix.row(j) != 0) + rho(j));
   }

   // === Main routine === //
   B_out = B_init;
   double hv;
   // double lambda_max;
   // arma::vec bic_array(nlambda - 1);
   arma::vec tmp(p);
   for (int v = 0; v < T_outer; v++)
   {
     Rcpp::checkUserInterrupt(); // checking interruption
     // ===== bandwidth ===== //
     // hv = std::sqrt(s * std::max(d1, d2)*std::log(d1 + d2) / N) +
     //   std::pow(s, -0.5) * std::pow((c0 * s * s * std::max(d1, d2)*std::log(d1 + d2) / n), (v + 1) * 1.0 / 2);
     hv = std::sqrt(s * (d1 + d2) * std::log(N) / N) +
          std::pow(s, -0.5) * std::pow((c0 * s * s * (d1 + d2) * std::log(N) / n),
                                       (v + 1) * 1.0 / 2);
     if (!quiet)
     {
       Rcpp::Rcout << "Outer iteration: v =" << v << ", hv =" << hv << "\n";
       // Rcpp::Rcout << "The density estimates at each node are:\n";
     }

     // ===== construct the pseudo-responses ===== //
     for (int j = 0; j < m; j++)
     {
       idx = calN_j_cpp(n, j);
       E = y(idx) - X.rows(idx) * B_out.col(j);
       fhat = arma::as_scalar(f0); // arma::as_scalar(kernel(E, hv * arma::ones(K, 1)));
       if (!quiet)
       {
         Rcpp::Rcout << fhat << std::endl;
       }
       yt(idx) = X.rows(idx) * B_out.col(j) - 1 / fhat * ((E <= 0) - tau);
     }

     // ===== Tune parameter ===== //
     // arma::vec lambda_array(nlambda);
     // lambda_max_array(v) = arma::norm(arma::reshape(X.t() * yt / N, d1, d2));
     // int nlambda = 100; // 选择lambda的数量
     // double lambda_min = 0.005;
     // double lambda_max = 0.5;

     // arma::vec lambda_array = arma::exp(arma::linspace(std::log10(lambda_min),
     //                                                   std::log10(lambda_max), nlambda));


     // lambda_max = std::min( //新注释掉的
     //   arma::norm(arma::reshape(X.t() * yt / N, d1, d2)),
     //   lambda_max);
     // lambda_max = .5;
     //
     // lambda_array = arma::exp(arma::linspace(std::log(1),
     //                                         std::log(lambda_factor), nlambda));
     // lambda_array.shed_row(0);
     // lambda_array *= lambda_max;


     if (!quiet)
     {
       Rcpp::Rcout << "lambda = " << lambda << std::endl;
     }
     // bic_array = arma::zeros(nlambda - 1, 1);
     // arma::vec shat_array = arma::zeros(nlambda - 1, 1);

     // double lambda;
     // lambda = lambda_array(ilambda);

     B_inner = B_out;
     P_beta = arma::zeros(p, m);
     for (int t = 0; t < T_inner; t++)
     {
       B_inner_old = B_inner;
       for (int j = 0; j < m; j++)
       {
         idx = calN_j_cpp(n, j);
         arma::uvec neighbors_j = arma::find(adjacency_matrix.row(j));
         P_beta.col(j) = P_beta.col(j) + tau_penalty_beta *
           arma::sum(arma::repmat(B_inner_old.col(j),
                                  1, neighbors_j.n_elem)
                       - B_inner_old.cols(neighbors_j), 1);

         tmp = omega(j) * (rho(j) * B_inner_old.col(j) - 1.0 / m / n * arma::trans(X.rows(idx)) * (X.rows(idx) * B_inner_old.col(j) - yt(idx)) -
                           P_beta.col(j) + tau_penalty_beta * arma::sum(arma::repmat(B_inner_old.col(j), 1, neighbors_j.n_elem) + B_inner_old.cols(neighbors_j), 1));
         // Update B
         // B_inner.col(j) = soft_thresholding_cpp(tmp, lambda * omega(j));
         B_inner.col(j) = SVT_vec_cpp(tmp, d1, d2, lambda * omega(j));
       }
       errors_inner(t, v) = std::sqrt(std::pow(arma::norm(B_inner - arma::repmat(betaT, 1, m), "fro"), 2) / m);
       if (!quiet)
       {
         Rcpp::Rcout << errors_inner(t, v) << "\t";
       }
     }
     // bic
     int shat = 0;

     for (int j = 0; j < m; j++)
     {
       idx = calN_j_cpp(n, j);
       bicj += quantile_lossCPP(y(idx) - X.rows(idx) * B_inner.col(j), tau);
       shat += int(arma::rank(arma::reshape(B_inner.col(j), d1, d2)));
       // shat = std::max(shat, int(arma::as_scalar(arma::accu(B_inner.col(j) != 0))));
       // shat += int(arma::as_scalar(arma::accu(B_inner.col(j) != 0)));
     }
     shat = 1.0 * shat / m;
     bicj = bicj / N + std::log(N) * std::log(std::log(N)) / N * std::max(shat, 1);

     B_out = B_inner;
   }
   // Rcpp::Named("history") = Rcpp::List::create(Rcpp::Named("errors_inner") = arma::vectorise(errors_inner),
   //             Rcpp::Named("errors_outer") = errors_inner.row(T_inner - 1),
   //             Rcpp::Named("shat_mat") = shat_mat,
   //             Rcpp::Named("lambda_max_array") = lambda_max_array,
   //             Rcpp::Named("bic_mat") = bic_mat));

   return Rcpp::List::create(Rcpp::Named("B") = B_out,
                             Rcpp::Named("BIC") = bicj,
                             Rcpp::Named("history") = Rcpp::List::create(Rcpp::Named("errors_inner") = arma::vectorise(errors_inner))
                               );

 }
