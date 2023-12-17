#ifndef CQR_H
#define CQR_H
#include <iostream>
#include <cmath>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat kernel(arma::mat E, arma::vec h);

double mod(arma::vec v);

double quantile_lossCPP(arma::vec x, double tau);

double cqr_loss_cpp(arma::mat X, arma::vec y, arma::vec beta, arma::vec alpha, arma::vec& tau);

arma::uvec calN_j_cpp(int n, int j);

arma::vec pmax_arma(arma::vec x, double bound);

arma::vec soft_thresholding_cpp(arma::vec x, double t);

arma::mat SVT_cpp(arma::mat &B, double lambda);

arma::vec SVT_vec_cpp(arma::vec &beta, int d1, int d2, double lambda);

#endif
