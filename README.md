
% -------------------------------------------------------------------------


#RobustLRL: Fast and Robust Low-Rank Learning over Networks: A Decentralized Matrix Quantile Regression Approach

*RobustLRL* is an R package for Fast and Robust Low-Rank Learning over Networks: A Decentralized Matrix Quantile Regression Approach. 

All files for reproducing the numerical results in the paper 
"Fast and Robust Low-Rank Learning over Networks: A Decentralized Matrix Quantile Regression Approach" are relegated to simulations subdirectory. We assume the root folder is the same folder the README.md resides in. 

All methods are implemented with R (Version 3.6.3) and conducted on a Linux server with 64-core Intel(R) Xeon(R) Gold 256 CPU (2.0GHz) and 62.6 GB RAM. We implement the distributed algorithms in a fully synchronized distributed setting.

% -------------------------------------------------------------------------




# Install the package

For installation, we recommend to unzip the `tar.gz` file first and then use `devtools::install()` to install the package, which can make sure to install all the depends. Make sure R package `igraph`, `pracma`, `Rcpp`, `MASS`, `ggplot2`, `ggpubr`, `doRNG`, `foreach`, `doFuture`, `parallel`, `xtable` have all been properly imported.


#  The following shows how to get the simulation results 

To facilitate your reproducible run, please run

```shell
bash main.sh
```
Note that please after installing the *RobustLRL* package, run the code again！

  Figure  1 by sim_iteration.R
  Table 1  by sim_heavy_tailed.R
  Table 2  by sim_heterogeneity.R & sim_heterogeneity_deSVQR.R
  Table 3  by sim_heterogeneity.R & sim_heterogeneity_deSVQR.R
  Table 4  by sim_number_of_nodes.R & sim_number_of_nodes_deSVQR.R
  Table 5  by sim_probability_of_connection.R

In Supplementary Material:

  Figure H.1 by sim_lambda.R
  Table H.1  by sim_local_samplesize.R
  Table H.2  by sim_localBIC.R


% -------------------------------------------------------------------------

This package is maintained by Canyi Chen. If you have any questions or find any buggs please contact [canyic@umich.edu]. 

Thank you.



