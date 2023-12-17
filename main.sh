#!/usr/bin/env bash


R CMD INSTALL decentralizedCQR_1.0.tar.gz
R CMD INSTALL FHDCQR_0.2.0.tar.gz
# R CMD INSTALL FHDQR_0.2.1.tar.gz



Rscript ./Sim/sim_iteration.R
Rscript ./Sim/sim_heavy_tailed.R
Rscript ./Sim/sim_heterogeneity.R
Rscript ./Sim/sim_heterogeneity_deSVQR.R
Rscript ./Sim/sim_number_of_nodes.R
Rscript ./Sim/sim_number_of_nodes_deSVQR.R
Rscript ./Sim/sim_probability_of_connection.R
Rscript ./Sim/sim_lambda.R
Rscript ./Sim/sim_local_samplesize.R
Rscript ./Sim/sim_localBIC.R





