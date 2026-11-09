This repository contains the R code and the data required to reproduce the simulation studies and real-data application presented in the paper: “Clustering with Compositional Data in the Presence of Spatio-Temporal Structures”, by Daniele Tancini and Francesco Bartolucci.

In particular, the repository contains two folders:

**Simulation_study** contains the R code required to reproduce the estimation results for the three synthetic datasets, as well as the two additional analyses presented in Appendix B. In particular:
- Exchange_simulation_1_Dirichlet_250_edges.R reproduces the estimation results for the Dirichlet model using the first synthetic dataset;
- Exchange_simulation_2_Dirichlet_250_edges.R reproduces the estimation results for the Dirichlet model using the second synthetic dataset;
- Exchange_simulation_3_Dirichlet_250_edges.R reproduces the estimation results for the Dirichlet model using the third synthetic dataset;
- Exchange_simulation_1_Gauss_250_edges.R reproduces the estimation results for the Gaussian model using the first synthetic dataset;
- Exchange_simulation_2_Gauss_250_edges.R reproduces the estimation results for the Gaussian model using the second synthetic dataset;
- Exchange_simulation_3_Gauss_250_edges.R reproduces the estimation results for the Gaussian model using the third synthetic dataset;
- Omega_sampler_general_model_speed.R contains the function implementing the auxiliary process used in the approximate exchange algorithm;
- Sim_1_Dir_k1.R, Sim_1_Dir_k2.R, Sim_1_Dir_k3.R, Sim_1_Gau_k1.R, Sim_1_Gau_k2.R, Sim_1_Gau_k3.R, Mod_sel_Dir.R, and Mod_sel_Gau.R reproduce the results for model selection comparison;
- Traceplot.R, Diagnostic.R, Posterior_summary.R, and Model_selection.R contain the code for producing the additional results reported in the paper.

**Application** contains the R code required to reproduce the estimation results for the application, as well as the data. In particular:
- K1_apl_sum_exch.R reproduces the estimation results for the Dirichlet model when k=1;
- K2_apl_sum_exch.R reproduces the estimation results for the Dirichlet model when k=2;
- K3_apl_sum_exch.R reproduces the estimation results for the Dirichlet model when k=3;
- K4_apl_sum_exch.R reproduces the estimation results for the Dirichlet model when k=4;
- K1_apl_sum_exch_gauss.R reproduces the estimation results for the Gaussian model when k=1;
- K2_apl_sum_exch_gauss.R reproduces the estimation results for the Gaussian model when k=2;
- K3_apl_sum_exch_gauss.R reproduces the estimation results for the Gaussian model when k=3;
- K4_apl_sum_exch_gauss.R reproduces the estimation results for the Gaussian model when k=4;
- Omega_sampler_general_model_speed.R contains the function implementing the auxiliary process used in the approximate exchange algorithm;
- Mod_sel_Dir_apl.R, Mod_sel_Gaus_apl.R, Plot_applicazione.R, and Analysis_results_apl.R contain the code for producing the additional results reported in the paper;
- Dati_provinciali.xlsx is the dataset and ProvCM01012025_g_WGS84 are the files for graphical representations.
