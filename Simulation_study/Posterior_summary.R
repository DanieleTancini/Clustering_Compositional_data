library(mcmcse)
#simulation 1
load("Simulation_1_Dirichlet_250_edges.RData")

mean(Betas_es[1,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[2,1,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[2,1,5001:10000])

load("Simulation_1_Gauss_250_edges.RData")

mean(Betas_es[1,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[2,1,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[2,1,5001:10000])

#simulation 2
load("Simulation_2_Dirichlet_250_edges.RData")

mean(Betas_es[1,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[2,1,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[2,1,5001:10000])

load("Simulation_2_Gauss_250_edges.RData")

mean(Betas_es[1,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[2,1,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[2,1,5001:10000])

#simulation 3
load("Simulation_3_Dirichlet_250_edges.RData")

mean(Betas_es[1,5001:10000])
mean(Betas_es[2,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Gammas_es[1,3,5001:10000])
mean(Gammas_es[2,3,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[1,3,5001:10000])
mean(Delta_es[2,1,5001:10000])
mean(Delta_es[2,3,5001:10000])
mean(Delta_es[3,1,5001:10000])
mean(Delta_es[3,2,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Betas_es[2,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Gammas_es[1,3,5001:10000])
sd(Gammas_es[2,3,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[1,3,5001:10000])
sd(Delta_es[2,1,5001:10000])
sd(Delta_es[2,3,5001:10000])
sd(Delta_es[3,1,5001:10000])
sd(Delta_es[3,2,5001:10000])


load("Simulation_3_Gauss_250_edges.RData")

mean(Betas_es[1,5001:10000])
mean(Betas_es[2,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Gammas_es[1,3,5001:10000])
mean(Gammas_es[2,3,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[1,3,5001:10000])
mean(Delta_es[2,1,5001:10000])
mean(Delta_es[2,3,5001:10000])
mean(Delta_es[3,1,5001:10000])
mean(Delta_es[3,2,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Betas_es[2,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Gammas_es[1,3,5001:10000])
sd(Gammas_es[2,3,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[1,3,5001:10000])
sd(Delta_es[2,1,5001:10000])
sd(Delta_es[2,3,5001:10000])
sd(Delta_es[3,1,5001:10000])
sd(Delta_es[3,2,5001:10000])

#------------------------
#100 edges

library(mcmcse)
#simulation 1
load("Exch_1.RData")

mean(Betas_es[1,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[2,1,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[2,1,5001:10000])

load("Exch_1_Gauss.RData")

mean(Betas_es[1,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[2,1,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[2,1,5001:10000])

#simulation 2
load("Exchange_2_Dirichlet.RData")

mean(Betas_es[1,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[2,1,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[2,1,5001:10000])

load("Exchange_2_Gauss.RData")

mean(Betas_es[1,5001:10000])
mean(Gammas_es[1,2,5001:10000])
mean(Delta_es[1,2,5001:10000])
mean(Delta_es[2,1,5001:10000])
sd(Betas_es[1,5001:10000])
sd(Gammas_es[1,2,5001:10000])
sd(Delta_es[1,2,5001:10000])
sd(Delta_es[2,1,5001:10000])





