#Diagnostic
load("Simulation_1_Dirichlet_250_edges.RData")
library(coda)
Z_a = 2.576
Res_l = matrix(NA,dim(alpha_par)[1],dim(alpha_par)[2])
for(u in 1:(dim(alpha_par)[2])){
  for(d in 1:(dim(alpha_par)[1])){
    Res_l[d,u] = (abs(geweke.diag(alpha_par[d,u,5001:10000])$z) <= Z_a)*1
  }
}

Res_b = rep(NA,k-1)
for(u in 1:(k-1)){
    Res_b[u] = (abs(geweke.diag(Betas_es[u,5001:10000])$z) <= Z_a)*1
}

Res_g = matrix(NA,k-1,k-1)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_g[u,v-u] = (abs(geweke.diag(Gammas_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

Res_d = matrix(NA,k,k)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

for(v in 1:(k-1)){
  for(u in (v+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

sum(c(Res_d[!is.na(Res_d)],Res_g,Res_b,Res_l)==1)/length(c(Res_d[!is.na(Res_d)],Res_g,Res_b,Res_l))

mean(c(accal,accbs,accd[accd!=0],accgs[accgs!=0])/10000)

#--------------------

#Diagnostic
load("Simulation_1_Gauss_250_edges.RData")
library(coda)
Z_a = 2.576
Res_m = matrix(NA,dim(mu)[1],dim(mu)[2])
for(u in 1:(dim(mu)[2])){
  for(d in 1:(dim(mu)[1])){
    Res_m[d,u] = (abs(geweke.diag(mu[d,u,5001:10000])$z) <= Z_a)*1
  }
}

Res_s = array(NA,dim = c(dim(sigma)[1],dim(sigma)[1],dim(sigma)[3]))
for(u in 1:(dim(sigma)[3])){
  for(d in 1:(dim(sigma)[1]-1)){
    for(dd in (d+1):(dim(sigma)[1])){
      Res_s[d,dd,u] = (abs(geweke.diag(sigma[d,dd,u,5001:10000])$z) <= Z_a)*1
    }
  }
}

Res_b = rep(NA,k-1)
for(u in 1:(k-1)){
  Res_b[u] = (abs(geweke.diag(Betas_es[u,5001:10000])$z) <= Z_a)*1
}

Res_g = matrix(NA,k-1,k-1)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_g[u,v-u] = (abs(geweke.diag(Gammas_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

Res_d = matrix(NA,k,k)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

for(v in 1:(k-1)){
  for(u in (v+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}


sum(c(Res_d[!is.na(Res_d)],Res_g,Res_b,Res_m,Res_s[!is.na(Res_s)])==1)/length(c(Res_d[!is.na(Res_d)],Res_g,Res_b,Res_m,Res_s[!is.na(Res_s)]))

mean(c(accbs,accd[accd!=0],accgs[accgs!=0])/10000)


#--------------

#Diagnostic
load("Simulation_2_Dirichlet_250_edges.RData")
library(coda)
Z_a = 2.576
Res_l = matrix(NA,dim(alpha_par)[1],dim(alpha_par)[2])
for(u in 1:(dim(alpha_par)[2])){
  for(d in 1:(dim(alpha_par)[1])){
    Res_l[d,u] = (abs(geweke.diag(alpha_par[d,u,5001:10000])$z) <= Z_a)*1
  }
}

Res_b = rep(NA,k-1)
for(u in 1:(k-1)){
  Res_b[u] = (abs(geweke.diag(Betas_es[u,5001:10000])$z) <= Z_a)*1
}

Res_g = matrix(NA,k-1,k-1)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_g[u,v-u] = (abs(geweke.diag(Gammas_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

Res_d = matrix(NA,k,k)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

for(v in 1:(k-1)){
  for(u in (v+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

sum(c(Res_d[!is.na(Res_d)],Res_g,Res_b,Res_l)==1)/length(c(Res_d[!is.na(Res_d)],Res_g,Res_b,Res_l))

mean(c(accal,accbs,accd[accd!=0],accgs[accgs!=0])/10000)

#--------------------

#Diagnostic
load("Simulation_2_Gauss_250_edges.RData")
library(coda)
Z_a = 2.576
Res_m = matrix(NA,dim(mu)[1],dim(mu)[2])
for(u in 1:(dim(mu)[2])){
  for(d in 1:(dim(mu)[1])){
    Res_m[d,u] = (abs(geweke.diag(mu[d,u,5001:10000])$z) <= Z_a)*1
  }
}

Res_s = array(NA,dim = c(dim(sigma)[1],dim(sigma)[1],dim(sigma)[3]))
for(u in 1:(dim(sigma)[3])){
  for(d in 1:(dim(sigma)[1]-1)){
    for(dd in (d+1):(dim(sigma)[1])){
      Res_s[d,dd,u] = (abs(geweke.diag(sigma[d,dd,u,5001:10000])$z) <= Z_a)*1
    }
  }
}

Res_b = rep(NA,k-1)
for(u in 1:(k-1)){
  Res_b[u] = (abs(geweke.diag(Betas_es[u,5001:10000])$z) <= Z_a)*1
}

Res_g = matrix(NA,k-1,k-1)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_g[u,v-u] = (abs(geweke.diag(Gammas_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

Res_d = matrix(NA,k,k)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

for(v in 1:(k-1)){
  for(u in (v+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}


sum(c(Res_d[!is.na(Res_d)],Res_g,Res_b,Res_m,Res_s[!is.na(Res_s)])==1)/length(c(Res_d[!is.na(Res_d)],Res_g,Res_b,Res_m,Res_s[!is.na(Res_s)]))

mean(c(accbs,accd[accd!=0],accgs[accgs!=0])/10000)


#--------------

#Diagnostic
load("Simulation_3_Dirichlet_250_edges.RData")
library(coda)
Z_a = 2.576
Res_l = matrix(NA,dim(alpha_par)[1],dim(alpha_par)[2])
for(u in 1:(dim(alpha_par)[2])){
  for(d in 1:(dim(alpha_par)[1])){
    Res_l[d,u] = (abs(geweke.diag(alpha_par[d,u,5001:10000])$z) <= Z_a)*1
  }
}

Res_b = rep(NA,k-1)
for(u in 1:(k-1)){
  Res_b[u] = (abs(geweke.diag(Betas_es[u,5001:10000])$z) <= Z_a)*1
}

Res_g = matrix(NA,k,k)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_g[u,v] = (abs(geweke.diag(Gammas_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

Res_d = matrix(NA,k,k)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

for(v in 1:(k-1)){
  for(u in (v+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

sum(c(Res_d[!is.na(Res_d)],Res_g[!is.na(Res_g)],Res_b,Res_l)==1)/length(c(Res_d[!is.na(Res_d)],Res_g[!is.na(Res_g)],Res_b,Res_l))

mean(c(accal,accbs,accd[accd!=0],accgs[accgs!=0])/10000)

#--------------------

#Diagnostic
load("Simulation_3_Gauss_250_edges.RData")
library(coda)
Z_a = 2.576
Res_m = matrix(NA,dim(mu)[1],dim(mu)[2])
for(u in 1:(dim(mu)[2])){
  for(d in 1:(dim(mu)[1])){
    Res_m[d,u] = (abs(geweke.diag(mu[d,u,5001:10000])$z) <= Z_a)*1
  }
}

Res_s = array(NA,dim = c(dim(sigma)[1],dim(sigma)[1],dim(sigma)[3]))
for(u in 1:(dim(sigma)[3])){
  for(d in 1:(dim(sigma)[1]-1)){
    for(dd in (d+1):(dim(sigma)[1])){
      Res_s[d,dd,u] = (abs(geweke.diag(sigma[d,dd,u,5001:10000])$z) <= Z_a)*1
    }
  }
}

Res_b = rep(NA,k-1)
for(u in 1:(k-1)){
  Res_b[u] = (abs(geweke.diag(Betas_es[u,5001:10000])$z) <= Z_a)*1
}

Res_g = matrix(NA,k,k)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_g[u,v] = (abs(geweke.diag(Gammas_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

Res_d = matrix(NA,k,k)
for(u in 1:(k-1)){
  for(v in (u+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}

for(v in 1:(k-1)){
  for(u in (v+1):k){
    Res_d[u,v] = (abs(geweke.diag(Delta_es[u,v,5001:10000])$z) <= Z_a)*1
  }
}


sum(c(Res_d[!is.na(Res_d)],Res_g[!is.na(Res_g)],Res_b,Res_m,Res_s[!is.na(Res_s)])==1)/length(c(Res_d[!is.na(Res_d)],Res_g[!is.na(Res_g)],Res_b,Res_m,Res_s[!is.na(Res_s)]))

mean(c(accbs,accd[accd!=0],accgs[accgs!=0])/10000)




