#Packages
library(LaplacesDemon)
library(igraph)
load("Simulation_2_Gauss_250_edges.RData")

source("Omega_sampler_general_model_speed.R")

dim(y)


x1 = exp(y[1,])/(1+exp(y[1,])+exp(y[2,]))
x2 = exp(y[2,])/(1+exp(y[1,])+exp(y[2,]))
x3 = 1/(1+exp(y[1,])+exp(y[2,]))
y_ = rbind(x1,x2,x3)
y= y_
Y = array(y, dim = c(dim(y)[1],n,Ti))

#Initial settings
n = 100
Ti = 5
k = 2
U = 1:k

#Parameters
Gammas_ = matrix(c(0,1/2,1/2,0), 2,2)
Betas_ = c(1/2,0)
Delta_ = matrix(c(0,-1,-1,0), 2,2)

#Packages
#K=2 
library(LaplacesDemon)
library(igraph)

logbetas = function(Betan, Beta, k, j, u){
  I = as.numeric(table(factor(u, 1:k))[j])
  (Betan[j] - Beta[j])*I
}

loggammas = function(Gamman, Gamma, n, Ti, k, j, l, u, C){
  I = matrix(0,k,k)
  for(t in 1:Ti){
    I = I + table(factor(matrix(u[,t], n,n),1:k),factor(t(matrix(u[,t],n,n)*t(C)),1:k))
  }
  (Gamman[j,l] - Gamma[j,l])*(I[j,l]+I[l,j])
}

logdelta = function(Deltan, Delta, Ti, k, j, l, u){
  I = matrix(0,k,k)
  for(t in 2:Ti){
    I = I + table(factor(u[,t-1],1:k),factor(u[,t],1:k))
  }
  (Deltan[j,l] - Delta[j,l])*I[j,l]
}

#Prior
log_prior_theta = function(x, mb, sb){
  dnorm(x, mb, sb, log = TRUE)
}

log_prior_alpha = function(x, mb, sb){
  sum(dgamma(x, mb, sb, log = TRUE))
}

#Hyperparameter
mb = 0
sb = 1
sb_b = 1

gam_a = 1
gam_b = 0.1

#Iterations
R = 10000

R_set = rep(2, R)

#Initial settings
set.seed(1234)
betas = Betas = Betasn = c(rnorm(k-1, 0, 1),0)
delta = Delta = Deltan = matrix(rnorm(k*k, 0, 1), k, k)
gammas = Gammas = Gammasn = matrix(rnorm(k*k, 0, 1), k, k)

alpha_par_iniz = matrix(rgamma(k*dim(y)[1], 1, 1), dim(y)[1], k)

diag(delta) = 0
diag(Delta) = 0
diag(Deltan) = 0

diag(gammas) = 0
diag(Gammas) = 0
diag(Gammasn) = 0

#Initial U
U_iniz = sample(1:k,n*Ti,replace = TRUE)
u = matrix(U_iniz,n,Ti)
U_ex = array(NA, dim=c(n,Ti,R))
U_ex[,,1] = u

#Cont
alpha_par = array(0, dim =c(dim(y)[1], k, R))
dens = matrix(0, k, n)
Betas_es = matrix(0, k, R)
Delta_es = array(0, dim = c(k,k,R))
Gammas_es = array(0, dim = c(k,k,R))

#Init
alpha_par[,,1] = alpha_par_iniz
Betas_es[,1] = betas
Gammas_es[,,1] = gammas
Delta_es[,,1] = delta

#Acceptance
accbs = rep(0, k-1)
accd = matrix(0, k, k)
accgs = matrix(0, k, k)
accal = rep(0, k)

#Sigma RW
sbs_rw = c(0.8)
sgs_rw = matrix(0.8,2,2)
diag(sgs_rw) = 0
sd_rw = matrix(0.8,2,2)
diag(sd_rw) = 0
alp_rw = c(0.1,0.1)

#Adaptive
stop_ = 0.5*R
lambda_bs = matrix(1,k-1,stop_)
lambda_gs = lambda_d = array(1, dim = c(k,k,stop_))
mu_bs = matrix(0,k-1,stop_)
sigma_bs = matrix(0,k-1,stop_)
sigma_bs[,1] = sbs_rw
mu_gs = mu_d = array(0,dim = c(k,k,stop_))
sigma_gs = sigma_d = array(0,dim = c(k,k,stop_))
sigma_gs[,,1] = sgs_rw
sigma_d[,,1] = sd_rw
alph_opt = 0.4
Cost = 0.001

set.seed(1234)
for(r in 2:R){
  ga = Cost/r
  RR = R_set[r]
  #---------------------
  for(j in 1:k){
    yy = t(y)[c(u)==j,]
    if(sum(c(u)==j)==0){
      log_alpha_par_new = log(alpha_par[,j,r-1]) + rnorm(dim(alpha_par)[1],0,alp_rw[j]) 
      alpha_par_new = exp(log_alpha_par_new)
      ratio_al = log_prior_alpha(alpha_par_new,gam_a,gam_b) - log_prior_alpha(alpha_par[,j,r-1],gam_a,gam_b) + sum(log(alpha_par_new) - log(alpha_par[,j,r-1]))
      alpha_al = min(1, exp(ratio_al))
      if(runif(1)<= alpha_al){
        alpha_par[,j,r] = alpha_par_new
        accal[j] = accal[j] +1
      }else{
        alpha_par[,j,r] = alpha_par[,j,r-1]
      }
    }else{
      log_alpha_par_new = log(alpha_par[,j,r-1]) + rnorm(dim(alpha_par)[1],0,alp_rw[j]) 
      alpha_par_new = exp(log_alpha_par_new)
      ratio_al = sum(ddirichlet(yy,alpha_par_new,log = TRUE)) - sum(ddirichlet(yy,alpha_par[,j,r-1],log = TRUE)) + log_prior_alpha(alpha_par_new,gam_a,gam_b) - log_prior_alpha(alpha_par[,j,r-1],gam_a,gam_b) + sum(log(alpha_par_new) - log(alpha_par[,j,r-1]))
      alpha_al = min(1, exp(ratio_al))
      if(runif(1)<= alpha_al){
        alpha_par[,j,r] = alpha_par_new
        accal[j] = accal[j] +1
      }else{
        alpha_par[,j,r] = alpha_par[,j,r-1]
      }
    }
  }
  
  #Identif constraint
  norm_e = alpha_par[1,,r]
  ind = order(norm_e)
  if(any(ind!=(1:k))){
    alpha_par[,,r] = alpha_par[,ind,r]
  }
  
  #Betas
  for(j in 1:(k-1)){
    Betasn[j] = Betas[j] + rnorm(1, 0, sbs_rw[j]) 
    
    omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Betas = Betasn, Delta = Delta, Gammas = Gammas, iter = RR, Com = Com)
    omega = omega$u
    
    #Pseudo betas
    ratio_b = log_prior_theta(Betasn[j], mb, sb_b) - log_prior_theta(Betas[j], mb, sb_b) + logbetas(Betan =  Betasn, Beta = Betas, k = k, j = j, u = u) + logbetas(Betan =  Betas, Beta = Betasn, k = k, j = j, u = omega)
    
    alpha_b = min(1, exp(ratio_b))
    if(runif(1)<= alpha_b){
      Betas[j] = Betasn[j]
      accbs[j] = accbs[j] +1
    }
    Betasn = Betas
    Betas_es[j,r] = Betas[j]
    if(r <= stop_){
      lambda_bs[j,r] = exp(log(lambda_bs[j,r-1]) + ga*(alpha_b - alph_opt))
      mu_bs[j,r] = mu_bs[j,r-1] + ga*(Betas[j] - mu_bs[j,r-1])
      sigma_bs[j,r] = sigma_bs[j,r-1] + ga*((Betas[j] - mu_bs[j,r-1])^2 - sigma_bs[j,r-1])
      sbs_rw[j] = lambda_bs[j,r]*sigma_bs[j,r]
    }
  }
  
  #Gammas
  for(j in 1:(k-1)){
    for(l in (j+1):k){
      Gammasn[j,l] = Gammas[j,l] + rnorm(1, 0, sgs_rw[j,l])
      Gammasn[l,j] = Gammasn[j,l]
      
      omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Betas = Betas, Delta = Delta, Gammas = Gammasn, iter = RR, Com = Com)
      omega = omega$u
      
      ratio_g = log_prior_theta(Gammasn[j,l], mb, sb) - log_prior_theta(Gammas[j,l], mb, sb) + loggammas(Gammasn, Gammas, n, Ti, k, j, l, u = u, C = C) + loggammas(Gammas, Gammasn, n, Ti, k, j, l, u = omega, C = C)
      
      alpha_g = min(1, exp(ratio_g))
      if(runif(1)<= alpha_g){
        Gammas[j,l] = Gammasn[j,l]
        Gammas[l,j] = Gammasn[l,j]
        accgs[j,l] = accgs[j,l] +1
      }
      Gammasn = Gammas
      Gammas_es[j,l,r] = Gammas[j,l]
      Gammas_es[l,j,r] = Gammas[l,j]
      if(r <= stop_){
        lambda_gs[j,l,r] = exp(log(lambda_gs[j,l,r-1]) + ga*(alpha_g - alph_opt))
        mu_gs[j,l,r] = mu_gs[j,l,r-1] + ga*(Gammas[j,l] - mu_gs[j,l,r-1])
        sigma_gs[j,l,r] = sigma_gs[j,l,r-1] + ga*((Gammas[j,l] - mu_gs[j,l,r-1])^2 - sigma_gs[j,l,r-1])
        sgs_rw[j,l] = lambda_gs[j,l,r]*sigma_gs[j,l,r]
      }
    }
  }
  
  #Delta
  for(j in 1:(k-1)){
    for(l in (j+1):k){
      Deltan[j,l] = Delta[j,l] + rnorm(1, 0, sd_rw[j,l])
      omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Betas = Betas, Delta = Deltan, Gammas = Gammas, iter = RR, Com = Com)
      omega = omega$u
      
      ratio_d = log_prior_theta(Deltan[j,l], mb, sb) - log_prior_theta(Delta[j,l], mb, sb) + logdelta(Deltan, Delta, Ti, k, j, l, u = u) + logdelta(Delta, Deltan, Ti, k, j, l, u = omega)
      
      alpha_d = min(1, exp(ratio_d))
      if(runif(1)<= alpha_d){
        Delta[j,l] = Deltan[j,l]
        accd[j,l] = accd[j,l] +1
      }
      Deltan = Delta
      Delta_es[j,l,r] = Delta[j,l]
      if(r <= stop_){
        lambda_d[j,l,r] = exp(log(lambda_d[j,l,r-1]) + ga*(alpha_d - alph_opt))
        mu_d[j,l,r] = mu_d[j,l,r-1] + ga*(Delta[j,l] - mu_d[j,l,r-1])
        sigma_d[j,l,r] = sigma_d[j,l,r-1] + ga*((Delta[j,l] - mu_d[j,l,r-1])^2 - sigma_d[j,l,r-1])
        sd_rw[j,l] = lambda_d[j,l,r]*sigma_d[j,l,r]
      }
    }
  }
  
  for(l in 1:(k-1)){
    for(j in (l+1):k){
      Deltan[j,l] = Delta[j,l] + rnorm(1, 0, sd_rw[j,l])
      
      omega = Omega_sim(n = n, Ti = Ti, k = k, u = u, Betas = Betas, Delta = Deltan, Gammas = Gammas, iter = RR, Com = Com)
      omega = omega$u
      
      ratio_d = log_prior_theta(Deltan[j,l], mb, sb) - log_prior_theta(Delta[j,l], mb, sb) + logdelta(Deltan, Delta, Ti, k, j, l, u = u) + logdelta(Delta, Deltan, Ti, k, j, l, u = omega)
      
      alpha_d = min(1, exp(ratio_d))
      if(runif(1)<= alpha_d){
        Delta[j,l] = Deltan[j,l]
        accd[j,l] = accd[j,l] +1
      }
      Deltan = Delta
      Delta_es[j,l,r] = Delta[j,l]
      if(r <= stop_){
        lambda_d[j,l,r] = exp(log(lambda_d[j,l,r-1]) + ga*(alpha_d - alph_opt))
        mu_d[j,l,r] = mu_d[j,l,r-1] + ga*(Delta[j,l] - mu_d[j,l,r-1])
        sigma_d[j,l,r] = sigma_d[j,l,r-1] + ga*((Delta[j,l] - mu_d[j,l,r-1])^2 - sigma_d[j,l,r-1])
        sd_rw[j,l] = lambda_d[j,l,r]*sigma_d[j,l,r]
      }
    }
  }
  
  #Full U
  res = matrix(NA, n, k)
  res2 = matrix(NA, n, k)
  res3 = matrix(NA, n, k)
  res4 = matrix(NA, n, k)
  fin = efin = efinb = efinc = eq_a = eq_b = eq_c = matrix(NA, n, k)
  
  for(t in 1:Ti){
    if(t==1){
      for(i in 1:k) dens[i,] = ddirichlet(t(Y[,,t]), alpha_par[,i,r], log = TRUE)
      for(l in 1:n){
        neigh = which(Com[l,]==1)
        counts = table(factor(u[neigh,t], levels=1:k))
        res[l,] = Gammas %*% as.numeric(counts)
        res2[l,] = Delta[1:k,u[l,2]]
        res3[l,] = Betas[1:k]
        fin[l,] = res[l,]+res2[l,]+res3[l,]
        efin[l,] = exp(fin[l,]+dens[,l])
        eq_a[l,] = efin[l,]/sum(efin[l,])
        u[l,t] = sample(U, 1, prob = eq_a[l,], replace = TRUE)
      }
    }else if(t == Ti){
      for(i in 1:k) dens[i,] = ddirichlet(t(Y[,,t]), alpha_par[,i,r], log = TRUE)
      for(l in 1:n){
        neigh = which(Com[l,]==1)
        counts = table(factor(u[neigh,t], levels=1:k))
        res[l,] = Gammas %*% as.numeric(counts)
        res4[l,] = Delta[u[l,Ti-1],1:k]
        res3[l,] = Betas[1:k]
        fin[l,] = res[l,]+res3[l,]+res4[l,]
        efinb[l,] = exp(fin[l,]+dens[,l])
        eq_b[l,] = efinb[l,]/sum(efinb[l,])
        u[l,t] = sample(U, 1, prob = eq_b[l,], replace = TRUE)
      }
    }else{
      for(i in 1:k) dens[i,] = ddirichlet(t(Y[,,t]), alpha_par[,i,r], log = TRUE)
      for(l in 1:n){
        neigh = which(Com[l,]==1)
        counts = table(factor(u[neigh,t], levels=1:k))
        res[l,] = Gammas %*% as.numeric(counts)
        res2[l,] = Delta[1:k,u[l,t+1]]
        res4[l,] = Delta[u[l,t-1],1:k]
        res3[l,] = Betas[1:k]
        fin[l,] = res[l,]+res2[l,]+res3[l,]+res4[l,]
        efinc[l,] = exp(fin[l,]+dens[,l])
        eq_c[l,] = efinc[l,]/sum(efinc[l,])
        u[l,t] = sample(U, 1, prob = eq_c[l,], replace = TRUE)
      }
    }
  }
  U_ex[,,r] = u 
}

U_est = matrix(NA,n,Ti)
for(t in 1:Ti){
  for(i in 1:n){
    U_est[i,t] = which.max(table(factor(U_ex[i,t,],levels = 1:2)))
  }
}

save.image(file = "Simulation_2_Dirichlet_250_edges.RData")

