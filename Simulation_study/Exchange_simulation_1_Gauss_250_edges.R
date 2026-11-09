#carica sim 1
library(LaplacesDemon)
library(igraph)
library(mvtnorm)
library(bvartools)
source("Omega_sampler_general_model_speed.R")

load("Simulation_1_Dirichlet_250_edges.RData")

#Initial settings
n = 100
Ti = 5
k = 2
U = 1:k

dim(Y)
Z = array(NA, dim = c(dim(Y)[1]-1,n,Ti))

for(i in 1:n){
  for(t in 1:Ti){
    Z[,i,t] = log(Y[-4,i,t]/Y[4,i,t])
  }
}     

z = matrix(NA,dim(Y)[1]-1,n*Ti)
for(i in 1:(n*Ti)){
  z[,i] = log(y[-4,i]/y[4,i])
}

y = z
d= dim(Z)[1]

yf = array(y, dim = c(d,n,Ti))

#Parameters
Gammas_ = matrix(c(0,1/2,1/2,0), 2,2)
Betas_ = c(1/2,0)
Delta_ = matrix(c(0,-1,-1,0), 2,2)

mm = rep(0,d)
Vm = diag(100, d)

half_alpha = ceiling((d+1)/2) + 1
alpha = half_alpha*2
S = matrix(0,d,d)
diag(S) = alpha

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

mu_iniz = matrix(rnorm(d*k, 0, 100), d, k)
sigma_iniz = array(rinvwishart(alpha,S), dim = c(d, d, k))

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
mu = array(0, dim = c(d, k, R))
sigma = array(0, dim = c(d, d, k, R))
dens = matrix(0, k, n)
Betas_es = matrix(0, k, R)
Delta_es = array(0, dim = c(k,k,R))
Gammas_es = array(0, dim = c(k,k,R))

#Init
mu[,,1] = mu_iniz
sigma[,,,1] = sigma_iniz
Betas_es[,1] = betas
Gammas_es[,,1] = gammas
Delta_es[,,1] = delta

#Full par
solvesigmap1 = solve(Vm)
solvesigmap1mup1 = solve(Vm)%*%mm

id_con = 1

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
  for(v in 1:k){
    ns = ifelse(as.vector(u)==v,1,0)
    yn = y*matrix(rep(ns,d), nrow = d, byrow = TRUE)
    mu[,v,r] = rmvnorm(1,solve(solvesigmap1+sum(ns)*solve(sigma[,,v,r-1]))%*%(solvesigmap1mup1 + solve(sigma[,,v,r-1])%*%rowSums(yn)),solve((solvesigmap1+sum(ns)*solve(sigma[,,v,r-1]))))
    sum = colSums(yn)
    pos = ifelse(sum==0,1,0)
    posiz = which(pos==1,arr.ind = TRUE)
    ynn = as.matrix(yn[,-c(posiz)])
    if(dim(ynn)[2]==0){
      Sn = S
    }else{
      Sa = array((ynn-mu[,v,r]), dim = c(d,1,ncol(ynn)))
      St = matrix(rowSums(apply(Sa, 3, function(x) x%*%t(x))), nrow = d)
      Sn = S + St
    }
    an = alpha + sum(ns)
    sigma[,,v,r] = rinvwishart(an,Sn)
  }
  
  #Identif constraint
  norm_e = mu[id_con,,r]
  ind = order(norm_e)
  if(any(ind!=(1:k))){
    mu[,,r] =  mu[,ind,r]
    sigma[,,,r] = sigma[,,ind,r]
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
      for(i in 1:k) dens[i,] = loglik_normal(yf[,,t]-mu[,i,r],sigma[,,i,r])
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
      for(i in 1:k) dens[i,] = loglik_normal(yf[,,t]-mu[,i,r],sigma[,,i,r])
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
      for(i in 1:k) dens[i,] = loglik_normal(yf[,,t]-mu[,i,r],sigma[,,i,r])
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

save.image(file = "Simulation_1_Gauss_250_edges.RData")
