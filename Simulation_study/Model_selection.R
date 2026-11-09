#Packages
library(mvtnorm)
library(LaplacesDemon)

#---------------------------------------------

#100 edges

#---------------------------------------------
#Scenario 1
#---------------------------------------------

load("Exch_1_Gauss.RData")
y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y3_ = exp(y[3,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y4_ = 1/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y_ = rbind(y1_,y2_,y3_,y4_)

loglik_z = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_z[i + (t-1)*n,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}

WAIC_z_sim_1 = WAIC(loglik_z)

#-----------------------------------
load("Exch_1.RData")
loglik_y = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_y[i + (t-1)*n,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}

WAIC_y_sim_1 = WAIC(loglik_y)

#---------------------------------------------
#Scenario 2
#---------------------------------------------

load("Exchange_2_Gauss.RData")
y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,]))
y3_ = 1/(1+exp(y[1,])+exp(y[2,]))
y_ = rbind(y1_,y2_,y3_)

loglik_z = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_z[i + (t-1)*n,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}

WAIC_z_sim_2 = WAIC(loglik_z)

#-----------------------------------
load("Exchange_2_Dirichlet.RData")
loglik_y = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_y[i + (t-1)*n,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}

WAIC_y_sim_2 = WAIC(loglik_y)


U_est = matrix(NA, n, Ti)
for(i in 1:n){
  for(t in 1:Ti){
    U_est[i,t] = which.max(table(factor(U_ex[i,t,], levels = 1:2)))
  }
}

sum(U_est==u_real$U_est[,,10000])/(n*Ti)

#---------------------------------------------

#250 edges

#---------------------------------------------
#Scenario 1
#---------------------------------------------

load("Simulation_1_Gauss_250_edges.RData")
y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y3_ = exp(y[3,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y4_ = 1/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y_ = rbind(y1_,y2_,y3_,y4_)

loglik_z = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_z[i + (t-1)*n,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}

WAIC_z_sim_1 = WAIC(loglik_z)

#-----------------------------------
load("Simulation_1_Dirichlet_250_edges.RData")
loglik_y = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_y[i + (t-1)*n,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}

WAIC_y_sim_1 = WAIC(loglik_y)


U_est = matrix(NA, n, Ti)
for(i in 1:n){
  for(t in 1:Ti){
    U_est[i,t] = which.max(table(factor(U_ex[i,t,], levels = 1:2)))
  }
}

sum(U_est==u_real$U_est[,,10000])/(n*Ti)



#---------------------------------------------
#Scenario 2
#---------------------------------------------

load("Simulation_2_Gauss_250_edges.RData")
y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,]))
y3_ = 1/(1+exp(y[1,])+exp(y[2,]))
y_ = rbind(y1_,y2_,y3_)

loglik_z = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_z[i + (t-1)*n,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}

WAIC_z_sim_2 = WAIC(loglik_z)

#-----------------------------------
load("Simulation_2_Dirichlet_250_edges.RData")
loglik_y = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_y[i + (t-1)*n,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}

WAIC_y_sim_2 = WAIC(loglik_y)


U_est = matrix(NA, n, Ti)
for(i in 1:n){
  for(t in 1:Ti){
    U_est[i,t] = which.max(table(factor(U_ex[i,t,], levels = 1:2)))
  }
}

sum(U_est==u_real$U_est[,,10000])/(n*Ti)



#---------------------------------------------
#Scenario 3
#---------------------------------------------

load("Simulation_3_Gauss_250_edges.RData")
y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y3_ = exp(y[3,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y4_ = 1/(1+exp(y[1,])+exp(y[2,])+exp(y[3,]))
y_ = rbind(y1_,y2_,y3_,y4_)

loglik_z = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_z[i + (t-1)*n,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}

WAIC_z_sim_3 = WAIC(loglik_z)

U_est = matrix(NA, n, Ti)
for(i in 1:n){
  for(t in 1:Ti){
    U_est[i,t] = which.max(table(factor(U_ex[i,t,], levels = 1:3)))
  }
}

sum(U_est==u_real$U_est[,,10000])/(n*Ti)

#-----------------------------------
load("Simulation_3_Dirichlet_250_edges.RData")
loglik_y = matrix(NA,n*Ti,5000)

for(t in 1:Ti){
  for(i in 1:n){
    for(r in 1:5000){
      loglik_y[i + (t-1)*n,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}

WAIC_y_sim_3 = WAIC(loglik_y)


U_est = matrix(NA, n, Ti)
for(i in 1:n){
  for(t in 1:Ti){
    U_est[i,t] = which.max(table(factor(U_ex[i,t,], levels = 1:3)))
  }
}

sum(U_est==u_real$U_est[,,10000])/(n*Ti)

