library("mvtnorm")
library("LaplacesDemon")
WAIC1 = WAIC2 = WAIC3 = WAIC4 = 0

load("K1_apl_sum_exch.RData")

#K1
bb = 5000
loglike1 = matrix(NA, n*Ti, 10000-bb)
for(r in 1:(10000-bb)){
  for(i in 1:n){
    for(t in 1:Ti){
      loglike1[t+(i-1)*Ti,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}
  
WAIC1 = WAIC(loglike1)$WAIC
  
#K2
load("K2_apl_sum_exch.RData")
bb = 5000
loglike2 = matrix(NA, n*Ti, 10000-bb)
for(r in 1:(10000-bb)){
  for(i in 1:n){
    for(t in 1:Ti){
      loglike2[t+(i-1)*Ti,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}
  
WAIC2 = WAIC(loglike2)$WAIC
  
#K3
load("K3_apl_sum_exch.RData")
bb = 5000
loglike3 = matrix(NA, n*Ti, 10000-bb)
for(r in 1:(10000-bb)){
  for(i in 1:n){
    for(t in 1:Ti){
      loglike3[t+(i-1)*Ti,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}
  
WAIC3 = WAIC(loglike3)$WAIC

#K3
load("K4_apl_sum_exch.RData")
bb = 5000
loglike4 = matrix(NA, n*Ti, 10000-bb)
for(r in 1:(10000-bb)){
  for(i in 1:n){
    for(t in 1:Ti){
      loglike4[t+(i-1)*Ti,r] = ddirichlet(y[,i + (t-1)*n], alpha_par[,U_ex[i,t,r+5000],r+5000], log = TRUE) 
    }
  }
}

WAIC4 = WAIC(loglike4)$WAIC

save.image(paste0("Mod_sel_dir_apl.Rdata"))



