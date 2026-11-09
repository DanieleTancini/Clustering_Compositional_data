library("mvtnorm")
library("LaplacesDemon")
WAIC1 = WAIC2 = WAIC3 = WAIC4 = 0
load("K1_apl_sum_exch_gauss.RData")

y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y3_ = exp(y[3,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y4_ = exp(y[4,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y5_ = 1/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y_ = rbind(y1_,y2_,y3_,y4_,y5_)
#K1
bb = 5000
loglike1 = matrix(NA, n*Ti, 10000-bb)
for(r in 1:(10000-bb)){
  for(i in 1:n){
    for(t in 1:Ti){
      loglike1[t+(i-1)*Ti,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}
  
WAIC1 = WAIC(loglike1)$WAIC
  
#K2
load("K2_apl_sum_exch_gauss.RData")

y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y3_ = exp(y[3,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y4_ = exp(y[4,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y5_ = 1/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y_ = rbind(y1_,y2_,y3_,y4_,y5_)
bb = 5000
loglike2 = matrix(NA, n*Ti, 10000-bb)
for(r in 1:(10000-bb)){
  for(i in 1:n){
    for(t in 1:Ti){
      loglike2[t+(i-1)*Ti,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}
  
WAIC2 = WAIC(loglike2)$WAIC
  
load("K3_apl_sum_exch_gauss.RData")

y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y3_ = exp(y[3,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y4_ = exp(y[4,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y5_ = 1/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y_ = rbind(y1_,y2_,y3_,y4_,y5_)
bb = 5000
loglike3 = matrix(NA, n*Ti, 10000-bb)
for(r in 1:(10000-bb)){
  for(i in 1:n){
    for(t in 1:Ti){
      loglike3[t+(i-1)*Ti,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}
WAIC3 = WAIC(loglike3)$WAIC

load("K4_apl_sum_exch_gauss.RData")

y1_ = exp(y[1,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y2_ = exp(y[2,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y3_ = exp(y[3,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y4_ = exp(y[4,])/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y5_ = 1/(1+exp(y[1,])+exp(y[2,])+exp(y[3,])+exp(y[4,]))
y_ = rbind(y1_,y2_,y3_,y4_,y5_)
bb = 5000
loglike4 = matrix(NA, n*Ti, 10000-bb)
for(r in 1:(10000-bb)){
  for(i in 1:n){
    for(t in 1:Ti){
      loglike4[t+(i-1)*Ti,r] = dmvnorm(y[,i + (t-1)*n], mu[,U_ex[i,t,r+5000],r+5000], sigma[,,U_ex[i,t,r+5000],r+5000], log = TRUE) - sum(log(y_[,i + (t-1)*n]))
    }
  }
}
WAIC4 = WAIC(loglike4)$WAIC


save.image(paste0("Mod_sel_gaus_apl.Rdata"))



