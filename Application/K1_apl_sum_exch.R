#Packages
library(readxl)
library(sf)
library(spdep)
library(giscoR)
library(maps)
library(igraph)
library(LaplacesDemon)

Dataset_compositional = read_excel("Dati_provinciali.xlsx")

Data = Dataset_compositional[,10:15]
Territorio = Dataset_compositional$TERRITORIO[1:107]
Y = array(NA, dim = c(6,107,4))

for(t in 1:4){
  for(i in 1:107){
    Y[,i,t] = as.numeric(Data[107*(t-1)+i,])
  }
}

dim(Y)
Y[3, , ] = Y[3,,]+Y[6,,]
Y = Y[-6,,]

# Read shapefile
prov = st_read("ProvCM01012025_g_WGS84.shp")

# Check data
names(prov)
plot(st_geometry(prov))

nb = poly2nb(prov, queen = TRUE)
W = nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(W) = prov$DEN_UTS
colnames(W) = prov$DEN_UTS

#Uniforma i nomi di W a quelli di Territorio
nm = rownames(W)

nm[nm == "Aosta"] = "Valle d'Aosta"
nm[nm == "Forli'-Cesena"] = "Forlì"
nm[nm == "Massa Carrara"] = "Massa-Carrara"
nm[nm == "Pesaro e Urbino"] = "Pesaro-Urbino"
nm[nm == "Reggio di Calabria"] = "Reggio Calabria"
nm[nm == "Reggio nell'Emilia"] = "Reggio Emilia"
nm[nm == "Verbano-Cusio-Ossola"] = "Verbania"

rownames(W) = nm
colnames(W) = nm

#--------

order_ind = order(Territorio)
Territorio[order_ind]
Y_fin = Y[,order_ind,]
Y = Y[,order_ind,]
W = W[order(rownames(W)), order(rownames(W))]

#check
Territorio[order_ind] == rownames(W)

y = matrix(Y, nrow = dim(Y)[1], ncol = dim(Y)[2]*dim(Y)[3])
Com = W
C = Com
C[lower.tri(C)] = 0

#Initial settings
n = 107
Ti = 4
k = 1
U = 1:k

#functions

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

alpha_par_iniz = matrix(rgamma(k*dim(y)[1], 1, 1), dim(y)[1], k)

#Initial U
U_ex = array(1, dim=c(n,Ti,R))
u = U_ex[,,1]

#Cont
alpha_par = array(0, dim =c(dim(y)[1], k, R))

#Init
alpha_par[,,1] = alpha_par_iniz

accal = 0

#Sigma RW
alp_rw = c(0.05)

set.seed(1234)
for(r in 2:R){
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
  
}

save.image(file = "K1_apl_sum_exch.RData")

