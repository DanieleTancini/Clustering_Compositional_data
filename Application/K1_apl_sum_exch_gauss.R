#Packages
library(readxl)
library(sf)
library(spdep)
library(giscoR)
library(maps)
library(igraph)
library(LaplacesDemon)
library(mvtnorm)
library(bvartools)

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

dim(Y)
Z = array(NA, dim = c(dim(Y)[1]-1,n,Ti))

for(i in 1:n){
  for(t in 1:Ti){
    Z[,i,t] = log(Y[-5,i,t]/Y[5,i,t])
  }
}     

z = matrix(NA,dim(Y)[1]-1,n*Ti)
for(i in 1:(n*Ti)){
  z[,i] = log(y[-5,i]/y[5,i])
}

y = z
d= dim(Z)[1]

yf = array(y, dim = c(d,n,Ti))

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

mu_iniz = matrix(rnorm(d*k, 0, 100), d, k)
sigma_iniz = array(rinvwishart(alpha,S), dim = c(d, d, k))

#Initial U
U_ex = array(1, dim=c(n,Ti,R))
u = U_ex[,,1]

#Cont
mu = array(0, dim = c(d, k, R))
sigma = array(0, dim = c(d, d, k, R))

#Init
mu[,,1] = mu_iniz
sigma[,,,1] = sigma_iniz

#Full par
solvesigmap1 = solve(Vm)
solvesigmap1mup1 = solve(Vm)%*%mm

id_con = 1

set.seed(1234)
for(r in 2:R){
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
  
  
}

save.image(file = "K1_apl_sum_exch_gauss.RData")



