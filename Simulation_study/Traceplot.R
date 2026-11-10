load("Simulation_1_Dirichlet_250_edges.RData")

plot(alpha_par[1,1,], type = "l")
plot(alpha_par[2,1,], type = "l")
plot(alpha_par[3,1,], type = "l")
plot(alpha_par[4,1,], type = "l")

plot(alpha_par[1,2,], type = "l")
plot(alpha_par[2,2,], type = "l")
plot(alpha_par[3,2,], type = "l")
plot(alpha_par[4,2,], type = "l")

quartz(width = 14, height = 8)
par(mfrow = c(2, 4), mar = c(6, 6, 5, 3),
    mgp = c(3.5, 1.5, 0), cex.main = 2.5,  cex.lab  = 2.5,  cex.axis = 2.5)

plot(alpha_par[1,1,], type = "l",
     main = expression(italic("\u03B1"[11]^"")), 
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(alpha_par[2,1,], type = "l",
     main = expression(italic("\u03B1"[12]^"")), 
     xlab = "Iteration", ylab = "", lwd = 1.5)

plot(alpha_par[3,1,], type = "l",
     main = expression(italic("\u03B1"[13]^"")), 
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(alpha_par[4,1,], type = "l",
     main = expression(italic("\u03B1"[14]^"")), 
     xlab = "Iteration", ylab = "", lwd = 1.5)

plot(alpha_par[1,2,], type = "l",
     main = expression(italic("\u03B1"[21]^"")), 
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(alpha_par[2,2,], type = "l",
     main = expression(italic("\u03B1"[22]^"")), 
     xlab = "Iteration", ylab = "", lwd = 1.5)

plot(alpha_par[3,2,], type = "l",
     main = expression(italic("\u03B1"[23]^"")), 
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(alpha_par[4,2,], type = "l",
     main = expression(italic("\u03B1"[24]^"")), 
     xlab = "Iteration", ylab = "", lwd = 1.5)

#---------------------------------------------------

load("Simulation_2_Gauss_250_edges.RData")


quartz(width = 14, height = 8)
par(mfrow = c(2, 4), mar = c(6, 6, 5, 3),
    mgp = c(3.5, 1.5, 0), cex.main = 2.5,  cex.lab  = 2.5,  cex.axis = 2.5)

plot(mu[1,1,], type = "l",
     main = expression(italic("\u03BC"[11]^"")),  ylim = c(-4,4),
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(mu[2,1,], type = "l",
     main = expression(italic("\u03BC"[12]^"")), ylim = c(-4,4),
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(sigma[1,1,1,], type = "l",
     main = expression(italic("\u03C3"[111]^"2")), 
     xlab = "Iteration", ylab = "", lwd = 1.5)

plot(sigma[2,2,1,], type = "l",
     main = expression(italic("\u03C3"[122]^"2")), 
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(mu[1,2,], type = "l",
     main = expression(italic("\u03BC"[21]^"")), ylim = c(-4,4),
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(mu[2,2,], type = "l",
     main = expression(italic("\u03BC"[22]^"")), 
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

plot(sigma[1,1,2,], type = "l",
     main = expression(italic("\u03C3"[211]^"")), 
     xlab = "Iteration", ylab = "", lwd = 1.5)

plot(sigma[2,2,2,], type = "l",
     main = expression(italic("\u03C3"[222]^"")), 
     xlab = "Iteration", ylab = "Value", lwd = 1.5)

