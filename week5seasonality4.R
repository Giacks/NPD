#Adding detritus to the system and change units to nutrients

library(deSolve)
library(spam)
library(viridis)
library(viridisLite)
library(fields)

# Create an empty list to add the parameters
param <- list()

# Define parameters 
param$up <- 0.5 # sinking rate of phyto (m/d)
param$ud <- 15 # Settling velocity detritus (m/d)
param$DV <- 5*10^-5*3600*24 # Diffusivity (m^2/d)
param$dz <- 4 # grid spacing (m)
param$z <- seq(param$dz/2, 100, by=param$dz) # depth (m)
param$n <- length(param$z) # number of grid cells
param$kc <- 0.05 # light absorbed by phyto plankton and detritus (m2/mmol N)
param$kw <- 0.0375 # light absorbed by water (1/m)
param$I0 <- 350 #350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
param$alp <- 0.1/(param$I0/200) # alpha slope of PI curve (1/uEinstein*d)
param$gmax <- 1.5 # g max (1/d)
param$HN <- 0.3 # Half saturation constant for nutrients (mmol /m^3)
param$m <- 0.03 # specific loss rate or mortality (1/d)
param$gamma <- 0.1 # grazing mortality m3 (mmol N)-1 d-1
param$eps <- 0.1 # remineralization coefficient (d-1)
param$Nb <- 30 # nutrient content at the bottom (mmol/m^3) 

# Define initial conditions - concentration of phyto plankton
PND <- c(rep(10,length(param$z)), rep(30,length(param$z)), rep(0, length(param$z)))

# Create a function
CalLightPDseason <- function(t, P, D, param) {
  
  season <- 1+cos(pi+(2*pi/365)*t+22*pi/365)
  
  Q <- param$kc * param$dz * ((cumsum(P) - P/2)+(cumsum(D) - D/2))
  
  I <- param$I0 * exp(-param$kw * param$z - Q)*season
  
  return(I)
}


# Create function that creates the differential equation at each step
FuncPNDseason <- function(t, PND, param) {
  P <- PND[1:param$n]
  N <- PND[(1+param$n):(2*param$n)]
  D <- PND[(2*param$n+1):(3*param$n)]
  
  # Phytoplankton flux
  Jap <- rep(0,param$n+1)
  Jdp <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jap[i] <- param$up * P[i-1]
    Jdp[i] <- -param$DV * (P[i] - P[i-1]) / param$dz
  }
  
  # Advective flux boundary
  Jap[1] = 0
  Jap[param$n+1] = 0
  
  # Diffusive flux boundary
  Jdp[1] = 0
  Jdp[param$n+1] = 0
  
  Jp = Jap + Jdp
  
  # Nutrient flux
  Jdn <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jdn[i] <- -param$DV * (N[i] - N[i-1]) / param$dz
  }
  
  
  # Diffusive flux boundary
  Jdn[1] = 0
  Jdn[param$n+1] = -param$DV*(param$Nb-N[param$n])/param$dz
  
  Jn = Jdn
  
  # Detritus flux
  Jad <- rep(0,param$n+1)
  Jdd <- rep(0,param$n+1)
  
  for (i in 2:param$n){
    Jad[i] <- param$ud * D[i-1]
    Jdd[i] <- -param$DV * (D[i] - D[i-1]) / param$dz
  }
  
  # Advective flux boundary
  Jad[1] = 0
  Jad[param$n+1] = param$ud*D[param$n]
  
  # Diffusive flux boundary
  Jdd[1] = 0
  Jdd[param$n+1] = 0
  
  Jd = Jad + Jdd
  
  I <- CalLightPDseason(t, P, D, param)
  
  g <- param$gmax * pmin(param$alp*I / sqrt(param$gmax^2 + (param$alp*I)^2),N/(N+param$HN))
  
  dPdT <- - (Jp[2:(param$n+1)] - Jp[1:param$n]) / param$dz + (g - param$m) * P - param$gamma*P*P
  
  dNdT <- - g*P + param$eps*D - (Jn[2:(param$n+1)] - Jn[1:param$n]) / param$dz
  
  dDdT <- - (Jd[2:(param$n+1)] - Jd[1:param$n]) / param$dz + param$m*P + param$gamma*P*P - param$eps*D
  
  return(list(c(dPdT, dNdT, dDdT)))
}

# Define time step
time <- seq(0,365*5, by=1)
timeplot <- seq(365,1096, by=1)

# Solve our differential equation
start.time <- Sys.time()
res = ode(PND, time, FuncPNDseason, param)
end.time <- Sys.time()
time.taken <- end.time-start.time
print(paste("Run time=",as.character(round(time.taken,1)),"s"))
res
dim(res)

Pplot <- res[365:1096,2:(param$n+1)]
Nplot <- res[365:1096,(param$n+2):(2*param$n+1)]
Dplot <- res[365:1096,(2*param$n+2):(3*param$n+1)]

Pplot <- res[,2:(param$n+1)]
Nplot <- res[,(param$n+2):(2*param$n+1)]
Dplot <- res[,(2*param$n+2):(3*param$n+1)]

dim(Nplot)
nrow(Pplot)
nrow(Dplot)

#Plot

#With fields package
phytolab <- expression("Phytoplankton (mmol N m"^{-3}~")")
nlab <- expression("Nutrients (mmol N m"^{-3}~")")
dlab <- expression("Detritus (mmol N m"^{-3}~")")

par(mfrow=c(1,3))
par(mar = c(4, 3, 3, 1), mgp = c(2, 0.5, 0))
image.plot(timeplot/365, param$z, Pplot, ylim=rev(range(param$z)), ylab="Depth (m)", xlab="Years",
           main=phytolab, col=topo.colors(100)[20:100], cex.axis = 1.3, cex.lab = 1.3)
box()
par(mar = c(4, 3, 3, 1), mgp = c(2, 0.5, 0))
image.plot(timeplot/365, param$z, Nplot, ylim=rev(range(param$z)), ylab="", xlab="Years",
           main=nlab, col=rev(rainbow(300, start = 0.1, end = 0.6)), cex.axis = 1.3, cex.lab = 1.3)
box()
par(mar = c(4, 3, 3, 1), mgp = c(2, 0.5, 0))
image.plot(timeplot/365, param$z, Dplot, ylim=rev(range(param$z)), ylab="", xlab="Years",
           main=dlab, col=rev(rainbow(300, start = 0.1, end = 0.6)), cex.axis = 1.3, cex.lab = 1.3)
box()
#970 420

#vertical Profiles of Phytoplankton at the summer solstice
par(mfrow=c(1,4))
par(mar = c(5, 3, 4, 1), mgp = c(2, 0.5, 0))
plot(Pplot[365/2,], param$z, ylim = rev(range(param$z)), type="l", lwd=3,col = "darkolivegreen", ylab="Depth (m)", xlab=phytolab, main = "Year 1 (summer solstice)", cex.axis = 1.3, cex.lab = 1.3)
plot(Pplot[(365+365/2),], param$z, ylim = rev(range(param$z)), type="l", lwd=3,col = "darkolivegreen", ylab="", xlab=phytolab, main = "Year 2 (summer solstice)", yaxt = "n", cex.axis = 1.3, cex.lab = 1.3)
plot(Pplot[(2*365+365/2),], param$z, ylim = rev(range(param$z)), type="l", lwd=3,col = "darkolivegreen", ylab="", xlab=phytolab, main = "Year 3 (summer solstice)", yaxt = "n", cex.axis = 1.3, cex.lab = 1.3)
plot(Pplot[(3*365+365/2),], param$z, ylim = rev(range(param$z)), type="l", lwd=3,col = "darkolivegreen", ylab="", xlab=phytolab, main = "Year 4 (summer solstice)", yaxt = "n", cex.axis = 1.3, cex.lab = 1.3)

#970 420
","#998ec3","#f1a340","#e66101"

#plot P and N at equilibrium
#winter solstice is the 21st of dec and corresponds to the day 2*365-11
dphylab <- expression("Concentration (mmol N m"^{-3}~")")
nutleg <- expression("Nutrients"*10^-2)
nutleg1 <- expression("Nutrients"*10^-1)
q1 <- expression(paste(1^st," Quarter (year 3)"))
q3 <- expression(paste(3^rd," Quarter (year 3)"))
wsol <- expression("Winter solstice (year 2)")
ssol <- expression("Summer solstice (year 3)")
jan1 <- expression(paste(1^st," of January (year 3)"))

par(mfrow=c(1,2))
par(mar = c(5, 3, 4, 0.5)) 
plot(Pplot[2*365-11,], param$z, ylim = rev(range(param$z)), xlim=c(0,0.4), type='l', lwd=3, col = "#03AC13", xlab = dphylab, ylab = "Depth (m)", main=wsol)
points(Dplot[2*365-11,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "#e66101", xlab = "", ylab = "")
points(Nplot[2*365-11,]/100, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "#998ec3", xlab = "", ylab = "")
legend("topright", legend = c("Phytoplankton", "Detritus", nutleg), col = c("#03AC13", "#e66101","#998ec3"), lty = 1, lwd=3, cex = 0.8, bty = "n")

par(mar = c(5, 3, 4, 0.5)) 
plot(Pplot[2*365-11+365/4,], param$z, ylim = rev(range(param$z)), xlim=c(0,3), type='l', lwd=3, col = "#03AC13", xlab = dphylab, ylab = "", main=q1)
points(Dplot[2*365-11+365/4,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "#e66101", xlab = "", ylab = "")
points(Nplot[2*365-11+365/4,]/10, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "#998ec3", xlab = "", ylab = "")
legend("topright", legend = c("Phytoplankton", "Detritus", nutleg1), col = c("#03AC13", "#e66101","#998ec3"), lty = 1, lwd=3, cex = 0.8, bty = "n")
#800 420

par(mfrow=c(1,2))
par(mar = c(5, 3, 4, 0.5)) 
plot(Pplot[2*365-11+365*2/4,], param$z, ylim = rev(range(param$z)), xlim=c(0,3), type='l', lwd=3, col = "#03AC13", xlab = dphylab, ylab = "Depth (m)", main=ssol)
points(Dplot[2*365-11+365*2/4,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "#e66101", xlab = "", ylab = "")
points(Nplot[2*365-11+365*2/4,]/10, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "#998ec3", xlab = "", ylab = "")
legend("topright", legend = c("Phytoplankton", "Detritus", nutleg1), col = c("#03AC13", "#e66101","#998ec3"), lty = 1, lwd=3, cex = 0.8, bty = "n")

par(mar = c(5, 3, 4, 0.5)) 
plot(Pplot[2*365-11+365*3/4,], param$z, ylim = rev(range(param$z)), xlim=c(0,3), type='l', lwd=3, col = "#03AC13", xlab = dphylab, ylab = "", main=q3)
points(Dplot[2*365-11+365*3/4,], param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "#e66101", xlab = "", ylab = "")
points(Nplot[2*365-11+365*3/4,]/10, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "#998ec3", xlab = "", ylab = "")
legend("topright", legend = c("Phytoplankton", "Detritus", nutleg1), col = c("#03AC13", "#e66101","#998ec3"), lty = 1, lwd=3, cex = 0.8, bty = "n")


#Light and nutrient limitation
#1
par(mfrow=c(1,2))
fiqeresp <- expression(frac(alpha*I[eq], sqrt(g[max]^2+(alpha*I[eq])^2)))
fnqeresp <- expression(frac(N[eq], N[eq]+H[N]))
Ieqws <- CalLightPDseason(719, Pplot[719,], Dplot[719,], param)

fiws <- param$alp*Ieqws / sqrt(param$gmax^2 + (param$alp*Ieqws)^2)
fnws <- Nplot[719,] / (Nplot[719,] + param$HN)

par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fiws, param$z, ylim = rev(range(param$z)), xlim=c(0,1), type='l', lwd=3, col = "gold", xlab="Numerical value (unitless)", ylab="Depth (m)", main = wsol)
points(fnws, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred")


fiqeresp <- expression(frac(alpha*I[eq], sqrt(g[max]^2+(alpha*I[eq])^2)))
fnqeresp <- expression(frac(N[eq], N[eq]+H[N]))
Ieq <- CalLightPDseason(810, Pplot[810,], Dplot[810,], param)

fi <- param$alp*Ieq / sqrt(param$gmax^2 + (param$alp*Ieq)^2)
fn <- Nplot[810,] / (Nplot[810,] + param$HN)

par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fi, param$z, ylim = rev(range(param$z)), xlim=c(0,1), type='l', lwd=3, col = "gold", xlab="Numerical value (unitless)", ylab="", main = q1)
points(fn, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred")

legend("top", legend = c(fiqeresp, fnqeresp), col = c("gold", "darkred"), lwd=3, lty = 1, cex = 0.8, bty = "n")


#2
par(mfrow=c(1,2))
fiqeresp <- expression(frac(alpha*I[eq], sqrt(g[max]^2+(alpha*I[eq])^2)))
fnqeresp <- expression(frac(N[eq], N[eq]+H[N]))
Ieqws <- CalLightPDseason(902, Pplot[902,], Dplot[902,], param)

fiws <- param$alp*Ieqws / sqrt(param$gmax^2 + (param$alp*Ieqws)^2)
fnws <- Nplot[902,] / (Nplot[902,] + param$HN)

par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fiws, param$z, ylim = rev(range(param$z)), xlim=c(0,1), type='l', lwd=3, col = "gold", xlab="Numerical value (unitless)", ylab="Depth (m)", main = ssol)
points(fnws, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred")


fiqeresp <- expression(frac(alpha*I[eq], sqrt(g[max]^2+(alpha*I[eq])^2)))
fnqeresp <- expression(frac(N[eq], N[eq]+H[N]))
Ieq <- CalLightPDseason(993, Pplot[810,], Dplot[993,], param)

fi <- param$alp*Ieq / sqrt(param$gmax^2 + (param$alp*Ieq)^2)
fn <- Nplot[993,] / (Nplot[993,] + param$HN)

par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fi, param$z, ylim = rev(range(param$z)), xlim=c(0,1), type='l', lwd=3, col = "gold", xlab="Numerical value (unitless)", ylab="", main = q3)
points(fn, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred")

#1st of jan
par(mfrow=c(1,1))
fiqeresp <- expression(frac(alpha*I[eq], sqrt(g[max]^2+(alpha*I[eq])^2)))
fnqeresp <- expression(frac(N[eq], N[eq]+H[N]))
Ieqws <- CalLightPDseason(730, Pplot[730,], Dplot[730,], param)

fiws <- param$alp*Ieqws / sqrt(param$gmax^2 + (param$alp*Ieqws)^2)
fnws <- Nplot[730,] / (Nplot[730,] + param$HN)

par(mar = c(5, 4, 4, 4) + 0.3) 
plot(fiws, param$z, ylim = rev(range(param$z)), xlim=c(0,1), type='l', lwd=3, col = "gold", xlab="Numerical value (unitless)", ylab="Depth (m)", main = jan1)
points(fnws, param$z, ylim = rev(range(param$z)), type='l', lwd=3, col = "darkred")

