##3 years
dzzzsensitivity <- function(dzzz,day, colour, legends){
  for (i in 1:length(dzzz)){ 
    # Create an empty list to add the parameters
    param <- list()
    
    # Define parameters 
    param$up <- 0.5 # sinking rate of phyto (m/d)
    param$ud <- 15 # Settling velocity detritus (m/d)
    param$DV <- 4.32 # Diffusivity (m^2/d)
    param$dz <- dzzz[i] # grid spacing (m)
    param$z <- seq(param$dz/2, 100, by=param$dz) # depth (m)
    param$n <- length(param$z) # number of grid cells
    param$kc <- 0.05 # light absorbed by phyto plankton and detritus (m2/mmol N)
    param$kw <- 0.0375 # light absorbed by water (1/m)
    param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
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
    time <- seq(0,4*365, by=1)
    
    # Solve our differential equation
    res = ode(PND, time, FuncPNDseason, param)
    
    Pplot <- res[,2:(param$n+1)]
    Nplot <- res[,(param$n+2):(2*param$n+1)]
    Dplot <- res[,(2*param$n+2):(3*param$n+1)]
    
    phytolab <- expression("Phytoplankton (mmol N m"^{-3}~")")
    
    if (i==1){
      par(mfrow=c(1,1))
      par(mar = c(5, 4, 4, 2))
      plot(Pplot[day,], param$z, ylim = rev(c(40,100)), type="l", lwd=3, col = colfunc[i], xlim=c(0,3),
           main = "Year 3 (summer solstice)", ylab = "Depth (m)", xlab = phytolab ,cex = 1)
    }
    points(Pplot[day,], param$z, ylim = rev(c(40,100)), type="l", lwd=3, col = colfunc[i],cex = 1)
    legend("topright",60,legend=legends,col=colfunc,lty=1, lwd=3,cex = 1, bty = "n")
    
  }
}
dzzz <- c(1,2,2.5,4,5,6)
colfunc1 <- colorRampPalette(c("white", "darkgreen"))
colfunc <- colfunc1(10)[2:8]
legends <- c("dz =1 m", "dz = 2 m", "dz = 2.5 m", "dz = 4 m", "dz = 5 m", "dz = 6 m")
day <- 3*365+365/2+0.5 # Which day to look at (The converged solution start at around day 100)

dzzzsensitivity(dzzz,day,colfunc, legends)  

#600 420




##7 years
dzzzsensitivity <- function(dzzz,day, colour, legends){
  for (i in 1:length(dzzz)){ 
    # Create an empty list to add the parameters
    param <- list()
    
    # Define parameters 
    param$up <- 0.5 # sinking rate of phyto (m/d)
    param$ud <- 15 # Settling velocity detritus (m/d)
    param$DV <- 4.32 # Diffusivity (m^2/d)
    param$dz <- dzzz[i] # grid spacing (m)
    param$z <- seq(param$dz/2, 100, by=param$dz) # depth (m)
    param$n <- length(param$z) # number of grid cells
    param$kc <- 0.05 # light absorbed by phyto plankton and detritus (m2/mmol N)
    param$kw <- 0.0375 # light absorbed by water (1/m)
    param$I0 <- 350 # incoming light radiation at surface (micro mol photons / (m2*s); uEinstein)
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
    time <- seq(0,8*365, by=1)
    
    # Solve our differential equation
    res = ode(PND, time, FuncPNDseason, param)
    
    Pplot <- res[,2:(param$n+1)]
    Nplot <- res[,(param$n+2):(2*param$n+1)]
    Dplot <- res[,(2*param$n+2):(3*param$n+1)]
    
    phytolab <- expression("Phytoplankton (mmol N m"^{-3}~")")
    
    if (i==1){
      par(mfrow=c(1,1))
      par(mar = c(5, 4, 4, 2))
      plot(Pplot[day,], param$z, ylim = rev(c(40,100)), type="l", lwd=3, col = colfunc[i], xlim=c(0,3),
           main = "Year 7 (summer solstice)", ylab = "Depth (m)", xlab = phytolab ,cex = 1)
    }
    points(Pplot[day,], param$z, ylim = rev(c(40,100)), type="l", lwd=3, col = colfunc[i],cex = 1)
    legend("topright",60,legend=legends,col=colfunc,lty=1, lwd=3,cex = 1, bty = "n")
    
  }
}
dzzz <- c(1,2,2.5,4,5,6)
colfunc1 <- colorRampPalette(c("white", "darkgreen"))
colfunc <- colfunc1(10)[2:8]
legends <- c("dz =1 m", "dz = 2 m", "dz = 2.5 m", "dz = 4 m", "dz = 5 m", "dz = 6 m")
day <- 7*365+365/2+0.5 # Which day to look at (The converged solution start at around day 100)

dzzzsensitivity(dzzz,day,colfunc, legends)  

#600 420


