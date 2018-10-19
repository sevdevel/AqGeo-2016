###############################################################################
# Simulation file 
# CSFe model 
# Authors: Filip Meysman - Sebastiaan van de Velde
###############################################################################

# Source file containing model function
source("CH2O CO2 O2 SO4 H2S FeS FeOH3 SVDV.R")

# Source file containing plotting info function
source("plotting info FeOOH.r")

# Number of simulations

model <- CSFe.model

#=============================================================================
# Dynamic simulation 2 (Db=0)
#=============================================================================

load("04 CSFe Steady Db~xL=1 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "04 CSFe Transient Db~xL=0"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$Db.0 <- 0
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","200 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 



#=============================================================================
# Dynamic simulation 2 (Db=0)
#=============================================================================

load("04 CSFe Transient Db~xL=0 1.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Transient Db~xL=0"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$Db.0 <- 0
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","200 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Steady simulation 2 (Db=2)
#=============================================================================

load("04 CSFe Transient Db~xL=0 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "04 CSFe Steady Db~xL=0"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$Db.0 <- 0
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Steady simulation 2 (Db=2)
#=============================================================================

load("00 CSFe Steady k.CSFO=494 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Steady Db~xL=2"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$Db.0 <- 2
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Steady simulation 2 (Db=1)
#=============================================================================

load("04 CSFe Steady Db~xL=2 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Steady Db~xL=1"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$Db.0 <- 1
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Dynamic simulation 2 (Db=0.5)
#=============================================================================

load("04 CSFe Steady Db~xL=1 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "04 CSFe Transient Db~xL=0.5"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$Db.0 <- 1/2
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","200 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Steady simulation 2 (Db=0.5)
#=============================================================================

load("04 CSFe Transient Db~xL=0.5 1.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Steady Db~xL=0.5"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$Db.0 <- 1/2
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Steady simulation 2 (Db=0.25)
#=============================================================================

load("04 CSFe Steady Db~xL=0.25 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Steady Db~xL=0.25"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$k.CSFO <- 494
PL$Db.0 <- 1/4
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Steady simulation 2 (Db=1/10)
#=============================================================================

load("04 CSFe Steady Db~xL=0.1 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Steady Db~xL=0.1"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$k.CSFO <- 494
PL$Db.0 <- 1/10
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Dynamic simulation 2 (Db=6)
#=============================================================================

load("04 CSFe Steady Db~xL=0.1 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "04 CSFe Transient Db~xL=0"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$k.CSFO <- 494
PL$Db.0 <- 0
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","200 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Steady simulation 2 (Db=0)
#=============================================================================

load("02 CSFe Steady Db=5 x.L=0 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Steady Db~xL=0"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$k.CSFO <- 494
#PL$Db.0 <- 0
#PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 


#=============================================================================
# Dynamic simulation 2 (Db=6)
#=============================================================================

load("01 CSFe Steady Db=6 x.L=5 3.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "04 CSFe Transient Db~xL=6"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$k.CSFO <- 494
PL$Db.0 <- 6
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,10000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","200 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Steady simulation 2 (Db=6)
#=============================================================================

load("04 CSFe Transient Db~xL=6 1.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Steady Db~xL=6"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$k.CSFO <- 494
PL$Db.0 <- 6
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Dynamic simulation 2 (Db=8)
#=============================================================================

load("04 CSFe Steady Db~xL=6 2.Rdata")

# Initialisation simulation type 

sim.info$index <- 1
sim.info$code <- "04 CSFe Transient Db~xL=8"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

PL$simulation.type <- "time.dependent"

PL$k.CSFO <- 494
PL$Db.0 <- 8
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Sequence of time points where output is needed
#-------------------------------------------------------------------------------

sim.info$time.seq <- c(0,1,6,12,24,48,24*10,1000000*24*365.25)/(24*365.25) 
sim.info$time.schedule <- c("0 hr","1 hr","6 hr","12 hr","1 d","2 d","10 d","200 yr")
sim.info$N.out <- length(sim.info$time.seq)

#-------------------------------------------------------------------------------
# Dynamic simulation
#-------------------------------------------------------------------------------

# Initialise simulation progress variables (no need to change)
sim.info$tc <- 0  # simulation time counter
sim.info$fc <- 0  # function call counter
sim.info$progress <- initialise.simulation.progress(PL$var.names) # matrix

# Initialise simulation progress plot (no need to change)
dev.set(2)
par(mfrow=c(1,1))
RC <- model(t=0,state=sim.info$SV.init,parameters=PL,summary.call=FALSE)[[1]]
RCC.t <- RCC(RC, sim.info$SV.init, SV.ref = rep(1E-03,PL$N.var))
RCC.max <- round(log10(sqrt(sum(RCC.t^2)/length(RCC.t)))) + 1
plot(seq(-10,3,length.out=10),seq(-5,RCC.max,length.out=10),main="progress to steady state",xlab="simulation time (log10 yr)",ylab="steady state index",type="l",col="white")
abline(h=-3,col="blue",lwd=2)

# Actual dynamic simulation 
sim.info$run.time <- print(system.time({
  out <- ode.1D(y=sim.info$SV, times = sim.info$time.seq, func = model, parms = PL, nspec = PL$N.var, method = "vode", bandwidth = 2, maxsteps = 200000, atol=1E-11, rtol=1E-8, verbose=TRUE)
}))["elapsed"]
print(sim.info$run.time)

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 

#=============================================================================
# Steady simulation 2 (Db=8)
#=============================================================================

load("04 CSFe Transient Db~xL=8 1.Rdata")

# Initialisation simulation type 

sim.info$index <- 2
sim.info$code <- "04 CSFe Steady Db~xL=8"
sim.info$name <- paste(sim.info$code,sim.info$index)

# Initialisation parameter list
sim.start <- simulation.list[[sim.info$N.out]]
PL <- sim.start$PL

# Adapt parameter list

#PL$simulation.type <- "time.dependent"
PL$simulation.type <- "direct.steady.state"

PL$k.CSFO <- 494
PL$Db.0 <- 8
PL$x.L <- 1+9*(1-exp(-PL$Db.0/3))

PL <- initialise.parameters(PL)

# Initialise set of state variables
PL <- initialise.state.variables(PL)

# Create and initialize object structures
SL <- initialise.species.arrays(PL,SL) 
OL <- c(SL,RL) # object list (OL)

# Initialisation of matrix SV of state variables
SV <- crude.initial.state(PL) # Crude initialise of matrix SV
head(SV)

# Initialisation of matrix SV with values from start sim
for (spec.id in PL$var.names)  SV[,spec.id] <- sim.start$SL[[spec.id]]$C
head(SV)

# Store initial state in sim.info
sim.info$SV.init <- SV
remove(SV)

#-------------------------------------------------------------------------------
# Steady simulation
#-------------------------------------------------------------------------------

#Actual steady state simulation

sim.info$run.time <- system.time({
  out <- steady.1D(y=sim.info$SV, func=model, parms=PL, names=PL$var.names, nspec=PL$N.var, positive=TRUE) #
  steady.state.reached <- attributes(out)$steady
  if (steady.state.reached) {SV <- out$y} else stop ("steady state not reached")
})["elapsed"]
print(sim.info$run.time)

#prepare matrix out

out <- matrix(nrow = 2, ncol = (PL$N.var*PL$N+1), data = 0)
out[1,] <- c(0,as.vector(sim.info$SV))
out[2,] <- c(1,as.vector(SV))
simulation.list <- setup.simulation.list(out,PL,model)
sim.info$N.out <- 2

# Packaging the simulation output in one data structure
# function "setup.simulation.list" from package "diagenesis"
simulation.list <- setup.simulation.list(out,PL,model)

#-------------------------------------------------------------------------------
# Plotting preparation
#-------------------------------------------------------------------------------

# Screen summary of last time point 
screenprint.simulation.summary(simulation.list[[sim.info$N.out]],F.fac=(10/365.25))

# initialise.plotting.info 
plot.info <- setup.plot.info(SL,RL,depth.unit="[cm]",spec.unit="[mM]",reac.unit="[umol cm-3 yr-1]")
plot.info <- initialise.plotting.info(plot.info,PL)

# Setting the line color for each simulation 
palette(rainbow(sim.info$N.out))
line.color <- rainbow(sim.info$N.out)
line.color[sim.info$N.out] <-  "black"
for (i in 1:sim.info$N.out) simulation.list[[i]]$sim.color <- line.color[i]

x11()
par(mfrow=plot.info$spec.mfrow,
    mar=plot.info$spec.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="species",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

x11()
par(mfrow=plot.info$reac.mfrow,
    mar=plot.info$reac.mar,
    lwd=2,col="black",cex.main=par("cex.lab"))

plot.simulation.list(x=simulation.list[1:sim.info$N.out],what="reactions",plot.info)
legend("bottomright", legend = sim.info$time.schedule[1:sim.info$N.out], col = line.color[1:sim.info$N.out],lwd="2")

#-------------------------------------------------------------------------------
# Save the simulation output list
#-------------------------------------------------------------------------------

# sim.info: stores info how simulation is performed (initial conditions, etc)
# plot.info: stores info how simulation output is plotted 
# simulation.list: stores output of simulation

save(sim.info, plot.info, simulation.list, file = paste(sim.info$name,".Rdata",sep="")) 
