#######################
## Testing for Chaos
##########################


# load packages
require(Chaos01)
require(tseriesChaos)


# Boom_1sp_Model.R: Boom Life History population model
# 1spFunctions.R: bespoke functions
source("Scripts/Boom_1sp_Model.R")
source("Scripts/1spFunctions.R")

# Simulation time in days (~100 years)
Time <- c(1:36500)

# Day-of-year index repeated for 100 years
Date <- rep(1:365, times = 100)

# Actual calendar dates spanning 2022–2121
Day <- seq(as.Date("2022-01-01"), as.Date("2121-12-31"), by="days")
# find and remove leap days
find_leap = function(x){
  day(x) == 29 & month(x) == 2 
}
Day <- Day[which(find_leap(Day) == F)]

# set Temperature to mean temperature of simulated temperature curve (Appendix S1: Equations S4)
Temperature <-13.95

# Combine time, date, and temperature into a data frame
temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")

# Convert to timestep-based temperature object
temp <- TimestepTemperature(temp)

# Keep only timestep and temperature columns
temp <- temp[c(1,3)]

# no disturbance, no hydropeaking
discharge <- rep(0.1, times = length(temp$dts))

# Run single-species model across fecundity values and test for chaos
chaos1 <- vector()
fecs <- seq(1200, 100000, by = 2000)
for (fec in 1:length(fecs)){
  out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = fecs[fec])
  out <- mean.data.frame(out, 250, 2)
  chaos <- testChaos01(out$mean.abund)
  chaos1[fec] <- chaos
}

# compile data
chaostestdf <- as.data.frame(cbind(fecs, chaos1))

# create an example of what chaotic dynamics look like 
# we will use a F3 = 5200
outchaos <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = 5200 )
outchaos2 <- mean.data.frame(outchaos, 250, 2)
length(unique(outchaos2$mean.abund))
chaosdf <- as.data.frame(cbind(outchaos2$mean.abund[1:2359], outchaos2$mean.abund[2:2360]))

# create an example of what stable dynamics look like 
# we will use a F3 = 1200
outstable <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = 1200)
outstable2 <- mean.data.frame(outstable, 250, 2)
length(unique(outstable2$mean.abund))
stabledf <- as.data.frame(cbind(outstable2$mean.abund[1:2359], outstable2$mean.abund[2:2360]))

# try with Lyapunov Exponent

le <- vector()
fecs <- seq(1200, 100000, by = 2000)
for (fec in 1:length(fecs)){
  out1 <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = fecs[fec])
  out2 <- mean.data.frame(out1, 250, 2)
  out3 <- as.vector(out2$mean.abund)
  lyapun <- lyap_k(out3, m=1, d=2, s=200, t=40, ref=1700, k=2, eps=4)
  le[fec] <- max(lyapun)
}

# same results as chaos 0-1 test 

