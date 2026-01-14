####################
# Code for density independent runs
####################

library(lubridate)

# read in life history strategy models and bespoke functions
source("Scripts/Boom_1sp_Model.R")
source("Scripts/Fast_1sp_Model.R")
source("Scripts/Moderate_1sp_Model.R")
source("Scripts/Slow_1sp_Model.R")
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

# Generate temperature time series set at 14 C
Temperature <-  rep(14, times = length(Time))

# Combine time, date, and temperature into a data frame
temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]

# no disturbance, no hydropeaking
discharge <- rep(0.1, time = length(temp$dts))

# Moderate LH 
A_out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dens.dep = F)
A_out <- mean.data.frame(A_out, burnin = 250, iteration = 2) # clean burn in 
# combine into dataframe with label
A_dens.ind <- cbind.data.frame(temp$dts[250:last(A_out$timesteps)], A_out, rep("A", length(A_out$mean.abund)))

# Boom LH
B_out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dens.dep = F)
B_out <- mean.data.frame(B_out, burnin = 250, iteration = 2)# clean burn in 
# combine into dataframe with label
B_dens.ind <- cbind.data.frame(temp$dts[250:last(B_out$timesteps)], B_out, rep("B", length(B_out$mean.abund)))

# Fast LH
C_out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dens.dep = F)
C_out <- mean.data.frame(C_out, burnin = 250, iteration = 2)# clean burn in 
# combine into dataframe with label
C_dens.ind <- cbind.data.frame(temp$dts[250:last(C_out$timesteps)], C_out, rep("C", length(C_out$mean.abund)))

# Slow LH
D_out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dens.dep = F)
D_out <- mean.data.frame(D_out, burnin = 250, iteration = 2)# clean burn in 
# combine into dataframe with label
D_dens.ind <- cbind.data.frame(temp$dts[250:last(C_out$timesteps)], D_out, rep("D", length(B_out$mean.abund)))

# make sure column names are the same
colnames(A_dens.ind)<- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(B_dens.ind) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(C_dens.ind) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(D_dens.ind) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")

# combine into large dataframe, make sure dates are in date format
dens.ind <- rbind(A_dens.ind, B_dens.ind, C_dens.ind, D_dens.ind)
dens.ind$Date <- as.Date(dens.ind$Date)
#subset to a 3 year period for easier graphing
dens.ind <- subset(dens.ind, Date >= "2030-01-01" & Date <= "2033-12-31")
