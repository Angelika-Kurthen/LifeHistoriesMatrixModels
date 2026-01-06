###########################
## Slow Life History Degree Day toggle
##########################
library(lubridate)

# Slow_1sp_Model.R: Slow Life History population model
# 1spFunctions.R: bespoke functions
# NegExpSurv.R: negative exponential survival function for disturbance
source("Scripts/D_1sp_Model.R")
source("Scripts/1spFunctions.R")
source("Scripts/NegExpSurv.R")

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

# Generate seasonal temperature time series
Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 13.956243

# Combine time, date, and temperature into a data frame
temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")'
'
# Convert to timestep-based temperature object
temp <- TimestepTemperature(temp)

# Keep only timestep and temperature columns
temp <- temp[c(1,3)]

Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)
# no disturbance, no hydropeaking
discharge <- rep(0.1, time = length(temp$dts))

# Degree-day requirement sensitivity range ±10% around baseline degree-day requirement
dd_seq <- seq(1500*0.9, 1500*1.1, length.out = 21)


# Run single-species model across DD values
dd_means <- vector()
for (dd in 1:length(dd_seq)){
  print(dd)
  out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dds = dd_seq[dd])
  # Calculate mean abundance after burn-in
  means.list.D <- mean.data.frame(out, burnin = 250, iteration = 2) 
  dd_means[dd] <- mean(means.list.D$mean.abund)
}

# Assemble results 
ddd_df <- as.data.frame(cbind(dd_seq, dd_means, rep("D", length(dd_means))))
ddd_df$dd_seq <- as.numeric(ddd_df$dd_seq)
ddd_df$dd_means <- as.numeric(ddd_df$dd_means)

# Fit linear model: relative abundance vs DD requirement
ddd_lm <- lm((dd_means/10000) ~ dd_seq, data = ddd_df)
