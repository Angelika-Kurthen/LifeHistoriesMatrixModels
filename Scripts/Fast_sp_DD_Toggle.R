#############################
# Fast Life History degree day toggle
############################
library(lubridate)

# source datasets
# Fast_1sp_Model.R: Fast Life History population model
# 1spFunctions.R: bespoke functions
# NegExpSurv.R: negative exponential survival function for disturbance
source("Scripts/Fast_1sp_Model.R")
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
colnames(temp) <- c("Time", "Date", "Temperature")

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
dd_seq <- seq(900*0.9, 900*1.1, length.out = 21)

# Run single-species model across DD values
dd_means <- vector()
for (dd in 1:length(dd_seq)){
  print(dd)
  out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dds = dd_seq[dd])
  # Calculate mean abundance after burn-in
  means.list.C <- mean.data.frame(out, burnin = 250, iteration = 2) 
  dd_means[dd] <- mean(means.list.C$mean.abund)
}

# Assemble results 
cdd_df <- as.data.frame(cbind(dd_seq, dd_means, rep("C", length(dd_means))))
cdd_df$dd_seq <- as.numeric(cdd_df$dd_seq)
cdd_df$dd_means <- as.numeric(cdd_df$dd_means)

# Fit linear model: relative abundance vs DD requirement
cdd_lm <- lm((dd_means/10000) ~ dd_seq, data = cdd_df)