#############################
# Fast Life History fecundity toggle
############################
library(lubridate)

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

# Combine time, date, and temperature into a data frame
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]

Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)

# no disturbance, no hydropeaking
discharge <- rep(0.1, time = length(temp$dts))

# Degree-day requirement sensitivity range ±10% around baseline fecundity 
fec_seq <- seq(0.9*500, 1.1*500, length.out = 21)

# Run single-species model across fecundity values
fec_means <- vector()
for (fec in 1:length(fec_seq)){
print(fec)
out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = fec_seq[fec])
means.list.C <- mean.data.frame(out, burnin = 250, iteration = 99) 
fec_means[fec] <- mean(means.list.C$mean.abund)
}

# Assemble results 
c_fec_df <- as.data.frame(cbind(fec_seq, fec_means, rep("C", times = length(fec_means))))
c_fec_df$fec_seq <- as.numeric(c_fec_df$fec_seq)
c_fec_df$fec_means <- as.numeric(c_fec_df$fec_means)

# Fit linear model: relative abundance vs fecundity
cfec_lm <- lm(fec_means ~ fec_seq, data = c_fec_df)

#clear means
rm(fec_means)