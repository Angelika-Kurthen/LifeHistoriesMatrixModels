#############################
# Boom Life History degree day toggle
############################
library(lubridate)

# source datasets
# Boom_1sp_Model.R: Boom Life History population model
# 1spFunctions.R: bespoke functions
# NegExpSurv.R: negative exponential survival function for disturbance
source("Scripts/Boom_1sp_Model.R")
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

# Degree-day requirement sensitivity range ±10% around baseline degree-day requirement (500)
dd_seq <- seq(500*0.9, 500*1.1, length.out = 21)


# Run Bsp (itermediate response to different scenarios) with no HFE, no hydropeaking (so no disturbance), temperature of North American temperate stream between 3 and 18 C, and fecundities from between 20 and 2000 eggs per brood. 
dd_means <- vector()
for (dd in 1:length(dd_seq)){
  print(dd)
  out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), dds = dd_seq[dd])
  means.list.B <- mean.data.frame(out, burnin = 250, iteration = 2) 
  dd_means[dd] <- mean(means.list.B$mean.abund)
}

bdd_df <- as.data.frame(cbind(dd_seq, dd_means, rep("B", length(dd_means))))
bdd_df$dd_seq <- as.numeric(bdd_df$dd_seq)
bdd_df$dd_means <- as.numeric(bdd_df$dd_means)
bdd_lm <- lm((dd_means/10000)~dd_seq, data = bdd_df)
# bdd <- ggplot(data = dd_df, mapping = aes(x = dd_seq, y = dd_means/10000))+
#   geom_point(size = 1, col = "#228833")+
#   xlab("Degree Day Requirement")+
#   ylab("B sp Abundance Relative to K")+
#   theme_bw()

