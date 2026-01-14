######################
# Multivoltinism Moderate Life History
######################

# source Moderate LH strategy model + bespoke functions
source("Scripts/Moderate_1sp_Model.R")
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

# Generate seasonal temperature time series
Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 13.956243

# Combine time, date, and temperature into a data frame
temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]

# no disturbance, no hydropeaking
discharge <- rep(0.1, time = length(temp$dts))

# different temp regimes (-2 and +6)
runs <- c(-2, 6)

# for adults, loop through temperature regimes, extract adult abundance in each timestep
for (i in 1:length(runs)){
  temp$Temperature <- temp$Temperature + runs[i]
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "3")
  m <- rowMeans(out)
  m <- cbind.data.frame(temp$dts, m[-1], rep(i, length(temp$dts)))
  assign(paste0("m",i), m)
  temp$Temperature <- temp$Temperature - runs[i]
}

# combine the different temperature regimes into one dataframe
mlist <- rbind(m1, m2)
# rename columns
colnames(mlist) <- c("Date", "Abund", "MeanTemp")
# subset one year 
A.oneyear <- mlist[which(mlist$Date >= "2035-01-01" & mlist$Date <= "2035-12-31"), ]
# label peaks
# peaks <- oneyear[c(19, 43, 71, 87, 98), ]
# # create tibble for where we want arrows to point at peaks
# arrows <- tibble(
#   x1 = peaks$Date,
#   x2 = peaks$Date,
#   y1 = c(1650, 675, 1800, 475, 1350), 
#   y2 = c(1550, 575, 1700, 375, 1250)
# )
# 
# # assign colors for arrows
# arrowcols <- c("#4477AA", "#EE6677", "#228833","#CCBB44","#CCBB44")
# 
# # make sure all dates are in date format
# arrows$x1 <- as.POSIXct(arrows$x1, format = "%Y-%m-%d")
# arrows$x1 <- as.Date(arrows$x1)
# arrows$x2 <- as.POSIXct(arrows$x2, fomrat = "%Y-%m_%d")
# arrows$x2 <- as.Date(arrows$x2)
# A.oneyear$Date <- as.Date(A.oneyear$Date)
# A.oneyear$Strategy <- rep("Moderate", times = length(A.oneyear$Date))