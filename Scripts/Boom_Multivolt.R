######################
# Multivoltinism Boom Life History
######################

# Source Boom life history strategy model and bespoke functions
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

# Generate seasonal temperature time series
Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 13.956243

# Combine time, date, and temperature into a data frame
temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]

#no discharge, no hydropeaking
discharge <- rep(0.1, time = length(temp$dts))

# different temp regimes (-2 and +6)
runs <- c(-2, 6)

# for adults, loop through temperature regimes, extract adult abundance in each timestep
for (i in 1:length(runs)){
  temp$Temperature <- temp$Temperature + runs[i]
  out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "3")
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
B.oneyear <- mlist[which(mlist$Date >= "2035-01-01" & mlist$Date <= "2035-12-31"), ]
B.oneyear$Date <- as.Date(B.oneyear$Date)
B.oneyear$Strategy <- rep("Boom", times = length(B.oneyear$Date))