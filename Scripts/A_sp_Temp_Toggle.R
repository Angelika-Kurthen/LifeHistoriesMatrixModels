#############################
# Moderate Life History temperature toggle
############################

library(lubridate)

# Moderate_1sp_Model.R: Moderate Life History population model
# 1spFunctions.R: bespoke functions
# NegExpSurv.R: negative exponential survival function for disturbance
source("Scripts/Moderate_1sp_Model.R")
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


###################################
# How does abundance, size, and 
# biomass respond to temperature regime
###################################

# Combine time, date, and temperature into a data frame
temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]

Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)

# no disturbance, no hydropeaking
discharge <- rep(0.1, time = length(temp$dts))

# vectors for data to go into 
temp_regime <- vector() # mean temperature of temperature regime
temp_means <- vector() # mean abundance 
temp_seq <- seq(-10, 10, by = 1) # sequence to modify temperature regime
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te] # modify temperature +-10 C 
  temp_regime[te] <- mean(temp$Temperature) # get mean temperature
  # run model
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  # reset temperature regime
  temp$Temperature <- temp$Temperature - temp_seq[te] 
  # clean up burn in 
  means.list.A <- mean.data.frame(out, burnin = 250, iteration = 2) 
  # get mean 
  temp_means[te] <- mean(means.list.A$mean.abund)
}

# create dataframe with output
a_temp_adjust_df <- as.data.frame(cbind(temp_regime, temp_means, rep("A", times = length(temp_means))))
a_temp_adjust_df$temp_regime <- as.numeric(a_temp_adjust_df$temp_regime)
a_temp_adjust_df$temp_means <- as.numeric(a_temp_adjust_df$temp_means)

# make a vector for what means are 0 or not
binary <- as.integer(a_temp_adjust_df$temp_means != 0)

# vectors for data to go into 
size_means <- vector() 
stage3s_means <- vector()
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te] # modify temperature +-10 C 
  temp_regime[te] <- mean(temp$Temperature) # mean temperature of temperature regime
  # run model
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "size")
  temp$Temperature <- temp$Temperature - temp_seq[te]   # reset temperature regime
  size_means[te] <- colMeans(out) # get mean size of all stages
  size_means[te] <- sum(0.0094*(size_means[te])^2.754)  # multiply relative size (which is also biologically plausible) by Benke et al 1999 Table 2 a and b params (M(mg) = aL^b) 
  stage3s_means[te] <- mean(out[,3,1]) # extract only stage 3 mean size
}

# calculate stage 3 biomass
stage3s_means <-  (0.0094*(stage3s_means)^2.754)
stage3s_means <- stage3s_means * binary  # figure out which are associated with non-zero abundances
stage3s_means[which(stage3s_means == 0)] <- NA # set those to NA (no abundance, no biomass)

# combine these into a dataframe and make sure numbers are numeric
a_size_df <- as.data.frame(cbind(temp_regime, size_means, stage3s_means, rep("A", times = length(temp_means))))
a_size_df$stage3s_means <- as.numeric(a_size_df$stage3s_means)
a_size_df$temp_regime <- as.numeric(a_size_df$temp_regime)
a_size_df$size_means <- as.numeric(a_size_df$size_means)

# we can plot
atemp <- ggplot(data = a_temp_adjust_df, mapping = aes(x = temp_seq, y = temp_means))+
  geom_line(size= 1, col = "#66CCEE")+
  xlab("Degree C Change")+
  ylab("Sp A Abundance ")+
  theme_bw()


################################
# How does abundance, size, and biomass 
# repsond to seasonal timing of disturbance
###############################

# -----------------------------
# For a winter time disturbance 
# -----------------------------
# vectors for data to go into
temp_regime <- vector() # mean temperature of temperature regime
temp_means <- vector() # mean abundance 
temp_seq <- seq(-10, 10, by = 1) # sequence to modify temperature regime
short <- vector() # mean abundance shortly post disturbance (6 time steps or 3 months)
discharge[259] <- 1 # set a bankfull discharge disturbance in wintertime (Nov 24) 
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te]  # modify temperature +-10 C 
  temp_regime[te] <- mean(temp$Temperature) # mean temperature
  # run model
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  temp$Temperature <- temp$Temperature - temp_seq[te]   # reset temperature regime
  means.list.A<- mean.data.frame(out, burnin = 250, iteration = 2)   # get means without burn in 
  short[te] <- mean(means.list.A$mean.abund[10:16]) # get means right after disturbance 
}

# compile into dataframe
winter <- as.data.frame(cbind(temp_regime, short, log(short)))

size_means <- vector() # vector for individual biomass data to go into
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te] # modify temperature +-10 C 
  temp_regime[te] <- mean(temp$Temperature) # mean temperature
  # run model
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "size")
  temp$Temperature <- temp$Temperature - temp_seq[te] # reset temperature regime
  size_means[te] <- colMeans(out) # get mean size of all stages
  size_means[te] <- sum(0.0094*(size_means[te])^2.754)  # multiply relative size (which is also biologically plausible) by Benke et al 1999 Table 2 a and b params (M(mg) = aL^b) 
}

# compile into dataframe
winter_size_means <- as.data.frame(cbind(temp_regime, size_means))

#---------------------------------
# For a summer time disturbance
#---------------------------------
temp_regime <- vector() # mean temperature of temperature regime
temp_means <- vector()  # mean abundance 
temp_seq <- seq(-10, 10, by = 1) # sequence to modify temperature regime
short <- vector()# mean abundance shortly post disturbance (6 time steps or 3 months)
discharge[272] <- 1 # set a bankfull discharge disturbance in summertime (May 25) 
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te]  # modify temperature +-10 C 
  temp_regime[te] <- mean(temp$Temperature) # mean temperature
  #run model
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
  temp$Temperature <- temp$Temperature - temp_seq[te] # reset temperature
  means.list.A<- mean.data.frame(out, burnin = 250, iteration = 2) # remove burn in 
  short[te] <- mean(means.list.A$mean.abund[23:29]) # get means right after disturbance
}

#compile into dataframe
summer <- as.data.frame(cbind(temp_regime,short,log(short)))


size_means <- vector() # vector for individual biomass data to go into
for (te in 1:length(temp_seq)){
  temp$Temperature <- temp$Temperature + temp_seq[te] # modify temperature +-10 C 
  temp_regime[te] <- mean(temp$Temperature) # mean temperature
  # run model
  out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 1, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "size")
  temp$Temperature <- temp$Temperature - temp_seq[te] # reset temperature
  size_means[te] <- colMeans(out) # get mean size of all stages
  size_means[te] <- sum(0.0094*(size_means[te])^2.754)  # multiply relative size (which is also biologically plausible) by Benke et al 1999 Table 2 a and b params (M(mg) = aL^b) 
}

# compile into dataframe
summer_size_means <- as.data.frame(cbind(temp_regime, size_means))

# bind together, 1 = winter 2 = summer
temp_dist_a <- bind_rows(winter, summer, .id = "season") # dataframe for winter and summer abundances

# now we want to calculate the difference between mean abundances, so wake summer - winter
deltatemp_a <- as.data.frame(cbind(rep(3, times = length(temp_regime)),temp_regime, summer[,3]-winter[,3]))

sizes <- rbind(winter_size_means, summer_size_means) # dataframe for winter and summer individual biomass
temp_size_a <- bind_rows(winter_size_means, summer_size_means, .id = "season") # dataframe for winter and summer individual biomass

# standing biomass is abundance * individual biomass, to get difference in standing biomass
# between winter and summer, we multiply to get standing biomas, then substract winter from summer
deltasize_a <- as.data.frame(cbind(rep(3, times = length(temp_regime)), temp_regime, (summer_size_means[,2]*summer[,2])-(winter_size_means[,2]*winter[2])))

# make a column for standing biomass produced just after a disturbance
temp_size_a <- mutate(.data = temp_size_a, size_means = temp_dist_a$short * sizes$size_means )

# now we want to compile abundances for summer, winter, and difference between summer and winter 
deltatemp_a <- setNames(deltatemp_a, names(temp_dist_a[c(1,2,4)])) #need to rename columns to match 
temp_dist_a <- rbind(temp_dist_a[c(1,2,4)], deltatemp_a) # combine

# now we want to compile standing biomass for summer, winter, and difference between summer and winter 
deltasize_a <- setNames(deltasize_a, names(temp_size_a))#need to rename columns to match 
temp_size_a <- rbind(temp_size_a, deltasize_a) # combine
