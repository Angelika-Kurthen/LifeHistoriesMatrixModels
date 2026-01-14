#############################
## Code to run model with and without disturbance
## No temperature changes
############################

#load library
library(lubridate)

# load in model for 4 life histories and bespoke functions 
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

#For each life history strategy, run for 100 years

# Moderate
A_out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
A_out <- mean.data.frame(A_out, burnin = 250, iteration = 2) # clean burn in 
# combine with timesteps, and label
A_annual <- cbind.data.frame(temp$dts[250:last(A_out$timesteps)], A_out, rep("A", length(A_out$mean.abund)))


#Boom 
B_out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
B_out <- mean.data.frame(B_out, burnin = 250, iteration = 2)# clean burn in 
# combine with timesteps, and label
B_annual <- cbind.data.frame(temp$dts[250:last(B_out$timesteps)], B_out, rep("B", length(B_out$mean.abund)))

# Fast
C_out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
C_out <- mean.data.frame(C_out, burnin = 250, iteration = 2)# clean burn in 
# combine with timesteps, and label
C_annual <- cbind.data.frame(temp$dts[250:last(C_out$timesteps)], C_out, rep("C", length(C_out$mean.abund)))


#Slow
D_out <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
D_out <- mean.data.frame(D_out, burnin = 250, iteration = 2)# clean burn in 
# combine with timesteps, and label
D_annual <- cbind.data.frame(temp$dts[250:last(C_out$timesteps)], D_out, rep("D", length(B_out$mean.abund)))

#match column names
colnames(A_annual) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(B_annual) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(C_annual) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(D_annual) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")

#combine into large dataframe and subset 
annual <- rbind(A_annual, B_annual, C_annual, D_annual)
annual$Date <- as.Date(annual$Date)
annual <- subset(annual, Date >= "2035-01-01" & Date <= "2035-12-31")

# ggplot(data = annual, aes(x = Date, y  =Abundance, color = Taxa))+
#   geom_point(size = 1, alpha = 0.5)+
#   geom_line(size = 1)+
#   xlab("Month")+
#   ylab("Abundance")+
#   scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

##----------------------------------------------------------------------
# Run a comparison between no disturbance and a small pulse disturbance
##----------------------------------------------------------------------
# threeyear = three years of model run, no disturbance, 2030 - 2033
# add small spate in 2035
# pulse = 1 year of model run, include the spate 

# select dates for spate (mid 2035)
selected_date <- temp$dts[temp$dts >= as.Date("2035-05-01") & temp$dts <= as.Date("2035-06-15")]
# set to 27% bankfull discharge
discharge[match(selected_date, temp$dts)] <- 0.27

#For each life history strategy, run

# Moderate
A_pout <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
A_pout <- mean.data.frame(A_pout, burnin = 250, iteration = 2)# clean burn in 
# combine with timesteps, and label
A_pulse <-  cbind.data.frame(temp$dts[250:last(A_pout$timesteps)], A_pout, rep("A", length(A_pout$mean.abund)))


#Boom
B_pout <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
B_pout <- mean.data.frame(B_pout, burnin = 250, iteration = 2)# clean burn in 
# combine with timesteps, and label
B_pulse <-  cbind.data.frame(temp$dts[250:last(B_pout$timesteps)], B_pout, rep("B", length(B_pout$mean.abund)))


#Fast
C_pout <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
C_pout <- mean.data.frame(C_pout, burnin = 250, iteration = 2)# clean burn in 
# combine with timesteps, and label
C_pulse <-  cbind.data.frame(temp$dts[250:last(C_pout$timesteps)], C_pout, rep("C", length(C_pout$mean.abund)))

#Slow
D_pout <- Dmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
D_pout <- mean.data.frame(D_pout, burnin = 250, iteration = 2)# clean burn in 
# combine with timesteps, and label
D_pulse <-  cbind.data.frame(temp$dts[250:last(D_pout$timesteps)], D_pout, rep("D", length(D_pout$mean.abund)))

#match column names
colnames(A_pulse) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(B_pulse) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(C_pulse) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
colnames(D_pulse) <- c("Date", "timesteps", "Abundance", "sd", "se", "Taxa")
# combine
pulse <- rbind(A_pulse, B_pulse, C_pulse, D_pulse)
# make sure dates are dates
pulse$Date <- as.Date(pulse$Date)

# 3 years, no disturbance
threeyear <- subset(pulse, Date >= "2030-01-01" & Date <= "2033-12-31")

# 1 year, spate in the summer
pulse <- subset(pulse, Date >= "2035-01-01" & Date <= "2035-12-31")
# 
# ggplot(data = pulse, aes(x = Date, y  =Abundance/10000, color = Taxa))+
#   geom_point(size = 1, alpha = 0.5)+
#   geom_line(size = 1)+
#   ylim(c(0,60))+
#   xlab("Month")+
#   ylab("Relativized Abundance")+
#   geom_vline(xintercept = as.numeric(as.Date("2035-03-13")),
#              color = "black",
#              lwd = 1,
#              linetype = "dotted") +
#   scale_color_manual(name = "Taxa", values=c("#66CCEE", "#228833", "#EE6677", "#AA3377"))+
#   scale_x_date(date_labels="%B", date_breaks  ="1 month")+
#   theme_bw()+
#   theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, angle = 45, size = 12.5),
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))
