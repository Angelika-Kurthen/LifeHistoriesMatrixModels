######################
# Multivoltinism Fast Life History
######################
library(pracma)
# load Fast life history strategy model + bespoke functions
source("Scripts/Fast_1sp_Model.R")

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
  out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), stage_output = "3")
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
C.oneyear <- mlist[which(mlist$Date >= "2035-01-01" & mlist$Date <= "2035-12-31"), ]
C.oneyear$Date <- as.Date(C.oneyear$Date)
C.oneyear$Strategy <- rep("Fast", times = length(C.oneyear$Date))

# go over a series of temperatures 
seqs <- seq(-10, 10, by = 1)

for (i in 1:length(seqs)){
  temp$Temperature <- temp$Temperature + seqs[i]
  out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, 
                Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, 
                peakeach = length(temp$Temperature), stage_output = "3")
  m <- rowMeans(out)
  m <- cbind.data.frame(temp$dts, m[-1], rep(seqs[i], length(temp$dts)))
  assign(paste0("m", seqs[i]), m)
  temp$Temperature <- temp$Temperature - seqs[i]
}

seq_names <- paste0("m", seqs)
mseq_c <- do.call(rbind, mget(seq_names))
colnames(mseq_c) <- c("Date", "Abund", "MeanTemp")
mseq_c$Date <- as.Date(mseq_c$Date)
mseq_c$Strategy <- "Fast"

# subset one year 
mseq_c_1 <- mseq_c[which(mseq_c$Date >= "2035-01-01" & mseq_c$Date <= "2035-12-31"), ]

# The simulation starts 2022-01-01, so burn-in ends ~2032
mseq_c <- mseq_c %>%
  filter(Date >= "2032-01-01") %>%
  mutate(Year = year(Date),
         log_abund = log(Abund))

# Peak Detection Function
detect_peaks <- function(abund_vector, dates, min_height = NULL, min_distance = 6) {
  peaks_mat <- findpeaks(abund_vector,
                         minpeakheight  = min_height,
                         minpeakdistance = min_distance)
  
  if (is.null(peaks_mat)) return(data.frame(date      = as.Date(character()),
                                            abund     = numeric(),
                                            peak_idx  = integer()))
  data.frame(
    date     = dates[peaks_mat[, 2]],
    abund    = peaks_mat[, 1],
    peak_idx = peaks_mat[, 2]
  )
}

MIN_HEIGHT_LOG  <- log(50)
MIN_DISTANCE_TS <- 6

# Detect peaks per Strategy x MeanTemp x Year
peak_results_c <- mseq_c %>%
  group_by(Strategy, MeanTemp, Year) %>%
  group_modify(~{
    detect_peaks(
      abund_vector = .x$log_abund,
      dates        = .x$Date,
      min_height   = MIN_HEIGHT_LOG,
      min_distance = MIN_DISTANCE_TS
    )
  }) %>%
  ungroup()

# Count peaks per Strategy x MeanTemp x Year
peak_counts_by_year_c <- peak_results_c %>%
  dplyr::count(Strategy, MeanTemp, Year, name = "n_peaks")

# summarize to mean peaks per Strategy x MeanTemp
peak_counts_c <- peak_counts_by_year_c %>%
  group_by(Strategy, MeanTemp) %>%
  dplyr::summarise(
    mean_peaks = mean(n_peaks),
    sd_peaks   = sd(n_peaks),
    .groups    = "drop"
  )%>%
  complete(
    Strategy,
    MeanTemp = seqs
  ) %>%
  mutate(
    mean_peaks = ifelse(is.na(mean_peaks), 0, mean_peaks),
    sd_peaks   = ifelse(is.na(sd_peaks), 0, sd_peaks)
  )



peak_results_c_1 <- mseq_c_1 %>%
  mutate(log_abund = log(Abund)) %>%
  group_by(Strategy, MeanTemp) %>%
  group_modify(~{
    peaks <- detect_peaks(
      abund_vector  = .x$log_abund,
      dates         = .x$Date,
      min_height    = MIN_HEIGHT_LOG,
      min_distance  = MIN_DISTANCE_TS
    )
    peaks
  }) %>%
  ungroup()

peak_counts_c_1 <- peak_results_c_1 %>%
  dplyr::count(Strategy, MeanTemp, name = "n_peaks")%>%
  complete(
    Strategy,
    MeanTemp = seqs
  ) %>%
  mutate(
    mean_peaks = ifelse(is.na(n_peaks), 0, n_peaks)
  )


