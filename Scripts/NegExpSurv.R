##############################
# Code to calculate negative exponential fit of survivorships
##############################

# load libraries
library(minpack.lm)
library(readxl)

# Function to fit a negative exponential survival–flow relationship  
flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  df <- rbind(df, c(Qmin, 1))  # make sure we specify that we have 100% survival at Qmin (0.25)
  nls.fit <- nlsLM(formula = (y ~ k* exp(-h*x)), data = df, start = c(k = 1, h = 1)) # fit to negative exponential function (see Equation 12)
  return(nls.fit)
}

# Function to generate a smooth survival curve from fitted parameters
flow.surv.rate <- function(h, k, max, min, interval, Qmin) {
  Q <- seq(min, max, by = interval)   # Generate a sequence of flow magnitudes
  surv <- k*(exp(-h*Q))   # Calculate survival at each flow magnitude
  surv.df <- as.data.frame(cbind(Q, surv)) 
  surv.df$surv[which(surv.df$Q <= Qmin)] <- 1   # Enforce full survival below the minimum flow threshold
  return(surv.df)
}

# High flow sensitivity scenario
# Assumed flow magnitudes and mortality rates
mag <- c(0.3, 0.5, 0.75, 1, 3)
mort <- c(0.3, 0.5, 0.7, 0.75, 0.9)

# Fit the negative exponential survival model
high <- flow.surv.fit(mag, mort, 0.25)
# Generate a smooth survival curve using fitted parameters
high.df <- flow.surv.rate(high$m$getPars()[2], high$m$getPars()[1],2, 0.001, 0.001, 0.25 )
# Add response category label
high.df <- as.data.frame(cbind(high.df, rep("High Response", times = length(high.df$Q))))
colnames(high.df) <- c("Q", "Survival", "Response") # Standardize column names



# Medium flow sensitivity scenario
# Assumed flow magnitudes and mortality rates
mag <- c(0.75, 1, 3)
mort <- c(0.6, 0.65, 0.75)
# Fit the negative exponential survival model
med <- flow.surv.fit(mag, mort, 0.25)
# Generate a smooth survival curve using fitted parameters
med.df <- flow.surv.rate(med$m$getPars()[2], med$m$getPars()[1],2, 0.001, 0.001, 0.25 )
med.df <- as.data.frame(cbind(med.df, rep("Medium Response", times = length(med.df$Q))))
colnames(med.df) <- c("Q", "Survival", "Response") # Standardize column names



# Low flow sensitivity scenario
# Assumed flow magnitudes and mortality rates
mag <- c(0.75, 1, 2, 3)
mort <- c(0.4, 0.45, 0.55, 0.6)
# Fit the negative exponential survival model
low <- flow.surv.fit(mag, mort, 0.25)
# Generate a smooth survival curve using fitted parameters
low.df <- flow.surv.rate(low$m$getPars()[2], low$m$getPars()[1],2, 0.001, 0.001, 0.25 )
low.df <- as.data.frame(cbind(low.df, rep("Low Response", times = length(low.df$Q))))
colnames(low.df) <- c("Q", "Survival", "Response") # Standardize column names

# Combine all response curves into a single data frame
df <- as.data.frame(rbind(high.df, med.df, low.df))



# Load empirical Hydropsyche temperature survival data
HYOSSurvRates <- read_excel("Data/VitalRates.xlsx", sheet = "Hydropsyche Survival Rates ")
HYOSSurvRates <- as.data.frame(HYOSSurvRates)

# Temperature-dependent survival function 
TempSurv <- function(n){
  if (n <= 0.5){   # Enforce zero survival below a minimum temperature threshold
    a <- 0
  }else{ # Calculate survival using a negbinom distribution and scale it to match observed max survival
    a <-  dnbinom(as.integer(-n + 29), size = 2.8835371 , prob = 0.1932115)*(max(HYOSSurvRates$Survival)/max(dnbinom(as.integer(-HYOSSurvRates$Temperature + 32), size =2.8835371, prob = 0.1932115 )))
  }
  return((a))
}
