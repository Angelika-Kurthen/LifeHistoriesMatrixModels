##########################
# Moderate Life History model
###########################

library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)

# source bespoke functions
source("Scripts/1spFunctions.R")

Cmodel <- function(flow.data, temp.data, baselineK, disturbanceK, Qmin, extinct, iteration, peaklist = NULL, peakeach = NULL, fecundity = 500, dds = 900, stage_output = "all", dens.dep = T){
  
  # set up model
  source("Scripts/NegExpSurv.R")
  Q <- as.numeric(flow.data)
  temps <- temp.data
  
  degreedays <- as.data.frame(cbind(temps$dts, temps$Temperature * 14))
  colnames(degreedays) <- c("dts", "DegreeDay")
  degreedays$DegreeDay[degreedays$DegreeDay<0] <- 0
  degreedays$dts <- as.Date(degreedays$dts, origin = "1970-01-01")
  
  # need to make ramped increasing hydropeaking index 
  hp <- c(rep(peaklist, each = peakeach))
  
  # specify iterations
  iterations <- iteration
  
  # baseline K in the absence of disturbance
  Kb <- as.numeric(baselineK)
  # max K after a big disturbance
  Kd <- as.numeric(disturbanceK)
  
  # specify baseline transition probabilities for each species at mean temps
  G1 =  0.125#move to Stage2 (subimago)
  G2 =  0.075 #move to Stage3 (adult)
  P1 =  0.75#stay in Stage1 (larvae)
  P2 =  0.75#stay in Stage2 (subimago)
  
  # want to run this for one year, in 14 day timesteps 
  timestep <- seq(2, (length(temps$Temperature) + 1), by = 1)

  # create array to put the total N of all species into
  Total.N <- array(0,
                   dim  <-c((length(timestep) +1 ), iterations),
                   dimnames <- list(1:(length(timestep) + 1), 1:iterations))
  
  # create list of arrays w/ abundance data for each spp
  output.N.list <- array(0,
                    
                    dim = c(length(timestep) + 1, 3, iterations),
                    dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations)
  )
  
  sizelist <- array(0,
                    
                    dim = c(length(timestep) + 1, 3, iterations),
                    dimnames = list(1:(length(timestep)+1), c("S1", "S2", "S3"), 1:iterations)
  )
  
  Qmin <- Qmin
  a <- 0.01
  g <- 0.075
  h <- med$m$getPars()[2]  
  k <- med$m$getPars()[1] 
  
  iterlist <- c(1:iterations)
  extinction <- extinct
  
  #-------------------------
  # Outer Loop of Iterations
  #--------------------------
  for (iter in c(1:iterations)) {
         K = Kb # need to reset K for each iteration
    
    output.N.list[1,1:3, iter]<- c(5000, 3000, 100)
    
    # we often want to look at different parameter values after we run code, so we create some lists
    # list to input Ks
    Klist <- vector()
    Klist[1] <- 10000
    
    # list to imput flow morts
    flowmortlist <- vector()
    
    Flist <- vector()
    
    emergetime <- vector()

    TempSurvival <- vector()
    for(c in temps$Temperature){
      
      b <- TempSurv(c)
      
      TempSurvival <- append(TempSurvival, b)
    }
    #-------------------------
    # Inner Loop of Timesteps
    #-------------------------
    for (t in timestep) {
      
      #----------------------------------------------------------
      # Calculate how many timesteps emerging adults have matured
      emergetime <- append(emergetime, back.count.degreedays(t, dds, degreedays)) # value from Sweeney et al 2017
      #---------------------------------------------------------
      # Calculate fecundity per adult
      
      # we start by pulling fecundities from normal distribution
      # assuming 50 50 sex ration, 0.22 of egg masses 'dissapearred', and 0.2 desiccation because of rock drying
      F3 = fecundity  * hydropeaking.mortality(0.8, 0.99, h = hp[t-1])
      
      # we can also relate fecundities to body mass.
      # in order to iterate through different fecundities
      # emergetimes for our temp regime are between 3 and 9 
      # create a lm for that data, with +10% and -10% of fecundity
      x <- c(5,13)
      y <- c(fecundity*0.9, fecundity*1.1)
      mod <- lm(y~x)
      
      if (t > 19) { # will be erased in burn
        size <- emergetime[t-1]
        sizes <- c(mean(c(1, emergetime[t-1]/2), na.rm = T),
                   mean(c(emergetime[t-1]/2, emergetime[t-1]), na.rm = T),
                   emergetime[t-1])
        sizelist[t, 1:3, iter] <- sizes 
        F3 <- ((size*mod$coefficients[2])+mod$coefficients[1]) * hydropeaking.mortality(0.8, 0.99, h = hp[t-1])
      }
      #--------------------------------------------------
      # Calculate the disturbance magnitude-K relationship
      # Sets to 0 if below the Qmin
      Qf <- Qf.Function(Q[t-1], Qmin, a)
      
      #-------------------------------------------------------------------
      # Calculate K carrying capacity immediately following the disturbance
      K0 <- K + ((Kd-K)*Qf)
      
      # Calculate final K for timestep, including relationship between K and time since disturbance
      K <- post.dist.K(K0, Kb, g, t, Q, Qmin)
      
      Klist <- append(Klist, K)
      #---------------------------------------------
      # Calculate effect of density dependence on fecundity 
      if (dens.dep == T){
      # Logistic via Rogosch et al. Fish Model
      # no immediate egg mortality incorporated
      F3 <- Logistic.Dens.Dependence(F3, K, Total.N[t-1, iter])
      # 
      # add F_BAET to list
      Flist <- append(Flist, F3)}
      
      if (dens.dep == F){
        F3 <- F3
      }
      #-----------------------------------------------
      # Calculate new transition probabilities based on temperature
      # This is the growth v development tradeoff
      
      # development measures
      # in this function, we assume that if below the min temp threshold (9) no maturation occurs (slow maturation, large growth)
      # if above the max temp threshold (15), no one remains more than 1 timestep in each stage (fast maturation, small growth)
      if (is.na(emergetime[t-1]) == F) {
      P1 <- (1-(1/((emergetime[t-1])/2))) *TempSurvival[t-1]
      P2 <- (1-(1/((emergetime[t-1])/2)))*TempSurvival[t-1]
      G1 <- (0.9/((emergetime[t-1])/2))*TempSurvival[t-1]
      G2 <- (0.7/((emergetime[t-1])/2))*TempSurvival[t-1]
    }

      
      if (is.na(emergetime[t-1]) == T) {
        G1 <- (0.9/((-0.72 * temps$Temperature[t-1]) + 19.54))*TempSurvival[t-1]
        P1 <- (1-(1/((-0.72 * temps$Temperature[t-1]) + 19.54)))*TempSurvival[t-1]
        G2 <- (0.7/((-0.72 * temps$Temperature[t-1]) + 19.54))*TempSurvival[t-1]
        P2 <- (1-(1/((-0.72 * temps$Temperature[t-1]) + 19.54)))*TempSurvival[t-1]
      }
      
      if (G1 > 1) G1 <- 1
      if (G1 < 0) G1 <- 0
      if (G2 > 1) G2 <- 1
      if (G2 < 0) G2 <- 0
      if (P1 > 1) P1 <- 1
      if (P2 < 0) P2 <- 0
      
      #-----------------------------------------------
      # Create Lefkovitch Matrix
      
      T1 <- c(P1, 0, F3)
      T2 <- c(G1, P2, 0)
      T3 <- c(0, G2, 0) 
      
      A <- rbind( T1, T2, T3)
      
      #--------------------------------------
      # Calculate abundances for each stage
      
      output.N.list[t, 1:3, iter] <- A %*% output.N.list[t-1, 1:3, iter] 
      
      #------------------------------------------
      #Calculate immediate mortality due to flows
      # mortality due to flooding follows N0 = Nz*e^-hQ
      
      #s1
      output.N.list[t, 1, iter] <- flood.mortality(output.N.list[t, 1, iter], k, h, Q[t-1], Qmin)
      #s2Q
      output.N.list[t,2,iter] <- flood.mortality(output.N.list[t,2,iter], k, h, Q[t-1], Qmin)
      
      output.N.list[t,3,iter] <- flood.mortality(output.N.list[t,3,iter], k, h, Q[t-1], Qmin)
      
      flowmortlist <- append(flowmortlist, flood.mortality(1, k, h, Q[t-1], Qmin))
      
      #------------------------------------------------------
      # if any values infinity, turn to highest integer, since model breaks and considered Inf non numeric 
      if( any(is.infinite(output.N.list[t,,iter]))== T){
        output.N.list[t, 1:3, iter] <- 5.6e+307 # max value div by 3
        Total.N[t, iter] <- 1.7e+308 ## max integer value in R, if we allow it to be INF, causes errors
      }
      else {
        # check extinction threshold and if below set to 0
        Total.N[t,iter] <- sum(output.N.list[t,,iter])
        if (Total.N[t, iter] < extinction){
          output.N.list[t,,iter] <- 0
          Total.N[t, iter] <- 0}
      }
    } #-------------------------
    # End Inner Loop  
    #------------------------- 
  } #----------------------
  # End Outer Loop
  #----------------------
  if (stage_output == "all"){
    return(output.N.list[ , 1:3, ])
  }
  if (stage_output == "3"){
    return(output.N.list[ , 3, ])
  }
  
  if (stage_output == "size"){
    return(sizelist)
  }
}

