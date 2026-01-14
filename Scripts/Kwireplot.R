################################
## Wireplot for Q, K, and tau
################################
# Code taken and tweaked from McMullen et al., 2017. 

## Minimum threshold of what is considered a flood
Qmin <- 0.25 
## Half saturation constant
a <- 0.01
## Rate that K returns to pre-disturbance level
g <- 0.1

## Maximum flood size to run model to
Qmax = 2

## Timesteps to run model out to - days
t = 55

dfK <- data.frame(Q = seq(0, Qmax, by = 0.002))
for(i in 1:length(dfK$Q)){
  if (dfK$Q[i] < Qmin) {
    dfK$Q[i] <- 0
  } else {
    dfK$Q[i] <- (dfK$Q[i] - Qmin)/(a + dfK$Q[i]- Qmin)
  }}

KQT <- data.frame(Q = numeric(((Qmax/0.002) + 1) * (t + 1)))
## Filling in Q
KQT$Q <- rep(seq(0, Qmax, by =0.002), each = t + 1)

for(i in 1:length(KQT$Q)){
  if (KQT$Q[i] < Qmin) {
    KQT$Qf[i] <- 0
  } else {
    KQT$Qf[i] <- (KQT$Q[i] - Qmin)/(a + KQT$Q[i]- Qmin)
  }}


for (i in 1:length(KQT$Q)){
  # Calculate K arrying capacity immediately following the disturbance
  KQT$K0[i]  <- 10000 + ((40000-10000)*KQT$Qf[i])
}

KQT$t <- rep(seq(0, t), (Qmax/0.002) + 1)

for (i in 1:length(KQT$Q)){
  KQT$K[i] <- 10000 + (KQT$K0[i] - 10000) * exp(-g * KQT$t[i])}
