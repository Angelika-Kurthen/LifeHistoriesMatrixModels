# Modeling the temperature-mediated link between demography and biomass in ectotherms
### Time-varying matrix population models for macroinvertebrates in disturbed environments
This repository contains data, scripts, and outputs associated with the manuscript: *Modeling the temperature-mediated link between demography and biomass in ectotherms* (working title).

We present a stage-structured population model that links temperature-mediated life-history tradeoffs to both abundance and biomass under changing temperature regimes and seasonal disturbance. Our model reveals emergent properties, such as population cycles and cohort structuring, typical of ectotherm populations. Under warming temperature scenarios, population size increased, but standing population biomass did not due to declining individual biomass. The seasonal timing of large pulse disturbances also interacted with temperature, altering post-disturbance recovery dynamics. 

## Repository Structure

```
/LifeHistoriesMatrixModels
│
├── Data/ # all data used in the modeling
│
├── Output/ # output produced in the modeling │   
│
├── Scripts/ # all scripts that can be run in R
│
├── Results/ 
│   └── FigGen.R # R script to produce figures in the manuscript
│
└── README.md   # This file
```

## Reproducibility Instructions

1. Clone the repository

```
git clone https://github.com/Angelika-Kurthen/LifeHistoriesMatrixModels.git
cd LifeHistoriesMatrixModels
```

2. Install required R packages
A complete list is found at the top of each script.

3. Run ```FigGen.R```, which sources all other scripts and data

For questions about the modeling workflow, data, or manuscript, please contact:

Angelika Kurthen

Email: akurthen@berkeley.edu
