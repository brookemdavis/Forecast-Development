# For each stock, fit Larken, Ricker, Recursive Bayes 
  # 1) current forecast priors
  # 2) mimicking same priors in TMB
  # 3) With no priors

# How far back to go? Maybe 10 years retrospective?
# For each get alpha, Smax, and forecast value for next year to compare]

source("Code/Functions.R")
library(dplyr)
library(tidyr)

# Read in Current SR data from 2020 forecast
SRdat <- read.csv("SRDATA2020.csv")