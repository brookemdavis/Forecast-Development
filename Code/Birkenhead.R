# Aug 31 2020

# Birkenhead analysis for Karen Rickards

library(dplyr)
library(tidyr)
library(TMB)
library(ggplot2)

#read in data

SRdat <- read.csv("DataIn/SRDATA2020.csv")
AllData <- SRdat %>% filter(PopID == 10)


plot(AllData$yr, AllData$eff*1000000, type = "l", 
     ylab = "Effective Female Spawners", xlab = "Brood Year")


# but these numbers are eff, not total esc, which I think was used to get BM's

# fit Ricker to eff and put those BMs on top

Data <- AllData %>% filter(is.na(rec) == F)

#compile("Code/TMB/Single_Stock_Ricker.cpp")
dyn.load(dynlib("Code/TMB/Single_Stock_Ricker"))

Scale <- 1

data <- list() #data inputs
data$S <- Data$eff/Scale 
data$logR <- log(Data$rec/Scale) 
data$logA_mean <- 1
data$Sig_Gam_Dist <- 1
data$logSmax_mean <- -2
# set Bayes to 0, switch to one if using tmbstan
data$Bayes <- 0

# set up starting values
param <- list()
param$logA <- 1
data$logR <- log(Data$rec/Scale)
data$Scale <- Scale
data$Priors <- 0
param$logSmax <- log(as.numeric(quantile(Data$eff, 0.8)/Scale))

data$logA_sig <- 0.1
data$logSmax_sig <- 0.1/Scale
param$logSigma <- -2

# Fit model
obj <- MakeADFun(data, param, DLL="Single_Stock_Ricker", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))

# pull out estimates from ML fit
# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

# not sure why sgen calc not working in tmb, calculate outside

SMSY <- All_Ests$Estimate[All_Ests$Param == "SMSY"]*1000000
alpha <- All_Ests$Estimate[All_Ests$Param == "A"]
beta <- 1/(All_Ests$Estimate[All_Ests$Param == "Smax"]*1000000)

# solve for Sgen
ObjectiveSGen <- function( S, SpMSY, alpha, beta ) {
  # Recruits to get to SMSY (Holt et al. 2009)
  R <- alpha * S * exp( -beta * S )
  # Residual difference between ln(SMSY) and ln(R)
  delta <- log( SpMSY ) - log( R )
  # Calculate the negative log-likelihood
  negLL <- -1 * dnorm( x=delta, mean=0, sd=1, log=TRUE )
  # Return the negative log-likelihood (this is the value to minimize)
  return( negLL )
}  # End ObjectiveSGen function

opt <- optimize( f=ObjectiveSGen, interval=c(0, SMSY), SpMSY=SMSY,
                 alpha=alpha, beta=beta )

# Get SGen from the optimized output (i.e., minimum neg log-like)
SGen <- opt$minimum


# Sgen is tiny becayse productivity is very high

# what is last gen average?

RecentGen_avg <- AllData %>% filter(yr %in% 2014:2017) %>%
                          summarise(Eff_Avg = mean(eff)*1000000)
# about 17000 -- well below smsy, but above very low Sgen

pdf("Figures/Birkenhead.pdf")
plot(AllData$yr, AllData$eff*1000000, type = "l", 
     ylab = "Effective Female Spawners", xlab = "Brood Year")
abline(h = SMSY, col = "darkgreen")
abline(h=SGen, col = "red")
dev.off()

# make plot to show age-at-recruitment

Rec_By_Age <- AllData %>%  mutate("3" = rec3/rec, "4" = rec4/rec, "5" = rec5/rec) %>% 
                           filter(is.na("rec")==F) %>%
                          select(yr, "3", "4", "5") %>% 
                         pivot_longer(cols =  -c(yr),  names_to = "Age",  values_to = "Prop")
                          

ggplot(Rec_By_Age ,aes(x = yr, y = Prop, fill = Age)) + 
  geom_bar(position = "fill", stat = "identity") +
  theme_bw()



