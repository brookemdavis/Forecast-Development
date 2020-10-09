# Test/Compare SR models using different software/priors

library(dplyr)
library(ggplot2)
library(TMB)
library(R2jags)
source("Code/Functions.R")


##==================================================================================

#Get TMB model ready
# only need to compile if changed model
#dyn.unload(dynlib("Code/TMB/Single_Stock_Larkin"))
#compile("Code/TMB/Single_Stock_Larkin.cpp") # can't seem to compile if have TMBstan loaded
dyn.load(dynlib("Code/TMB/Single_Stock_Larkin"))
# simulate some Larkin model data

SimData <- Sim_Larkin_SR_Data( leng=70, age=4, Sig_Larkin = 0.2, true_a = 4, 
                               true_b=c(1/10000, 0.00006, 0.00004, 0.00002),
                                hr_min = 0.2, hr_max = 0.8, lnorm_corr = F, autoCorr = F, rho=NA)

# get too perfect of cycles at beginning -- do "burn in" of 20 years -- do after get fits

SimDataDF <- data.frame(S = round(SimData$S), R = (SimData$R), Year = 1:(length(SimData$S)))
  
# Create DF to store true and fitted values
DataDF <- SimDataDF[, c("S", "R", "Year")]
DataDF$Fit <- NA
# How to do this for Larkin?
for(i in 20:70){
     DataDF$Fit[i] <- SimData$true_a * DataDF$S[i] * exp(-SimData$true_b[1]*DataDF$S[i] - 
                            SimData$true_b[2]*DataDF$S[i-1] - 
                            SimData$true_b[3]*DataDF$S[i-2] - 
                            SimData$true_b[4]*DataDF$S[i-3] )
}


DataDF$Mod <- "True"
DataDF$CI_low <- DataDF$CI_up  <-  DataDF$Pred <- DataDF$Pred_low <- DataDF$Pred_up <- DataDF$Fit  

DataDF <- DataDF %>% filter(Year > 20)
  
ggplot(data = DataDF, aes(x=S, y=Fit)) +
  geom_line(size = 1.5, col = "red") +
  geom_point(aes(x=S, y=R), col = "black") 

ggplot(data = DataDF, aes(x=Year, y=R)) +
  geom_line(size = 1.5) +
  geom_line(aes(x=Year, y=Fit), col = "red")


# run using TMB, without prior

TMB_No_Prior <- RunLarkin(Data = DataDF, 
                          Fitting_SW = "TMB", 
                          Priors = F, Name = "TMB_No_Prior")

# In forecast version all priors for betas are 0, 0.001


# Now run with priors same as "true" values
TMB_True_Prior <- RunLarkin(Data = DataDF,  Fitting_SW = "TMB", 
                            Priors = T, Name = "TMB_True_Prior",
                            logA_mean = log(SimData$true_a), logA_sig = 0.5, # priors on logAlpha
                            Sig_Gam_Dist = 0.001, # inverse gamma shape and scale param
                            B_means = SimData$true_b, B_sigs = SimData$true_b/4) 

# Run same model as mcmc using tmbstan
library(tmbstan)
tmbstan_True_Prior <- RunLarkin(Data = DataDF,  Fitting_SW = "tmbstan", 
                                Priors = T,  Name = "tmbstan_True_Prior",
                                logA_mean = log(SimData$true_a), logA_sig = 0.5, # priors on logAlpha
                                Sig_Gam_Dist = 0.001, # inverse gamma shape and scale param
                                B_means = SimData$true_b, B_sigs = SimData$true_b/4) 

# Fit using JAGS
library(R2jags)
JAGS_True_Prior <- RunLarkin(Data = DataDF,  Fitting_SW = "JAGS", 
                             Priors = T,  Name = "JAGS_True_Prior",
                             logA_mean = log(SimData$true_a), logA_sig = 0.5, # priors on logAlpha
                             Sig_Gam_Dist = 0.001, # inverse gamma shape and scale param
                             B_means = SimData$true_b, B_sigs = SimData$true_b/4) 

# Fit using Stan
Stan_True_Prior <- RunLarkin(Data = DataDF,  Fitting_SW = "Stan", 
                             Priors = T,   Name = "Stan_True_Prior",
                             logA_mean = log(SimData$true_a), logA_sig = 0.5, # priors on logAlpha
                             Sig_Gam_Dist = 0.001, # inverse gamma shape and scale param
                             B_means = SimData$true_b, B_sigs = SimData$true_b/4) 


All_Ests <- bind_rows(DataDF, TMB_No_Prior$Fit, TMB_True_Prior$Fit, tmbstan_True_Prior$Fit,
                      JAGS_True_Prior$Fit, Stan_True_Prior$Fit)

# Now plot all to compare 
ggplot(data = All_Ests, aes(x=S, y=Fit, ymin = CI_low, ymax = CI_up, col = Mod, fill= Mod)) +
  geom_line(size = 1.5) +
  geom_ribbon( alpha = 0.1) +
  geom_point(aes(x=S, y=R), col = "black") +
  geom_ribbon(aes(x=S, y=Pred, ymin = Pred_low, ymax = Pred_up, fill= Mod), 
              alpha = 0.05) +
  theme_bw()

#===========================================================================
# Also look at posteriors of alpha and beta from each bayesian model?

# compile all posts into DF
Bayes_Mods <- list("tmbstan_True_Prior" = tmbstan_True_Prior, 
                   "JAGS_True_Prior" = JAGS_True_Prior, 
                   "Stan_True_Prior" = Stan_True_Prior)

Posts <- data.frame(A = numeric(), Mod = character())#, B = numeric())

for(i in 1:length(Bayes_Mods)){
  New_Rows <- data.frame(A = Bayes_Mods[[i]]$A_Post, Mod = names(Bayes_Mods)[[i]])#, B = Bayes_Mods[[i]]$B_Post)
  Posts <- bind_rows(Posts, New_Rows)
}



# Also Add two tmb estimates
A_No_Prior <- TMB_No_Prior$Ests %>% filter(Param == "A") %>% pull(Estimate)
A_True_Prior <- TMB_True_Prior$Ests %>% filter(Param == "A") %>% pull(Estimate)


ggplot(Posts, aes(A, stat(density), col = Mod)) +
  geom_freqpoly() +
  geom_vline(xintercept = A_No_Prior, col = "black") +
  geom_vline(xintercept = A_True_Prior, col = "darkgrey") +
  geom_vline(xintercept = SimData$true_a, col = "black", linetype = "dashed") +
  theme_bw()

# JAGS looks very different -- but closer to true.

