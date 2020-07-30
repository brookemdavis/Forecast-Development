# run power model for fennel, maybe Nadina? to get idea of possible parameters

library(dplyr)
library(ggplot2)
library(TMB)
library(tmbstan)
library(R2jags)
source("Code/Functions.R")


##==================================================================================

#Get TMB model ready
# only need to compile if changed model
# dyn.unload(dynlib("Code/TMB/Single_Stock_Power"))
# TMB::compile("Code/TMB/Single_Stock_Power.cpp") # can't seem to compile if have TMBstan loaded
dyn.load(dynlib("Code/TMB/Single_Stock_Power"))
# simulate some power model data

SimData <- Sim_Power_SR_Data(leng=50, age=4, Sig_Mod = 0.2, true_a = 6, true_b=0.5,
                              hr_min = 0.2, hr_max = 0.8, lnorm_corr = F, max_pop = 100000)

  
SimDataDF <- data.frame(S = round(SimData$S), R = (SimData$R), Year = 1:length(SimData$S))
  
# Create DF to store true and fitted values
DataDF <- SimDataDF[, c("S", "R", "Year")]
DataDF$Fit <- exp(SimData$true_a + SimData$true_b*log(DataDF$S))
DataDF$Mod <- "True"
DataDF$CI_low <- DataDF$CI_up  <-  DataDF$Pred <- DataDF$Pred_low <- DataDF$Pred_up <- DataDF$Fit  
  

ggplot(data = DataDF, aes(x=S, y=Fit)) +
  geom_line(size = 1.5, col = "red") +
  geom_point(aes(x=S, y=R), col = "black") 

# run using TMB, without prior

TMB_No_Prior <- RunPower(Data = SimDataDF, 
                          Fitting_SW = "TMB", 
                          Priors = F, Name = "TMB_No_Prior")

# Now run with priors same as "true" values
TMB_True_Prior <- RunPower(Data = SimDataDF,  Fitting_SW = "TMB", 
                            Priors = T, Name = "TMB_True_Prior",
                            A_mean = SimData$true_a,
                            A_sig = SimData$true_a/4, 
                            B_mean = SimData$true_b, 
                            B_sig = SimData$true_b/4) 

# Run same model as mcmc using tmbstan
library(tmbstan)
tmbstan_True_Prior <- RunPower(Data = SimDataDF,  Fitting_SW = "tmbstan", 
                                Priors = T,  Name = "tmbstan_True_Prior",
                                A_mean = SimData$true_a,
                                A_sig =SimData$true_a/4, 
                                B_mean = SimData$true_b, 
                                B_sig = SimData$true_b/4) 

# Fit using JAGS
library(R2jags)
JAGS_True_Prior <- RunPower(Data = SimDataDF,  Fitting_SW = "JAGS", 
                             Priors = T,  Name = "JAGS_True_Prior",
                            A_mean = SimData$true_a,
                            A_sig = SimData$true_a/4, 
                            B_mean = SimData$true_b, 
                            B_sig = SimData$true_b/4) 

# Fit using Stan
Stan_True_Prior <- RunPower(Data = SimDataDF,  Fitting_SW = "Stan", 
                             Priors = T,   Name = "Stan_True_Prior",
                            A_mean = SimData$true_a,
                            A_sig = SimData$true_a/4, 
                            B_mean = SimData$true_b, 
                            B_sig = SimData$true_b/4) 


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

Posts <- data.frame(A_Scaled = numeric(), B = numeric(), Mod = character())

for(i in 1:length(Bayes_Mods)){
  New_Rows <- data.frame(A_Scaled = Bayes_Mods[[i]]$A_Post, B = Bayes_Mods[[i]]$B_Post, Mod = names(Bayes_Mods)[[i]])
  Posts <- bind_rows(Posts, New_Rows)
}

#  A is messed up by scale, need to put back into same scale as true
# Scale will be same across all models
Scale <- TMB_No_Prior$Scale
Posts <- Posts %>% mutate(A = A_Scaled + (1-B)*log(Scale))

# Also Add two tmb estimates
A_No_Prior_Scaled <- TMB_No_Prior$Ests %>% filter(Param == "A") %>% pull(Estimate)
A_True_Prior_Scaled <- TMB_True_Prior$Ests %>% filter(Param == "A") %>% pull(Estimate)
# Also Add two tmb estimates
B_No_Prior <- TMB_No_Prior$Ests %>% filter(Param == "B") %>% pull(Estimate) 
B_True_Prior <- TMB_True_Prior$Ests %>% filter(Param == "B") %>% pull(Estimate) 

# Fix these too
A_No_Prior <- A_No_Prior_Scaled + (1-B_No_Prior)*log(Scale)
A_True_Prior <- A_True_Prior_Scaled + (1-B_True_Prior)*log(Scale)

ggplot(Posts, aes(A, stat(density), col = Mod)) +
  geom_freqpoly() +
  geom_vline(xintercept = A_No_Prior, col = "black") +
  geom_vline(xintercept = A_True_Prior, col = "darkgrey") +
  geom_vline(xintercept = SimData$true_a, col = "black", linetype = "dashed") +
  theme_bw()

# look good, priors have small effect
# don't return exact true value though

ggplot(Posts, aes(B, stat(density), col = Mod)) +
  geom_freqpoly() +
  geom_vline(xintercept = B_No_Prior, col = "black") +
  geom_vline(xintercept = B_True_Prior, col = "darkgrey") +
  geom_vline(xintercept = SimData$true_b, col = "black", linetype = "dashed") +
  theme_bw()


