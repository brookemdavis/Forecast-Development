# Look at lognormal correction

library(dplyr)
library(ggplot2)
library(TMB)
#library(tmbstan)
source("Code/Functions.R")

#================================================================================

# only need to compile if changed model
#dyn.unload(dynlib("Code/TMB/Single_Stock_Ricker"))
#TMB::compile("Code/TMB/Single_Stock_Ricker.cpp") # can't seem to compile if have TMBstan loaded
dyn.load(dynlib("Code/TMB/Single_Stock_Ricker"))

# Simulate some basic Ricker data - don't use lognormal correction

set.seed(1004)
ntrials <- 1#100
RicPars <- matrix(NA, ntrials, 4)
Fitting_SW <- "tmbstan"#"TMB"

for (i in 1:ntrials) {
  print(cat(i," of ", ntrials, "ntrials"))
  SimData <- Sim_Ricker_SR_Data(leng=50, age=4, Sig_Ricker = 0.8, true_a = 1.93, true_b=1/159744, #true_a = rnorm(1,5, 2)
                                hr_min = 0.25, hr_max = 0.35, lnorm_corr = T)
  # The code below has the expanded alpha value (Hilborn and Walters 1992), and gives the same simulated data as above with lnorm_corr= F. 
  # Note, Hilborn and  Walaters 1992 show corrections to alpha and SREP, but SMAX (1/b) remains constant
  # SimData <- Sim_Ricker_SR_Data(leng=50, age=4, Sig_Ricker = 0.8, true_a = exp(log(5)+0.8^2/2), true_b=1/100000, #true_a = rnorm(1,5, 2)
  #                              hr_min = 0.2, hr_max = 0.8, lnorm_corr = T)
  
  SimDataDF <- data.frame(S = round(SimData$S), R = (SimData$R), Year = 1:length(SimData$S))
  
  #ggplot(SimDataDF, aes(x=S, y=R)) + geom_point() + coord_fixed()
  
  # Create DF to store true and fitted values
  DataDF <- SimDataDF[, c("S", "R", "Year")]
  DataDF$Fit <- SimData$true_a * DataDF$S * exp( -SimData$true_b * DataDF$S ) 
  DataDF$Mod <- "True"
  DataDF$CI_low <- DataDF$CI_up  <-  DataDF$Pred <- DataDF$Pred_low <- DataDF$Pred_up <- DataDF$Fit
  
  # run using TMB, without prior
  # Switch BiasCorr between T and F and plot to see impact of bias correction in LL
  TMB_No_Prior <- RunRicker(Data = SimDataDF, 
                            Fitting_SW = Fitting_SW, #"tmbstan",#"TMB", 
                            Priors = T, BiasCorr=F, Name = "TMB_No_Prior")
  
  if( Fitting_SW == "tmbstan") {
    Ests_quant <- apply ( as.data.frame(TMB_No_Prior$Ests) , 2, quantile, probs = c(0.025, 0.5, 0.975) )
    logA <- as.data.frame(Ests_quant)$logA[2] # median
    logSmax <- as.data.frame(Ests_quant)$logSmax[2] # median
    logSigma <- as.data.frame(Ests_quant)$logSigma[2] # median
    RicPars[i,] <- c(logA, logSmax, logSigma, TMB_No_Prior$Scale)  
  }
  if( Fitting_SW == "TMB") {
    logA <- TMB_No_Prior$Ests %>% filter(Param=="logA") %>% pull(Estimate)
    logSmax <- TMB_No_Prior$Ests %>% filter(Param=="logSmax") %>% pull(Estimate)
    logSigma <- TMB_No_Prior$Ests %>% filter(Param=="logSigma") %>% pull(Estimate)
    RicPars[i,] <- c(logA, logSmax, logSigma, TMB_No_Prior$Scale)  
  }
  
}


# Set up data frame of Ricker paramters for all MC trials

colnames(RicPars) <- c("logA", "logSmax", "logSigma", "Scale")
RicPars <- as.data.frame(RicPars)
RicPars <- RicPars %>% mutate(Smax = exp(logSmax)*Scale) 

# Tune simulations by removing all MC trials with productivity less than 1 (logA < 0)

RicPars <- RicPars %>% filter(logA >= 0)
ntrialsTuned <- length (RicPars$logA)


# Calculate predicted recruitments for all MC trials along a range of spawner abundances
S <- NA
Rpred <- matrix(NA, 2000, ntrialsTuned)

for (i in 1:ntrialsTuned) {
  for (j in 1:2000){ #Assuming true smax ~ 100,000, and plot extends to 2 x Smax
    S[j] <- j*100
    Rpred[j,i] <-  exp( RicPars$logA[i] ) * S[j] * exp( -S[j] / RicPars$Smax[i] )
  }
}

# Calculate 5th, 50th, and 95th percentiles of the distribution of recruitment along
# the range of spawner abundances

Rpred_dist <- apply(Rpred, 1, quantile, probs = c(0.025, 0.5, 0.975))
Rpred_mean <- apply(Rpred, 1, mean) 


# Create a dataframe of these percentiles
EstAgg <- data.frame(S=S, Fit=Rpred_dist[2,], CI_low=Rpred_dist[1,], CI_up = Rpred_dist[3,])
EstAgg$Mod <- "TMB_No_Prior"

# I get the same answer when I plot mean of the disribiton of recruitments instead of the median
# EstAgg <- data.frame(S=S, Fit=Rpred_mean, CI_low=Rpred_dist[1,], CI_up = Rpred_dist[3,])


if(ntrials ==1){
  All_Ests <- bind_rows(DataDF, TMB_No_Prior[[1]])
}

if(ntrials > 1){
  All_Ests <- bind_rows(DataDF, EstAgg)
}

# Now plot all to compare 
ggplot(data = All_Ests, aes(x=S, y=Fit, ymin = CI_low, ymax = CI_up, col = Mod, fill= Mod)) +
  geom_line(size = 1.5) +
  geom_ribbon( alpha = 0.1) +
  geom_point(aes(x=S, y=R), col = "black") +
  #geom_ribbon(aes(x=S, y=Pred, ymin = Pred_low, ymax = Pred_up, fill= Mod), 
  #            alpha = 0.05) +
  scale_color_discrete(name = "", labels = c("Estimated curve w/o\nbias correction in LL", "OM")) + 
  guides(fill = FALSE)  +
  xlab("Spawners") + 
  ylab("Recruitment") +
  theme_bw(base_size=16)


# now take estimated values and simulate using those, with and without correction

# now re-simulate data from these estimates
Ests <- TMB_No_Prior[[2]]

SimData_NoCorr <- Sim_Ricker_SR_Data(leng=50, age=4, Sig_Ricker = Ests$Estimate[Ests$Param=="sigma"], 
                               true_a = Ests$Estimate[Ests$Param=="A"], 
                               true_b = 1/(Ests$Estimate[Ests$Param=="Smax"]*TMB_No_Prior[[3]]),
                               hr_min = 0.2, hr_max = 0.8, lnorm_corr = F)

NoCorr_DF <- data.frame(S = round(SimData_NoCorr$S), R = (SimData_NoCorr$R), 
                        Year = 1:length(SimData_NoCorr$S), Mod = "NoCorr_True")


SimData_Corr <- Sim_Ricker_SR_Data(leng=50, age=4, Sig_Ricker = Ests$Estimate[Ests$Param=="sigma"], 
                             true_a = Ests$Estimate[Ests$Param=="A"], 
                             true_b = 1/(Ests$Estimate[Ests$Param=="Smax"]*TMB_No_Prior[[3]]),
                             hr_min = 0.2, hr_max = 0.8, lnorm_corr = T)

Corr_DF <-  data.frame(S = round(SimData_Corr$S), R = (SimData_Corr$R), 
                       Year = 1:length(SimData_Corr$S), Mod = "Corr_True")

Data_All <- bind_rows(DataDF, NoCorr_DF, Corr_DF)

ggplot(Data_All, aes(x=S, y=R, col = Mod)) + geom_point() + coord_fixed()


# Now fit to both and compare
# run using TMB, without prior, and estimating median curve (without bias correction in LL)
TMB_No_Corr <- RunRicker(Data = SimData_NoCorr,
                          Fitting_SW = "TMB",
                          Priors = F, BiasCorr = F, Name = "TMB_No_Corr")
TMB_Corr <- RunRicker(Data = SimData_Corr, 
                         Fitting_SW = "TMB", 
                         Priors = F, BiasCorr = F, Name = "TMB_Corr")

Data_All <- Data_All %>% filter (Mod=="True") #Removed simulated data as plot was too busy; keep curves estimated from simulated data
All_Ests <- bind_rows(Data_All, TMB_No_Corr[[1]], TMB_Corr[[1]])

 ggplot(data = All_Ests, aes(x=S, y=Fit, ymin = CI_low, ymax = CI_up, col=Mod, fill= Mod)) +
   geom_line(size = 1.5) +
   geom_ribbon( alpha = 0.1) +
   geom_point(aes(x=S, y=R, col = Mod) ) + #c(rep("black",150)) ) ) +
   #geom_ribbon(aes(x=S, y=Pred, ymin = Pred_low, ymax = Pred_up, fill= Mod), 
  #             alpha = 0.05) +
   #scale_color_discrete(name = "", labels = c("Projected with bias correction", "Projected w/o bias correction", "True: mean")) + 
   #scale_fill_discrete(name = "", labels = c("Projected with bias correction", "Projected w/o bias correction", "True: mean")) + 
   #scale_color_manual(values = RColorBrewer::brewer.pal(3, "Set2"), name = "", labels = c("Projected with bias correction", "Projected w/o bias correction", "True: mean")) + 
   #scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Set2"), name = "", labels = c("Projected with bias correction", "Projected w/o bias correction", "True: mean")) + 
   scale_color_manual(values = c( "#FC8D62", "#66C2A5", "#8DA0CB"), name = "", labels = c("Projected with bias correction", "Projected w/o bias correction", "True: mean")) + 
   scale_fill_manual(values =c( "#FC8D62", "#66C2A5", "#8DA0CB"), name = "", labels = c("Projected with bias correction", "Projected w/o bias correction", "True: mean")) + 
   #scale_color_manual(values = c( "#FC8D62", "#66C2A5", "black"), name = "", labels = c("Projected with bias correction", "Projected w/o bias correction", "True: mean")) + 
   #scale_fill_manual(values =c( "#FC8D62", "#66C2A5", "black"), name = "", labels = c("Projected with bias correction", "Projected w/o bias correction", "True: mean")) + 
   xlab("Spawners") + 
   ylab("Recruitment") +
   xlim(0,250000) + 
   theme_bw(base_size=16)
 
c( "#FC8D62", "#66C2A5", "#8DA0CB")
 
 # do this 100 times and see what difference in estimtes or trajectories are?
 