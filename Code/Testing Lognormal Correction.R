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

set.seed(100)#(1004)
ntrials <- 1000
RicPars <- matrix(NA, ntrials, 4)

for (i in 1:ntrials) {
  SimData <- Sim_Ricker_SR_Data(leng=50, age=4, Sig_Ricker = 0.8, true_a = 5, true_b=1/100000, #true_a = rnorm(1,5, 2)
                                hr_min = 0.2, hr_max = 0.8, lnorm_corr = T)
  
  SimDataDF <- data.frame(S = round(SimData$S), R = (SimData$R), Year = 1:length(SimData$S))
  
  #ggplot(SimDataDF, aes(x=S, y=R)) + geom_point() + coord_fixed()
  
  # Create DF to store true and fitted values
  DataDF <- SimDataDF[, c("S", "R", "Year")]
  DataDF$Fit <- SimData$true_a * DataDF$S * exp( -SimData$true_b * DataDF$S ) #* exp(Sig_Ricker^2/2)
  DataDF$Mod <- "True"
  DataDF$CI_low <- DataDF$CI_up  <-  DataDF$Pred <- DataDF$Pred_low <- DataDF$Pred_up <- DataDF$Fit
  
  # run using TMB, without prior
  # Switch BiasCorr between T and F and plot to see impact of bias correction in LL
  TMB_No_Prior <- RunRicker(Data = SimDataDF, 
                            Fitting_SW = "TMB", 
                            Priors = F, BiasCorr=T, Name = "TMB_No_Prior")
  #RicPars[i,] <-  TMB_No_Prior$Ests$Estimate[1:3] #logA, logSigma, logSMAX
  logA <- TMB_No_Prior$Ests %>% filter(Param=="logA") %>% pull(Estimate)
  logSmax <- TMB_No_Prior$Ests %>% filter(Param=="logSmax") %>% pull(Estimate)
  logSigma <- TMB_No_Prior$Ests %>% filter(Param=="logSigma") %>% pull(Estimate)
  RicPars[i,] <- c(logA, logSmax, logSigma, TMB_No_Prior$Scale)
}



colnames(RicPars) <- c("logA", "logSmax", "logSigma", "Scale")
RicPars <- as.data.frame(RicPars)

# Remove all MC trials with productivity less than 1 (logA < 0)
RicPars <- RicPars %>% filter(logA >= 0)
# Remove all MC trials with Smax > 100x the true Smax (1 million). I don't think this constratint is needed
RicPars <- RicPars %>% mutate(Smax = exp(logSmax)*Scale) #%>% filter(Smax <= 1000000)
ntrialsTuned <- length (RicPars$logA)


S <- NA
Rpred <- matrix(NA, 1500, ntrialsTuned)

for (i in 1:ntrialsTuned) {
  for (j in 1:1500){ #Assuming true smax ~ 100,000, and plot extends to 1.5 x Smax
    S[j] <- j*100
    Rpred[j,i] <-  exp( RicPars$logA[i] ) * S[j] * exp( -S[j] / RicPars$Smax[i] )
  }
}

Rpred_dist <- apply(Rpred, 1, quantile, probs = c(0.025, 0.5, 0.975))
EstAgg <- data.frame(S=S, Fit=Rpred_dist[2,], CI_low=Rpred_dist[1,], CI_up = Rpred_dist[3,])
EstAgg$Mod <- "TMB_No_Prior"


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
  scale_color_discrete(name = "", labels = c("Estimated curve w/o\nbias correction in LL", "True")) + 
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
 