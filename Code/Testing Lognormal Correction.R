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
ntrials <- 1000#10000#1#100
RicPars <- matrix(NA, ntrials, 4)
RicPars_nBC <- matrix(NA, ntrials, 4) #Estimated Ricker pars with NO bias correction
RicPars_wBC <- matrix(NA, ntrials, 4) #Estimated Ricker pars with bias correction
Fitting_SW <- "TMB"#"tmbstan"#"TMB"
Rdist <- matrix(NA, ntrials,4)
outTrueR <- list()
#true_UMSY <- ( 0.5 * log(true_a) - 0.07*log(true_a)^2 )
#true_SMSY <- ( 0.5 * log(1.93) - 0.07*log(1.93)^2 ) = 0.2984
# true_a <- 1.93
# sig_Ricker <-0.76
# hr_min <- 0.25
# hr_max <- 0.35

for (i in 1:ntrials) {
  #print(cat(i," of ", ntrials, "ntrials"))
  # Assuming true parameters are median-unbiased so lnorm_corr = F in initial data generation. 
    # If we instead assume true parmaeters are mean-unbiased and lnorm_corr=T in initital data generation,
    # results in terms of relative errors in distribution of recruits are similar but ave "true" recruits are lower
  SimData <- Sim_Ricker_SR_Data(leng=30, age=4, Sig_Ricker = 0.9, true_a = 1.93, true_b=1/159744, #true_a = rnorm(1,5, 2)
                                hr_min = 0.25, hr_max = 0.35, EscPolicy = F, constEsc = 75000, lnorm_corr = F)
  # The code below has the expanded alpha value (Hilborn and Walters 1992), and gives the same simulated data as above with lnorm_corr= F. 
  # Note, Hilborn and  Walters 1992 show corrections to alpha and SREP, but SMAX (1/b) remains constant
  # SimData <- Sim_Ricker_SR_Data(leng=50, age=4, Sig_Ricker = 0.8, true_a = exp(log(5)+0.8^2/2), true_b=1/100000, #true_a = rnorm(1,5, 2)
  #                              hr_min = 0.2, hr_max = 0.8, lnorm_corr = T)
  
  # SimDataDF <- data.frame(S = round(SimData$S), R = (SimData$R), Year = 1:length(SimData$S))
  SimDataDF <- data.frame(S = SimData$S, R = (SimData$R), Year = 1:length(SimData$S))
  
  outTrueR[[i]] <- data.frame(i=i, R=SimData$R) 
  #ggplot(SimDataDF, aes(x=S, y=R)) + geom_point() + coord_fixed()

  # Create DF to store true and fitted values
  DataDF <- SimDataDF[, c("S", "R", "Year")]
  DataDF$Fit <- SimData$true_a * DataDF$S * exp( -SimData$true_b * DataDF$S ) 
  DataDF$Fitadj <- SimData$true_a * DataDF$S * exp( -SimData$true_b * DataDF$S ) * exp( SimData$sigma^2 /2) 
  DataDF$Mod <- "True"
  DataDF$CI_low <- DataDF$CI_up  <-  DataDF$Pred <- DataDF$Pred_low <- DataDF$Pred_up <- DataDF$Fit
  
  # What is the distribution of simulated recruits
  resid <- data.frame(resid = (SimDataDF$R - DataDF$Fit), residadj = (SimDataDF$R - DataDF$Fitadj), 
                      year = c(1:length(SimDataDF$R) ) )
  #ggplot(resid, aes(year, resid) ) + geom_point() + geom_hline(yintercept=0)
  Rdist[i,1] <- mean(resid$resid)
  Rdist[i,2] <- median(resid$resid)
  
  Rdist[i,3] <- mean(resid$residadj)
  Rdist[i,4] <- median(resid$residadj)

    # # run using TMB, without prior
  # # Switch BiasCorr between T and F and plot to see impact of bias correction in LL
  TMB_nBC <- RunRicker(Data = SimDataDF,
                            Fitting_SW = Fitting_SW, #"tmbstan",#"TMB",
                            Priors = F, BiasCorr=F, Name = "TMB_nBC")
  TMB_wBC <- RunRicker(Data = SimDataDF,
                       Fitting_SW = Fitting_SW, #"tmbstan",#"TMB",
                       Priors = F, BiasCorr=T, Name = "TMB_wBC")
  
  if( Fitting_SW == "tmbstan") {
    Ests_quant_nBC <- apply ( as.data.frame(TMB_nBC$Ests) , 2, quantile, probs = c(0.025, 0.5, 0.975) )
    logA_nBC <- as.data.frame(Ests_quant_nBC)$logA[2] # median
    logSmax_nBC <- as.data.frame(Ests_quant_nBC)$logSmax[2] # median
    logSigma_nBC <- as.data.frame(Ests_quant_nBC)$logSigma[2] # median
    RicPars_nBC[i,] <- c(logA_nBC, logSmax_nBC, logSigma_nBC, TMB_nBC$Scale)
    
    Ests_quant_wBC <- apply ( as.data.frame(TMB_wBC$Ests) , 2, quantile, probs = c(0.025, 0.5, 0.975) )
    logA_wBC <- as.data.frame(Ests_quant_wBC)$logA[2] # median
    logSmax_wBC <- as.data.frame(Ests_quant_wBC)$logSmax[2] # median
    logSigma_wBC <- as.data.frame(Ests_quant_wBC)$logSigma[2] # median
    RicPars_wBC[i,] <- c(logA_wBC, logSmax_wBC, logSigma_wBC, TMB_wBC$Scale)
    
  }
  if( Fitting_SW == "TMB") {
    logA_nBC <- TMB_nBC$Ests %>% filter(Param=="logA") %>% pull(Estimate)
    logSmax_nBC <- TMB_nBC$Ests %>% filter(Param=="logSmax") %>% pull(Estimate)
    logSigma_nBC <- TMB_nBC$Ests %>% filter(Param=="logSigma") %>% pull(Estimate)
    RicPars_nBC[i,] <- c(logA_nBC, logSmax_nBC, logSigma_nBC, TMB_nBC$Scale)
    
    logA_wBC <- TMB_wBC$Ests %>% filter(Param=="logA") %>% pull(Estimate)
    logSmax_wBC <- TMB_wBC$Ests %>% filter(Param=="logSmax") %>% pull(Estimate)
    logSigma_wBC <- TMB_wBC$Ests %>% filter(Param=="logSigma") %>% pull(Estimate)
    RicPars_wBC[i,] <- c(logA_wBC, logSmax_wBC, logSigma_wBC, TMB_wBC$Scale)
    
  }
  
}


# Set up data frame of Ricker paramters for all MC trials

colnames(RicPars_nBC) <- c("logA", "logSmax", "logSigma", "Scale")
RicPars_nBC <- as.data.frame(RicPars_nBC)
RicPars_nBC <- RicPars_nBC %>% mutate(Smax = exp(logSmax)*Scale) 
colnames(RicPars_wBC) <- c("logA", "logSmax", "logSigma", "Scale")
RicPars_wBC <- as.data.frame(RicPars_wBC)
RicPars_wBC <- RicPars_wBC %>% mutate(Smax = exp(logSmax)*Scale) 


# Plot distribution of simulated recruits (in response to Murdoch's comments)
RdistDF <- data.frame( Metric = as.factor(c(rep("mean", ntrials), rep("median", ntrials))), 
                       Resids = c(Rdist[,1], Rdist[,2]) , 
                       ResidsAdj = c(Rdist[,3], Rdist[,4]))
ggplot(RdistDF, aes(x=Metric, y=Resids, fill = Metric)) + 
  geom_boxplot() + 
  annotate("text", x=1.45, y=mean(Rdist[,1]), label = round(mean(Rdist[,1],0)) ) + 
  annotate("text", x=2.46, y=mean(Rdist[,2]), label = round(mean(Rdist[,2],0)) ) +
  geom_hline(yintercept=0)

ggplot(RdistDF, aes(x=Metric, y=ResidsAdj, fill = Metric)) + 
  geom_boxplot() + 
  annotate("text", x=1.45, y=mean(Rdist[,3]), label = round(mean(Rdist[,3],0)) ) + 
  annotate("text", x=2.46, y=mean(Rdist[,4]), label = round(mean(Rdist[,4],0)) ) +
  geom_hline(yintercept=0)

cat("mean=", mean(Rdist[,1]) )
cat("median=", mean(Rdist[,2]) )


# Tune simulations by removing all MC trials with productivity less than 1 (logA < 0)

 # RicPars_nBC <- RicPars_nBC %>% filter(logA >= 0)
 # ntrialsTuned <- length (RicPars_nBC$logA)
 # RicPars_wBC <- RicPars_wBC %>% filter(logA >= 0)
 # ntrialsTuned <- length (RicPars_wBC$logA)


# Calculate predicted recruitments for all MC trials along a range of spawner abundances for plots
S <- NA
Rpred_nBC <- matrix(NA, 2000, ntrials)#Tuned)
Rpred_wBC <- matrix(NA, 2000, ntrials)#Tuned)

for (i in 1:ntrials) {#Tuned) {
  for (j in 1:2000){ #Assuming true smax ~ 100,000, and plot extends to 2 x Smax
    S[j] <- j*100
    Rpred_nBC[j,i] <-  exp( RicPars_nBC$logA[i] ) * S[j] * exp( -S[j] / RicPars_nBC$Smax[i] )
    Rpred_wBC[j,i] <-  exp( RicPars_wBC$logA[i] ) * S[j] * exp( -S[j] / RicPars_wBC$Smax[i] )
  }
}

# Calculate 5th, 50th, and 95th percentiles of the distribution of recruitment along
# the range of spawner abundances for plots

Rpred_nBC_dist <- apply(Rpred_nBC, 1, quantile, probs = c(0.025, 0.5, 0.975))
Rpred_nBC_mean <- apply(Rpred_nBC, 1, mean) 
Rpred_wBC_dist <- apply(Rpred_wBC, 1, quantile, probs = c(0.025, 0.5, 0.975))
Rpred_wBC_mean <- apply(Rpred_wBC, 1, mean) 


# Create a dataframe of these percentiles
EstAgg_nBC <- data.frame(S=S, Fit=Rpred_nBC_dist[2,], CI_low=Rpred_nBC_dist[1,], CI_up = Rpred_nBC_dist[3,])
EstAgg_nBC$Mod <- "TMB_nBC"
EstAgg_wBC <- data.frame(S=S, Fit=Rpred_wBC_dist[2,], CI_low=Rpred_wBC_dist[1,], CI_up = Rpred_wBC_dist[3,])
EstAgg_wBC$Mod <- "TMB_wBC"

# I get the same answer when I plot mean of the disribiton of recruitments instead of the median
# EstAgg <- data.frame(S=S, Fit=Rpred_mean, CI_low=Rpred_dist[1,], CI_up = Rpred_dist[3,])


if(ntrials ==1){
  All_Ests <- bind_rows(DataDF, TMB[[1]])
}

if(ntrials > 1){
  All_Ests <- bind_rows(DataDF, EstAgg_nBC)#EstAgg_wBC
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


# For a single iteration: now take estimated values and simulate using those, with and without correction

Ests <- TMB_nBC[[2]]

SimData_NoCorr <- Sim_Ricker_SR_Data(leng=30, age=4, Sig_Ricker = Ests$Estimate[Ests$Param=="sigma"], 
                               true_a = Ests$Estimate[Ests$Param=="A"], 
                               true_b = 1/(Ests$Estimate[Ests$Param=="Smax"]*TMB_nBC[[3]]),
                               hr_min = 0.25, hr_max = 0.35, lnorm_corr = F)

NoCorr_DF <- data.frame(S = round(SimData_NoCorr$S), R = (SimData_NoCorr$R), 
                        Year = 1:length(SimData_NoCorr$S), Mod = "NoCorr_True")

Ests <- TMB_wBC[[2]]

SimData_Corr <- Sim_Ricker_SR_Data(leng=30, age=4, Sig_Ricker = Ests$Estimate[Ests$Param=="sigma"], 
                             true_a = Ests$Estimate[Ests$Param=="A"], 
                             true_b = 1/(Ests$Estimate[Ests$Param=="Smax"]*TMB_wBC[[3]]),
                             hr_min = 0.25, hr_max = 0.35, lnorm_corr = T)

Corr_DF <-  data.frame(S = round(SimData_Corr$S), R = (SimData_Corr$R), 
                       Year = 1:length(SimData_Corr$S), Mod = "Corr_True")

Data_All <- bind_rows(DataDF, NoCorr_DF, Corr_DF)

ggplot(Data_All, aes(x=S, y=R, col = Mod)) + geom_point() + coord_fixed()

Data_All


#iterate generation of new data!
out <- list()
outSp <- list()

for (i in 1:ntrials){
  #set.seed(2422)
  SimData_nBC <- Sim_Ricker_SR_Data(leng=30, age=4, Sig_Ricker = exp(RicPars_nBC$logSigma[i]), 
                                    true_a = exp(RicPars_nBC$logA[i]), 
                                    true_b = 1/(RicPars_nBC$Smax[i]),
                                    hr_min = 0.25, hr_max = 0.35, lnorm_corr = F)
  
  
  #set.seed(2422)
  
  SimData_wBC <- Sim_Ricker_SR_Data(leng=30, age=4, Sig_Ricker = exp(RicPars_wBC$logSigma[i]), 
                                    true_a = exp(RicPars_wBC$logA[i]), 
                                    true_b = 1/(RicPars_wBC$Smax[i]),
                                    hr_min = 0.25, hr_max = 0.35, lnorm_corr = T)
  
  
  #set.seed(2422)
  SimData_wrongBC <- Sim_Ricker_SR_Data(leng=30, age=4, Sig_Ricker = exp(RicPars_nBC$logSigma[i]), 
                                    true_a = exp(RicPars_nBC$logA[i]), 
                                    true_b = 1/(RicPars_nBC$Smax[i]),
                                    hr_min = 0.25, hr_max = 0.35, lnorm_corr = T)
  
  SimData_HilBC <- Sim_Ricker_SR_Data(leng=30, age=4, Sig_Ricker = exp(RicPars_nBC$logSigma[i]), 
                                        true_a = exp(RicPars_nBC$logA[i] + exp(RicPars_nBC$logSigma[i])^2/2), 
                                        true_b = 1/(RicPars_nBC$Smax[i]),
                                        hr_min = 0.25, hr_max = 0.35, lnorm_corr = T)
  
  out[[i]] <- data.frame(i = i, SimData_nBC = SimData_nBC$R, SimData_wBC = SimData_wBC$R, 
                         SimData_wrongBC = SimData_wrongBC$R, SimData_HilBC= SimData_HilBC$R)
  outSp[[i]] <- data.frame(i = i, SimData_nBC = SimData_nBC$S, SimData_wBC = SimData_wBC$S, 
                         SimData_wrongBC = SimData_wrongBC$S, SimData_HilBC= SimData_HilBC$S)
  
}
out <- bind_rows(out)

out <- reshape2::melt(out, id.vars = "i") # or better, tidyr::pivot_longer

#Get spanwers from the simulated data
outSp <- bind_rows(outSp)
outSp <- reshape2::melt(outSp, id.vars = "i") # or better, tidyr::pivot_longer
# Get fitted R values from simulated spawners
get_Rhat <- function(a, b, S, sigma) (a*S*exp(-b*S))
# True Ricker values for Harrison: Sig_Ricker = 0.9, true_a = 1.93, true_b=1/159744,
Rhat <- get_Rhat(a=1.93, b=1/159744, outSp$value, 0.9)
outRhat <- outSp %>% mutate(Rhat = Rhat) %>% select(i, variable, Rhat)

resid <- data.frame(cbind(out, Rhat))  %>% mutate(resid=(value - Rhat)) %>% select(i, variable, resid)
#resid <- out %>% left_join(outRhat, by=c("i", "variable"))
resid <- group_by(resid, variable, i) %>%
  summarise(mean = mean(resid), median = median(resid), .groups = "drop")
resid <- reshape2::melt(resid, id.vars = c("i", "variable"), variable.name = "summary")
resid <- resid %>% filter(summary %in% c("mean", "median")) %>% mutate( residuals_medPars = round(as.numeric(value),0) )

result <- group_by(out, variable, i) %>%
  summarise(mean = mean(value), median = median(value), .groups = "drop")

result <- reshape2::melt(result, id.vars = c("i", "variable"), variable.name = "summary")

outTrueR <- bind_rows(outTrueR)
outTrueR <- reshape2::melt(outTrueR, id.vars = "i") # or better, tidyr::pivot_longer
resultTrueR <- group_by(outTrueR, variable, i) %>%
  summarise(mean = mean(value), median = median(value), .groups = "drop")
resultTrueR <- reshape2::melt(resultTrueR, id.vars = c("i", "variable"), variable.name = "summary")
resultTrueR <- resultTrueR %>%  filter(summary == c("mean", "median")) %>% 
  mutate( trueR= round(as.numeric(value),0) ) %>% select(i, summary, trueR)

# What is the mean and median of the true R over all MC trials. Plot this on plot of distribution of simulated Recruits
mean(outTrueR$value)
median(outTrueR$value)
lines <- tibble(summary = c("mean", "median"),
                value = c(mean(outTrueR$value), median(outTrueR$value)))

# These are the distribution of mean and mean simulated recruits over MC trials
result <- result %>% 
  filter(summary == c("mean", "median")) %>% mutate( value= round(as.numeric(value),0) ) 

# This is the distribution of relative errors over MC trials of the mean/median simulated recruits for each MC trials from its MC-specific true mean/median recruitment
# For each MC trial, I estimated the mean & median of the true recruits and mean & median of simulated recruits, and calcuated RE as the ratio of sim/true
rel_error <- result %>% left_join(resultTrueR, by=c("i", "summary") )  %>% mutate(RelativeError= value/trueR)


# Plot distribution of simulated Recruits
result %>%
  ggplot(aes(variable, value)) +
  #geom_point() +
  geom_jitter(alpha = 0.5, width=0.2) +
  ylim (0,300000) + 
  facet_wrap(vars(summary), ncol = 1) +
  geom_hline(data = lines, mapping = aes(yintercept = value))

# What is the relative error in simulated Rs from the overall mean and median of true Recruits?
group_by(result, summary, variable) %>%
  summarise(
    rel_error_mean = mean(as.numeric(value) / mean(outR$value)),
    rel_error_median = mean(as.numeric(value) / median(outR$value))) %>%
  knitr::kable(digits = 2)

# Plot residuals, but only for median of residiuals becuase we assumed the original parameters were median-unbiased
# Could switch to look at residuals from mean-unbiased parameters if we sitch the intial data generation to mean unbiased
# This plot is less clear than rel error plot
resid %>% filter(summary=="median") %>% #filter(variable=="SimData_nBC") %>%
  ggplot(aes(variable, residuals_medPars)) +
  geom_jitter(alpha = 0.5, width=0.2) +
  ylim (-2000,75000) + 
  facet_wrap(vars(summary), ncol = 1) +
  geom_hline(mapping = aes(yintercept = 0))


# Plot distribution or relative errors
REplot <- rel_error %>%
  ggplot(aes(variable, RelativeError)) +
  geom_jitter(alpha = 0.5, width=0.2) +
  ylim (0,5) + 
  facet_wrap(vars(summary), ncol = 1) +
  geom_hline(mapping = aes(yintercept = 1))

ggsave("DataOut/RelErrorRec.png", REplot)

group_by(rel_error, summary, variable) %>%
  summarise(
    RelativeError = median(RelativeError)) %>%
  knitr::kable(digits = 2)

# Now fit to both and compare (this code is for only 1 iteration)
# run using TMB, without prior, and estimating median curve (without bias correction in LL)
TMB_No_Corr <- RunRicker(Data = SimData_nBC,
                          Fitting_SW = "TMB",
                          Priors = F, BiasCorr = F, Name = "TMB_No_Corr")
TMB_Corr <- RunRicker(Data = SimData_wBC, 
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
 