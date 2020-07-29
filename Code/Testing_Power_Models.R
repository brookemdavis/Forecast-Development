# run power model for fennel, maybe Nadina? to get idea of possible parameters

library(dplyr)
library(ggplot2)
library(TMB)
library(tmbstan)
library(R2jags)
source("Code/Functions.R")


#==================================================================================
# simulate some power model data

SimData <- Sim_Power_SR_Data(leng=50, age=4, Sig_Mod = 0.3, true_a = 6, true_b=0.5,
                              hr_min = 0.2, hr_max = 0.8, lnorm_corr = F, max_pop = 100000)

  
SimDataDF <- data.frame(S = round(SimData$S), R = (SimData$R), Year = 1:length(SimData$S))
  
ggplot(SimDataDF, aes(x=S, y=R)) + geom_point() + coord_fixed()
  
# Create DF to store true and fitted values
DataDF <- SimDataDF[, c("S", "R", "Year")]
DataDF$Fit <- SimData$true_a * DataDF$S * exp( -SimData$true_b * DataDF$S )
DataDF$Mod <- "True"
DataDF$CI_low <- DataDF$CI_up  <-  DataDF$Pred <- DataDF$Pred_low <- DataDF$Pred_up <- DataDF$Fit  
  
Data <- SimDataDF
# Set up data and starting value lists list to go into model
data <- list() #data inputs
data$R_Obs <- Data$R 
data$S <- Data$S
data$alpha_mean <- 5
data$beta_mean <-  0
data$alpha_tau <- data$beta_tau <- 0.01
data$N <- dim(Data)[1]
data$Sig_Gam_Dist <- 0.001
# set up starting values

param <- list()
param$alpha <- param$beta <- 0
param$tau <- 1

init_vals <- list(param, param, param) # come back, this could be better

JagsFit <- jags(data, inits = init_vals, model.file = Power.model.MCMC, 
                n.chains =3, n.iter=10000, n.burnin = 4000, n.thin = 3, 
                parameters.to.save = c("alpha","beta","sigma","R_Fit", "R_Pred"))

#alpha = 6, beta = 0.5

# Turn into Data Frame
All_Ests <- data.frame(JagsFit$BUGSoutput$summary)
All_Ests$Param <- row.names(All_Ests)

R_Ests_Jags <- All_Ests[grepl("R_Fit", All_Ests$Param),  ]

FitsDF <- data.frame(S = data$S, R = data$R, Fit = R_Ests_Jags$X50. , 
                     Year = 1:dim(R_Ests_Jags)[1], 
                     CI_up = R_Ests_Jags$X97.5.,
                     CI_low = R_Ests_Jags$X2.5.)
                    
ggplot(data = FitsDF, aes(x=S, y=Fit, ymin = CI_low, ymax = CI_up)) +
  geom_line(size = 1.5) +
  geom_ribbon( alpha = 0.1) +
  geom_point(aes(x=S, y=R), col = "black") +
  theme_bw()


# now TMB
# only need to compile if changed model
dyn.unload(dynlib("Code/TMB/Single_Stock_Power"))
TMB::compile("Code/TMB/Single_Stock_Power.cpp") # can't seem to compile if have TMBstan loaded
dyn.load(dynlib("Code/TMB/Single_Stock_Power"))

Scale <- 10^(floor(log(mean(Data$R), 10)))

# change param values

TMBdata <- list() #data inputs
TMBdata$logR <- log(Data$R/Scale) 
TMBdata$S <- Data$S/Scale
TMBdata$A_mean <- 5
TMBdata$B_mean <-  0
TMBdata$A_sig <- TMBdata$B_sig <- sqrt(1/data$alpha_tau)
TMBdata$Sig_Gam_Dist <- 0.001
TMBdata$Scale <- Scale
TMBdata$Priors <- 0
TMBdata$Bayes <- 0

# set up starting values

param <- list()
param$logA <- param$logB <- 0
param$logSigma <- 1

obj <- MakeADFun(TMBdata, param, DLL="Single_Stock_Power", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))

# pull out estimates from ML fit
# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

# pull out fitted values
R_Ests <- All_Ests[grepl("R_Fit", All_Ests$Param),  ] [, -3 ]
names(R_Ests) <- c("R_Fit", "StdErr")

# create new rows with fitted values
FitsDF_TMB <- data.frame(S = Data$S, R = NA, Fit = R_Ests$R_Fit, Year = 1:dim(R_Ests)[1], 
                     Mod = "TMB_No_Prior")
FitsDF_TMB$CI_low <- R_Ests$R_Fit - 1.96*R_Ests$StdErr
FitsDF_TMB$CI_up <- R_Ests$R_Fit + 1.96*R_Ests$StdErr

# create prediction interval using simulate
R_Preds <- matrix(nrow = 1000, ncol = 50)
for(i in 1:1000){
  R_Preds[i, ] <- obj$simulate()$R_Pred
}
R_Pred_Summ <- apply(R_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))

FitsDF_TMB$Pred <- R_Pred_Summ[2,]
FitsDF_TMB$Pred_low <- R_Pred_Summ[1,]
FitsDF_TMB$Pred_up<- R_Pred_Summ[3,]


All_Ests <- bind_rows(FitsDF, FitsDF_TMB)
ggplot(data = All_Ests, aes(x=S, y=Fit, ymin = CI_low, ymax = CI_up, col = Mod, fill= Mod)) +
  geom_line(size = 1.5) +
  geom_ribbon( alpha = 0.1) +
  geom_point(aes(x=S, y=R), col = "black") +
  theme_bw()
