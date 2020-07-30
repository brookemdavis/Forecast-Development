#single population recruits per spawner model - made to be modular to be able to incorperate interrpopulation variability properly

library(dplyr)
library(TMB)
library(ggplot2)

# read in data and make sure it matches
SRdat_raw <- read.csv("DataIn/Chilliwack_updaterec.csv")

SRdat_Start2001 <- SRdat_raw %>% filter(Year >= 2001, is.na(rec) == F)

dyn.unload(dynlib("Code/TMB/Single_Stock_Ricker"))
TMB::compile("Code/TMB/Single_Stock_Ricker.cpp")
dyn.load(dynlib("Code/TMB/Single_Stock_Ricker"))
# want to scale down all data so close to 0
Scale <- 10000

Dat <- SRdat_Start2001 %>% filter(is.na(rec)==F)
data <- list()
data$S <- Dat$Escape/Scale
data$logR <- log(Dat$rec/Scale)
data$Priors <- 0
# need to have these in Data even though aren't used
data$logA_mean <- 0
data$logA_sig <- 0
data$Sig_Gam_Dist <- 0
data$Smax_mean <- 0
data$Smax_sig <- 0
data$Scale <- Scale

param <- list()
param$logA <- 1
param$Smax <- as.numeric(quantile(Dat$Escape, 0.8, na.rm = T)/Scale)
param$logSigma <- -2

# Now Fit TMB model
obj <- MakeADFun(data, param, DLL="Single_Stock_Ricker", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))

# pull out estimates from ML fit
# Create Table of outputs
All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

# pull out predictions
R_Ests <- All_Ests[grepl("R_Pred", All_Ests$Param),  ] [, -3 ]
names(R_Ests) <- c("R_Pred", "StdErr")

# add to simDataDF
FitsDF <- data.frame(S = SRdat_Start2001$Escape, R = SRdat_Start2001$rec, 
                     Fit = R_Ests$R_Pred, Year = 1:dim(R_Ests)[1])
FitsDF$CI_low <- FitsDF$Fit - 1.96*R_Ests$StdErr
FitsDF$CI_up <- FitsDF$Fit + 1.96*R_Ests$StdErr

FitsDF <- FitsDF %>% arrange(S)


ggplot(data = FitsDF, aes(x=S, y=Fit, ymin = CI_low, ymax = CI_up)) +
  geom_line() +
  geom_ribbon( alpha = 0.1) +
  geom_point(aes(x=S, y=R), col = "black") +
  theme_bw()

# extremely uncertain -- but stan model seems to have problem with negative values