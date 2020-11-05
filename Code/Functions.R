# Misc Functions for testing model fits


# Simulate a single stock of Ricker SR data =============================================================================================
Sim_Ricker_SR_Data <- function( leng=20, age=4, Sig_Ricker = 0.2, true_a = 3, true_b=1/5000,
                          hr_min = 0.2, hr_max = 0.8, lnorm_corr = F, autoCorr = F, rho=NA){
  
  # Estimate Sgen 
  SRep<- log(true_a) / true_b
  # Calculate spawners at MSY (approximation; Hilborn and Walters 1992)
  SMSY <- SRep * ( 0.5 - 0.07*log(true_a) )

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
  
  # opt <- optimize( f=ObjectiveSGen, interval=c(0, SMSY), SpMSY=SMSY,
  #                  alpha=true_a, beta=true_b )
  # 
  # # Get SGen from the optimized output (i.e., minimum neg log-like)
  # SGen <- opt$minimum
  
  # initiate population somwhere between 100 and Smax
  init <- round(runif(4, 100, 1/true_b))
  hr_vec <- runif(leng+age, hr_min, hr_max)
  
  esc<-rep(NA,leng+age)
  esc[1:age] <- init
  catch<-rep(NA,leng+age)
  rec<-rep(NA,leng+age)
  eps <- rep(NA,leng+age)
  
  for(i in (age+1):(leng+age)){
    # don't let esc or rec get lower than 100
    #rec[i] <- max (rlnorm(1, log(true_a) + log(esc[i-4]) - true_b*esc[i-4] , Sig_Ricker), 100)
    R_mean <- true_a*esc[i-4] * exp(-true_b*esc[i-4])
    # random recruitment residual
    if(autoCorr == F){
       eps[i] <- rnorm(1, 0, Sig_Ricker)
       Sig <- Sig_Ricker
    } else {
      if(i == age+1){
        # if first year just simulate random resid
        eps[i] <- rnorm(1, 0, Sig_Ricker)
        Sig <- Sig_Ricker
      } else {
        # in subsequent years start from previous e
        eps[i] <- rnorm(1, rho*eps[i-1], Sig_Ricker*sqrt(1-rho^2))
        Sig <- Sig_Ricker*sqrt(1-rho^2)
      }
    }
    if(lnorm_corr == F){ 
      rec[i] <- max(R_mean*exp(eps[i]), 100) 
    } else {
      rec[i] <- max(R_mean*exp(eps[i]-0.5*Sig^2), 100)
    }
    # want to find way to make sure pop gets knocked down every once in a while
    if(sum(esc[(i-4):(i-1)] > SMSY) == 4) { hr_vec[i] <- hr_max}
    esc[i] <- max((1-hr_vec[i])*rec[i], 100)
    catch[i] <- hr_vec[i]*rec[i]
  }

  
  catch_obs<-catch[(age+1):(leng+age)]
  esc_obs<-esc[(age+1):(leng+age)]
  
  #esc_obs<-pmax(1,rnorm(length(esc_obs), esc_obs, esc_obs*esc_cv))
  
  output<-list(true_a=true_a, true_b=true_b, sigma=Sig_Ricker, 
               R=rec[(age+1):(leng+age)],
               S=esc[1:leng], 
               catch=catch[(age+1):(leng+age)], 
               catch_obs=catch_obs, esc_obs=esc_obs)#, SGen = SGen)
  return(output)
  
}

#=============================================================================

Sim_Power_SR_Data <- function(leng=20, age=4, Sig_Mod = 0.2, true_a = 6, true_b=0.5,
                              hr_min = 0.2, hr_max = 0.8, lnorm_corr = F, max_pop = 100000){
  # initiate population somwhere between 100 and Smax
  init <- round(runif(4, 100, max_pop))
  hr_vec <- runif(leng+age, hr_min, hr_max)
  
  esc<-rep(NA,leng+age)
  esc[1:age] <- init
  catch<-rep(NA,leng+age)
  rec<-rep(NA,leng+age)
  
  for(i in (age+1):(leng+age)){
    # don't let esc or rec get lower than 100
    R_mean <- exp(true_a + true_b*log(esc[i-4]))
    # random recruitment residual
    e <- rnorm(1, 0, Sig_Mod)
    if(lnorm_corr == F){ 
      rec[i] <- max(R_mean*exp(e), 100) 
    } else {
      rec[i] <- max(R_mean*exp(e-0.5*Mod^2), 100)
    }
    # want to find way to make sure pop gets knocked down every once in a while
    if(sum(esc[(i-4):(i-1)] > max_pop) == 4) { hr_vec[i] <- hr_max}
    esc[i] <- max((1-hr_vec[i])*rec[i], 100)
    catch[i] <- hr_vec[i]*rec[i]
  }
  
  catch_obs<-catch[(age+1):(leng+age)]
  esc_obs<-esc[(age+1):(leng+age)]
  
  
  output<-list(true_a=true_a, true_b=true_b, sigma=Sig_Mod, 
               R=rec[(age+1):(leng+age)],
               S=esc[1:leng], 
               catch=catch[(age+1):(leng+age)], 
               catch_obs=catch_obs, esc_obs=esc_obs)
  return(output)
  
}



#============================================================================
# Simulate a single stock of Larkin SR data
Sim_Larkin_SR_Data <- function( leng=20, age=4, Sig_Larkin = 0.2, true_a = 3, true_b=c(1/10000, 0.00006, 0.00004, 0.00002),
                                hr_min = 0.2, hr_max = 0.8, lnorm_corr = F, autoCorr = F, rho=NA){
  
  
  # need extra years of lead-in for density-dependent affects
  LeadIn <- age+length(true_b)-1
  
  # initiate population somwhere between 100 and Smax_0
  init <- round(runif(LeadIn, 100, 1/true_b[1]))
  # Try to initiate with dominant cycle in 1st and 4th year
  init[c(1, 5)] <- init[c(1,5)]*2
  
  # Initiate other vectors
  esc<-rep(NA,leng + LeadIn)
  esc[1:(LeadIn)] <- init
  catch<-rep(NA, leng +LeadIn)
  rec<-rep(NA, leng +LeadIn)
  eps <- rep(NA, leng +LeadIn)
  hr_vec <- runif(leng+LeadIn, hr_min, hr_max)
  
  for(i in (LeadIn+1):(leng+LeadIn)){
    # don't let esc or rec get lower than 100
    R_mean <- true_a * esc[i-4] * exp(-true_b[1]*esc[i-4]-
                                     true_b[2]*esc[i-5] - 
                                     true_b[3]*esc[i-6] - 
                                     true_b[4]*esc[i-7])
    # random recruitment residual
    if(autoCorr == F){
      eps[i] <- rnorm(1, 0, Sig_Larkin)
      Sig <- Sig_Larkin
    } else {
      if(i == age+1){
        # if first year just simulate random resid
        eps[i] <- rnorm(1, 0, Sig_Larkin)
        Sig <- Sig_Larkin
      } else {
        # in subsequent years start from previous e
        eps[i] <- rnorm(1, rho*eps[i-1], Sig_Larkin*sqrt(1-rho^2))
        Sig <- Sig_Larkin*sqrt(1-rho^2)
      }
    }
    if(lnorm_corr == F){ 
      rec[i] <- max(R_mean*exp(eps[i]), 100) 
    } else {
      rec[i] <- max(R_mean*exp(eps[i]-0.5*Sig^2), 100)
    }
    # want to find way to make sure pop gets knocked down every once in a while
    # hard with cyclic -- look at 4 years ago and 8 years ago instead
    if(i>=8){
      if(sum(esc[c(i-4,i-8)] > 1/true_b[1]) == 2) { hr_vec[i] <- hr_max}
    }
    esc[i] <- max((1-hr_vec[i])*rec[i], 100)
    catch[i] <- hr_vec[i]*rec[i]
  }
  
  
  catch_obs <- catch[(LeadIn+1):(leng+LeadIn)]
  esc_obs <- esc[(age):(leng+age-1)]
  
  
  output<-list(true_a=true_a, true_b=true_b, sigma=Sig_Larkin, 
               R=rec[(LeadIn+1):(leng+LeadIn)],
               S=esc_obs, 
               catch=catch[(LeadIn+1):(leng+LeadIn)], 
               catch_obs=catch_obs, esc_obs=esc_obs)
  
  
  return(output)
  
}
#===========================================================================
# Generic function to Run Ricker model using TMB, stan, TMBstan, jags
# returns estimates and prediction intervals

RunRicker <- function(Data, 
                      Fitting_SW = "TMB", # which model fitting software to use? TMB, Stan, tmbstan, JAGS
                      Priors = T, # if fitting with TMB, have option not use priors
                      BiasCorr = F, #Should a bias correction be included in the LL?
                      Name = "Test", # Name to put in Mod column
                      logA_mean = 0,logA_sig = 0, # priors on logAlpha
                      Sig_Gam_Dist = 0.001, # inverse gamma shape and scale param
                      Smax_mean = 0, Smax_sig = 0) { # priors on capacity, Smax
 
  # Want to scale down obs to make models more stable
  # only really required for TMB, but will use for all models
    Scale <- 10^(floor(log(mean(Data$R), 10)))

  # Set up data and starting value lists list to go into model
  # these are the same regardless of software
  data <- list() #data inputs
  data$S <- Data$S/Scale 
  data$logA_mean <- logA_mean
  data$Sig_Gam_Dist <- Sig_Gam_Dist
  data$logSmax_mean <- log(Smax_mean/Scale)
  # set Bayes to 0, switch to one if using tmbstan
  data$Bayes <- 0
  # Add a bias correction in Ricker estimation
  data$BiasCorr <- as.numeric(BiasCorr)
  # set up starting values
  param <- list()
  param$logA <- 1
  
  
  # if using TMB input logR, rather than R_Obs, need Scale, Prior indicator
  if(Fitting_SW %in% c("TMB", "tmbstan")){
    data$logR <- log(Data$R/Scale)
    data$Scale <- Scale
    data$Priors <- as.numeric(Priors)
    param$logSmax <- log(as.numeric(quantile(Data$S, 0.8)/Scale))
    if(Fitting_SW =="tmbstan") data$Bayes <- 1
  }
  
  # if using TMB, stan variance is defined using sd
  if(Fitting_SW %in% c("TMB", "Stan", "tmbstan")){
    data$logA_sig <- logA_sig
    data$logSmax_sig <- Smax_sig/Scale
    param$logSigma <- -2
  }
  
  # if jags or stan need R_Obs rather than logR, and N
  if(Fitting_SW %in% c("JAGS", "Stan")){
    data$R_Obs <- Data$R/Scale
    data$N <- dim(Data)[1]
    param$Smax <- as.numeric(quantile(Data$S, 0.8)/Scale)
  }
  
  # if using jags, variance is defined using precision (1/variance)
  if(Fitting_SW == "JAGS"){
    data$logA_tau <- 1/(logA_sig)^2
    data$logSmax_tau <- 1/(Smax_sig)^2
    param$tau <- 1/exp(-2)^2
  }

  #==========
  # TMB Fit
  #==========
   if( Fitting_SW == "TMB") {
    # Now Fit TMB model
    obj <- MakeADFun(data, param, DLL="Single_Stock_Ricker", silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
    
    # pull out estimates from ML fit
    # Create Table of outputs
    All_Ests <- data.frame(summary(sdreport(obj)))
    All_Ests$Param <- row.names(All_Ests)
    
    # pull out fitted values
    R_Ests <- All_Ests[grepl("R_Fit", All_Ests$Param),  ] [, -3 ]
    names(R_Ests) <- c("R_Fit", "StdErr")
    
    # create new rows with fitted values
    FitsDF <- data.frame(S = Data$S, R = NA, Fit = R_Ests$R_Fit, Year = 1:dim(R_Ests)[1], 
                           Mod = Name)
    FitsDF$CI_low <- R_Ests$R_Fit - 1.96*R_Ests$StdErr
    FitsDF$CI_up <- R_Ests$R_Fit + 1.96*R_Ests$StdErr
    
    # create prediction interval using simulate
    R_Preds <- matrix(nrow = 1000, ncol = 50)
    for(i in 1:1000){
      R_Preds[i, ] <- obj$simulate()$R_Pred
    }
    R_Pred_Summ <- apply(R_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    FitsDF$Pred <- R_Pred_Summ[2,]
    FitsDF$Pred_low <- R_Pred_Summ[1,]
    FitsDF$Pred_up<- R_Pred_Summ[3,]
    
   } # end TMB fit
  
  #====================
  # Fit using tmbstan
  #====================
  if(Fitting_SW == "tmbstan"){
    # set up and fit TMB model
    obj <- MakeADFun(data, param, DLL="Single_Stock_Ricker", silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
    # now fit as mcmc
    fitmcmc <- tmbstan(obj, chains=3, iter=10000, init=list(opt$par), 
                       control = list(adapt_delta = 0.95))
    # pull out posterior vals
    All_Ests <- as.matrix(fitmcmc)
    
    # pull out R_Fit values
    R_Fit_Med <- obj$report(All_Ests[1,-ncol(All_Ests)])$R_Fit         
    R_Fit <- matrix(NA, nrow=nrow(All_Ests), ncol = length(R_Fit_Med))
    for(i in 1:nrow(All_Ests)){
      r <- obj$report(All_Ests[i,-ncol(All_Ests)])
      R_Fit[i,] <- r$R_Fit
    }
    
    # now for each column get median, quantiles and add to DF
    R_Fit_Summ <- apply(R_Fit, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    # Do same simulate() routine to get prediction intervals
    R_Preds <- matrix(nrow = 1000, ncol = 50)
    for(i in 1:1000){
      R_Preds[i, ] <- obj$simulate()$R_Pred
    }
    R_Pred_Summ <- apply(R_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    
    FitsDF <- data.frame(S = Data$S, R = Data$R, Fit = R_Fit_Summ[2,], 
                               Year = 1:dim(R_Pred_Summ)[2],   Mod = Name,
                               CI_up = R_Fit_Summ[1,],
                               CI_low = R_Fit_Summ[3,],
                               Pred = R_Pred_Summ[2,],
                               Pred_low = R_Pred_Summ[1,],
                               Pred_up = R_Pred_Summ[3,])
    
    # prepare posteriors for return
    A_Post <- exp(All_Ests[, c("logA")])
    Smax_Post <- exp(All_Ests[, c("logSmax")])*Scale
  } #end tmbstan fit

  #===========
  # JAGS Fit
  #===========
  if(Fitting_SW == "JAGS"){

    init_vals <- list(param, param, param) # come back, this could be better
    
    JagsFit <- jags(data, inits = init_vals, model.file = Ricker.model.MCMC, 
                    n.chains =3, n.iter=10000, n.burnin = 4000, n.thin = 3, 
                    parameters.to.save = c("R_Fit", "R_Pred", "logA", "Smax", "sigma"))
    
    # Turn into Data Frame
    All_Ests <- data.frame(JagsFit$BUGSoutput$summary)
    All_Ests$Param <- row.names(All_Ests)
    
    R_Ests_Jags <- All_Ests[grepl("R_Fit", All_Ests$Param),  ]
    R_Preds_Jags <- All_Ests[grepl("R_Pred", All_Ests$Param),  ]
    
    FitsDF <- data.frame(S = Data$S, R = Data$R, Fit = R_Ests_Jags$X50. * Scale, 
                              Year = 1:dim(R_Ests_Jags)[1],   Mod = Name,
                              CI_up = R_Ests_Jags$X97.5. * Scale,
                              CI_low = R_Ests_Jags$X2.5. * Scale,
                              Pred = R_Preds_Jags$X50. * Scale,
                              Pred_low = R_Preds_Jags$X2.5. * Scale,
                              Pred_up = R_Preds_Jags$X97.5. * Scale)
    
    # Prep A and Smax posteriors for outputs
    A_Post <- exp(JagsFit$BUGSoutput$sims.list$logA)
    Smax_Post <- JagsFit$BUGSoutput$sims.list$Smax * Scale
  } # end JAGS fit
  
  #===========
  # Stan Fit
  #===========
  
  if(Fitting_SW == "Stan"){
  
  #run stan model
  stan_fit <- stan(file = 'Code/STAN/Single_Stock_Ricker.stan', data = data, iter = 10000, 
                   chains = 3,  control = list(adapt_delta = 0.95))
  
  All_Ests <- data.frame(summary(stan_fit)$summary)
  All_Ests$Param <- row.names(All_Ests)
  
  R_Fits_Stan <- All_Ests[grepl("R_Fit", All_Ests$Param),  ]
  R_Preds_Stan <- All_Ests[grepl("R_Pred", All_Ests$Param),  ]
  
  FitsDF <- data.frame(S = Data$S, R = Data$R, Fit = R_Fits_Stan$X50. * Scale, 
                            Year = 1:dim(R_Fits_Stan)[1],   Mod = Name,
                            CI_up = R_Fits_Stan$X97.5. * Scale,
                            CI_low = R_Fits_Stan$X2.5. * Scale,
                            Pred = R_Preds_Stan$X50. * Scale,
                            Pred_up = R_Preds_Stan$X97.5. * Scale,
                            Pred_low = R_Preds_Stan$X2.5. * Scale)
  
  # get A and Smax posteriors
  fit_values <- extract(stan_fit)
  A_Post <- exp(fit_values$logA)
  Smax_Post <- fit_values$Smax * Scale
  
  } # end stan fit
  
# Return fit and predicted values  
  out <- list()
  out$Fits <- FitsDF
  out$Ests <- All_Ests
  out$Scale <- Scale
  # If Bayesian model return variable posteriors
  if(Fitting_SW %in% c("tmbstan", "JAGS", "Stan")){
    out$A_Post <- A_Post
    out$Smax_Post <- Smax_Post
  }
  out
  
}
#========================================================================================

# Now same function but for power model
RunPower <- function(Data, 
                      Fitting_SW = "TMB", # which model fitting software to use? TMB, Stan, tmbstan, JAGS
                      Priors = T, # if fitting with TMB, have option not use priors
                      Name = "Test", # Name to put in Mod column
                      A_mean = 0, A_sig = 0, # priors on Alpha
                      B_mean = 0, B_sig = 0, # priors on beta
                      Sig_Gam_Dist = 0.001 ) {  # inverse gamma shape and scale param
  
  # Want to scale down obs to make models more stable
  # only really required for TMB, but will use for all models
  Scale <- 10^(floor(log(mean(Data$R), 10)))
  
  # Set up data and starting value lists list to go into model
  # these are the same regardless of software
  data <- list() #data inputs
  data$S <- Data$S/Scale 
  data$A_mean <- A_mean
  data$B_mean <- B_mean
  data$Sig_Gam_Dist <- Sig_Gam_Dist
  # set Bayes to 0, switch to one if using tmbstan
  data$Bayes <- 0
  
  # set up starting values
  param <- list()
  param$logA <- param$logB <- 0
  param$logSigma <- 1
  
  # if using TMB input logR, rather than R_Obs, need Scale, Prior indicator
  if(Fitting_SW %in% c("TMB", "tmbstan")){
    data$logR <- log(Data$R/Scale)
    data$Scale <- Scale
    data$Priors <- as.numeric(Priors)
    if(Fitting_SW =="tmbstan"){ data$Bayes <- 1}
  }
  
  # if using TMB, stan variance is defined using sd
  if(Fitting_SW %in% c("TMB", "Stan", "tmbstan")){
    data$A_sig <- A_sig
    data$B_sig <- B_sig
  }
  
  # if jags or stan need R_Obs rather than logR, and N
  if(Fitting_SW %in% c("JAGS", "Stan")){
    data$R_Obs <- Data$R/Scale
    data$N <- dim(Data)[1]
  }
  
  # if using jags, variance is defined using precision (1/variance)
  if(Fitting_SW == "JAGS"){
    data$A_tau <- 1/(A_sig)^2
    data$B_tau <- 1/(B_sig)^2
    param$tau <- 1/exp(-2)^2
  }
  
  #==========
  # TMB Fit
  #==========
  if( Fitting_SW == "TMB") {
    # Now Fit TMB model
    obj <- MakeADFun(data, param, DLL="Single_Stock_Power", silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
    
    # pull out estimates from ML fit
    # Create Table of outputs
    All_Ests <- data.frame(summary(sdreport(obj)))
    All_Ests$Param <- row.names(All_Ests)
    
    # pull out fitted values
    R_Ests <- All_Ests[grepl("R_Fit", All_Ests$Param),  ] [, -3 ]
    names(R_Ests) <- c("R_Fit", "StdErr")
    
    # create new rows with fitted values
    FitsDF <- data.frame(S = Data$S, R = NA, Fit = R_Ests$R_Fit, Year = 1:dim(R_Ests)[1], 
                         Mod = Name)
    FitsDF$CI_low <- R_Ests$R_Fit - 1.96*R_Ests$StdErr
    FitsDF$CI_up <- R_Ests$R_Fit + 1.96*R_Ests$StdErr
    
    # create prediction interval using simulate
    R_Preds <- matrix(nrow = 1000, ncol = 50)
    for(i in 1:1000){
      R_Preds[i, ] <- obj$simulate()$R_Pred
    }
    R_Pred_Summ <- apply(R_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    FitsDF$Pred <- R_Pred_Summ[2,]
    FitsDF$Pred_low <- R_Pred_Summ[1,]
    FitsDF$Pred_up<- R_Pred_Summ[3,]
    
  } # end TMB fit
  
  #====================
  # Fit using tmbstan
  #====================
  if(Fitting_SW == "tmbstan"){
    # set up and fit TMB model
    obj <- MakeADFun(data, param, DLL="Single_Stock_Power", silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
    # now fit as mcmc
    fitmcmc <- tmbstan(obj, chains=3, iter=10000, init=list(opt$par), 
                       control = list(adapt_delta = 0.95))
    # pull out posterior vals
    All_Ests<-as.matrix(fitmcmc)
    
    # pull out R_Fit values
    R_Fit_Med <- obj$report(All_Ests[1,-ncol(All_Ests)])$R_Fit         
    R_Fit <- matrix(NA, nrow=nrow(All_Ests), ncol = length(R_Fit_Med))
    for(i in 1:nrow(All_Ests)){
      r <- obj$report(All_Ests[i,-ncol(All_Ests)])
      R_Fit[i,] <- r$R_Fit
    }
    
    # now for each column get median, quantiles and add to DF
    R_Fit_Summ <- apply(R_Fit, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    # Do same simulate() routine to get prediction intervals
    R_Preds <- matrix(nrow = 1000, ncol = 50)
    for(i in 1:1000){
      R_Preds[i, ] <- obj$simulate()$R_Pred
    }
    R_Pred_Summ <- apply(R_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    FitsDF <- data.frame(S = Data$S, R = Data$R, Fit = R_Fit_Summ[2,], 
                         Year = 1:dim(R_Pred_Summ)[2],   Mod = Name,
                         CI_up = R_Fit_Summ[1,],
                         CI_low = R_Fit_Summ[3,],
                         Pred = R_Pred_Summ[2,],
                         Pred_low = R_Pred_Summ[1,],
                         Pred_up = R_Pred_Summ[3,])
    
    # prepare posteriors for return
    A_Post <- exp(All_Ests[, c("logA")])
    B_Post <- exp(All_Ests[, c("logB")])
    
  } #end tmbstan fit
  
  #===========
  # JAGS Fit
  #===========
  if(Fitting_SW == "JAGS"){
    
    init_vals <- list(param, param, param) # come back, this could be better
    
    JagsFit <- jags(data, inits = init_vals, model.file = Power.model.MCMC, 
                    n.chains =3, n.iter=10000, n.burnin = 4000, n.thin = 3, 
                    parameters.to.save = c("R_Fit", "R_Pred", "A", "B", "sigma"))
    
    # Turn into Data Frame
    All_Ests <- data.frame(JagsFit$BUGSoutput$summary)
    All_Ests$Param <- row.names(All_Ests)
    
    R_Ests_Jags <- All_Ests[grepl("R_Fit", All_Ests$Param),  ]
    R_Preds_Jags <- All_Ests[grepl("R_Pred", All_Ests$Param),  ]
    
    FitsDF <- data.frame(S = Data$S, R = Data$R, Fit = R_Ests_Jags$X50. * Scale, 
                         Year = 1:dim(R_Ests_Jags)[1],   Mod = Name,
                         CI_up = R_Ests_Jags$X97.5. * Scale,
                         CI_low = R_Ests_Jags$X2.5. * Scale,
                         Pred = R_Preds_Jags$X50. * Scale,
                         Pred_low = R_Preds_Jags$X2.5. * Scale,
                         Pred_up = R_Preds_Jags$X97.5. * Scale)
    
    # Prep A and Smax posteriors for outputs
    A_Post <- JagsFit$BUGSoutput$sims.list$A
    B_Post <- JagsFit$BUGSoutput$sims.list$B
  } # end JAGS fit
  
  #===========
  # Stan Fit
  #===========
  
  if(Fitting_SW == "Stan"){
    
    #run stan model
    stan_fit <- stan(file = 'Code/STAN/Single_Stock_Power.stan', data = data, iter = 10000, 
                     chains = 3,  control = list(adapt_delta = 0.95))
    
    All_Ests <- data.frame(summary(stan_fit)$summary)
    All_Ests$Param <- row.names(All_Ests)
    
    R_Fits_Stan <- All_Ests[grepl("R_Fit", All_Ests$Param),  ]
    R_Preds_Stan <- All_Ests[grepl("R_Pred", All_Ests$Param),  ]
    
    FitsDF <- data.frame(S = Data$S, R = Data$R, Fit = R_Fits_Stan$X50. * Scale, 
                         Year = 1:dim(R_Fits_Stan)[1],   Mod = Name,
                         CI_up = R_Fits_Stan$X97.5. * Scale,
                         CI_low = R_Fits_Stan$X2.5. * Scale,
                         Pred = R_Preds_Stan$X50. * Scale,
                         Pred_up = R_Preds_Stan$X97.5. * Scale,
                         Pred_low = R_Preds_Stan$X2.5. * Scale)
    
    # get A and Smax posteriors
    fit_values <- extract(stan_fit)
    A_Post <- fit_values$A
    B_Post <- fit_values$B

  } # end stan fit
  
  # Return fit and predicted values  
  out <- list()
  out$Fit <- FitsDF
  out$Ests <- All_Ests
  out$Scale <- Scale
  # If Bayesian model return variable posteriors
  if(Fitting_SW %in% c("tmbstan", "JAGS", "Stan")){
    out$A_Post <- A_Post
    out$B_Post <- B_Post
  }
  out
  
}

#========================================================================================
# Same function but for Larkin model
#======================================================================================

RunLarkin <-  function(Data, 
                       Fitting_SW = "TMB", # which model fitting software to use? TMB, Stan, tmbstan, JAGS
                       Priors = T, # if fitting with TMB, have option not use priors
                       Name = "Test", # Name to put in Mod column
                       logA_mean = 0,logA_sig = 1, # priors on logAlpha
                       Sig_Gam_Dist = 0.001, # inverse gamma shape and scale param
                       B_means = rep(0,4), B_sigs = rep(0.001, 4)) { # priors on Bs
  
 
  
  # Want to scale down obs to make models more stable
  # only really required for TMB, but will use for all models
  Scale <- 10^(floor(log(mean(Data$R), 10)))
  
  # Set up data and starting value lists list to go into model
  # these are the same regardless of software
  data <- list() #data inputs
  data$S <- Data$S/Scale 
  data$logA_mean <- logA_mean
  data$Sig_Gam_Dist <- Sig_Gam_Dist
  data$B_means <- B_means*Scale
  # set Bayes to 0, switch to one if using tmbstan
  data$Bayes <- 0
  
  # set up starting values
  param <- list()
  param$logA <- 1
  
  
  # if using TMB input logR, rather than R_Obs, need Scale, Prior indicator
  if(Fitting_SW %in% c("TMB", "tmbstan")){
    data$logR <- log(Data$R/Scale)
    data$Scale <- Scale
    data$Priors <- as.numeric(Priors)
    B0 = 1/as.numeric(quantile(Data$S, 0.8)/Scale)
    param$logB <- log(c(B0, B0/2, B0/3, B0/4)) 
    if(Fitting_SW =="tmbstan") data$Bayes <- 1
  }
  
  # if using TMB/stan variance is defined using sd
  if(Fitting_SW %in% c("TMB", "Stan", "tmbstan")){
    data$logA_sig <- logA_sig
    data$B_sigs <- B_sigs*Scale
    param$logSigma <- -2
  }
  
  # if jags or stan need R_Obs rather than logR, and N
  if(Fitting_SW %in% c("JAGS", "Stan")){
    data$R_Obs <- Data$R/Scale
    data$N <- dim(Data)[1]
    B0 = 1/as.numeric(quantile(Data$S, 0.8)/Scale)
    param$beta0 <- B0
    param$beta1 <- B0/2
    param$beta2 <- B0/3
    param$beta3 <- B0/4
    data$Bayes <- NULL
  }
  
  # if using jags, variance is defined using precision (1/variance)
  if(Fitting_SW == "JAGS"){
    data$logA_tau <- 1/(logA_sig)^2
    data$B_taus <- 1/(B_sigs)^2
    param$tau <- 1/exp(-2)^2
  }
  
  #==========
  # TMB Fit
  #==========
  if( Fitting_SW == "TMB") {
    
    # Now Fit TMB model
    obj <- MakeADFun(data, param, DLL="Single_Stock_Larkin", silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
    
    # pull out estimates from ML fit
    # Create Table of outputs
    All_Ests <- data.frame(summary(sdreport(obj)))
    All_Ests$Param <- row.names(All_Ests)
    
    # pull out fitted values
    R_Ests_Raw <- All_Ests[grepl("R_Fit", All_Ests$Param),  ] [, -3 ]
    
    # add 3 NAs to R_Ests at beginning so aligns with data
    R_Ests <- data.frame(R_Fit = c(rep(NA, (length(B_means)-1)), R_Ests_Raw$Estimate),
                         StdErr = c(rep(NA, (length(B_means)-1)), R_Ests_Raw$Std..Error) )
    
    
    # create new rows with fitted values
    FitsDF <- data.frame(S = Data$S, R = NA, Fit = R_Ests$R_Fit, Year = 1:dim(R_Ests)[1], 
                         Mod = Name)
    FitsDF$CI_low <- R_Ests$R_Fit - 1.96*R_Ests$StdErr
    FitsDF$CI_up <- R_Ests$R_Fit + 1.96*R_Ests$StdErr
    
    # create prediction interval using simulate
    R_Preds <- matrix(nrow = 1000, ncol = (dim(Data)[1]-length(B_means)+1) )
    for(i in 1:1000){
      R_Preds[i, ] <- obj$simulate()$R_Pred
    }
    R_Pred_Summ <- apply(R_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    # also need to add NAs to beginning of these vectors
    FitsDF$Pred <- c(rep(NA, (length(B_means)-1)), R_Pred_Summ[2,])
    FitsDF$Pred_low <- c(rep(NA, (length(B_means)-1)),R_Pred_Summ[1,] )
    FitsDF$Pred_up<- c(rep(NA, (length(B_means)-1)),R_Pred_Summ[3,])
    
  } # end TMB fit
  
  #====================
  # Fit using tmbstan
  #====================
  if(Fitting_SW == "tmbstan"){
    # set up and fit TMB model
    obj <- MakeADFun(data, param, DLL="Single_Stock_Larkin", silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
    # now fit as mcmc
    fitmcmc <- tmbstan(obj, chains=3, iter=10000, init=list(opt$par), 
                       control = list(adapt_delta = 0.95))
    # pull out posterior vals
    All_Ests <- as.matrix(fitmcmc)
    
  
    # pull out R_Fit values
    R_Fit_Med <- obj$report(All_Ests[1,-ncol(All_Ests)])$R_Fit   
    # why is this 67 long instead of 47 long?
    R_Fit <- matrix(NA, nrow=nrow(All_Ests), ncol = length(R_Fit_Med))
    for(i in 1:nrow(All_Ests)){
      r <- obj$report(All_Ests[i,-ncol(All_Ests)])
      R_Fit[i,] <- r$R_Fit
    }
    
    
    # now for each column get median, quantiles and add to DF
    R_Fit_Summ <- apply(R_Fit, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    # Do same simulate() routine to get prediction intervals
    R_Preds <- matrix(nrow = 1000, ncol = (dim(Data)[1]-length(B_means)+1))
    for(i in 1:1000){
      R_Preds[i, ] <- obj$simulate()$R_Pred
    }
    R_Pred_Summ <- apply(R_Preds, 2, quantile, probs = c(0.025, 0.5, 0.975))
    
    FitsDF <- data.frame(S = Data$S, R = Data$R, 
                         Fit = c(rep(NA, length(B_means)-1), R_Fit_Summ[2,]), 
                         Year = 1:length(Data$S),   Mod = Name,
                         CI_up = c(rep(NA, length(B_means)-1), R_Fit_Summ[1,]),
                         CI_low = c(rep(NA, length(B_means)-1), R_Fit_Summ[3,]),
                         Pred = c(rep(NA, length(B_means)-1), R_Pred_Summ[2,]),
                         Pred_low = c(rep(NA, length(B_means)-1), R_Pred_Summ[1,]),
                         Pred_up = c(rep(NA, length(B_means)-1), R_Pred_Summ[3,]) )
    
    # prepare posteriors for return
    A_Post <- exp(All_Ests[, c("logA")])

  } #end tmbstan fit
  
  #===========
  # JAGS Fit
  #===========
  if(Fitting_SW == "JAGS"){
    
    
    init_vals <- list(param, param, param) # come back, this could be better
    
    JagsFit <- jags(data, inits = init_vals, model.file = Larkin.model.MCMC, 
                    n.chains =3, n.iter=10000, n.burnin = 4000, n.thin = 3, 
                    parameters.to.save = c("R_Fit", "R_Pred", "logA", 
                                    "beta0", "beta1", "beta2", "beta3", "sigma"))
    
    # Turn into Data Frame
    All_Ests <- data.frame(JagsFit$BUGSoutput$summary)
    All_Ests$Param <- row.names(All_Ests)
    
    R_Ests_Jags <- All_Ests[grepl("R_Fit", All_Ests$Param),  ]
    R_Preds_Jags <- All_Ests[grepl("R_Pred", All_Ests$Param),  ]
    
    FitsDF <- data.frame(S = Data$S, R = Data$R, 
                         Fit = c(rep(NA, length(B_means)-1), R_Ests_Jags$X50. * Scale), 
                         Year = 1:dim(Data)[1],   Mod = Name,
                         CI_up = c(rep(NA, length(B_means)-1), R_Ests_Jags$X97.5. * Scale),
                         CI_low = c(rep(NA, length(B_means)-1), R_Ests_Jags$X2.5. * Scale),
                         Pred = c(rep(NA, length(B_means)-1), R_Preds_Jags$X50. * Scale),
                         Pred_low = c(rep(NA, length(B_means)-1), R_Preds_Jags$X2.5. * Scale),
                         Pred_up = c(rep(NA, length(B_means)-1), R_Preds_Jags$X97.5. * Scale))
    
    # Prep A and Smax posteriors for outputs
    A_Post <- exp(JagsFit$BUGSoutput$sims.list$logA)
    #Smax_Post <- JagsFit$BUGSoutput$sims.list$Smax * Scale
  } # end JAGS fit
  
  #===========
  # Stan Fit
  #===========
  
  if(Fitting_SW == "Stan"){
    
    #run stan model
    stan_fit <- stan(file = 'Code/STAN/Single_Stock_Larkin.stan', data = data, iter = 10000, 
                     chains = 3,  control = list(adapt_delta = 0.95))
    
    All_Ests <- data.frame(summary(stan_fit)$summary)
    All_Ests$Param <- row.names(All_Ests)
    
    R_Fits_Stan <- All_Ests[grepl("R_Fit", All_Ests$Param),  ]
    R_Preds_Stan <- All_Ests[grepl("R_Pred", All_Ests$Param),  ]
    
    FitsDF <- data.frame(S = Data$S, R = Data$R, 
                         Fit =c(rep(NA, length(B_means)-1),  R_Fits_Stan$X50. * Scale), 
                         Year = 1:dim(Data)[1],   Mod = Name,
                         CI_up = c(rep(NA, length(B_means)-1), R_Fits_Stan$X97.5. * Scale),
                         CI_low = c(rep(NA, length(B_means)-1), R_Fits_Stan$X2.5. * Scale),
                         Pred = c(rep(NA, length(B_means)-1), R_Preds_Stan$X50. * Scale),
                         Pred_up = c(rep(NA, length(B_means)-1), R_Preds_Stan$X97.5. * Scale),
                         Pred_low = c(rep(NA, length(B_means)-1), R_Preds_Stan$X2.5. * Scale))
    
    # get A and Smax posteriors
    fit_values <- extract(stan_fit)
    A_Post <- exp(fit_values$logA)
    #Smax_Post <- fit_values$Smax * Scale
    
  } # end stan fit
  
  # Return fit and predicted values  
  out <- list()
  out$Fits <- FitsDF
  out$Ests <- All_Ests
  out$Scale <- Scale
  # If Bayesian model return variable posteriors
  if(Fitting_SW %in% c("tmbstan", "JAGS", "Stan")){
    out$A_Post <- A_Post
  }
  out
  
}

#========================================================================================

# Jags Models

Ricker.model.MCMC <- function(){
  for (i in 1:N) {                       #loop over N sample points
    R_Obs[i] ~ dlnorm(logR_Fit[i], tau)          #likelihood -> predicted value for NA in data set
    logR_Fit[i] <-  logA - beta * S[i] + log(S[i])               # calc log(R) - fitted values  
    R_Fit[i] <- exp(logR_Fit[i])
    R_Pred[i] ~ dlnorm(logR_Fit[i],tau)
  }
  
  logA ~ dnorm(logA_mean, logA_tau)       
  beta <-1/Smax					   # prior for beta
  Smax ~ dlnorm(logSmax_mean, logSmax_tau)       			   # prior for beta 
  tau ~ dgamma(Sig_Gam_Dist,Sig_Gam_Dist)                    #prior for precision parameter
  sigma <- 1/sqrt(tau) 		   	
  
}

# Power model from forecast model
Power.model.MCMC <- function(){
  for (i in 1:N) {                             # loop over N sample points
    R_Obs[i] ~ dlnorm(logR_Fit[i], tau)          # likelihood -> predicted value for NA in data set
    logR_Fit[i] <- A + B * log(S[i])       # power model
    R_Fit[i] <- exp(logR_Fit[i])
    R_Pred[i] ~ dlnorm(logR_Fit[i],tau)
  }
  
  A ~ dnorm(A_mean, A_tau)             # prior for alpha
  B ~ dnorm(B_mean, B_tau)                # prior for beta
  tau ~ dgamma(Sig_Gam_Dist,Sig_Gam_Dist)    # prior for precision parameter
  sigma <- 1/sqrt(tau)   		                  	
  
}


# Larkin Model
Larkin.model.MCMC <- function(){
  for(i in 4:N) {               		#loop over N sample points
    logR_Fit[i-3] <- logA + log(S[i])-beta0*S[i]-beta1*S[i-1]-beta2*S[i-2]-beta3*S[i-3]
    R_Obs[i] ~ dlnorm(logR_Fit[i-3], tau)	 # likelihood
    R_Pred[i-3] ~ dlnorm(logR_Fit[i-3],tau)
    R_Fit[i-3] <- exp(logR_Fit[i-3])
  }
  
  logA ~ dnorm(logA_mean,logA_tau)     	# prior for alpha
  beta0 ~ dnorm(B_means[1],B_taus[1])	# prior for beta0  
  beta1 ~ dnorm(B_means[2],B_taus[2])# prior for beta1   
  beta2 ~ dnorm(B_means[3],B_taus[3])# prior for beta2   
  beta3 ~ dnorm(B_means[4],B_taus[4])	# prior for beta3  
  tau ~ dgamma(Sig_Gam_Dist, Sig_Gam_Dist)		# prior for precision parameter
  sigma <- 1/sqrt(tau)		# transformation of precision to sd
  
}





  
  

