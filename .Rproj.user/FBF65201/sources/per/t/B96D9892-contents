
#########################################################################
# SR MODEL SPECIFICATIONS FOR FRASER SOCKEYE FORECASTING
#########################################################################

# This file reads in all the functions and 
# at the end creates a nested list object with the 
# following elements for each model:
# MCMC model form:  BUGS/JAGS model as a function
# R model nodelist : lists of nodes to trck for different settings
# R model form - single inputs : R function to plug the SR parameters into
# R model form - array inputs : R function to plug the SR parameters into (pars in large array form)
# Inits values: default inits values as a list (same for all stocks, but can override later, not used in JAGS!)
# List of parameters and priors: 2 sepafrate lists
# Data Object Check: R function to use to check the data list object 
# merge all the pieces for the model


# BUGS/JAGS model forms and data file list object formats
# are adapted from code by Cass, Michielsens, Grant, MacDonald and others.

# BDavis 5/6/2020
# Models in this file:
#1) RICKER MODEL  WITH FIXED PRIORS (Rickerfp)
#2) BASIC RICKER MODEL WITH USER-SPECIFIED PRIORS (Ricker) -- note process error prior not user specified
#3) POWER MODEL (Power)
#4) CYCLE LINE RICKER MODEL (RickerCycleModel) - confused why need this? Two versions?
#5) RICKER MODEL WITH CoVARIATE (ENVIRONMENTAL VARIABLES) (RickerCov) -- looks like can only accommodate 1 variable
#6) LARKIN MODEL (Larkin)
#7) KALMAN FILTER (Kalman)
#8) SMOLT-JACK (SmoltJack)
#9) POWER_JUV COVARIATE MODEL (PowerCov)
#10) SIBLING MODEL (Sibling) -- 2 versions, does not work with current FC Code



#### ========================================== RICKER MODEL  WITH FIXED PRIORS (FP)============================== ####

# MCMC model form 

Rickerfp.model.MCMC <- function(){
    for (i in 1:N) {                       #loop over N sample points
      R_Obs[i] ~ dlnorm(Y[i],tau_R)          #likelihood -> predicted value for NA in data set
      Y[i] <- RS[i] +log(S[i])               # calc log(R) - fitted values
      RS[i] <- alpha - beta * S[i]           # ricker model 
      Rep[i] ~ dlnorm(Y[i],tau_R)
   }

    Rec_track_BY1<- Rep[N-track-3]  
    Rec_track_BY2<- Rep[N-track-2]  
    Rec_track_BY3<- Rep[N-track-1]  
    Rec_track_BY4<- Rep[N-track]
        
    alpha ~ dnorm(0,0.0001)        #prior for alpha (actually ln.alpha!) -> fix? 
	  #could change to dnorm(p.alpha,tau_alpha), but how to feed into r2winbugs?
    beta <-1/C					   # prior for beta
    C~ dlnorm(1,0.1)       			   # prior for beta -> could change to dlnorm(p.beta, tau_beta)
    tau_R ~ dgamma(0.001,0.001)                    #prior for precision parameter
    sigma <- 1/sqrt(tau_R) 			# changed based on Fleishman and Evenson     	

}


# List of nodes that may be tracked under various commands
Rickerfp.model.nodelist <- list(All=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta","sigma", "Rep","deviance"),
                              RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta","sigma","deviance"),
                              RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","deviance"))

# Do we need one form for single par and another to work with multi-stock arrays and all MCMC outputs in one go?
# OR can we use this with apply()?

# R model form - single inputs

Rickerfp.model.R <- function(spn,alpha,beta,sigma,ln.alpha=FALSE,predicted=FALSE){
	if(ln.alpha){alpha <- ln(alpha)}
	rs <- exp(alpha-beta*spn)
	if(predicted){rs <- rs + exp(rnorm(1,rec,sigma))} # need to test this step
	rec <- rs *spn
	rec
}

# R model form - array inputs

Rickerfp.model.R.array  <- function(Spn,alpha.pars,beta.pars,out="quants"){
# Spn is an array of dim: years or intervals by pop
  # alpha.pars,beta.pars are arrays of dim: pops by par.est (bootstrap or MCMC estimates)
  # output -> if out = quants then for each Spn level, calculate R based on each alternative par set, then take quantiles of R
  #           if out = samples then just output fitted R
  
  vals <- dimnames(Spn)[[1]]  # years or index of Spn interval 
  pops <- dimnames(alpha.pars)[[1]]
  reps <- dimnames(alpha.pars)[[2]]
  
  spn.arr <- array(NA,dim= c(length(vals),length(pops),length(reps)),dimnames= list(vals,pops,reps ))
  rec.arr<- alpha.arr <- beta.arr <- spn.arr
  
  
  # NEED TO FIX THIS BECAUSE SLOW!!!!!!!!!
  #populate each slice of the array NOTE: should be able to this in one step, but keeps having bug, so do brute force way for now  
  for(repl in reps){
  	for(pop in pops){               # need this second loop in case the order of pop is not the same!!!	
  		
  		#print(repl);print(pop)
  		#print(Spn[1:2,1:2])
  		spn.arr[,pop,repl] <- as.matrix(Spn[,pop])
  	}
  } 
  
  for(val in vals){
  	alpha.arr[val,,] <- as.matrix(alpha.pars) 
  	beta.arr[val,,] <- as.matrix(beta.pars) 
  } 
       
  # apply ricker calculation  
  rec.arr <- exp(log(alpha.arr)-beta.arr*spn.arr)*spn.arr
  
  # INCLUDE OPTION TO ADD RANDOM ERROR BASED ON SIGMA PAR
  # BUILD IN MULTIPLE SIGMA SAMPLES BY FEEDING IN REPLICATES OF PAR DATA SET
  
  # calculate percentiles across repl values
  if(out=="quants"){out.arr <- apply(rec.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)}
  if(out=="samples"){out.arr=rec.arr}
  
  out.arr

} # end Ricker.model.R.array 



# Inits values

Rickerfp.model.inits<-list(list(tau_R=3, C=1),list(tau_R=7, C=2))

# List of parameters and priors
Rickerfp.par.list<-c("S","R_Obs","N","track")
Rickerfp.prior.list<-NULL  # priors are hardwired in the model form above?


# Data Object Check
Rickerfp.data.check <- function(x,par.list,prior.list){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  chk.list <- c(par.list,prior.list)  
  chk <- !(chk.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",chk.list[chk]))}

}



# merge all the pieces for the Ricker Model
RickerfpModel.List <- list(MCMC.Model=Rickerfp.model.MCMC, nodelist=Rickerfp.model.nodelist, R.Model=Rickerfp.model.R, 
							R.Model.array=Rickerfp.model.R.array, inits=Rickerfp.model.inits, Data.Check=Rickerfp.data.check,
							parlist = Rickerfp.par.list, priorlist = Rickerfp.prior.list)
							
							
							
											
							
							
							
							

# ------------------------------------------- END RICKER FP MODEL --------------------------------------------------------------# 


#### ====================================== BASIC RICKER MODEL WITH USER-SPECIFIED PRIORS ============================== ####
# MCMC model form 

Ricker.model.MCMC <- function(){
    for (i in 1:N) {                       #loop over N sample points
      R_Obs[i] ~ dlnorm(Y[i],tau_R)          #likelihood -> predicted value for NA in data set
      Y[i] <- RS[i] +log(S[i])               # calc log(R) - fitted values
      RS[i] <- alpha - beta * S[i]           # ricker model 
      Rep[i] ~ dlnorm(Y[i],tau_R)
    }

    Rec_track_BY1<- Rep[N-track-3]  
    Rec_track_BY2<- Rep[N-track-2]  
    Rec_track_BY3<- Rep[N-track-1]  
    Rec_track_BY4<- Rep[N-track]
    
    alpha ~ dnorm(p.alpha,tau_alpha)        #prior for alpha (actually ln.alpha!) -> fix? 
    beta <-1/C					   # prior for beta
    C~ dlnorm(p.beta, tau_beta)       			   # prior for beta 
    tau_R ~ dgamma(0.001,0.001)                    #prior for precision parameter
    sigma <- 1/sqrt(tau_R) 			# changed based on Fleishman and Evenson     	

}

# Do we need one form for single par and another to work with multi-stock arrays and all MCMC outputs in one go?
# OR can we use this with apply()?


# R model nodelist 
# List of nodes that may be tracked under various commands
Ricker.model.nodelist <- list(All=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta","sigma", "Rep","deviance"),
                              RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta","sigma","deviance"),
                              RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","deviance"))

# Do we need one form for single par and another to work with multi-stock arrays and all MCMC outputs in one go?
# OR can we use this with apply()?



# R model form - single inputs

Ricker.model.R <- function(spn,alpha,beta,sigma,ln.alpha=FALSE,predicted=FALSE){
	if(ln.alpha){alpha <- ln(alpha)}
	rs <- exp(alpha-beta*spn)
	if(predicted){rs <- rs + exp(rnorm(1,rec,sigma))} # need to test this step
	rec <- rs *spn
	rec
}

# R model form - array inputs 
Ricker.model.R.array  <- function(Spn,alpha.pars,beta.pars,out="quants"){
# Spn is an array of dim: years or intervals by pop
  # alpha.pars,beta.pars are arrays of dim: pops by par.est (bootstrap or MCMC estimates)
  # output -> if out = quants then for each Spn level, calculate R based on each alternative par set, then take quantiles of R
  #           if out = samples then just output fitted R
  
  vals <- dimnames(Spn)[[1]]  # years or index of Spn interval 
  pops <- dimnames(alpha.pars)[[1]]
  reps <- dimnames(alpha.pars)[[2]]
  
  spn.arr <- array(NA,dim= c(length(vals),length(pops),length(reps)),dimnames= list(vals,pops,reps ))
  rec.arr<- alpha.arr <- beta.arr <- spn.arr
  
  
  # NEED TO FIX THIS BECAUSE SLOW!!!!!!!!!
  #populate each slice of the array NOTE: should be able to this in one step, but keeps having bug, so do brute force way for now  
  for(repl in reps){
  	for(pop in pops){               # need this second loop in case the order of pop is not the same!!!	
  		
  		#print(repl);print(pop)
  		#print(Spn[1:2,1:2])
  		spn.arr[,pop,repl] <- as.matrix(Spn[,pop])
  	}
  } 
  
  for(val in vals){
  	alpha.arr[val,,] <- as.matrix(alpha.pars) 
  	beta.arr[val,,] <- as.matrix(beta.pars) 
  } 
       
  # apply ricker calculation  
  rec.arr <- exp(log(alpha.arr)-beta.arr*spn.arr)*spn.arr
  
  # INCLUDE OPTION TO ADD RANDOM ERROR BASED ON SIGMA PAR
  # BUILD IN MULTIPLE SIGMA SAMPLES BY FEEDING IN REPLICATES OF PAR DATA SET
  
  # calculate percentiles across repl values
  if(out=="quants"){out.arr <- apply(rec.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)}
  if(out=="samples"){out.arr=rec.arr}
  
  out.arr

} # end Ricker.model.R.array 

   
# Inits values
Ricker.model.inits<-list(list(tau_R=3, C=1),list(tau_R=7, C=2))


# List of parameters and priors
Ricker.par.list<-c("S","R_Obs","N","track")
Ricker.prior.list<-c("p.alpha","tau_alpha","p.beta","tau_beta")


# Data Object Check
Ricker.data.check <- function(x,par.list,prior.list){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  chk.list <- c(par.list,prior.list) 
  print("Required Inputs:"); print(chk.list)
  chk <- !(chk.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",chk.list[chk]))}

}


# merge all the pieces for the Ricker Model
RickerModel.List <- list(MCMC.Model=Ricker.model.MCMC, nodelist=Ricker.model.nodelist, R.Model=Ricker.model.R, 
							R.Model.array=Ricker.model.R.array, inits=Ricker.model.inits, Data.Check=Ricker.data.check,
							parlist = Ricker.par.list, priorlist = Ricker.prior.list)

							

# ----------------------------------------- END RICKER MODEL ----------------------------------------------------------#


#### ========================================== POWER MODEL  =========================================== ####



# NOTE: SHOULD CONSIDER CHANGING Y_log to Y here, or making the opposite change in the Ricker models
# in both the Ricker and power models the Y variable is in log-transform, and labelling should be consistent across



# MCMC model form 

Power.model.MCMC <- function(){
  for (i in 1:N) {                             # loop over N sample points
    R_Obs[i] ~ dlnorm(Y_log[i],tau_R)          # likelihood -> predicted value for NA in data set
    Y_log[i] <- alpha - beta * log(S[i])       # power model
    Rep[i] ~ dlnorm(Y_log[i],tau_R)
  }
  
  Rec_track_BY1<- Rep[N-track-3]  
  Rec_track_BY2<- Rep[N-track-2]  
  Rec_track_BY3<- Rep[N-track-1]  
  Rec_track_BY4<- Rep[N-track]
  
  alpha ~ dnorm(p.alpha,tau_alpha)             # prior for alpha
  beta ~ dnorm(p.beta,tau_beta)                # prior for beta
  tau_R ~ dgamma(0.001,0.001)                  # prior for precision parameter
  #sigma <- 1/pow(tau_R,2)                     # old sigma calc
  sigma <- 1/sqrt(tau_R)   		               # changed based on Fleishman and Evenson     	
  
}


# List of nodes that may be tracked under various commands
Power.model.nodelist <- list(All=c("alpha","beta","sigma", "Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","Rep","deviance"),
                              RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta","sigma","deviance"),
                              RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","deviance"))

							  
							  
# R model form - single inputs

Power.model.R <- function(spn,alpha,beta,sigma,predicted=FALSE){
  rs <- exp(alpha-beta*log(spn))
  if(predicted){rs <- rs + exp(rnorm(1,rec,sigma))} # need to test this step
  rec <- rs 
  rec
}


# R model form - array inputs

Power.model.R.array  <- function(Spn,alpha.pars,beta.pars,out="quants"){
  # Spn is an array of dim: years or intervals by pop
  # alpha.pars,beta.pars are arrays of dim: pops by par.est (bootstrap or MCMC estimates)
  # output -> if out = quants then for each Spn level, calculate R based on each alternative par set, then take quantiles of R
  #           if out = samples then just output fitted R
  
  vals <- dimnames(Spn)[[1]]  # years or index of Spn interval 
  pops <- dimnames(alpha.pars)[[1]]
  reps <- dimnames(alpha.pars)[[2]]
  
  spn.arr <- array(NA,dim= c(length(vals),length(pops),length(reps)),dimnames= list(vals,pops,reps ))
  rec.arr<- alpha.arr <- beta.arr <- spn.arr
  
  
  # NEED TO FIX THIS BECAUSE SLOW!!!!!!!!!
  #populate each slice of the array NOTE: should be able to this in one step, but keeps having bug, so do brute force way for now  
  for(repl in reps){
    for(pop in pops){               # need this second loop in case the order of pop is not the same!!!	
      
      #print(repl);print(pop)
      #print(Spn[1:2,1:2])
      spn.arr[,pop,repl] <- as.matrix(Spn[,pop])
    }
  } 
  
  for(val in vals){
    alpha.arr[val,,] <- as.matrix(alpha.pars) 
    beta.arr[val,,] <- as.matrix(beta.pars) 
  } 
  
  # apply power calculation  
  rec.arr <- exp(alpha.arr-beta.arr*log(spn.arr))
  
  # INCLUDE OPTION TO ADD RANDOM ERROR BASED ON SIGMA PAR
  # BUILD IN MULTIPLE SIGMA SAMPLES BY FEEDING IN REPLICATES OF PAR DATA SET
  
  # calculate percentiles across repl values
  if(out=="quants"){out.arr <- apply(rec.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)}
  if(out=="samples"){out.arr=rec.arr}
  
  out.arr
  
} # end Power.model.R.array 


# Inits values

Power.model.inits<-list(list(tau_R=3),list(tau_R=7))



# List of parameters and priors

Power.par.list<-c("S","R_Obs","N","track")
Power.prior.list<-c("p.alpha","tau_alpha","p.beta","tau_beta")



# Data Object Check


Power.data.check <- function(x,par.list,prior.list){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  chk.list <- c(par.list,prior.list) 
  print("Required Inputs:"); print(chk.list)
  chk <- !(chk.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",chk.list[chk]))}

}




# merge all the pieces for the Power Model
PowerModel.List <- list(MCMC.Model=Power.model.MCMC, nodelist=Power.model.nodelist, R.Model=Power.model.R, 
						R.Model.array=Power.model.R.array, inits=Power.model.inits, Data.Check=Power.data.check,
							parlist = Power.par.list, priorlist = Power.prior.list)


# --------------------------------------- END POWER MODEL SECTION ------------------------------------------------------# 


#### ============================================== CYCLE LINE RICKER MODEL ==============================================####

###  NOTE: THERE ARE 2 VERSIONS OF THE RICKER CYCLE MODEL CODE BELOW
###  1) IDENTICAL TO BASIC RICKER MODEL, JUST THAT DATA IS SUBSET BASED ON fc.yr and avg.gen BEFORE CALLING THE MCMC
###  2) VERSION ORIGINALLY CREATED BY BRONWYN MACDONALD TO BE CONSISTENT WITH OLD MCMC CODE doing 2 brood years at once
###  => FURTHER DISCUSSION REQUIRED
#### IF go with option 1, don't need 2 copies, just call up the Ricker model
#### IF go  with option 2, need to update the model code here and fundamentally rethink how this is handled (given the new structure)
###  OPTION 1 works because not tracking Ret anymore, but doing Rec->Ret calcs after (i.e. 2 different brood years), BUT HAVE RUN INTO 
###  PROBLEMS WITH MERGING TEH TWO POSTERIORS (RESULTS NOT MATCHING UP WITH OLD CODE RESULTS)
###  => SWITCHING OVER TO OPTION 2


########################################
# To use Option 1, set this section to TRUE

if(FALSE){

# merge all the pieces for the Ricker Cycle Model using the same components as previously read in for the basic Ricker
RickerCycleModel.List <- list(MCMC.Model=Ricker.model.MCMC, nodelist=Ricker.model.nodelist, R.Model=Ricker.model.R, 
							R.Model.array=Ricker.model.R.array, inits=Ricker.model.inits, Data.Check=Ricker.data.check,
							parlist = Ricker.par.list, priorlist = Ricker.prior.list)

}  # end Ricker cycle model specified the same as basic Ricker



########################################
# To use Option 2, set this section to TRUE

if(TRUE){

RickerCycle.model.MCMC <- function(){
  
  for (i in 1:N4) {                       #loop over N sample points
    R4_Obs[i] ~ dlnorm(Y4[i],tau_R4)          #likelihood 
    Y4[i] <- RS4[i] +log(S4[i])               #calc log(R)
    RS4[i] <- alpha4 - beta4 * S4[i]           # ricker model 
    Rep4[i] ~ dlnorm(Y4[i],tau_R4)
    #year4[i]<-i    # ??? why need this? testing to just leave it out 
  }
  
  for (j in 1:N5) {                       #loop over N sample points
    R5_Obs[j] ~ dlnorm(Y5[j],tau_R5)          #likelihood 
    Y5[j] <- RS5[j] +log(S5[j])               #calc log(R)
    RS5[j] <- alpha5 - beta5 * S5[j]           # ricker model 
    Rep5[j] ~ dlnorm(Y5[j],tau_R5)
    #year5[j]<-j    # ??? why need this?    #  testing to just leave it out 
  }
  
	# note: think that track variable is obsolete altogether (i.e always set to 0), so tried not even building it in here)
	# however: it still needs to be in the data.obj for subsequent steps (for now, until removed throughout)
	# and WinBUGS/OpenBUGS crash if data.obj has a variable that's not used in the model....
	# So putting it in here as well for now
	
	# NOTE: THE BY lining up is turning out to be tricky -> NEED TO TRACE THROUGH CAREFULLY UNTIL TEH FINAL RETURN CALC
 

    Rec_track_BY1 <- Rep4[N4-track] * 0.00000001    # 7 yr olds => need this to be var else node not tracked
    Rec_track_BY2 <- Rep4[N4-track]  * 0.00000001   # 6 yr olds => need this to be var else node not tracked
    Rec_track_BY3 <- Rep5[N5-track]   # 5 yr olds       
    Rec_track_BY4 <- Rep4[N4-track]   # 4 yr olds 
  
  
  alpha4 ~ dnorm(p.alpha,tau_alpha)               #prior for alpha
  alpha5 ~ dnorm(p.alpha,tau_alpha)  
  beta4 <-1/C4
  beta5 <- 1/C5
  C4~ dlnorm(p.beta,tau_beta)           #prior for stock size at maximum mean recruitment
  C5~ dlnorm(p.beta,tau_beta)           #prior for stock size at maximum mean recruitment
  tau_R4 ~ dgamma(0.001,0.001)                    #prior for precision parameter
  tau_R5 ~ dgamma(0.001,0.001) 
  sigma4 <- 1/sqrt(tau_R4)              # changed based on Fleishman and Evenson       
  sigma5 <- 1/sqrt(tau_R5)              # changed based on Fleishman and Evenson
}


# List of nodes that may be tracked under various commands
RickerCycle.model.nodelist <- list(All=c("alpha4", "alpha5","beta4","beta5", "sigma4", "sigma5", 
                                         "Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","Rep4","Rep5","deviance"),
                             RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha4","beta4","sigma4", 
                                      "alpha5","beta5","sigma5","deviance"),
                             RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","deviance"))



# R model form - single inputs
RickerCycle.model.R <- function(spn4,spn5,alpha4,alpha5,beta4,beta5,sigma,ln.alpha=FALSE,predicted=FALSE){
  if(ln.alpha){
    alpha4 <- ln(alpha4)
    alpha5 <- ln(alpha5)
  }
  rs4 <- exp(alpha4-beta4*spn4)
  rs5 <- exp(alpha5-beta5*spn5)  
  
  if(predicted){
    rs4 <- rs4 + exp(rnorm(1,rec4,sigma))
    rs5 <- rs5 + exp(rnorm(1,rec5,sigma))
  } # need to test this step                                *** this step doesnt make sense, rec hasnt been calculated yet
  rec4 <- rs4 *spn4
  rec5 <- rs5*spn5
  list(rec4=rec4, rec5=rec5)
}


# R model form - array inputs

RickerCycle.model.R.array  <- function(Spn4,Spn5,alpha4.pars,alpha5.pars,beta4.pars,beta5.pars,out="quants"){  # ************** NOT TESTED
  # Spn is an array of dim: years or intervals by pop
  # alpha.pars,beta.pars are arrays of dim: pops by par.est (bootstrap or MCMC estimates)
  # output -> if out = quants then for each Spn level, calculate R based on each alternative par set, then take quantiles of R
  #           if out = samples then just output fitted R
  # env will likely be just one vector of values across years BUT g.pars will be array because g parameters will vary by stock and iteration
  
  vals <- dimnames(Spn4)[[1]]  # years or index of Spn interval 
  pops <- dimnames(alpha4.pars)[[1]]
  reps <- dimnames(alpha4.pars)[[2]]
  
  spn4.arr <- array(NA,dim= c(length(vals),length(pops),length(reps)),dimnames= list(vals,pops,reps ))
  rec4.arr <- rec5.arr <- alpha4.arr <- alpha5.arr <- beta4.arr <- beta5.arr <- spn5.arr <- spn4.arr
  
  
  # NEED TO FIX THIS BECAUSE SLOW!!!!!!!!!                                                        
  #populate each slice of the array NOTE: should be able to this in one step, but keeps having bug, so do brute force way for now  
  for(repl in reps){
    for(pop in pops){               # need this second loop in case the order of pop is not the same!!!  
      
      #print(repl);print(pop)
      #print(Spn[1:2,1:2])
      spn4.arr[,pop,repl] <- as.matrix(Spn4[,pop])
      spn5.arr[,pop,repl] <- as.matrix(Spn5[,pop])
    }
  } 
  
  for(val in vals){
    alpha4.arr[val,,] <- as.matrix(alpha4.pars);  alpha5.arr[val,,] <- as.matrix(alpha5.pars)  
    beta4.arr[val,,] <- as.matrix(beta4.pars);  beta5.arr[val,,] <- as.matrix(beta5.pars)
  } 
  
  # apply ricker calculation
  
  rec4.arr <- exp(log(alpha4.arr)-beta4.arr*spn4.arr)*spn4.arr
  rec5.arr <- exp(log(alpha5.arr)-beta5.arr*spn5.arr)*spn5.arr
  # INCLUDE OPTION TO ADD RANDOM ERROR BASED ON SIGMA PAR
  # BUILD IN MULTIPLE SIGMA SAMPLES BY FEEDING IN REPLICATES OF PAR DATA SET
  
  # calculate percentiles across repl values
  if(out=="quants"){out4.arr <- apply(rec4.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)
                    out5.arr <- apply(rec5.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)
                    out.arr=list(rec4.arr=out4.arr, rec5.arr=out5.arr)
  }
  if(out=="samples"){out.arr=list(rec4.arr=rec4.arr, rec5.arr=rec5.arr)}
  
  out.arr
  
} # end RickerCycle.model.R.array



# Inits values  
RickerCycle.model.inits<-list(list(tau_R4=3, tau_R5=3, C4=1, C5=1),list(tau_R4=7, tau_R5=7, C4=2, C5=2))


# List of parameters and priors
RickerCycle.par.list<-c("S4","S5","R4_Obs","R5_Obs","N4","N5","track")
RickerCycle.prior.list<-c("p.alpha","tau_alpha","p.beta","tau_beta")

# Data Object Check
RickerCycle.data.check <- function(x,par.list,prior.list){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  chk.list <- c(par.list,prior.list) 
  print("Required Inputs:"); print(chk.list)
  chk <- !(chk.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",chk.list[chk]))}

}

# merge  pieces for the Ricker Cycle  Model
RickerCycleModel.List <- list(MCMC.Model=RickerCycle.model.MCMC, nodelist=RickerCycle.model.nodelist, R.Model=RickerCycle.model.R, 
							R.ModelCycle.array=Ricker.model.R.array, inits=RickerCycle.model.inits, Data.Check=RickerCycle.data.check,
							parlist = RickerCycle.par.list, priorlist = RickerCycle.prior.list)
							
							
}  # end Ricker cycle model with 2 year classes at the same time






# ------------------------------------------- END RICKER CYCLE MODEL --------------------------------------------------# 


#### ================================== RICKER MODEL WITH CoVARIATE (ENVIRONMENTAL VARIABLES) ============================ ####

RickerCov.model.MCMC <- function(){
    for (i in 1:N) {                                    #loop over N sample points      
      R_Obs[i] ~ dlnorm(Y[i],tau_R)                   #likelihood 
      Y[i] <- RS[i] + log(S[i]) 
      RS[i] <- alpha - beta *S[i] + g * env[i]    #Ricker model  
      Rep[i] ~ dlnorm(Y[i],tau_R)
    }
    
    Rec_track_BY1<- Rep[N-track-3]  
    Rec_track_BY2<- Rep[N-track-2]  
    Rec_track_BY3<- Rep[N-track-1]  
    Rec_track_BY4<- Rep[N-track]
    
    g ~ dnorm(p.g,tau_g)     
    alpha ~ dnorm(p.alpha,tau_alpha)               #prior for alpha
    beta <- 1/C					   # prior for beta
    C~ dlnorm(p.beta, tau_beta)       			   # prior for beta 
    tau_R ~ dgamma(0.001,0.001)                      #prior for precision parameter
    sigma <- 1/sqrt(tau_R)                            	# changed based on Fleishman and Evenson     		

}

	# 2016 version (matching "envrickermod.txt" from old code)
	#beta ~ dlnorm(p.beta,tau_beta)           #prior for stock size at maximum mean recruitment
	
	# 2017 testing (matching "srmodel.txt" from old code)
	#beta <- 1/C
	#C ~ dlnorm(p.beta,tau_beta)  
	#beta <-1/C
    #C~ dlnorm(1,0.1)   
	

	
# List of nodes that may be tracked under various commands
RickerCov.model.nodelist <- list(All=c("alpha","beta","g","sigma", "Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","Rep","deviance"),
                             RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta","g","sigma","deviance"),
                             RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","deviance"))


# R model form - single inputs

RickerCov.model.R <- function(spn,alpha,beta,g,env,sigma,ln.alpha=FALSE,predicted=FALSE){
  if(ln.alpha){alpha <- ln(alpha)}
  rs <- exp(alpha-beta*spn+g*env)
  if(predicted){rs <- rs + exp(rnorm(1,rec,sigma))} # need to test this step
  rec <- rs *spn
  rec
}


# R model form - array inputs

RickerCov.model.R.array  <- function(Spn,alpha.pars,beta.pars,g.pars,env,out="quants"){  # ************** NOT COMPLETED
  # Spn is an array of dim: years or intervals by pop
  # alpha.pars,beta.pars are arrays of dim: pops by par.est (bootstrap or MCMC estimates)
  # output -> if out = quants then for each Spn level, calculate R based on each alternative par set, then take quantiles of R
  #           if out = samples then just output fitted R
  # env will likely be just one vector of values across years BUT g.pars will be array because g parameters will vary by stock and iteration
  
  vals <- dimnames(Spn)[[1]]  # years or index of Spn interval 
  pops <- dimnames(alpha.pars)[[1]]
  reps <- dimnames(alpha.pars)[[2]]
  
  spn.arr <- array(NA,dim= c(length(vals),length(pops),length(reps)),dimnames= list(vals,pops,reps ))
  rec.arr<- alpha.arr <- beta.arr <- g.arr <- env.arr <- spn.arr
  
  
  # NEED TO FIX THIS BECAUSE SLOW!!!!!!!!!                                                        
  #populate each slice of the array NOTE: should be able to this in one step, but keeps having bug, so do brute force way for now  
  for(repl in reps){
    for(pop in pops){               # need this second loop in case the order of pop is not the same!!!	
      
      #print(repl);print(pop)
      #print(Spn[1:2,1:2])
      spn.arr[,pop,repl] <- as.matrix(Spn[,pop])
      env.arr[,pop,repl] <- env
    }
  } 
  
  for(val in vals){
    alpha.arr[val,,] <- as.matrix(alpha.pars) 
    beta.arr[val,,] <- as.matrix(beta.pars) 
    g.arr[val,,] <- as.matrix(g.pars)
  } 
  
  # apply ricker calculation
  
  rec.arr <- exp(log(alpha.arr)-beta.arr*spn.arr+g.arr*env.arr)*spn.arr
  
  # INCLUDE OPTION TO ADD RANDOM ERROR BASED ON SIGMA PAR
  # BUILD IN MULTIPLE SIGMA SAMPLES BY FEEDING IN REPLICATES OF PAR DATA SET
  
  # calculate percentiles across repl values
  if(out=="quants"){out.arr <- apply(rec.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)}
  if(out=="samples"){out.arr=rec.arr}
  
  out.arr
  
} # end RickerCov.model.R.array


# Inits values
# RickerCov.model.inits<-list(list(tau_R=3, beta=0.1),list(tau_R=7,beta=0.2)) # OLD VERSION: CRASHES OPNEBUGS AND JAGS (Can't initialize of beta is derived!)
RickerCov.model.inits <- list(list(tau_R=3, C=1),list(tau_R=7, C=2)) # instead use same as for basic ricker above

# List of parameters and priors
RickerCov.par.list<-c("S","R_Obs","N","env","track")
RickerCov.prior.list<-c("p.alpha","tau_alpha","p.beta","tau_beta","p.g","tau_g" )


# Data Object Check
RickerCov.data.check <- function(x,par.list,prior.list){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  chk.list <- c(par.list,prior.list) 
  print("Required Inputs:"); print(chk.list)
  chk <- !(chk.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",chk.list[chk]))}

}



# merge all the pieces for the Ricker COV Model
RickerCovModel.List <- list(MCMC.Model=RickerCov.model.MCMC, nodelist=RickerCov.model.nodelist, R.Model=RickerCov.model.R, 
							R.Model.array=RickerCov.model.R.array, inits=RickerCov.model.inits, Data.Check=RickerCov.data.check,
							parlist = RickerCov.par.list, priorlist = RickerCov.prior.list)

# ------------------------------------------ END EV RICKER MODEL --------------------------------------------------------#


#### =================================================== LARKIN MODEL ======================================================####

# MCMC model form

Larkin.model.MCMC <- function(){
    for(i in 4:N) {               		#loop over N sample points
        R_Obs[i] ~ dlnorm(Y[i],tau_R)		        #likelihood function  NOTE: Changed R[i] to Y[i] as in  Ricker (old version gave an error!)
        Y[i] <- RS_log[i] + log(S[i])						#NOTE: Changed R[i] to Y[i] as in Ricker(old version gave an error!)
        RS_log[i] <-alpha-beta0*S[i]-beta1*S[i-1]-beta2*S[i-2]-beta3*S[i-3]  
        Rep[i] ~ dlnorm(Y[i],tau_R)
    }
                                	
    Rec_track_BY1<- Rep[N-track-3]  
    Rec_track_BY2<- Rep[N-track-2]  
    Rec_track_BY3<- Rep[N-track-1]  
    Rec_track_BY4<- Rep[N-track]
    
    alpha ~ dnorm(p.alpha,tau_alpha)     	# prior for alpha
    beta0 ~ dnorm(p.beta0,tau_beta0)	# prior for beta0  # to constrain %_%I(0,)
    beta1 ~ dnorm(p.beta1,tau_beta1)# prior for beta1   # to constrain %_%I(0,)
    beta2 ~ dnorm(p.beta2,tau_beta2)# prior for beta2   # to constrain %_%I(0,)
    beta3 ~ dnorm(p.beta3,tau_beta3)	# prior for beta3  # to constrain %_%I(0,)
    tau_R ~ dgamma(0.001,0.001)		# prior for precision parameter
    sigma <- 1/sqrt(tau_R)		# transformation of precision to sd

}


# List of nodes that may be tracked under various commands
Larkin.model.nodelist <- list(All=c("alpha","beta0", "beta1", "beta2", "beta3","sigma", "Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","Rep","deviance"),
                                 RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta0", "beta1", "beta2", "beta3","sigma","deviance"),
                                 RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","deviance"))


# R model form - single inputs

Larkin.model.R <- function(spn0,spn1,spn2,spn3,alpha,beta0,beta1,beta2,beta3,sigma,ln.alpha=FALSE,predicted=FALSE){
  if(ln.alpha){alpha <- ln(alpha)}
  rs <- exp(alpha-beta0*spn0-beta1*spn1-beta2*spn2-beta3*spn3)
  if(predicted){rs <- rs + exp(rnorm(1,rec,sigma))} # need to test this step
  rec <- rs *spn
  rec
}


# R model form - array inputs ****************************************************************************************


Larkin.model.R.array  <- function(Spn0,Spn1,Spn2,Spn3,alpha.pars,beta0.pars,bea1.pars,beta2.pars,beta3,pars,out="quants"){  # ************** NOT COMPLETED
  # Spn is an array of dim: years or intervals by pop
  # alpha.pars,beta.pars are arrays of dim: pops by par.est (bootstrap or MCMC estimates)
  # output -> if out = quants then for each Spn level, calculate R based on each alternative par set, then take quantiles of R
  #           if out = samples then just output fitted R
  # env will likely be just one vector of values across years BUT g.pars will be array because g parameters will vary by stock and iteration
  
  vals <- dimnames(Spn0)[[1]]  # years or index of Spn interval 
  pops <- dimnames(alpha.pars)[[1]]
  reps <- dimnames(alpha.pars)[[2]]
  
  spn0.arr <- array(NA,dim= c(length(vals),length(pops),length(reps)),dimnames= list(vals,pops,reps ))
  rec.arr<- alpha.arr <- beta0.arr <- beta1.arr <- beta2.arr <- beta3.arr <- spn1.arr <- spn2.arr <- spn3.arr <- spn0.arr
  
  
  # NEED TO FIX THIS BECAUSE SLOW!!!!!!!!!                                                        
  #populate each slice of the array NOTE: should be able to this in one step, but keeps having bug, so do brute force way for now  
  for(repl in reps){
    for(pop in pops){               # need this second loop in case the order of pop is not the same!!!  
      
      #print(repl);print(pop)
      #print(Spn[1:2,1:2])
      spn0.arr[,pop,repl] <- as.matrix(Spn0[,pop])
      spn1.arr[,pop,repl] <- as.matrix(Spn1[,pop])
      spn2.arr[,pop,repl] <- as.matrix(Spn2[,pop])
      spn3.arr[,pop,repl] <- as.matrix(Spn3[,pop])
    }
  } 
  
  for(val in vals){
    alpha.arr[val,,] <- as.matrix(alpha.pars) 
    beta0.arr[val,,] <- as.matrix(beta0.pars) 
    beta1.arr[val,,] <- as.matrix(beta1.pars)
    beta2.arr[val,,] <- as.matrix(beta2.pars)
    beta3.arr[val,,] <- as.matrix(beta3.pars)
  } 
  
  # apply ricker calculation
  
  rec.arr <- exp(log(alpha.arr)-beta0.arr*spn0.arr-beta1.arr*spn1.arr-beta2.arr*spn2.arr-beta3.arr*spn3.arr)*spn0.arr
  
  # INCLUDE OPTION TO ADD RANDOM ERROR BASED ON SIGMA PAR
  # BUILD IN MULTIPLE SIGMA SAMPLES BY FEEDING IN REPLICATES OF PAR DATA SET
  
  # calculate percentiles across repl values
  if(out=="quants"){out.arr <- apply(rec.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)}
  if(out=="samples"){out.arr=rec.arr}
  
  out.arr
  
} # end Larkin.model.R.array



# Inits values

Larkin.model.inits<-list(list(tau_R=3),list(tau_R=7))


# List of parameters and priors
Larkin.par.list<-c("S","R_Obs","N","track")
Larkin.prior.list<-c("p.alpha","tau_alpha","p.beta0","p.beta1","p.beta2","p.beta3","tau_beta0","tau_beta1","tau_beta2","tau_beta3")

# Data Object Check

Larkin.data.check <- function(x,par.list,prior.list){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  
  chk.list <- c(par.list,prior.list) 
  print("Required Inputs:"); print(chk.list)
  chk <- !(chk.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",chk.list[chk]))}

  
}


# merge all the pieces for the Larkin Model

LarkinModel.List <- list(MCMC.Model=Larkin.model.MCMC, nodelist=Larkin.model.nodelist, R.Model=Larkin.model.R, R.Model.array=NA, inits=Larkin.model.inits,
                         Data.Check=Larkin.data.check,parlist = Larkin.par.list, priorlist = Larkin.prior.list)



# --------------------------------------END LARKIN MODEL --------------------------------------- #


#### ================================================ KALMAN FILTER MODEL ======================================================####

# MCMC model form

Kalman.model.MCMC <- function(){ 
    for (i in 1:N){
      R_Obs[i] ~ dlnorm(Y[i],tau_R)          #likelihood 
      Y[i] <- RS[i] +log(S[i])               #calc log(R)
      RS[i] <- alpha[i] - beta * S[i] + v[i]
      v[i] ~dnorm(0, tauv)  
      Rep[i] ~ dlnorm(Y[i],tau_R)
    }
    
    for (i in 2:N){
      alpha[i] <- alpha[i-1] + w[i]
      w[i]~ dnorm(0,tauw)
    }
    
    Rec_track_BY1<- Rep[N-track-3]  
    Rec_track_BY2<- Rep[N-track-2]  
    Rec_track_BY3<- Rep[N-track-1]  
    Rec_track_BY4<- Rep[N-track]
    
    alpha[1]~ dnorm(p.alpha,tau_alpha)
    beta ~ dlnorm(p.beta,tau_beta)        
    tau_R~ dgamma(0.01,0.001)           
    tauv ~ dgamma(0.01,0.001)
    varv<- 1/tauv
    tauw~ dgamma(0.01,0.001)
    varw<- 1/tauw
    sigma <- 1/sqrt(tau_R) 
}


# List of nodes that may be tracked under various commands
Kalman.model.nodelist <- list(All=c("alpha","beta","sigma", "Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","Rep","deviance"),
                              RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta","sigma","deviance"),
                              RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","deviance"))


# R model form - single inputs

  # LEAVING THIS OUT FOR NOW


# R model form - array inputs ****************************************************************************************

  #  LEAVING THIS OUT FOR NOW


# Inits values

Kalman.model.inits<-list(list(tau_R=3, tauw=2, tauv=2),list(tau_R=7, tauw=2, tauv=3))


# List of parameters and priors
Kalman.par.list<-c("S","R_Obs","N","track")
Kalman.prior.list<-c("p.alpha","tau_alpha","p.beta","tau_beta")



# Data Object Check

Kalman.data.check <- function(x,par.list,prior.list){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  chk.list <- c(par.list,prior.list) 
  print("Required Inputs:"); print(chk.list)
  chk <- !(chk.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",chk.list[chk]))}

}



# merge all the pieces for the Kalman Model

KalmanModel.List <- list(MCMC.Model=Kalman.model.MCMC, nodelist=Kalman.model.nodelist, R.Model=NA, R.Model.array=NA, inits=Kalman.model.inits,
                         Data.Check=Kalman.data.check,parlist = Kalman.par.list, priorlist = Kalman.prior.list)



# ---------------------------------------------- END KALMAN FILTER MODEL --------------------------------------------#


#### ======================================== SMOLT-JACK MODEL =======================================================####


SmoltJack.model.MCMC <- function(){ 
    # age3 datum used to update priors for parameters survival and Page in age 4 forecast yr 
    N_2age3 ~ dpois(m_2Page3)
    m_2Page3 <- juv2 * survival * Page3  
    
    # age3 datum used to update priors for parameters survival and Page in age 5 forecast yr 
    N_1age3 ~ dpois(m_1Page3)
    m_1Page3 <- survival * juv1 * Page3
    
    # priors for survival, Page3 and Page4
    survival  ~ dbeta(a1,b1)        #juvenile survival
    Page3 ~ dbeta(a2,b2)            #prortion age 3 (jacking rate) 
    Page4 ~ dbeta(a3,b3)            #proportion age 4
    
    # compute age 4&5 recruitment forecast distributions in forecast year
    N_2age45 ~ dpois(m_2Page45)  #forecast year recruits
    m_2Page45 <- survival * juv2 * (1 - Page3)
    
    # compute age 4 recruitment forecast distributions in forecast year
    N_2age4 ~ dpois(m_2Page4)                       #forecast year recruits
    m_2Page4 <- survival * juv2 * (1 - Page3) * Page4           #forecast age 4 
    
    # compute age 4&5 recruitment forecast distributions in forecast year - 1 to get age 5
    N_1age45 ~ dpois(m_1Page45)         #forecast year recruits
    m_1Page45 <- survival * juv1 * (1 - Page3)
    
    # compute age 5 recruitment forecast distributions in forecast year
    N_1age5 ~ dpois(m_1Page5)  #forecast year recruits
    m_1Page5 <- survival * juv1 * (1 - Page3) * (1-Page4)     #forecast age 5 
    
  
  Rall <- N_2age4 +  N_1age5
  
}


# List of nodes that may be tracked under various commands
SmoltJack.model.nodelist <- list(All=c("Rall","N_2age4","N_1age5", "survival", "Page3", "Page4"),
                              RecPar=c("Rall","N_2age4","N_1age5", "survival", "Page3", "Page4"),
                              RecOnly=c("Rall","N_2age4","N_1age5"))


# R model form - single inputs

# LEAVING THIS OUT FOR NOW


# R model form - array inputs ****************************************************************************************

#  LEAVING THIS OUT FOR NOW


# Inits values **************************************************************************** MADE THESE UP. ORIGINAL  MODEL HAS NO INITS

SmoltJack.model.inits<- list(list(survival=0.5,Page3=0.5,Page4=0.5),list(survival=0.1,Page3=0.1,Page4=0.1))


SmoltJack.data.check <- function(x){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  
  par.list<-c("juv1","juv2","N_2age3","a1","a2","a3","b1","b2","b3")
  chk <- !(par.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",par.list[chk]))}
  
}

# merge 4 pieces for the SmoltJack Model

SmoltJack.List <- list(MCMC.Model=SmoltJack.model.MCMC, nodelist=Kalman.model.nodelist, R.Model=NA, R.Model.array=NA, inits=SmoltJack.model.inits,
                         Data.Check=SmoltJack.data.check)

# -------------------------------------------------- END SMOLT JACK MODEL ------------------------------------------------------------ #


#### ========================================== POWER_JUV COVARIATE MODEL ================================================== ####

# MCMC model form 

PowerCov.model.MCMC <- function(){
  
  for (i in 1:N) {   
    R_Obs[i] ~ dlnorm(R_log[i],tau_R)                   #likelihood 
    R_log[i] <- alpha - beta * log(S[i]) + g * env[i]    # juv model 
    Rep[i] ~ dlnorm(R_log[i],tau_R)                #calculation of replicate data set
  }
  
  Rec_track_BY1<- Rep[N-track-3]  
  Rec_track_BY2<- Rep[N-track-2]  
  Rec_track_BY3<- Rep[N-track-1]  
  Rec_track_BY4<- Rep[N-track]
  
  g ~ dnorm(p.g,tau_g)     
  alpha ~ dnorm(p.alpha,tau_alpha)               #prior for alpha
  beta ~ dnorm(p.beta,tau_beta)                  #prior for beta
  tau_R ~ dgamma(0.001,0.001)                      #prior for precision parameter
  sigma <- 1/sqrt(tau_R)

}


# List of nodes that may be tracked under various commands
PowerCov.model.nodelist <- list(All=c("alpha","beta","g","sigma", "Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","Rep","deviance"),
                              RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","alpha","beta","g","sigma","deviance"),
                              RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","deviance"))


# Do we need one form for single par and another to work with multi-stock arrays and all MCMC outputs in one go?
# OR can we use this with apply()?
# R model form - single inputs

PowerCov.model.R <- function(spn,alpha,beta,g,env,sigma,ln.alpha=FALSE,predicted=FALSE){
  if(ln.alpha){alpha <- ln(alpha)}
  rs <- exp(alpha-beta*log(spn)+g*env)
  if(predicted){rs <- rs + exp(rnorm(1,rec,sigma))} # need to test this step
  rec <- rs
  rec
}


# R model form - array inputs

PowerCov.model.R.array  <- function(Spn,alpha.pars,beta.pars,g.pars,env,out="quants"){  # ************** NOT COMPLETED
  # Spn is an array of dim: years or intervals by pop
  # alpha.pars,beta.pars are arrays of dim: pops by par.est (bootstrap or MCMC estimates)
  # output -> if out = quants then for each Spn level, calculate R based on each alternative par set, then take quantiles of R
  #           if out = samples then just output fitted R
  # env will likely be just one vector of values across years BUT g.pars will be array because g parameters will vary by stock and iteration
  
  vals <- dimnames(Spn)[[1]]  # years or index of Spn interval 
  pops <- dimnames(alpha.pars)[[1]]
  reps <- dimnames(alpha.pars)[[2]]
  
  spn.arr <- array(NA,dim= c(length(vals),length(pops),length(reps)),dimnames= list(vals,pops,reps ))
  rec.arr<- alpha.arr <- beta.arr <- g.arr <- env.arr <- spn.arr
  
  
  # NEED TO FIX THIS BECAUSE SLOW!!!!!!!!!                                                        
  #populate each slice of the array NOTE: should be able to this in one step, but keeps having bug, so do brute force way for now  
  for(repl in reps){
    for(pop in pops){               # need this second loop in case the order of pop is not the same!!!  
      
      #print(repl);print(pop)
      #print(Spn[1:2,1:2])
      spn.arr[,pop,repl] <- as.matrix(Spn[,pop])
      env.arr[,pop,repl] <- env
    }
  } 
  
  for(val in vals){
    alpha.arr[val,,] <- as.matrix(alpha.pars) 
    beta.arr[val,,] <- as.matrix(beta.pars) 
    g.arr[val,,] <- as.matrix(g.pars)
  } 
  
  # apply ricker calculation
  
  rec.arr <- exp(log(alpha.arr)-beta.arr*log(spn.arr)+g.arr*env.arr)
  
  # INCLUDE OPTION TO ADD RANDOM ERROR BASED ON SIGMA PAR
  # BUILD IN MULTIPLE SIGMA SAMPLES BY FEEDING IN REPLICATES OF PAR DATA SET
  
  # calculate percentiles across repl values
  if(out=="quants"){out.arr <- apply(rec.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)}
  if(out=="samples"){out.arr=rec.arr}
  
  out.arr
  
} # end PowerCov.model.R.array


# Inits values
PowerCov.model.inits<-list(list(tau_R=3, beta=0.1),list(tau_R=7,beta=0.2))

# List of parameters and priors
 PowerCov.par.list<-c("S","R_Obs","N","env","track")
 PowerCov.prior.list<- c("p.alpha","tau_alpha","p.beta","tau_beta","p.g","tau_g" )


# Data Object Check
PowerCov.data.check <- function(x,par.list,prior.list){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  chk.list<-c(par.list, prior.list)
  print("Required Inputs:");print(chk.list)
  chk <- !(chk.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",par.list[chk]))}
  
}

# merge 4 pieces for the Ricker Model
PowerCovModel.List <- list(MCMC.Model=PowerCov.model.MCMC, nodelist=PowerCov.model.nodelist, R.Model=PowerCov.model.R, R.Model.array=PowerCov.model.R.array, 
						inits=PowerCov.model.inits, Data.Check=PowerCov.data.check, parlist=PowerCov.par.list,priorlist=PowerCov.prior.list)


#### ============================================== SIBLING MODEL ================================================ ####

### THIS SECTION CONTAINS 2 VERSIONS OF THE SIBLING MODEL
### OPTION 1 HAS BEEN MODIFIED TO REFLECT THE NAMING CONVENTIONS IN THE OTHER MODELS
### OPTION 2 IS THE ORIGINAL VERSION THAT RETAINS THE NAMING CONVENTION OF THE OLD MCMC CODE FOR THIS MODEL
# -> NEED TO DISCUSS
### NOTE OPTION 2 HAS NOT BEEN UPDATED TO WORK WITH THE NEW FC PACKAGE!!!!


########################
### SIBLING MODEL OPTION 1 -> WARNING: ASSUMES THAT SIBLING MODEL IS JUST A POWER MODEL WITH DIFFERENT PRIORS AND DIFFERENT INPUTS


if(TRUE){  # 

# MCMC model form 

# Use the pieces of the Power model previously read in
SiblingModel.List <- list(MCMC.Model=Power.model.MCMC, nodelist=Power.model.nodelist, R.Model=Power.model.R, 
						R.Model.array=Power.model.R.array, inits=Power.model.inits, Data.Check=Power.data.check,
							parlist = Power.par.list, priorlist = Power.prior.list)


} # END NEW VERSION OF SIBLING MODEL





########################
### SIBLING MODEL OPTION 2 -WARNING: NOT YET UPDATED

if(FALSE){  

# MCMC model form 

Sibling.model.MCMC <- function(){
  
  for(i in 1:N){
    ln_y[i] ~ dnorm(ln_pred_R[i], tau)  # likelihood function for L(data | parameters)
    ln_pred_R[i] <- intercept + slope * ln_x[i]
    Rep[i] ~ dnorm(ln_pred_R[i], tau)
  }
    
  #Priors
  slope ~ dnorm(0, 0.00000000001)
  intercept ~ dnorm(0, 0.00000000001)
  tau ~ dgamma(0.001, 0.001)
  sigma <- 1/sqrt(tau)
  
  Rec_track_BY4 <- exp(Rep[N-track])
  Rec_track_BY3 <- exp(Rep[N-track-1])
  Rec_track_BY2 <- exp(Rep[N-track-2])
  Rec_track_BY1 <- exp(Rep[N-track-3])
}


# List of nodes that may be tracked under various commands
Sibling.model.nodelist <- list(All=c("slope","intercept","sigma","Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","Rep"),
                                RecPar=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4","slope","intercept","sigma"),
                                RecOnly=c("Rec_track_BY1","Rec_track_BY2","Rec_track_BY3","Rec_track_BY4"))


# Do we need one form for single par and another to work with multi-stock arrays and all MCMC outputs in one go?
# OR can we use this with apply()?
# R model form - single inputs

Sibling.model.R <- function(ln_x,intercept,slope,sigma,predicted=FALSE){
  rec4 <- exp(alpha+beta*ln_x)
  if(predicted){rec4 <- rec4 + exp(rnorm(1,rec4,sigma))} # need to test this step
  rec <- rec4
  rec
}


# R model form - array inputs

Sibling.model.R.array  <- function(ln_x,intercept.pars,slope.pars,out="quants"){  # ************** NOT COMPLETED
  # Spn is an array of dim: years or intervals by pop
  # alpha.pars,beta.pars are arrays of dim: pops by par.est (bootstrap or MCMC estimates)
  # output -> if out = quants then for each Spn level, calculate R based on each alternative par set, then take quantiles of R
  #           if out = samples then just output fitted R
  # env will likely be just one vector of values across years BUT g.pars will be array because g parameters will vary by stock and iteration
  
  vals <- dimnames(ln_x)[[1]]  # years or index of Spn interval 
  pops <- dimnames(intercept.pars)[[1]]
  reps <- dimnames(intercept.pars)[[2]]
  
  ln_x.arr <- array(NA,dim= c(length(vals),length(pops),length(reps)),dimnames= list(vals,pops,reps ))
  rec4.arr<- intercept.arr <- slope.arr <- ln_x.arr
  
  
  # NEED TO FIX THIS BECAUSE SLOW!!!!!!!!!                                                        
  #populate each slice of the array NOTE: should be able to this in one step, but keeps having bug, so do brute force way for now  
  for(repl in reps){
    for(pop in pops){               # need this second loop in case the order of pop is not the same!!!  
      
      #print(repl);print(pop)
      #print(Spn[1:2,1:2])
      ln_x.arr[,pop,repl] <- as.matrix(ln_x[,pop])
    }
  } 
  
  for(val in vals){
    intercept.arr[val,,] <- as.matrix(intercept.pars) 
    slope.arr[val,,] <- as.matrix(slope.pars) 
  } 
  
  # apply ricker calculation
  
  rec4.arr <- exp(intercept.arr+slope.arr*ln_x.arr)
  
  # INCLUDE OPTION TO ADD RANDOM ERROR BASED ON SIGMA PAR
  # BUILD IN MULTIPLE SIGMA SAMPLES BY FEEDING IN REPLICATES OF PAR DATA SET
  
  # calculate percentiles across repl values
  if(out=="quants"){out.arr <- apply(rec4.arr,MARGIN=c(1,2),FUN=quantile,probs=seq(0.05,0.95,by=0.05),na.rm=TRUE)}
  if(out=="samples"){out.arr=rec4.arr}
  
  out.arr
  
} # end Sibling.model.R.array


# Inits values
Sibling.model.inits<-list(list(tau=100),list(tau=3))


# Data Object Check
Sibling.data.check <- function(x){
  #x is a list object with data 
  #for now this checks whether the required pieces are there
  #could be expanded to do quality check (e.g. vector lengths vs. N, juv data)
  
  par.list<-c("ln_x","ln_y","N")
  chk <- !(par.list %in% names(x)) # check which are not incl
  if(sum(chk)==0){print("All inputs available")}
  if(sum(chk)>0){print(paste("The following inputs are missing:",par.list[chk]))}
  
}

# merge 4 pieces for the Ricker Model
SiblingModel.List <- list(MCMC.Model=Sibling.model.MCMC, nodelist=Sibling.model.nodelist, R.Model=Sibling.model.R, R.Model.array=Sibling.model.R.array, inits=Sibling.model.inits, Data.Check=Sibling.data.check)

} # END OLD VERSION OF SIBLING MODEL




##### ********************************************  MERGING INTO SINGLE OBJECT ************************************* ####

# merge all
SR.Models.List <- list(Ricker.fp=RickerfpModel.List, Ricker=RickerModel.List, RickerCov=RickerCovModel.List, RickerCycle=RickerCycleModel.List,
                       Power=PowerModel.List, Larkin = LarkinModel.List, Kalman=KalmanModel.List, SmoltJack=SmoltJack.List, PowerCov=PowerCovModel.List,
                       Sibling=SiblingModel.List) # others to be added later







 
