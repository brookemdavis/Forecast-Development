library(tidyverse)
library(TMB)
source("Code/Functions.R")

# read in Sockeye data
SockDat <- read.csv("FraserSK/DataIn/Sockdat May 9 2016 including tlSPNplus2011rec(FRSSI).csv")
# read in Smax data to get CU/MU names
SmaxDat <- read.csv("FraserSK/DataIn/Smax Data.csv")
SockDat <- left_join(SockDat, SmaxDat, by="stk")

# Flag cyclic stocks
Cyclic_Stocks <- c("E. Stuart", "L. Stuart", "Quesnel", "L. Shuswap", "Seymour")
Cyclic_Data <- SockDat %>% filter(Stock_Name %in% Cyclic_Stocks) %>%
                  filter(is.na(ets) == F) %>% #& is.na(rec)==F) %>%
                  select(stk, Stock_Name, yr, rec, ets) %>% 
                  group_by(Stock_Name) %>% 
                  mutate( Cycle = 2015 - ((2015-yr) %% 4))
SockDat <- SockDat %>% left_join(Cyclic_Data) %>% select(-(cyc:rec5), -(TotalAdultSpawner:X), -(Start.Year:N))

# want to test all model forms to see how parameter estimates vary

# fit
#1) Base R model using lm
#2) Basic Ricker TMB
#3) Hier Ricker TMB
#4) Basic Ricker with LRP estimattion
#5) Hier Ricker with LRP estimation

# Cultus does not have ETS due to being lake spawners (can't estimate effective spawners)
# For CU assessment only Seymour used for SR modelling (Scotch shorter TS, and other issues)
# but if start excluding pops like this, will aggregate BM be as meaningful? SHould be flagged
# as main issue, becuase if agg. LRP's are going to be used "down river" they need to represent
# all abundance -- but is this even possible?
# also using ETS is a problem, since would not be counting ETS down river, would be counting
# total adults

# E. Stuart (which doesn't need agg. LRP), Quesnel, L. Stuart, Scotch-Seymour, 
# L. Shu determined cyclic in WSP (others show 
# some cyclic patterns too), which means basic Ricker might not be entirely suitable

# Besides Scotch, still getting ok SR relationships (excecpt maybe Harrison?) so will use
# for now 

# Hier model gives same alphas across stocks for S (4 CUs, including L Stu and Quesnel) and ES (7 CUs including seymour)
 # Maybe this isn't the model "not working" maybe this just says that alpha is so uncertain, that we get 
 # better fit overall if we just assume all same?

Stock_Tab <- unique(SockDat[, c("MU", "Stock_Name")]) %>% arrange(MU) 
MUs <- c("S", "ES", "L")

# There aren't isn't an MU that doesn't include a cyclic stock, so there isn't a 
# case study where we will be able to ignore without dealing with cyclic

# for now, will exclude scotch (since can't fit Ricker to, even base R gives negative B)

# Later will want to encorporate Larkin BMs?


# Basic Ricker model for all Fraser Sockeye stocks
# for cultus use TotalAdultSpawner since ets not available
SockDat$ets[SockDat$Stock_Name == "Cultus"] <- SockDat$TotalAdultSpawner[SockDat$Stock_Name == "Cultus"]
# remove one year with NA ets and rec, remove Scotch
SockDat <- SockDat %>% filter(is.na(ets) == F & is.na(rec) == F) %>% filter(Stock_Name != "Scotch")
SockDat$stk_num <-   group_indices(SockDat, stk) -1

# load all TMB models that will be used
#dyn.unload(dynlib("Code/Hier_Ricker"))
#compile("Code/Hier_Ricker.cpp")
dyn.load(dynlib("Code/Basic_Ricker"))
dyn.load(dynlib("Code/Hier_Ricker"))
dyn.load(dynlib("Code/Aggregate_LRPs_Hier"))
dyn.load(dynlib("Code/Aggregate_LRPs"))

# use TMB_Ricker and Run_Ricker_LRP model

Scale <- 1

All_Ests <- data.frame(Estimate = numeric(), Param = character(), Mod = character(),
                       stk = numeric(), Stock_Name = character(), MY = character())
Logistic_Data <-   data.frame( Mod = character(), MU=character(), Year = numeric(), 
                               N_CUs = numeric(), Agg_Abund = numeric())


for(mm in 1:length(MUs)){
  Data <- SockDat %>% filter(MU == MUs[mm])
  #Data <- Data %>% filter(!Stock_Name == "Seymour")
  Data$stk_num <- group_indices(Data, stk) -1
  Data$yr_num <- group_indices(Data, yr) - 1
  
  # Fit Basic and Hier Ricker
  All_Ests_Basic <- TMB_Ricker(Data, Mod = "Basic_Ricker",  Scale = 1)
  All_Ests_Hier <- TMB_Ricker(Data, Mod = "Hier_Ricker",  Scale = 1)
  
  # Fit Mods with LRPs
  Basic_LRP <- Run_Ricker_LRP(Data, Mod = "Aggregate_LRPs", Scale, Sim = F)
  Hier_LRP <- Run_Ricker_LRP(Data, Mod = "Aggregate_LRPs_Hier", Scale, Sim = F)
  All_Ests_Basic_LRP <- Basic_LRP[[1]]
  All_Ests_Hier_LRP <- Hier_LRP[[1]]
  
  # also pull out Logistic model data
  Logistic_Data <- bind_rows( Logistic_Data,  Basic_LRP[[2]], Hier_LRP[[2]]) 
  
  # now Hier without cyclic stocks
  Data <- Data %>% filter(is.na(Cycle))
  Data$stk_num <- group_indices(Data, stk) -1
  Data$yr_num <- group_indices(Data, yr) - 1
  
  # Hier model without cyclic stocks
  Hier_Ests2 <- TMB_Ricker(Data, Mod = "Hier_Ricker",  Scale = 1)
  Hier_Ests2$Mod <- "Hier_Non_Cyclic"
  
  All_Ests <- rbind(All_Ests, All_Ests_Basic, All_Ests_Hier, All_Ests_Basic_LRP, All_Ests_Hier_LRP, Hier_Ests2)
  
  # next fit basic Ricker using Base R
  Stks <- unique(Data$Stock_Name)
  Ests <- list()
  for(i in 1:length(Stks)){
    Ests[[i]] <- GetRickerParams(Stk = Stks[i], Data = SockDat)
    New_Rows <- data.frame(Estimate = unlist(Ests[[i]]), 
                           Param = names(unlist(Ests[i])), 
                           Mod = "R_Basic",
                           Stock_Name = Stks[i])
    # add stk number
    New_Rows <- left_join(New_Rows, unique(SockDat[, c("stk", "Stock_Name", "MU")]))
    All_Ests <- union(All_Ests, New_Rows)
  }
  
} # end MU loop

# also fit Basic Ricker for Estu -- not sure why TMB won't work
Ests<- GetRickerParams(Stk = "E. Stuart", Data = SockDat)
New_Rows <- data.frame(Estimate = unlist(Ests), 
                       Param = names(unlist(Ests)), 
                       Mod = "Basic_Ricker",
                       stk = NA,
                       Stock_Name = "E. Stuart",
                       MU = "EStu")
All_Ests <- rbind(All_Ests, New_Rows)
write.csv(All_Ests, "FraserSK/DataOut/All_Ests.csv", row.names = F)
write.csv(Logistic_Data,"FraserSK/DataOut/Logistic_Data.csv", row.names = F )
  

# For each stock plot S_R data and overlay curve, Sgen, and Smax 

# Compare basic R and TMB
Plot_2_SR_Curves(All_Ests) # all look theh same -- doesn't work

# Now compare Basic TMB mod with one where we est LRP
Plot_2_SR_Curves(All_Ests, Mod1 = "Aggregate_LRPs", Lab1 = "TMB Basic w/LRP",
              Mod2 = "Basic_Ricker", Lab2 = "TMB Basic Ricker" ,
              name = "Comp_Basic_LRP")
# no difference

# Same with hier model
Plot_2_SR_Curves(All_Ests, Mod1 = "Aggregate_LRPs_Hier", Lab1 = "TMB Hier w/LRP",
            Mod2 = "Hier_Ricker", Lab2 = "TMB Hier Ricker" ,
            name = "Comp_Hier_LRP")
# small differences, but not significant

# now basic with hier (no LRP)
Plot_2_SR_Curves(All_Ests, Mod1 = "Basic_Ricker", Lab1 = "TMB Ricker",
            Mod2 = "Hier_Ricker", Lab2 = "TMB Hier Ricker" ,
            name = "Comp_Hier_Basic_No_LRP")

# Now two LRP mods
Plot_2_SR_Curves(All_Ests, Mod1 = "Aggregate_LRPs", Lab1 = "TMB Basic w/LRP",
            Mod2 = "Aggregate_LRPs_Hier", Lab2 = "TMB Hier w/LRP" ,
            name = "Comp_Hier_Basic_With_LRP")

# Compare non-cyclic hier
Plot_2_SR_Curves(All_Ests, Mod1 = "Hier_Ricker", Lab1 = "TMB Hier",
                 Mod2 = "Hier_Non_Cyclic", Lab2 = "TMB Hier NC" ,
                 name = "Comp_Hier_Cyclic")

# now look at alpha shrinkage in hier vs non - hier

# for each MU make package of plots with SR relationships, Alphas,
# agg TS with BM's on top

MU_Summary_Plot(Data, MUs, Mod1 = "Aggregate_LRPs",  Mod2 = "Aggregate_LRPs_Hier", Lab1 = "Basic", Lab2 = "Hier.")

# all line up except sigma -- TMB consistently slightly smaller, is this a difference
# in waht each of these estimates mean? That's becasue it is residual standard error, 
# which is SSE/n-(1+k) instead of SSE/n (I think?)
# or is it because we are using lognormal and need correction factor?
# current model gives back simualted alpha, so those align

# these sigma values are much small than the ones I've been using!!

# These do not scale with S-R magnitude -- get same values if fit with non-scaled vals
# should increase ricker Sigma if want to be realistic! 0.73-1.34 (mean 0.98)

# now if we fit hierarchical models for the MUs with multiple stocks, 
#how different are our alpha estimates?

# Will have three examples, early summer, Late, Summer. Summer has same num of stocks
# over time, other two have variable TS length

# LOoks like my alpha values were low too, real vals 3.3-11 (mean = 6.6)

# not sure if should use cultus since only have TotalAdult Spawner, not ets
# -- also has lowest alpha, maybe should exclude, will still have 4 other pops
# apparently scotch and seymour are one CU, but from different MUs -- keep separated I guess?


# for later -- add trigger in TMB script so that there is a trigger if a BM is going to be manually input
# (maybe for a stock like Scotch, where SR curve can't be fit) or should use 25th percentile instead
# I think maybe Scotch & Seymour are supposed to be blended together anyways


#==============================================================================
# Larkin for cyclic stocks
# Want to estimate Sgen for both versions

# Can't get Sgen estimation to work -- plus want to flexibility on timeframe we are using

# Take Larkin parameters and estimate cycle-specific BMs using different subsets of data (full, 
  #last 5 gens (20 years), last 3 gens (12 years))
# then for each MU with a cyclic stock calculate both cycle-specific BMs, and just a single BM, and compare them
# look at what proportion of total escapement, on average, come from cyclic stock 
## may want to come up with guidline -- if xx% or more of total return come from cyclic stocks, use cycle-specific BMs
# do this using my estimates and Carrie's estimates

# also do same thing compare using hier and basic on non-cyclic
      # go back and fit hier models with only non-cyclic stocks (remove Seymour from ES, L Shu from late, L stu, Quesnel from S)
      # note that Estu is it's own MU -- shoudl automatically use cycle-specific BMs?

# Compare to not identifying cyclic stocks at all, as well

# first get table of Sgens

# Compile set of Sgens for non-cyclic stocks
All_Ests <- read.csv("FraserSK/DataOut/All_Ests.csv")
Logistic_Data <- read.csv("FraserSK/DataOut/Logistic_Data.csv")
NC_Sgens <- All_Ests %>% filter(Mod %in% c("Basic_Ricker", "Hier_Non_Cyclic"), Param == "Sgen") %>% 
               rename(Sgen = Estimate) %>% select(-Param)
# ADD Smsy
NC_Sgens <- All_Ests %>% filter(Mod %in% c("Basic_Ricker", "Hier_Non_Cyclic"), Param == "SMSY") %>% 
          rename(SMSY = Estimate) %>% select(-Param) %>% right_join(NC_Sgens)

# Now get cyclic Sgens
C_Sgens <- read.csv( "FraserSK/DataOut/Cyclic_Sgens.csv")

pdf("FraserSk/Figures/Sgens+SMSY.pdf")
par(mfrow= c(3,2), oma = c(2,3,2,2), mar = c(2,1,2,1))
for(ss in 1:length(Cyclic_Stocks)){
  Data <- SockDat %>% filter(Stock_Name == Cyclic_Stocks[ss])
  Sgens <- C_Sgens %>% filter(Stock_Name == as.character(Cyclic_Stocks[ss]))
  ymax = max(Data$ets[is.na(Data$ets) == F], Sgens$SMSY)
  plot(Data$yr, Data$ets, type = "o", lwd=1.5, pch=19, cex = 0.6, ann=F, ylim = c(0,ymax))
  lines(Sgens$Year, Sgens$Sgen, col = "red", type = "o", pch = 19, cex = 0.6)
  lines(Sgens$Year, Sgens$SMSY, col = "red", lty = 2)
  # add straight horizontal line for regular Ricker
  Ricker_Sgen <- NC_Sgens %>% filter(Stock_Name == as.character(Cyclic_Stocks[ss]))
  abline(h = Ricker_Sgen$Sgen, col = "blue", lwd = 1.2)
  abline(h = Ricker_Sgen$SMSY, col = "blue", lty=2)
  mtext(side = 3, text = Cyclic_Stocks[ss])
}
mtext(side = 1, text = "Year", outer = T)
mtext(side =2, text = "ETS (millions)", outer=T)
legend("topleft", col = c("black", "red", "blue", "black", "black"), 
       legend = c("ETS", "Larkin", "Ricker", "Sgen", "SMSY"), lty=c(1,1,1,1,2), bty="n")
dev.off()

# FOr Cyclic stocks, for each cycle line, Calculate BM
C_Sgens <- C_Sgens %>% mutate( Cycle = 2015 - ((2015-Year) %% 4)) %>%
           group_by(Stock_Name, Cycle, MU) %>% summarise(Sgen = median(Sgen)) %>%
           mutate(Mod = "Larkin")
# now join with non-cyclic
Sgens <- full_join(NC_Sgens, C_Sgens)

# now fit logistic  model to get agg. BM


# now compare with one gotten from assuming they were all Ricker/ All Hier

# Aggregate BMs to compare: All Ricker*, All Hier*, Ricker/Larkin, Hier Ricker/Larkin (*already have these)

# Add Sgen Column with Ricker/Larkin BMs
# for cyclic stocks use Larkin aligned with cycle  lines
# for non-cyclic use Basic_Ricker

Sgens_To_Join <- Sgens %>% filter(Mod %in% c("Basic_Ricker", "Larkin")) %>%
                           filter(!(Stock_Name %in% Cyclic_Stocks & Mod == "Basic_Ricker")) %>%
                           select(-stk)

LRP_Summ_Basic <- SockDat %>% filter(Stock_Name != "E. Stuart", Stock_Name != "Scotch", is.na(ets) == F) %>% 
                          left_join(Sgens_To_Join, by = c("Stock_Name", "Cycle", "MU")) %>% 
                          group_by(yr, MU) %>% summarise(Total_ETS = sum(ets),
                                                         CUs_Above_LRP = sum(ets > Sgen),
                                                         n_CUs = n(),
                                                         Prop_Above_LRP = sum(ets > Sgen)/n() ) %>%
                          group_by( MU) %>%  filter(n_CUs == max(n_CUs)) %>%
                          mutate( Cycle = 2015 - ((2015-yr) %% 4)) # need to add cycle back

# Same but hier model
Sgens_Hier <- Sgens %>% filter(Mod %in% c("Hier_Non_Cyclic", "Larkin")) %>%
  filter(!(Stock_Name %in% Cyclic_Stocks & Mod == "Basic_Ricker")) %>%
  select(-stk)

LRP_Summ_Hier <- SockDat %>% filter(Stock_Name != "E. Stuart",Stock_Name != "Scotch",  is.na(ets) == F) %>% 
  left_join(Sgens_Hier, by = c("Stock_Name", "Cycle", "MU")) %>% 
  group_by(yr, MU) %>% summarise(Total_ETS = sum(ets),
                                 CUs_Above_LRP = sum(ets > Sgen),
                                 n_CUs = n(),
                                 Prop_Above_LRP = sum(ets > Sgen)/n() ) %>%
  group_by(yr, MU) %>%  filter(n_CUs == max(n_CUs)) %>%
  mutate( Cycle = 2015 - ((2015-yr) %% 4)) # need to add cycle back

# Function to fit logistic model
Fit_Logistic <- function(Data, MU, byCycle = F, Cycle){
 LRP_Summ <- Data[Data$MU == MU, ]
 if(byCycle == T){
   LRP_Summ <- LRP_Summ[Data$Cycle == Cycle, ]
 }
 Prop_Mod = glm(cbind(CUs_Above_LRP, n_CUs - CUs_Above_LRP) ~ Total_ETS, family = binomial, data = LRP_Summ)
 LRP <- (log(0.95/(1-0.95)) - Prop_Mod$coefficients[1])/ Prop_Mod$coefficients[2]
 out <- list()
 out$Mod <- Prop_Mod
 out$LRP <- LRP
 out
}

MUs <- unique(LRP_Summ_Basic$MU)

# Now fit logistic model and compile results ----------------------------
#LRPs <- map_dbl( MUs, .f = Fit_Logistic, Data = LRP_Summ_Basic)
Logistic_Coeffs <- data.frame(MU = character(), Mod = character(), Cycle = numeric(), B_0 = numeric(), B_1 = numeric())
LRP_DF <- data.frame(MU = character(), LRP = numeric(), Cycle = numeric(), Mod = character())
Cycles <- unique(SockDat$Cycle[is.na(SockDat$Cycle) == F])

for(mm in 1:length(MUs)){
  out <- Fit_Logistic( Data = LRP_Summ_Basic, MU = MUs[mm])
  new.rows.coeffs <- data.frame( MU = MUs[mm], Mod =  "Basic_Larkin", Cycle = NA,
                                B_0 = as.numeric(out$Mod$coefficients[1]) ,
                                B_1 = as.numeric(out$Mod$coefficients[2]) )
  new.rows.LRP <- data.frame(MU = MUs[mm], LRP = as.numeric(out$LRP), Cycle = NA, Mod = "Basic_Larkin")
  LRP_DF <- rbind(LRP_DF, new.rows.LRP)
  Logistic_Coeffs <- rbind(Logistic_Coeffs, new.rows.coeffs)
  # Now fit cyclic
  LRPs_Cyclic <- NULL
  for(cc in 1:length(Cycles)){
      out <- Fit_Logistic( Data = LRP_Summ_Basic, MU = MUs[mm], byCycle = T, Cycle = Cycles[cc])
      LRPs_Cyclic <- c(LRPs_Cyclic, out$LRP)
      new.rows.coeffs <- data.frame(MU = MUs[mm], Mod = "Basic_Larkin", Cycle = Cycles[cc], 
                                    B_0 = as.numeric(out$Mod$coefficients[1]) ,
                                    B_1 = as.numeric(out$Mod$coefficients[2]) )
      Logistic_Coeffs <- rbind(Logistic_Coeffs, new.rows.coeffs)
  } # end cycle loop
  new.rows <- data.frame(MU = rep(MUs[mm], length(Cycles)), LRP = as.numeric(LRPs_Cyclic), Cycle = Cycles, Mod = "Basic_Larkin")
  LRP_DF <- rbind(LRP_DF, new.rows)
} # end MU loop

# same for hierarchical ---------------
for(mm in 1:length(MUs)){
  out <- Fit_Logistic( Data = LRP_Summ_Basic, MU = MUs[mm])
  new.rows.coeffs <- data.frame( MU = MUs[mm], Mod =  "Hier_Larkin", Cycle = NA,
                                 B_0 = as.numeric(out$Mod$coefficients[1]) ,
                                 B_1 = as.numeric(out$Mod$coefficients[2]) )
  new.rows.LRP <- data.frame(MU = MUs[mm], LRP = as.numeric(out$LRP), Cycle = NA, Mod = "Hier_Larkin")
  LRP_DF <- rbind(LRP_DF, new.rows.LRP)
  Logistic_Coeffs <- rbind(Logistic_Coeffs, new.rows.coeffs)
  # Now fit cyclic
  LRPs_Cyclic <- NULL
  for(cc in 1:length(Cycles)){
    out <- Fit_Logistic( Data = LRP_Summ_Basic, MU = MUs[mm], byCycle = T, Cycle = Cycles[cc])
    LRPs_Cyclic <- c(LRPs_Cyclic, out$LRP)
    new.rows.coeffs <- data.frame( MU = MUs[mm], Mod =  "Hier_Larkin", Cycle = Cycles[cc],
                                   B_0 = as.numeric(out$Mod$coefficients[1]) ,
                                   B_1 = as.numeric(out$Mod$coefficients[2]) )
    Logistic_Coeffs <- rbind(Logistic_Coeffs, new.rows.coeffs)
  } # end cycle loop
  new.rows <- data.frame(MU = rep(MUs[mm], length(Cycles)), LRP = as.numeric(LRPs_Cyclic), Cycle = Cycles, Mod = "Hier_Larkin")
  LRP_DF <- rbind(LRP_DF, new.rows)
} # end MU loop



# Add rows to Logistic_Data
Logistic_Data$Cycle <- rep(NA, dim(Logistic_Data)[1])
Logistic_Data_New <- data.frame(Mod = "Basic_Larkin", MU = LRP_Summ_Basic$MU, 
                                Year = LRP_Summ_Basic$yr, N_CUs = LRP_Summ_Basic$CUs_Above_LRP,
                                Agg_Abund = LRP_Summ_Basic$Total_ETS, Cycle = LRP_Summ_Basic$Cycle)
Logistic_Data <- full_join(Logistic_Data, Logistic_Data_New)


# same for logistic
Logistic_Data_Hier <- data.frame(Mod = "Hier_Larkin", MU = LRP_Summ_Hier$MU, 
                                Year = LRP_Summ_Hier$yr, N_CUs = LRP_Summ_Hier$CUs_Above_LRP,
                                Agg_Abund = LRP_Summ_Hier$Total_ETS, Cycle = LRP_Summ_Hier$Cycle)
Logistic_Data <- full_join(Logistic_Data, Logistic_Data_Hier)

# seems useful to add number of CUs per stock to each one
N_CUs <- LRP_Summ_Basic %>% group_by(MU) %>% summarise(Total_CUs = unique(n_CUs))
Logistic_Data <- Logistic_Data %>% left_join(N_CUs)

# Now add Aggregate_LRPs and Aggregate_LRPs_Hier to LRP DF to get ready to compare the 4
LRP_DF <- All_Ests %>% filter(Param == "Agg_BM", Mod %in% c("Aggregate_LRPs", "Aggregate_LRPs_Hier")) %>%
          mutate(LRP = Estimate) %>%
          select(LRP, Mod, MU) %>%
          full_join(LRP_DF)

# also add logistic coeffs from all_Ests
Logistic_Coeffs <- All_Ests %>% filter(Param %in% c("B_0", "B_1"), Mod %in% c("Aggregate_LRPs", "Aggregate_LRPs_Hier")) %>%
                    pivot_wider(values_from = Estimate, names_from = Param) %>%
                    select(Mod, MU, B_0, B_1) %>%
                    full_join(Logistic_Coeffs)
    
write.csv(Logistic_Data, "FraserSK/DataOut/Logistic_Data_Larkin.csv")         
write.csv(LRP_DF, "FraserSK/DataOut/LRP_DF.csv") 
write.csv(Logistic_Coeffs, "FraserSK/DataOut/Logistic_Coeffs.csv") 

Logistic_Data <- read.csv("FraserSK/DataOut/Logistic_Data_Larkin.csv")
LRP_DF <- read.csv("FraserSK/DataOut/LRP_DF.csv")
Logistic_Coeffs <- read.csv( "FraserSK/DataOut/Logistic_Coeffs.csv") 


# for each MU plot 4 different logistic model fits (on same scale) then plot agg. abund with 4 BMs
# Plot logistic Model fits
Mods_Vec <- unique(Logistic_Data$Mod)
cols <- c("red", "blue", "darkgreen", "purple")

Compare_Larkin <- function(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "S", name = "Larkin_Comp", ByCycle = F){

pdf(paste("FraserSK/Figures/", name, sep=""))
par(mfrow = c(3,2), mar = c(2,2,2,2), oma=c(2,1,2,1))

for(mm in 1:length(Mods_Vec)){

  if(ByCycle == T & Mods_Vec[mm] %in% c("Basic_Larkin", "Hier_Larkin")){
    Cycles <- unique(LRP_DF$Cycle[is.na(LRP_DF$Cycle) ==F])
  } else {
    Cycles <- 1
  }
  
  LRPs <- LRP_DF[LRP_DF$Mod == Mods_Vec[mm] & LRP_DF$MU == MU,]
  
for(cc in 1:length(Cycles)){
  if(ByCycle == T & Mods_Vec[mm] %in% c("Basic_Larkin", "Hier_Larkin")){
    Log1 <- Logistic_Data[Logistic_Data$Mod == Mods_Vec[mm] & Logistic_Data$MU == MU & 
                            is.na(Logistic_Data$Cycle)==F & Logistic_Data$Cycle == Cycles[cc],]
    Title <- paste(Mods_Vec[mm], Cycles[cc], "Cycle Line")
    Ests <- Logistic_Coeffs[Logistic_Coeffs$Mod == Mods_Vec[mm] &  Logistic_Coeffs$MU == MU & 
                              is.na(Logistic_Coeffs$Cycle)==F &
                             Logistic_Coeffs$Cycle == Cycles[cc],]
    LRP <- LRPs$LRP[is.na(LRPs$Cycle)==F & LRPs$Cycle == Cycles[cc]]
  } else {
    Log1 <- Logistic_Data[Logistic_Data$Mod == Mods_Vec[mm] & Logistic_Data$MU == MU,]
    Title <- Mods_Vec[mm]
    Ests <- Logistic_Coeffs[Logistic_Coeffs$Mod == Mods_Vec[mm] &  Logistic_Coeffs$MU == MU & is.na(Logistic_Coeffs$Cycle),]
    LRP <- LRPs$LRP[is.na(LRPs$Cycle)]
  }
 
  #x and y lims
  xmax = max(Log1$Agg_Abund, LRPs$LRP)
  plot(Log1$Agg_Abund, Log1$N_CUs/Log1$Total_CUs, pch=19, col = cols[mm], xlab = "Agg. Abund", ylab = "Proportion of Site > LRP",
       xlim = c(0, xmax), ylim = c(0,1), main = Title)

  # add logistic fits and LRP
  # xx values to predict from
  xx <- seq(0, xmax, by=0.05)

  #Pull out coeffs and sds
  yy1 <- inv_logit( Ests$B_0+ 
                      Ests$B_1*xx)
  lines(xx, yy1, col = cols[mm])

  # Now add LRPs
  abline(v = LRP, lty=2, col = cols[mm])
 } # end cycle loop
} # end mod loop

# Now plot TS with Agg LRP

ymax <- max(Logistic_Data$Agg_Abund[Logistic_Data$MU == MU], LRP_DF$LRP[LRP_DF$MU == MU])*1.25
plot(Log1$Year, Log1$Agg_Abund, type = "l", xlab = "Year", ylab = "Agg. ETS", 
     ylim = c(0, ymax), main = "Aggregate Abundance")
for(mm in 1:length(Mods_Vec)){
  # add LRPs
  
  if(ByCycle == T & Mods_Vec[mm] %in% c("Basic_Larkin", "Hier_Larkin")){
    LRPs <- LRP_DF[LRP_DF$Mod == Mods_Vec[mm] & LRP_DF$MU == MU & is.na(LRP_DF$Cycle)==F,]
    Agg_LRPs <-  Logistic_Data[Logistic_Data$Mod == Mods_Vec[mm] & Logistic_Data$MU == MU & 
                                 is.na(Logistic_Data$Cycle)==F,] %>% select( Year, Agg_Abund, Cycle) %>% full_join(LRPs)
    lines( Agg_LRPs$Year, Agg_LRPs$LRP, col = cols[mm])
  } else {
    LRPs <- LRP_DF[LRP_DF$Mod == Mods_Vec[mm] & LRP_DF$MU == MU,]
    abline(h = LRPs$LRP[is.na(LRPs$Cycle)], col = cols[mm])
  }
} # end mod loop
legend("topleft", col = cols[1:mm], legend = Mods_Vec, lty=1, bty="n", cex=0.6)
dev.off()
} # end function

Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "S", name = "Larkin_Comp_S.pdf")
Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "ES", name = "Larkin_Comp_ES.pdf")
Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "L", name = "Larkin_Comp_L.pdf")

Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "S", name = "Larkin_Comp_S_Cyclic.pdf", ByCycle = T)
Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "ES", name = "Larkin_Comp_ES_Cyclic.pdf", ByCycle = T)
Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "L", name = "Larkin_Comp_L_Cyclic.pdf", ByCycle = T)


# For two Larkin BMs plot fit for each cycle line. Then TS with line-specific BMs (same or diff plots?)

# now will have different fit for each cycle == plot on top of eachother?


Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "S", name = "Larkin_Comp_S.pdf")
Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "ES", name = "Larkin_Comp_ES.pdf")
Compare_Larkin(Logistic_Data, LRP_DF, Logistic_Coeffs, Mods_Vec, MU = "L", name = "Larkin_Comp_L.pdf")

#==========================================================================================================
# Compare My ref pts to LFRP's used in management -- and in years wheere MU is avoce LFRP, what proportion
# of CUs are above Sgen?

FRPs <- read.csv("FraserSK/DataIn/FR_SK_TAM_Rules.csv")

# For each MU, for 2007-2015, how many CUs were above Sgen in years met LFRP

# Use basic ricker (not hier for now) and Larkin

# attach FRPs to LRP summary DF, pick out years where above LFRP
LFRP_DF <- left_join(LRP_Summ_Basic, FRPs, by=c("yr" = "Year", "MU")) %>% 
         filter(is.na(LFRP) == F, Total_ETS*1000000 > LFRP) %>% arrange(MU) %>%
         mutate(ETS = Total_ETS*1000000) %>% select(-Total_ETS, -Cycle, -UFRP)

# same using UFRP
UFRP_DF <- left_join(LRP_Summ_Basic, FRPs, by=c("yr" = "Year", "MU")) %>% 
  filter(is.na(UFRP) == F, Total_ETS*1000000 > UFRP) %>% arrange(MU) %>%
  mutate(ETS = Total_ETS*1000000) %>% select(-Total_ETS, -Cycle, -LFRP)

# for each MU, plot proportion of CUs above Sgen above abundance
#LRPs from LRP_DF
LRPs <- LRP_DF %>% filter(Mod == "Basic_Larkin")

MUs <- unique(LRP_Summ_Basic$MU)
pdf("FraserSK/Figures/Compare_LRPs.pdf", width = 11, height = 9)
par(mfrow = c(2,2), oma=c(2,2,2,2), mar=c(2,3,1,3), mgp=c(2,0.75,0))
for(mm in 1:length(MUs)){
  Dat <- LRP_Summ_Basic %>% filter(MU == MUs[mm])
  LRP <- (LRPs %>% filter(MU == as.character(MUs[mm]) & is.na(Cycle)))$LRP
  ymax <- max(Dat$Total_ETS, LRP*1.25,1)
  plot(Dat$yr, Dat$Total_ETS, type = "o", pch=19, xlab = "Year", ylab = "Agg. ETS (millions)", 
       ylim = c(0, ymax), lwd = 2)
  abline(h = LRP, col = "blue", lwd=1.5)
  # now at LFRP
  LFRP <- FRPs %>% filter(MU == MUs[mm])
  lines(LFRP$Year, LFRP$LFRP/1000000, col = "red", lwd=1.5, type="o", pch=19)
  #plot below showing proportion of stocks above Sgen
  lines(Dat$yr, Dat$Prop_Above_LRP*round(ymax), col = "orange", type = "o", pch=19)
  roundmax <- round(ymax) 
  axis(4, at = c(0, roundmax/4, roundmax/2, roundmax*0.75, roundmax), labels = c(0,0.25, 0.5, 0.75, 1))
  mtext(side = 4, text = "Proportion of CUs > Sgen", line = 2)
  legend("topleft", col = c("black","blue", "red", "orange"), legend = c("ETS","Agg. LRP", "LFRP", "Prop. > Sgen"), lty=1, bty="n")
  mtext(side = 3, text = MUs[mm])
} # end mm loop
dev.off()

