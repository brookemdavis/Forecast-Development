#Update Table 2

library(dplyr)
library(zoo) # need rollapply for generational average
library(EnvStats)


# read in SRdat and CHilliwack data

Chill <- read.csv("DataIn/Chilliwack_updaterec.csv")
SRDATA_in <- read.csv("DataIn/SRDATA2020.csv")
Lookup <- read.csv("DataIn/FC_STockInfo.csv")

# Join with Lookup to get names
SRDat <- left_join(SRDATA_in, Lookup[, c("PopID", "PopNmL", "PopSeq")])


Chill$PopNmL <- "Chilliwack"
Chill$PopSeq <- 8.5 

# Join with Chilliwack
SRDat <- full_join(SRDat, Chill, by = c("rec4", "rec5", "rec", "PopNmL", "PopSeq", "yr" = "Year", "eff" = "eff_fem"))

# add one fish (0.000001) to all rec4 so that no zero survivals
SRDat$rec4 <- SRDat$rec4 + 0.000001

#============================================================================

# Now build function to make new tables for 2019 and 2020

# is 2005 still lowest year?

Surv_All <- SRDat %>% 
  mutate(Surv = rec4/eff)

Min_Year <- Surv_All %>%
  group_by(PopID, PopNmL) %>%
  summarise(Min_Year = yr[which.min(Surv)])

# Pull out 2013
Surv2013 <- Surv_All %>% filter(yr == 2013) %>% select(Surv)
Surv2015 <- Surv_All %>% filter(yr == 2015) %>% select(Surv)

write.csv(cbind(Surv2013, Surv2015), "DataOut/Surv2013_2015.csv")

# would have only had up to BY 2013 in 2018, have up to BY 2014 in 2020



Table2 <- function(Data = SRData, CompYear = 2005, RecentYear = 2013, Name = "2019Table2"){
  
  Surv_All <- Data %>% 
    mutate(Surv = rec4/eff)
  
  # Deal with special cases
  # Pitt is 5-year-old
  PittSurv <- Data %>% 
    filter(PopNmL == "Upper Pitt River") %>%
    mutate(Surv = rec5/eff) %>% 
    select(Surv)
  
  Surv_All[Surv_All$PopNmL == "Upper Pitt River", "Surv"] <- PittSurv
  
  
  # Harrison is total survival
  HarrSurv <- Data %>% 
    filter(PopNmL == "Harrison") %>%
    mutate(Surv = rec/eff) %>% 
    select(Surv)
 
  Surv_All[Surv_All$PopNmL == "Harrison", "Surv"] <- HarrSurv
  
  # Cultus survivals are by smolt
  CultusSurv <- Data %>% 
    filter(PopNmL == "Cultus") %>%
    mutate(Surv = rec4/juv) %>% 
    select(Surv)
  
  Surv_All[Surv_All$PopNmL == "Cultus", "Surv"] <- CultusSurv
  
  
  Geo_Mean <- Surv_All %>%
    group_by(PopNmL, PopSeq) %>% 
    summarise(Geo_Mean_Surv = geoMean(Surv[is.na(Surv) == F])) %>%
    arrange(PopSeq)
  
  
  SD <- Surv_All %>%
    group_by(PopNmL, PopSeq) %>% 
    summarise(SD_Surv = geoSD(Surv[is.na(Surv) == F])) %>%
    arrange(PopSeq)
  
  # get 4-year running avg then take max of those
  Surv_All <- Surv_All %>%
    group_by(PopNmL, PopSeq) %>% 
    mutate(Gen_Mean_Surv = rollapply(Surv, 4, geoMean, fill = NA, align="right"))  
  
  Max_Geo <- Surv_All %>%
    summarise(Peak_Surv = max(Gen_Mean_Surv, na.rm=T)) %>%
    arrange(PopSeq)
  
  # now get "comp year" value
  Comp_Surv <- Surv_All %>% filter(yr == CompYear) %>%
    select(PopNmL, PopSeq, Surv) %>%
    arrange(PopSeq)
  names(Comp_Surv)[3] <- paste("Surv", CompYear, sep="")
  
  # now get recent gen. geo. avg
  Recent <- Surv_All %>% filter(yr == RecentYear)  %>%
    select(PopNmL, PopSeq, Gen_Mean_Surv) %>%
    arrange(PopSeq)
  
  
  # recent for Pitt will need to be a year earlier
  Pitt_Recent <- Surv_All %>% filter(PopNmL == "Upper Pitt River" & yr == RecentYear-1)  
  
  Recent$Gen_Mean_Surv[Recent$PopNmL == "Upper Pitt River"] <- Pitt_Recent$Gen_Mean_Surv
  names(Recent)[3] <- paste("Surv", RecentYear, sep="")
  
 # Now put it all together 
  
  DFOut <- Geo_Mean %>% 
    left_join(Max_Geo) %>%
    left_join(Comp_Surv)%>%
    left_join(Recent) %>%
    left_join( SD) 
  
 write.csv(DFOut, paste(Name, ".csv", sep=""))
 
  
}

# for 2019 only use data up to that year
SRDat2019 <- SRDat %>% filter(yr <= 2014)
Table2(Data = SRDat2019, CompYear = 2005, RecentYear = 2014, Name = "2019Table2")
Table2(Data = SRDat, CompYear = 2005, RecentYear = 2015, Name = "2020Table2_Comp2005")
Table2(Data = SRDat, CompYear = 2013, RecentYear = 2015, Name = "2020Table2_Comp2013")
  

# See if can re-create original table

# regular stocks

Dat2012 <- SRDat %>% filter(yr <= 2012)

Surv_All <- Dat2012 %>% 
  mutate(Surv = rec4/eff)

Geo_Mean <- Surv_All %>%
  group_by(PopNmL, PopSeq) %>% 
  summarise(geoMean(Surv[is.na(Surv) == F])) %>%
  arrange(PopSeq)

# get 4-year running avg then take max of those
Max_Geo <- Surv_All %>%
  group_by(PopNmL, PopSeq) %>% 
  mutate(Gen_Mean_Surv = rollapply(Surv, 4, geoMean, fill = NA, align="right"))  %>%
  summarise(max(Gen_Mean_Surv, na.rm=T)) %>%
  arrange(PopSeq)

# CHilko is only one that is different -- not too worried


# Pitt is 5-year-old, Quesnel and Late Shu are cycle averages
Pitt <- Dat2012 %>% 
  filter(PopNmL == "Upper Pitt River") %>%
  mutate(Surv5 = rec5/eff) %>% 
  group_by(PopNmL, PopSeq) %>% 
  summarise(geoMean(Surv5[is.na(Surv5) == F])) 

# match!

# Harrison is total survival
Harr <- Dat2012 %>% 
  filter(PopNmL == "Harrison") %>%
  mutate(Surv = rec/eff) %>% 
  group_by(PopNmL, PopSeq) %>% 
  summarise(geoMean(Surv[is.na(Surv) == F]))
# Looks like there was an error in Harrison, and it was actually just age-4
# if it's total survival, it's much higher


# Cultus survivals are by smolt
Cultus <- Dat2012 %>% 
  filter(PopNmL == "Cultus") %>%
  mutate(Surv = rec4/juv) %>% 
  summarise(geoMean(Surv[is.na(Surv) == F]))
#match

# Quesnel and late Shu are cycle averages
# get cycle average surv, then take geo mean
QandL <- Dat2012 %>% 
  filter(PopNmL %in% c("Quesnel", "Late Shuswap")) %>%
  mutate(Surv = rec4/eff) %>% 
  group_by(PopID) %>% mutate(Gen_Mean_Surv = rollapply(Surv, 4, mean, fill = NA, align="right"))  %>%
  # 1956 is first year with completed cycle 4
  filter(yr > 1955, cyc == 4) %>%
  # now take average of cycle averages
  group_by(PopNmL, PopSeq) %>% 
  summarise(geoMean(Gen_Mean_Surv[is.na(Gen_Mean_Surv) == F]))

# looks like took arith. average for gen_mean
# can't get exact same answer but close


mutate(Surv = rec4/juv) %>% 
  summarise(geoMean(Surv[is.na(Surv) == F]))
  