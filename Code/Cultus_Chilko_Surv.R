library(dplyr)

# read in 

SRDATA_in <- read.csv("DataIn/SRDATA2020.csv")
Lookup <- read.csv("DataIn/FC_STockInfo.csv")

# Join with Lookup to get names
SRDat <- left_join(SRDATA_in, Lookup[, c("PopID", "PopNmL", "PopSeq")])

CCDat <- SRDat %>% filter(PopNmL %in% c("Cultus", "Chilko"))

CCDat <- CCDat %>% mutate(Surv = rec/juv) %>% mutate(Surv_Tot = (rec3+rec4+rec5)/juv)

write.csv(CCDat, "DataOut/Chilko_Cultus_Surv.csv")
