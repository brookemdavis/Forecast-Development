
SockDat <- read.csv("DataIn/SRDATA2020.csv")

StkInfo <- read.csv("DataIn/FC_StockInfo.csv")

library(dplyr)
SockDat_Names <- left_join(SockDat, StkInfo %>% select(c("PopID", "PopNmL", "PopNmS", "PopNmFC")) , by = "PopID")


write.csv(SockDat_Names, "DataOut/SRDATA2020_Names.csv")
