#------------------------------------#
#----Script to format observation----#
#----bias data from Mara Hyena-------# 
#----project session data------------#
#----Created by Matthew Farr---------#
#------------------------------------#

#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("~/IDM/DataFormatting")

#----------------#
#-Load libraries-#
#----------------#

library(tidyr)
library(dplyr)

#-----------#
#-Load data-#
#-----------#

talek <- tbl_df(read.csv(file = "RawData/ObsBiasTalek.csv")) %>% select(DATE, landmarkid, utme, utmn)
colnames(talek) <- c("date", "Name", "UTME", "UTMN")
serena <- tbl_df(read.csv(file = "RawData/ObsBiasSerena.csv")) %>% select(DATE, landmarkid, utme, utmn)
colnames(serena) <- c("date", "Name", "UTME", "UTMN")
hz <- tbl_df(read.csv(file = "RawData/ObsBiasHZ.csv")) %>% select(DATE, landmarkid, utme, utmn)
colnames(hz) <- c("date", "Name", "UTME", "UTMN")
ft <- tbl_df(read.csv(file = "RawData/ObsBiasFT.csv")) %>% select(DATE, landmarkid, utme, utmn)
colnames(ft) <- c("date", "Name", "UTME", "UTMN")
mr <- tbl_df(read.csv(file = "RawData/ObsBiasMR.csv")) %>% select(tmp1, tmp2, tmp3, tmp4)
colnames(mr) <- c("year", "month", "UTME", "UTMN") 

#----------------#
#-Load landmarks-#
#----------------#

narokL <- tbl_df(read.csv(file = "RawData/Landmarks/Narok_Landmarks.csv"))
serenaL <- tbl_df(read.csv(file = "RawData/Landmarks/Serena_Landmarks.csv"))

#------------------#
#-Filter landmarks-#
#------------------#

narokL <- narokL %>% select(NAME, UTME, UTMN)
colnames(narokL) <- c("Name", "UTME", "UTMN")
serenaL <- serenaL %>% select(Name, UTME, UTMN)
colnames(serenaL) <- c("Name", "UTME", "UTMN")

#----------------#
#-Remove strings-#
#----------------#

talek$Name <- gsub(" \\(.*", "", talek$Name)
serena$Name <- gsub(" \\(.*", "", serena$Name)
hz$Name <- gsub(" \\(.*", "", hz$Name)
ft$Name <- gsub(" \\(.*", "", ft$Name)

#-------------#
#-Replace NAs-#
#-------------#

talek[is.na(talek)] <- 0
serena[is.na(serena)] <- 0
hz[is.na(hz)] <- 0
ft[is.na(ft)] <- 0

#---------------------#
#-Remove problem data-#
#---------------------#

talek <- talek %>% 
  filter(UTME < 999999) %>% 
  filter(UTMN < 9999999) %>%
  filter(Name != "~location non standard")

serena <- serena %>% 
  filter(UTME < 999999) %>% 
  filter(UTMN < 9999999) %>%
  filter(Name != "")

hz <- hz %>%
  filter(UTME < 999999) %>% 
  filter(UTMN < 9999999)

ft <- ft %>%
  filter(UTME < 999999) %>% 
  filter(UTMN < 9999999)

#------------------------#
#-Replace missing values-#
#------------------------#

for(i in 1:length(talek$UTME)){
  if(talek[i,3] == 0)
    tryCatch({talek[i,3] <- narokL[which(narokL$Name == talek[i,]$Name),2]}, error = function(e){})
  if(talek[i,4] == 0)
    tryCatch({talek[i,4] <- narokL[which(narokL$Name == talek[i,]$Name),3]}, error = function(e){})
}

for(i in 1:length(serena$UTME)){
  if(serena[i,3] == 0)
    tryCatch({serena[i,3] <- serenaL[which(serenaL$Name == serena[i,]$Name),2]}, error = function(e){})
  if(serena[i,4] == 0)
    tryCatch({serena[i,4] <- serenaL[which(serenaL$Name == serena[i,]$Name),3]}, error = function(e){})
}

for(i in 1:length(hz$UTME)){
  if(hz[i,3] == 0)
    tryCatch({hz[i,3] <- serenaL[which(serenaL$Name == hz[i,]$Name),2]}, error = function(e){})
  if(hz[i,4] == 0)
    tryCatch({hz[i,4] <- serenaL[which(serenaL$Name == hz[i,]$Name),3]}, error = function(e){})
}

for(i in 1:length(ft$UTME)){
  if(ft[i,3] == 0)
    tryCatch({ft[i,3] <- narokL[which(narokL$Name == ft[i,]$Name),2]}, error = function(e){})
  if(ft[i,4] == 0)
    tryCatch({ft[i,4] <- narokL[which(narokL$Name == ft[i,]$Name),3]}, error = function(e){})
}

#---------------------#
#-Remove Missing Data-#
#---------------------#

talek <- talek %>%
  filter(UTME > 100000) %>%
  filter(UTMN > 1000000)

serena <- serena %>%
  filter(UTME > 100000) %>%
  filter(UTMN > 1000000)

hz <- hz %>%
  filter(UTME > 100000) %>%
  filter(UTMN > 1000000)

ft <- ft %>%
  filter(UTME > 100000) %>%
  filter(UTMN > 1000000)

#-----------#
#-Round UTM-#
#-----------#

talek$UTME <- round(talek$UTME)
talek$UTMN <- round(talek$UTMN)
serena$UTME <- round(serena$UTME)
serena$UTMN <- round(serena$UTMN)
hz$UTME <- round(hz$UTME)
hz$UTMN <- round(hz$UTMN)
ft$UTME <- round(ft$UTME)
ft$UTMN <- round(ft$UTMN)

#------------#
#-Split Date-#
#------------#

talek$date <- as.Date(talek$date, format = "%m/%d/%Y")
talek$year <- as.numeric(format(talek$date, format = "%Y"))
talek$month <- as.numeric(format(talek$date, format = "%m"))
serena$date <- as.Date(serena$date, format = "%m/%d/%Y")
serena$year <- as.numeric(format(serena$date, format = "%Y"))
serena$month <- as.numeric(format(serena$date, format = "%m"))
hz$date <- as.Date(hz$date, format = "%m/%d/%Y")
hz$year <- as.numeric(format(hz$date, format = "%Y"))
hz$month <- as.numeric(format(hz$date, format = "%m"))
ft$date <- as.Date(ft$date, format = "%m/%d/%Y")
ft$year <- as.numeric(format(ft$date, format = "%Y"))
ft$month <- as.numeric(format(ft$date, format = "%m"))

#---------------#
#-Remove fields-#
#---------------#

talek <- talek %>% select(year, month, UTME, UTMN)
colnames(talek) <- c("year", "month", "UTME", "UTMN")
serena <- serena %>% select(year, month, UTME, UTMN)
colnames(serena) <- c("year", "month", "UTME", "UTMN")
hz <- hz %>% select(year, month, UTME, UTMN)
colnames(hz) <- c("year", "month", "UTME", "UTMN")
ft <- ft %>% select(year, month, UTME, UTMN)
colnames(ft) <- c("year", "month", "UTME", "UTMN")

#-------------#
#-Bias output-#
#-------------#

bias <- rbind(talek, ft, mr, serena, hz)
bias <- bias %>% arrange(year, month) %>%
  filter(year > 2012 | month > 6)

#--------------#
#-Write to CSV-#
#--------------#

write.csv(bias, file = "RawData/Bias.csv")
