#--------------------------------------#
#----Code to extract bias data from----#
#----Mara River word documents---------#
#----Created by Matthew Farr-----------#
#--------------------------------------#

#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("~/IDM/DataFormatting/RawData/Notes")

#----------------#
#-Load libraries-#
#----------------#

library(textreadr)
library(dplyr)

#-----------#
#-Load data-#
#-----------#

#Landmark data for missing GPS data
lndmrk <- read.csv("~/ISDM/DataFormatting/RawData/Landmarks/Narok_Landmarks.csv")

#List of word document filenames
filenames <- c(list.files(path = "./Mara River 2012", pattern = "doc", full.names = TRUE), 
               list.files(path = "./Mara River 2013", pattern = "doc", full.names = TRUE), 
               list.files(path = "./Mara River 2014", pattern = "doc", full.names = TRUE))

#Import documents from July 2012 to March 2014
doc <- list()
for(i in 1:(length(filenames)-2)){
  doc[[i]] <- read_document(filenames[i])
}

#--------------#
#-Extract data-#
#--------------#

#Extract session data from documents
out <- list()
for(i in 1:length(doc)){
  out[[i]] <- grep("@", doc[[i]], perl = TRUE, value = TRUE)
}

#Create dataframe to store extracted data
MR <- data.frame(cbind(NA,NA,NA,NA))
colnames(MR) <- c("tmp1", "tmp2", "tmp3", "tmp4")

#Extract GPS coordinates, year, and month from documents
for(i in 1:length(out)){
  tmp1 <- as.numeric(substring(filenames[i], first = regexpr("201", filenames[i]), last = regexpr("[0-9]/", filenames[i])))
  tmp2 <- as.numeric(substring(filenames[i], first = regexpr("[0-9]/", filenames[i]) + 2, last = regexpr("[0-9]/", filenames[i]) + 3))
  for(j in 1:length(out[[i]])){
    tmp <- data.frame(tmp1, tmp2, NA, NA)
    colnames(tmp) <- c("tmp1", "tmp2", "tmp3", "tmp4")
    if(grepl("7[0-9][0-9][0-9][0-9][0-9]", out[[i]][j])){
       tmp3 <- as.numeric(substring(out[[i]][j], first = regexpr("7[0-9][0-9][0-9][0-9][0-9]", out[[i]][j]), last = (regexpr("7[0-9][0-9][0-9][0-9][0-9]", out[[i]][j]) + 6)))
       tmp4 <- as.numeric(substring(out[[i]][j], first = regexpr("98[0-9][0-9][0-9][0-9][0-9]", out[[i]][j]), last = (regexpr("98[0-9][0-9][0-9][0-9][0-9]", out[[i]][j]) + 6)))
       } else {
         tmp3 <- lndmrk[which(lndmrk[,3] == substring(out[[i]][j], first = (regexpr("@", out[[i]][j]) + 1))),5]
         tmp4 <- lndmrk[which(lndmrk[,3] == substring(out[[i]][j], first = (regexpr("@", out[[i]][j]) + 1))),6]
    }
    tryCatch({tmp <- data.frame(tmp1, tmp2, tmp3, tmp4)}, error = function(e){tmp <- data.frame()})
    MR <- rbind(MR, tmp)
    rm(tmp, tmp3, tmp4)
  }
}

#---------------------#
#-Remove Missing Data-#
#---------------------#

removal <- c(1,5,60,78,98,193,214,251,276,364,375,392,453,512,519,525,577,579)
MR <- MR[-removal,]
MR <- na.omit(MR)
rownames(MR) <- NULL

#--------------#
#-Write to CSV-#
#--------------#
write.csv(MR, file = "~/ISDM/DataFormatting/RawData/ObsBiasMR.csv")
