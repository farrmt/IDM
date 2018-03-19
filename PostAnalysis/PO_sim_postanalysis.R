#------------------------------------------#
#----Post analysis of simulation study.----#
#----Created by Matthew Farr---------------#
#------------------------------------------#

#----------------#
#-Load libraries-#
#----------------#

library(abind)

#-------------#
#-Merge files-#
#-------------#

#List of all filenames
filenames <- list.files(pattern = "output")

#Load first file
load(filenames[1])

#Initialize vector for all output
Out <- output$Out

#Harvest parameters from files and remove model runs with Rhat > 1.1
for(i in 2:length(filenames)){
  load(filenames[i])
  for(j in 1:length(output$Out[,1,1])){
    if(max(output$Out[j,4:5,], na.rm = TRUE) < 1.2)
      Out <- abind(Out, output$Out[j,,], along = 1)
  }
}

#----------------#
#-Bias of models-#
#----------------#

#Beta0
mean(abs(Out[,1,2]-Out[,2,2])) #DS
mean(abs(Out[,1,2]-Out[,3,2])) #DSalt

hist(Out[,1,2]-Out[,2,2], main = "PO no zeros", xlab = "Intercept Bias") #DS vs Truth
hist(Out[,1,2]-Out[,3,2], main = "PO zeros", xlab = "Intercept Bias") #DSalt vs Truth

#Beta1
mean(abs(Out[,1,3]-Out[,2,3])) #DS
mean(abs(Out[,1,3]-Out[,3,3])) #DSalt

hist(Out[,1,3]-Out[,2,3], main = "PO no zeros", xlab = "Effect Bias") #DS vs Truth
hist(Out[,1,3]-Out[,3,3], main = "PO zeros", xlab = "Effect Bias") #DSalt vs Truth

#Abundance
mean(abs(Out[,1,1]-Out[,2,1])) #DS
mean(abs(Out[,1,1]-Out[,3,1])) #DSalt

hist(Out[,1,1]-Out[,2,1], main = "PO no zeros", xlab = "Abundance Bias") #DS vs Truth
hist(Out[,1,1]-Out[,3,1], main = "PO zeros", xlab = "Abundance Bias") #DSalt vs Truth

#------------------------------#
#-Mean relative bias of models-#
#------------------------------#

#Beta0
mean(abs((Out[,1,2]-Out[,2,2])/Out[,1,2])) * 100 #DS
mean(abs((Out[,1,2]-Out[,3,2])/Out[,1,2])) * 100 #DSalt

#Beta1
mean(abs((Out[,1,3]-Out[,2,3])/Out[,1,3])) * 100 #DS
mean(abs((Out[,1,3]-Out[,3,3])/Out[,1,3])) * 100 #DSalt

#Abundance
mean(abs((Out[,1,1]-Out[,2,1])/Out[,1,1])) * 100 #DS
mean(abs((Out[,1,1]-Out[,3,1])/Out[,1,1])) * 100 #DSalt

