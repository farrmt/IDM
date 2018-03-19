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

#Time vector
Time <- output$Time

#Harvest parameters from files and remove model runs with Rhat > 1.1
for(i in 2:length(filenames)){
  load(filenames[i])
  for(j in 1:length(output$Out[,1,1])){
    if(max(output$Out[j,17:24,], na.rm = TRUE) < 1.2)
      Out <- abind(Out, output$Out[j,,], along = 1)
      Time <- c(Time, output$Time[j])
  }
}

#----------------#
#-Bias of models-#
#----------------#

#Beta0
Out[1,1,2]      #Truth
mean(Out[,2,2]) #PO
mean(Out[,3,2]) #ISDM 1% DS
mean(Out[,4,2]) #ISDM 5% DS
mean(Out[,5,2]) #ISDM 10% DS
mean(Out[,6,2]) #ISDM 15% DS
mean(Out[,7,2]) #ISDM 20% DS
mean(Out[,8,2]) #ISDM 25% DS
mean(Out[,9,2]) #ISDM 50% DS
mean(Out[,10,2]) #1% DS
mean(Out[,11,2]) #5% DS
mean(Out[,12,2]) #10% DS
mean(Out[,13,2]) #15% DS
mean(Out[,14,2]) #20% DS
mean(Out[,15,2]) #25% DS
mean(Out[,16,2]) #50% DS

mean(abs(Out[,1,2]-Out[,2,2])) #PO R
mean(abs(Out[,1,2]-Out[,3,2])) #PO S
mean(abs(Out[,1,2]-Out[,4,2])) #DS R
mean(abs(Out[,1,2]-Out[,5,2])) #DS S
mean(abs(Out[,1,2]-Out[,6,2])) #ISDM R
mean(abs(Out[,1,2]-Out[,7,2])) #ISDM S
mean(abs(Out[,1,2]-Out[,8,2])) #ISDM R DS S PO
mean(abs(Out[,1,2]-Out[,9,2])) #ISDM S DS R PO

hist(Out[,1,2]-Out[,2,2]) #PO R
hist(Out[,1,2]-Out[,3,2]) #PO S
hist(Out[,1,2]-Out[,4,2]) #DS R
hist(Out[,1,2]-Out[,5,2]) #DS S
hist(Out[,1,2]-Out[,6,2]) #ISDM R
hist(Out[,1,2]-Out[,7,2]) #ISDM S
hist(Out[,1,2]-Out[,8,2]) #ISDM R DS S PO
hist(Out[,1,2]-Out[,9,2]) #ISDM S DS R PO

#Beta1
Out[1,1,3]      #Truth
mean(Out[,2,3]) #PO
mean(Out[,3,3]) #ISDM 1% DS
mean(Out[,4,3]) #ISDM 5% DS
mean(Out[,5,3]) #ISDM 10% DS
mean(Out[,6,3]) #ISDM 15% DS
mean(Out[,7,3]) #ISDM 20% DS
mean(Out[,8,3]) #ISDM 25% DS
mean(Out[,9,3]) #ISDM 50% DS
mean(Out[,10,3]) #1% DS
mean(Out[,11,3]) #5% DS
mean(Out[,12,3]) #10% DS
mean(Out[,13,3]) #15% DS
mean(Out[,14,3]) #20% DS
mean(Out[,15,3]) #25% DS
mean(Out[,16,3]) #50% DS

mean(abs(Out[,1,3]-Out[,2,3])) #PO R
mean(abs(Out[,1,3]-Out[,3,3])) #PO S
mean(abs(Out[,1,3]-Out[,4,3])) #DS R
mean(abs(Out[,1,3]-Out[,5,3])) #DS S
mean(abs(Out[,1,3]-Out[,6,3])) #ISDM R
mean(abs(Out[,1,3]-Out[,7,3])) #ISDM S
mean(abs(Out[,1,3]-Out[,8,3])) #ISDM R DS S PO
mean(abs(Out[,1,3]-Out[,9,3])) #ISDM S DS R PO

hist(Out[,1,3]-Out[,2,3]) #PO R
hist(Out[,1,3]-Out[,3,3]) #PO S
hist(Out[,1,3]-Out[,4,3]) #DS R
hist(Out[,1,3]-Out[,5,3]) #DS S
hist(Out[,1,3]-Out[,6,3]) #ISDM R
hist(Out[,1,3]-Out[,7,3]) #ISDM S
hist(Out[,1,3]-Out[,8,3]) #ISDM R DS S PO
hist(Out[,1,3]-Out[,9,3]) #ISDM S DS R PO

#Abundance
Out[1,1,1]
mean(Out[,2,1])
mean(Out[,3,1])
mean(Out[,4,1])
mean(Out[,5,1])
mean(Out[,6,1])
mean(Out[,7,1])
mean(Out[,8,1])
mean(Out[,9,1])

mean(abs(Out[,1,1]-Out[,2,1]))
mean(abs(Out[,1,1]-Out[,3,1]))
mean(abs(Out[,1,1]-Out[,4,1]))
mean(abs(Out[,1,1]-Out[,5,1]))
mean(abs(Out[,1,1]-Out[,6,1]))
mean(abs(Out[,1,1]-Out[,7,1]))
mean(abs(Out[,1,1]-Out[,8,1]))
mean(abs(Out[,1,1]-Out[,9,1]))

hist(Out[,1,1]-Out[,2,1])
hist(Out[,1,1]-Out[,3,1])
hist(Out[,1,1]-Out[,4,1])
hist(Out[,1,1]-Out[,5,1])
hist(Out[,1,1]-Out[,6,1])
hist(Out[,1,1]-Out[,7,1])
hist(Out[,1,1]-Out[,8,1])
hist(Out[,1,1]-Out[,9,1])

#------------------------------#
#-Mean relative bias of models-#
#------------------------------#

#Beta0
mean(abs((Out[,1,2]-Out[,2,2])/Out[,1,2])) * 100 #DS
mean(abs((Out[,1,2]-Out[,3,2])/Out[,1,2])) * 100 #PO
mean(abs((Out[,1,2]-Out[,4,2])/Out[,1,2])) * 100 #ISDM

#Beta1
mean(abs((Out[,1,3]-Out[,2,3])/Out[,1,3])) * 100 #DS
mean(abs((Out[,1,3]-Out[,3,3])/Out[,1,3])) * 100 #PO
mean(abs((Out[,1,3]-Out[,4,3])/Out[,1,3])) * 100 #ISDM

#Abundance
mean(abs((Out[,1,1]-Out[,2,1])/Out[,1,1])) * 100 #DS
mean(abs((Out[,1,1]-Out[,3,1])/Out[,1,1])) * 100 #PO
mean(abs((Out[,1,1]-Out[,4,1])/Out[,1,1])) * 100 #ISDM

#---------#
#-Figures-#
#---------#

library(ggplot2)
library(ggthemes)
name <- c("0%", "1%", "5%", "10%", "15%", "20%", "25%", "50%", 
        "DS1%", "DS5%", "DS10%", "DS15%", "DS20%", "DS25%", "DS50%")
name <- as.character(name)
name <- factor(name, levels=unique(name))

type <- c("Presence Only", rep("Integrated", 7), rep("Distance Sampling", 7))
type <- as.character(type)
type <- factor(type, levels=unique(type))

y75 <- apply(Out[,2:16,], MARGIN = c(2,3), FUN = quantile, probs = 0.75, na.rm = TRUE)
y50 <- apply(Out[,2:16,], MARGIN = c(2,3), FUN = quantile, probs = 0.5, na.rm = TRUE)
y25 <- apply(Out[,2:16,], MARGIN = c(2,3), FUN = quantile, probs = 0.25, na.rm = TRUE)
ymax <-  ((y75 - y25) * 1.5) + y75
ymin <- y25 - ((y75 - y25) * 1.5)

ntot <- data.frame(name, type, ymax[,1], y75[,1], y50[,1], y25[,1], ymin[,1])
colnames(ntot) <- c("name", "type", "ymax", "y75", "y50", "y25", "ymin")

Fig1 <- ggplot(ntot[1:8,], aes(name)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", fill = "grey") +
  geom_hline(yintercept = Out[1,1,1], col = "red") +
  coord_cartesian(ylim = c(Out[1,1,1] + 2500, Out[1,1,1] - 10000)) +
  theme_few() +
  labs(y = "Abundance\n", x = "Coverage")

beta0 <- data.frame(name, type, ymax[,2], y75[,2], y50[,2], y25[,2], ymin[,2])
colnames(beta0) <- c("name", "type", "ymax", "y75", "y50", "y25", "ymin")

Fig2 <- ggplot(beta0[1:8,], aes(name)) + 
        geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                     stat = "identity", fill = "grey") +
        geom_hline(yintercept = Out[1,1,2], col = "red") +
        coord_cartesian(ylim = c(Out[1,1,2] + 1, Out[1,1,2] - 2)) +
        # coord_cartesian(ylim = c(Out[1,1,2] + 0.25, Out[1,1,2] - 0.25)) +
        theme_few() +
        theme(text = element_text(family = "Cambria", size = 24)) +
        labs(y = "Intercept\n", x = "\nCoverage")

FigFig2.5 <- ggplot(beta0[9:15,], aes(name)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", fill = "grey") +
  geom_hline(yintercept = Out[1,1,2], col = "red") +
  coord_cartesian(ylim = c(Out[1,1,2] + 0.25, Out[1,1,2] - 0.25)) +
  theme_few() +
  theme(text = element_text(family = "Cambria", size = 24)) +
  labs(y = "Intercept\n", x = "Coverage")

beta1 <- data.frame(name, type, ymax[,3], y75[,3], y50[,3], y25[,3], ymin[,3])
colnames(beta1) <- c("name", "type", "ymax", "y75", "y50", "y25", "ymin")

Fig3 <- ggplot(beta1[1:8,], aes(name)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", fill = "grey") +
  geom_hline(yintercept = Out[1,1,3], col = "red") +
  coord_cartesian(ylim = c(Out[1,1,3] + 0.125, Out[1,1,3] - 0.125)) +
  theme_few() +
  theme(text = element_text(family = "Cambria", size = 24)) +
  labs(y = "Effect\n", x = "\nCoverage")

Fig3.5 <- ggplot(beta1[9:15,], aes(name)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", fill = "grey") +
  geom_hline(yintercept = Out[1,1,3], col = "red") +
  coord_cartesian(ylim = c(Out[1,1,3] + 0.5, Out[1,1,3] - 0.5)) +
  theme_few() +
  labs(y = "Beta1\n", x = "Coverage")

Fig4 <- ggplot(ntot[c(1,5,12),], aes(type)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", fill = "grey") +
  geom_hline(yintercept = Out[1,1,1], col = "red") +
  coord_cartesian(ylim = c(Out[1,1,1] + 2500, Out[1,1,1] - 8000)) +
  theme_few() +
  theme(text = element_text(family = "Cambria", size = 24)) +
  labs(y = "Abundance\n", x = "")

Fig5 <- ggplot(beta0[c(1,6,13),], aes(type)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", fill = "grey") +
  geom_hline(yintercept = Out[1,1,2], col = "red") +
  coord_cartesian(ylim = c(Out[1,1,2] + 0.5, Out[1,1,2] - 2)) +
  theme_few() +
  theme(text = element_text(family = "Cambria", size = 24)) +
  labs(y = "Intercept\n", x = "")

Fig6 <- ggplot(beta1[c(1,6,13),], aes(type)) + 
  geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
               stat = "identity", fill = "grey") +
  geom_hline(yintercept = Out[1,1,3], col = "red") +
  coord_cartesian(ylim = c(Out[1,1,3] + 0.25, Out[1,1,3] - 0.25)) +
  theme_few() +
  theme(text = element_text(family = "Cambria", size = 24)) +
  labs(y = "Effect\n", x = "")
