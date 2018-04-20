#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("C:/Users/farrm/Documents/GitHub/ISDM/DataAnalysis/CaseStudy")

#------------------------------#
#-Libraries used in simulation-#
#------------------------------#

library(raster)
library(rgdal)
library(jagsUI)

#----------------#
#-Spatial Extent-#
#----------------#
SEstack <- raster("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/DSExtent.tif")
SEdata <- raster::as.matrix(SEstack)
add <- 1
coverage2 <- NULL
for(i in 1:length(SEdata[,1])){
  for(j in 1:length(SEdata[1,])){
    if(!is.na(SEdata[i,j])){
      SEdata[i,j] <- add
      add <- add + 1
    }
  }
  coverage2 <- c(coverage2, as.numeric(na.omit(rep(rep(SEdata[i,], each = 20, 20)))))
}

#----------------------#
#-Import Obs Bias data-#
#----------------------#

#Temporal extent/order
#te <- c(28,31,34,1,4,10,8,14,17,20,23,26,29,32,35,2,5,11,9,15,18)
te <- c(16,18,20,1,3,7,5,9,11,13,14,15,17,19,21,2,4,8,6,10,12)

dir <- "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/Bias_Data"
Biaslist <- list.files(path = dir, pattern = "tif$", full.names = TRUE)
Biasstack <- stack(Biaslist[te])
Biasdata <- raster::as.matrix(Biasstack)
Bias <- Biasdata[which(!is.na(Biasdata[,1])),]
Bias <- scale(Bias)

#----------------#
#-Import PO data-#
#----------------#

dir <- "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/PO_Data"
POlist <- list.files(path = dir, pattern = "tif$", full.names = TRUE)
POstack <- stack(POlist[te])
POdata <- raster::as.matrix(POstack)
PO <- POdata[which(!is.na(Biasdata[,1])),]
PO[is.na(PO)] <- 0

#----------------#
#-Import DS data-#
#----------------#

#Subregion of spatial extent that is covered by distance sampling (all DS pixels within PO pixel)
regionB <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/DSpixel.tif")
regionB <- raster::as.matrix(regionB)
#Number of DS pixels within PO pixels; 1 PO pixel is 20 x 20 50m^2 DS pixels; 235 PO pixels containing DS (20x20x235)
B <- length(which(!is.na(regionB)))

#Subregion of spatial extent that is covered by distance sampling (only DS pixels)
transect <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/Transect.tif")
transect <- raster::as.matrix(transect)
transect[transect == 128] <- NA
transect <- transect[which(!is.na(regionB))]
coverage1 <- which(!is.na(transect))

dir <- "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/DS_Data/"
DSlist <- list.files(path = dir, pattern = "tif$", full.names = TRUE)
DSstack <- stack(DSlist[te])
DS<- raster::as.matrix(DSstack)
DS <- DS[which(!is.na(regionB)),]
DS <- DS[coverage1,]
DS[is.na(DS)] <- 0




#----------------------------#
#-Variance covariance matrix-#
#----------------------------#

UTM <- as.matrix(xyFromCell(SEstack, which(SEdata[,1]==1)))
COV <- as.matrix(dist(UTM))

#---------#
#-Indices-#
#---------#

G <- length(PO[,1])
TT <- length(PO[1,])

#--------------#
#-Compile JAGS-#
#--------------#

data <- list(y.po = PO, w = Bias, G = G, TT = TT)

params <- c("alpha0", "alpha1", "tau", "beta1", "beta0")

inits <- function(){list(beta1 = runif(G, -1, 1), alpha1 = runif(1, -0.05, 0.05), 
                         beta0 = runif(1, -1, 1), alpha0 = runif(1, -9, -8))}

#-------------#
#-MCMC values-#
#-------------#

nb <- 20000
ni <- 29000
nt <- 3
nc <- 3
na <- 100

out2 <- jagsUI(data, inits, params, "PO2.txt", n.thin=nt, 
          n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

###
tst <- SEstack
values(tst)[which(values(tst)==0)] <- exp(out2$mean$beta1)
plot(tst)

which.max(out$mean$beta0)

expit <- function(eta) 
{
  1/(1+exp(-eta))
}

tst2 <- SEstack
values(tst2)[which(values(tst2)==0)] <- expit(out$mean$alpha0 + out$mean$alpha1*Bias[,1])
plot(tst2)
