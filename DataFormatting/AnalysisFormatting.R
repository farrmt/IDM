#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("~/IDM/DataFormatting")

#------------------------------#
#-Libraries used in simulation-#
#------------------------------#

library(raster)
library(rgdal)
library(igraph)
library(Matrix)

#----------------#
#-Spatial Extent-#
#----------------#

#Coverage by time


#Extent of distance sampling @ 50m*50m
scale3stack <- stack(c("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/tst5.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/tst6.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage1.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage2.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage3.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage4.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage5.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage6.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage7.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage1.3.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage2.3.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage3.3.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage4.3.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage5.3.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage6.3.tif",
                       "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage7.3.tif"))

scale3 <- raster::as.matrix(scale3stack)
scale4 <- scale3[!is.na(scale3[,2]),]
scale4 <- cbind(scale4, seq(1, 94000, 1))
scale4 <- scale4[order(scale4[,2]),]
#coverage2 <- scale4[,2]
scale4 <- cbind(scale4, seq(1, 94000, 1))

B1 <- B2 <- c1 <- c2 <- c3 <- c4 <- Bend <- Dend <- NULL
#fix na issue

for(i in 1:7){
  B1 <- scale4[,c(1, i+2, i+9, 17, 18)]
  B1 <- B1[order(B1[,3]),] #order by total sampled pixels
  c1 <- cbind(c1, B1[,3]) #1000m pixel ID
  c2 <- cbind(c2, B1[,5]) #50m pixel ID in order of 1000m pixels
  Bend <- c(Bend, max(which(!is.na(c1[,i])))) #Length of sampled pixels @ time t
  B2 <- B1[!is.na(B1[,1]),] #remove pixels not within 650m DS
  B2 <- B2[order(B2[,4]),] #order by upload ID
  B2 <- cbind(B2, seq(1, length(B2[,1]), 1)) #add upload ID for DS extent
  B2 <- B2[order(B2[,2]),] #order by DS pixels
  c3 <- cbind(c3, B2[,5]) #50m pixel ID in order of 1000m pixels
  c4 <- cbind(c4, B2[,6]) #50m pixel ID in order of upload
  Dend <- c(Dend, max(which(!is.na(B2[,2]))))
}

#Pixels for change of support during time t
# B <- scale4[,10:16]
# coverage2 <- scale4[,10:16]
# B[is.na(B)] <- 0
# for(i in 1:7){
#   for(j in 1:94000){
#     if(B[j,i]!=0){
#       B[j,i] <- scale4[j,18]
#     }
#   }
# }
# coverage2[,1] <- coverage2[order(B[,1]),1]
# B[,1] <- B[order(B[,1]),1]
# 
# #Pixels within distance sampling 650m
scale5 <- scale4[!is.na(scale4[,1]),]
scale5 <- scale5[order(scale5[,17]),]
# coverage1 <- scale5[,18]
# 
# #Pixels within distance sampling during time t
# DD <- scale5[,3:9]
# DD[is.na(DD)] <- 65727
# for(i in 1:7){
#   for(j in 1:length(coverage1)){
#     if(DD[j,i]!=65727){
#       DD[j,i] <- j
#     }
#   }
# }
# 
# DD <- rbind(DD, rep(65727, 7))

#Extent of opportunistic sampling @ 1000m*1000m
scale1stack <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/tst2.tif")
scale1 <- raster::as.matrix(scale1stack)
#Extent of distance sampling @ 1000m*1000m
scale2stack <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/tst7.tif",
                     "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage1.2.tif",
                     "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage2.2.tif",
                     "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage3.2.tif",
                     "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage4.2.tif",
                     "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage5.2.tif",
                     "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage6.2.tif",
                     "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/Coverage/Coverage7.2.tif")
scale2 <- raster::as.matrix(scale2stack)
#Remove pixels not within extent of opportunistic sampling
scale2 <- scale2[which(!is.na(scale1)),]
#Pixels within DS == 1
scale2[is.na(scale2)] <- 0
#Pixels not within DS == 1
scale <- scale2[,-1]
#Starting and ending values for change of support
cover <- c(4, 7, 4, 7, 4, 3, 7, 7, 7, NA, 7, 1, 6, 7, 5, 2, 7, 4, 3, 7, 4)
Gstart <- array(NA, dim = dim(scale))
Gend <- array(NA, dim = dim(scale))
start <- 1
end <- 400
for(i in 1:555){
  if(scale[i,7]==1){
    Gstart[i,7] <- start
    start <- start + 400
    Gend[i,7] <- end
    end <- end + 400
  }else{
    Gstart[i,7] <- c2[1,cover[7]]
    Gend[i,7] <- c2[1,cover[7]] + 1
  }
  for(j in 1:6){
    if(scale[i,j]==1){
      Gstart[i,j] <- Gstart[i,7]
      Gend[i,j] <- Gend[i,7]
    }else{
    Gstart[i,j] <- c2[1,j]
    Gend[i,j] <- c2[1,j] + 1
    }
  }
}

#----------------------#
#-Import Obs Bias data-#
#----------------------#

#Temporal extent/order
te <- c(16,18,20,1,3,7,5,9,11,13,14,15,17,19,21,2,4,8,6,10,12)

dir <- "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/Bias_Data"
Biaslist <- list.files(path = dir, pattern = "tif$", full.names = TRUE)
Biasstack <- stack(Biaslist[te])
Biasdata <- raster::as.matrix(Biasstack)
Bias1 <- Biasdata[which(!is.na(Biasdata[,1])),]
Bias <- scale(Bias1)

#----------------#
#-Import PO data-#
#----------------#

dir <- "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/PO_Data"
POlist <- list.files(path = dir, pattern = "tif$", full.names = TRUE)
POstack <- stack(POlist[te])
POdata <- raster::as.matrix(POstack)
PO <- POdata[which(!is.na(Biasdata[,1])),]
PO[is.na(PO)] <- 0

#---------#
#-Indices-#
#---------#

D <- length(c3[,1])
G <- length(PO[,1])
TT <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21)
scale <- scale

#----------------#
#-Import DS data-#
#----------------#

#Subregion of spatial extent that is covered by distance sampling (all DS pixels within PO pixel)
regionB <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/DSpixel.tif")
regionB <- raster::as.matrix(regionB)
#Number of DS pixels within PO pixels; 1 PO pixel is 20 x 20 50m^2 DS pixels; 235 PO pixels containing DS (20x20x235)
# B <- length(which(!is.na(regionB)))

dir <- "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/DS_Data/"
DSlist <- list.files(path = dir, pattern = "tif$", full.names = TRUE)
DSstack <- stack(DSlist[te])
DS<- raster::as.matrix(DSstack)
DS <- DS[which(!is.na(regionB)),]
DS <- DS[scale5[,17],]
DS[is.na(DS)] <- 0

y.ds <- array(NA, dim = dim(DS))
for(t in TT){
  for(d in 1:Dend[cover[t]]){
    y.ds[d,t] <- DS[c4[d, cover[t]],t]
  }
}

#Distance of each DS pixel
dststack <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/Distance.tif")
dst <- raster::as.matrix(dststack)
dst <- dst[which(!is.na(regionB))]
dst <- dst[scale5[,17]]

#PixelID for observed distances
# pixD <- NULL
# for(t in TT){
#   pixD <- c(pixD, which(y.ds[,t]>0))
# }
# 
# nobs <- length(pixD)

mdpt <- seq(25, 625, 50)
nD <- length(mdpt)
nobs <- 90
dclass <- c(4,1,1,3,1,1,13,1,2,1,5,2,2,3,1,1,1,3,1,1,1,5,1,1,
            4,4,2,2,3,1,3,1,3,2,2,1,3,1,2,1,3,1,2,3,3,3,2,1,
            1,1,1,1,2,1,1,1,2,4,3,4,1,1,5,2,2,1,1,1,1,1,1,2,
            1,1,2,1,1,2,1,3,2,1,1,2,2,3,1,2,4,2)

#------------#
#-Covariates-#
#------------#
regstack <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/Region.tif")
region <- raster::as.matrix(regstack)
region <- region[which(!is.na(scale1))]
region[is.na(region)] <- 0

borderstack <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/Border.tif")
border <- raster::as.matrix(borderstack)
border <- border[which(!is.na(border))]
border <- scale(border)

waterstack <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ArcGIS/BaseData/Water.tif")
water <- raster::as.matrix(waterstack)
water <- water[which(!is.na(water))]
water <- scale(water)

NDVInames <- list.files(path = "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/NDVI/", pattern = "tif$", full.names = TRUE)
NDVIstack <- stack(NDVInames[4:43])
values(NDVIstack) <- values(NDVIstack) * 0.0001
NDVI <- raster::as.matrix(NDVIstack)
NDVI <- NDVI[which(!is.na(NDVI[,1])),]
tmp <- matrix(NA, nrow = 555, ncol = 21)
for(i in 1:dim(NDVI)[1]){
tmp[i,] <- cbind(mean(NDVI[i,1:2]), mean(NDVI[i,3:4]), mean(NDVI[i,5:6]), mean(NDVI[i,7:8]), NDVI[i,9], mean(NDVI[i,10:11]), mean(NDVI[i,12:13]), mean(NDVI[i,14:15]), mean(NDVI[i,16:17]), mean(NDVI[i,18:19]),
              mean(NDVI[i,20:21]), mean(NDVI[i,22:23]), mean(NDVI[i,24:25]), mean(NDVI[i,26:27]), mean(NDVI[i,28:29]), NDVI[i,30], mean(NDVI[i,31:32]), mean(NDVI[i,33:34]), mean(NDVI[i,35:36]) ,mean(NDVI[i,37:38]),
              mean(NDVI[i,39:40]))
}
NDVI <- scale(tmp)

Lionnames <- list.files(path = "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/Lion_KDE/", pattern = "tif$", full.names = TRUE)
Lionstack <- stack(Lionnames[te])
Lion <- raster::as.matrix(Lionstack)
Lion <- Lion[which(!is.na(Lion[,1])),]
Lion <- Lion/(Bias1 + 0.1)
Lion <- scale(Lion)

LionSstack <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/Lion_KDE/KDEFULL.tif")
LionS <- raster::as.matrix(LionSstack)
LionS <- LionS[which(!is.na(LionS))]
BiasSstack <- stack("C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/Bias_Data/BiasFULL.tif")
BiasS <- raster::as.matrix(BiasSstack)
BiasS <- BiasS[which(!is.na(BiasS))]
LionS <- LionS/(BiasS + 0.1)
LionS <- scale(LionS)

#---------------#
#-CAR Matricies-#
#---------------#
load(file = "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/FormattedData/Adjacency.Rdata")
Wcar <- na.omit(W)
Wcar <- as_adjacency_matrix(graph.edgelist(as.matrix(Wcar)), sparse = TRUE)

Dcar <- diag(rowSums(Wcar))

# Ccar <- diag(sqrt(1/rowSums(Wcar))) %*% Wcar %*% diag(sqrt(1/rowSums(Wcar)))
# alpha <- max(eigen(C)$values) - 0.05

#---------#
#-Compile-#
#---------#

IDMdata <- list(D, Dend, Bend, Gstart, Gend, G, TT, c1, c2, c3, c4, scale, 
                 PO, Bias, y.ds, dst, mdpt, nD, dclass, nobs, cover, region, border, water, NDVI, Lion, LionS, Wcar, Dcar)
heads <- c("D", "Dend", "Bend", "Gstart", "Gend", "G", "TT", "c1", "c2", "c3", "c4", "scale", 
           "y.po", "Bias", "y.ds", "dst", "mdpt", "nD", "dclass", "nobs", "cover", "region", "border", "water", "NDVI", "Lion", "LionS", "Wcar", "Dcar")
IDMdata <- setNames(IDMdata, nm = heads)

#-------------#
#-Export data-#
#-------------#

save(ISDMdata, file = "C:/Users/farrm/Documents/GitHub/ISDM/DataFormatting/ISDM.Rdata")
