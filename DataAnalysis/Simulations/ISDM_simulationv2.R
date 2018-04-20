#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("C:/Users/farrm/Documents/GitHub/ISDM/DataAnalysis/Simulations")

#------------------------------#
#-Libraries used in simulation-#
#------------------------------#

library(mvtnorm)
library(jagsUI)
library(dplyr)
library(raster)

#------------------------------#
#-Functions used in simulation-#
#------------------------------#

#Euclidean distance (From library(AHMbook))
# e2dist <- function (x, y) 
# {
#   i <- sort(rep(1:nrow(y), nrow(x)))
#   dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
#   matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
# }

expit <- function(eta) 
{
  1/(1+exp(-eta))
}

logit <- function(pp) 
{ 
  log(pp) - log(1-pp) 
}

#--------------------------#
#-Create Sampling Region B-#
#--------------------------#

#Dimension of grid for Region B
W <- 100
#Pixel size
px <- 1
#Mid-point coordinates of first dimension
px.m <- seq(px/2, W - px/2, px)
#Expand mid-point to entire grid
gr <- expand.grid(px.m, px.m)

#-----------------------#
#-Draw parameter values-#
#-----------------------#

#Intercept parameter for intensity function
beta0 <- log(0.05) #Low Abundance (500 individuals; 0.01/pixel)
# beta0 <- log(0.5) #Medium Abundance (5,000 individuals; 1/pixel)
# beta0 <- log(5) #High Abundance (50,000 individuals; 10/pixel)

#Effect parameter of enivornment on intensity
# beta1 <- 0 #No effect on intensity (i.e., homogenous)
# beta1 <- 0.5 #Weak effect on intensity X 1.14 (14% increase)
beta1 <- 1 #Meadium effect on intensity X 1.66 (66% inrease)
# beta1 <- 1.25 #Strong effect on intensity X 2.19 (119% increase)

#Intercept parameter for prensence only (PO) detection
alpha0 <- logit(0.1) #Low detection (10%)
# alpha0 <- logit(0.5) #Medium detection (50%)
# alpha0 <- logit(0.99) #High detection (100%)

#Effect parameter of environment on PO detection
# alpha1 <- 0 #No effect on PO detection
alpha1 <- 2.25 #Weak effect on PO detection
# alpha1 <- 5 #Medium effect on PO detection
# alpha1 <- 10 #Strong effect on PO detection

#Scale parameter for DS detection
# sigma <- 0.5 #Low detection
sigma <- 1 #Medium detection
# sigma <- 2.5 #High detection

#-------------------------------------#
#-Draw environmental covariate values-#
#-------------------------------------#

#Variance-covariance matrix based on Euclidean distance
#V <- exp(-e2dist(gr, gr))
#Covariate values from a correlated multivariate normal (pg.534 AHM)
#x <- as.vector(t(chol(V))%*%rnorm(W^2))
#Load x covariate instead
load(file = "xcov.Rdata")

#Visualize covariate
par(mar=c(3,3,3,6))
image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20), axes = FALSE, xlab = "", ylab = "")

#--------------------------------#
#-Simulate number of individuals-#
#--------------------------------#

#Intensity function
intensity <- exp(beta0 + beta1 * x)
#Probability of each pixel based on intensity function
probs <- intensity/sum(intensity)
#Simulate number of individuals in Region B
N <- rpois(1, sum(intensity))
#Sample location, s, of simulated individuals
s <- sample(1:(W*W), N, replace=TRUE, prob=probs)
#Assign X and Y coordinates
u1 <- gr[s,1]
u2 <- gr[s,2]

#Visulaize individuals
points(u1, u2, col = 'black', pch = 20)

#-----------------------#
#-Create Transect Units-#
#-----------------------#

#Transect unit size
tu <- 10
#Mid-point coordinates of first dimension
tu.m <- seq(tu/2, W - tu/2, tu)
#Expand mid-point to entire grid
grt <- expand.grid(tu.m, tu.m)

#-----------------------#
#-Sample Transect Units-#
#-----------------------#

#Number of transect scenarios
nts <- 7
#Samples per scenario
ns <- c(50,25,20,15,10,5,1)
#Sampled transects for each scenario
tsamp <- array(NA, dim = c(50,nts))
#First scenario
tsamp[,1] <- sample(1:100, ns[1], replace = FALSE)
#Remanding scenarios
for(i in 2:nts){
tsamp[,i] <- c(sample(tsamp[1:ns[i-1],i-1], ns[i], replace = FALSE), rep(NA, ns[1] - ns[i]))
}

#----------------------------#
#-Simulate distance sampling-#
#----------------------------#

# #All possible transects
x.line <- rep(tu.m, each = W)
y.line <- rep(seq(0.5, W-0.5, 1), 10)
line <- cbind(x.line, y.line)
# #Number of points in transect line
# J <- length(x.line)
# #Distance array for all points on transect line
# d <- array(NA, dim = c(W*W, J))
# #Distance to nearest transect
# dist <- NULL
# #Simulate above quantities
# for(g in 1:(W*W)){
#   for(j in 1:J){
#     d[g,j] <- sqrt((gr[g,1] - x.line[j])^2 + (gr[g,2] - y.line[j])^2)
#   }
#   dist[g] <- min(d[g,])
# }

#Visualize transect lines
for(i in 1:nts){
  image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))
  points(u1, u2, col = 'black', pch = 20)
  for(j in 1:ns[i]){
    points(line[which((line[,1] > (grt[tsamp[j,i],1] - 5)) & 
                        ((grt[tsamp[j,i],1] + 5) > line[,1]) & 
                        (line[,2] > (grt[tsamp[j,i],2] - 5)) & 
                        ((grt[tsamp[j,i],2] + 5) > line[,2])),], col = 'grey', pch = 20)
  }
}


#Load distance to pixels for robust and sparse
load(file = "dist.Rdata")

#Sampled pixels and corresponding distances
fill <- seq(1,5100,by=100)
coverage <- array(NA, dim = c(5000,nts))
dst <- array(NA, dim = c(5000,nts))
for(i in 1:nts){
  for(j in 1:ns[i]){
    coverage[fill[j]:(fill[j+1]-1),i] <- which((gr[,1] > (grt[tsamp[j,i],1] - 5)) & 
                                           ((grt[tsamp[j,i],1] + 5) > gr[,1]) & 
                                           (gr[,2] > (grt[tsamp[j,i],2] - 5)) & 
                                           ((grt[tsamp[j,i],2] + 5) > gr[,2]))
    dst[fill[j]:(fill[j+1]-1),i] <- dist[coverage[fill[j]:(fill[j+1]-1),i]]
  }
}

#Individual (pixels) within region B
sds<- as.data.frame(s[s%in%coverage[,1]])
colnames(sds) <- "p"
#Detection probability for distance sampling
pds <- rep(NA, length(sds))
#Distance for each individual
dst50 <- as.data.frame(cbind(coverage[,1], dst[,1]))
colnames(dst50) <- c("p", "dst")
y50 <- inner_join(dst50, sds, by = "p")
dst.S <- y50$dst
#Individual presence/absence based on distance sampling
yds <- NULL
#Simulate distance sampling
pds <- exp(-dst.S * dst.S / (2 * sigma * sigma))
yds <- rbinom(length(sds[,1]), 1, pds)
y50$y <- yds
y50$uxds <- gr[y50[,1],1]
y50$uyds <- gr[y50[,1],2]
y50$dclass <- NA
y50 <- y50%>%arrange(p)

tmp2 <- as.data.frame(coverage[,2])
colnames(tmp2) <- c("p")
y25 <- inner_join(y50,tmp2,by="p")
y25 <- y25%>%arrange(p)

tmp2 <- as.data.frame(coverage[,3])
colnames(tmp2) <- c("p")
y20 <- inner_join(y50,tmp2,by="p")
y20 <- y20%>%arrange(p)

tmp2 <- as.data.frame(coverage[,4])
colnames(tmp2) <- c("p")
y15 <- inner_join(y50,tmp2,by="p")
y15 <- y15%>%arrange(p)

tmp2 <- as.data.frame(coverage[,5])
colnames(tmp2) <- c("p")
y10 <- inner_join(y50,tmp2,by="p")
y10 <- y10%>%arrange(p)

tmp2 <- as.data.frame(coverage[,6])
colnames(tmp2) <- c("p")
y5 <- inner_join(y50,tmp2,by="p")
y5 <- y5%>%arrange(p)

tmp2 <- as.data.frame(coverage[,7])
colnames(tmp2) <- c("p")
y1 <- inner_join(y50,tmp2,by="p")
y1 <- y1%>%arrange(p)

YD <- list(y50,y25,y20,y15,y10,y5,y1)
#Size of subregion B
B <- NULL
#Number of obesrvation per pixel of subregion
yds <- list()
#Distance class
for(i in 1:nts){
  B[i] <- length(na.omit(coverage[,i]))
  yds[[i]] <- as.data.frame(table(YD[[i]][YD[[i]][,3]==1,1]))
  yds[[i]][,1] <- as.numeric(as.character(yds[[i]][,1]))
  colnames(yds[[i]]) <- c("p", "freq")
  coverB <- as.data.frame(coverage[1:B[i],i])
  colnames(coverB) <- "p"
  yds[[i]] <- full_join(yds[[i]], coverB, by = "p")
  yds[[i]] <- yds[[i]]%>%arrange(p)
  yds[[i]]$dst <- dst[1:B[i],i]
}

yds <- rapply(yds, f=function(x) ifelse(is.na(x),0,x), how="replace" )

#---------------------------------------#
#-Draw covariate value for PO detection-#
#---------------------------------------#
#Environmental covariate on PO detection (Multivariate normal)
#Mean of x-dim distribution
mu.x <- 50
#Mean of y-dim distribution
mu.y <- 50
#Variance of x-dim distribution
sigmax <- 0.3*abs(W) #Robust PO
#Variance of y-dim distribution
sigmay <- 0.3*abs(W)
#Covariance of x-dim and y-dim distributions
rho.xy <- 0.1
mu <- c(mu.x, mu.y)
#Covariance matrix
Sigmaxy <- matrix(c(sigmax^2, rep(rho.xy*sigmax*sigmay, 2), sigmay^2), ncol=2)

w <- dmvnorm(gr, mean=mu, sigma=Sigmaxy)
w <- (w - mean(w))/sd(w)
#Visualize PO detection
image(rasterFromXYZ(cbind(gr,w)), col=topo.colors(20))

#----------------------------------#
#-Simulate opportunistic surveying-#
#----------------------------------#

#Detection probability of PO
ppo <- expit(alpha1*w + alpha0)

#Individuals detected in PO
ypo <- NULL
for(i in 1:N){
  ypo[i] <- rbinom(1, 1, ppo[s[i]])
}

#Pixel ID for true presence
pixpo <- s[ypo == 1]

#Coord of true presence
uxpo <- u1[ypo == 1]
uypo <- u2[ypo == 1]

# #Unobserved sampling / PO pixel ID
# error <- 0.8
# pixpo <- sample(pixpo, length(pixpo)*error, replace = FALSE)

#Number of PO detections per pixel
ypo <- as.data.frame(table(pixpo))
ypo$pixpo <- as.numeric(as.character(ypo$pixpo))

#Vector of pixels with PO detections
tmp <- rep(0, (W*W))
for(i in 1:length(ypo[,1])){
  tmp[ypo$pixpo[i]] <- ypo$Freq[i]
}

ypo <- tmp

#---------------#
#-Visualization-#
#---------------#

for(i in 1:nts){
  image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))
  points(u1, u2, col = 'black', pch = 20)
  for(j in 1:ns[i]){
    points(line[which((line[,1] > (grt[tsamp[j,i],1] - 5)) & 
                        ((grt[tsamp[j,i],1] + 5) > line[,1]) & 
                        (line[,2] > (grt[tsamp[j,i],2] - 5)) & 
                        ((grt[tsamp[j,i],2] + 5) > line[,2])),], col = 'grey', pch = 20)
  }
  points(YD[[i]][YD[[i]][,3]==1,4:5], col = 'orange', pch = 20)
  points(uxpo, uypo, col = "red", pch = 20)
}

#----------------------------------------------#
#-Compile BUGS data for each sampling scenario-#
#----------------------------------------------#

#Scenario 1: No DS
str(data1 <- list(x = x, G = (W*W),
                  w = w,  y.po = ypo))

#Scenario 2: 1 Transect
str(data2 <- list(x = x, G = (W*W), B = B[7], dst = yds[[7]]$dst,
                  y.ds = yds[[7]]$freq, coverage = yds[[7]]$p,
                  w = w, y.po = ypo))

#Scenario 3: 5 Transects
str(data3 <- list(x = x, G = (W*W), B = B[6], dst = yds[[6]]$dst,
                  y.ds = yds[[6]]$freq, coverage = yds[[6]]$p,
                  w = w, y.po = ypo))

#Scenario 4: 10 Transects
str(data4 <- list(x = x, G = (W*W), B = B[5], dst = yds[[5]]$dst,
                  y.ds = yds[[5]]$freq, coverage = yds[[5]]$p,
                  w = w, y.po = ypo))

#Scenario 5: 15 Transects
str(data5 <- list(x = x, G = (W*W), B = B[4], dst = yds[[4]]$dst,
                   y.ds = yds[[4]]$freq, coverage = yds[[4]]$p,
                   w = w, y.po = ypo))

#Scenario 6: 20 Transects
str(data6 <- list(x = x, G = (W*W), B = B[3], dst = yds[[3]]$dst,
                   y.ds = yds[[3]]$freq, coverage = yds[[3]]$p,
                   w = w, y.po = ypo))

#Scenario 7: 25 Transects
str(data7 <- list(x = x, G = (W*W), B = B[2], dst = yds[[2]]$dst,
                   y.ds = yds[[2]]$freq, coverage = yds[[2]]$p,
                   w = w, y.po = ypo))

#Scenario 8: 50 Transects
str(data8 <- list(x = x, G = (W*W), B = B[1], dst = yds[[1]]$dst,
                  y.ds = yds[[1]]$freq, coverage = yds[[1]]$p,
                  w = w, y.po = ypo))

#----------------#
#-Initial values-#
#----------------#

#Scenario 1: No DS
inits1 <- function(){list(beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25), 
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 2: 1 Transect
inits2 <- function(){list(N = yds[[7]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                           beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                           alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 3: Distance sampling robust
inits3 <- function(){list(N = yds[[6]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                           beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                           alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 4: Distance sampling robust
inits4 <- function(){list(N = yds[[5]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 5: ISDM robust DS & PO
inits5 <- function(){list(N = yds[[4]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 6: ISDM sparse DS & PO
inits6 <- function(){list(N = yds[[3]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 7: ISDM robust DS & sparse PO
inits7 <- function(){list(N = yds[[2]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 8: ISDM sparse DS & robust PO
inits8 <- function(){list(N = yds[[1]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#------------#
#-Parameters-#
#------------#

params <- c("Ntot", "sigma", "alpha0", "alpha1", "beta0", "beta1")

#-------------#
#-MCMC values-#
#-------------#

nb <- 2000
ni <- 8000
nt <- 2
nc <- 3
na <- 1000

start <- Sys.time()
#----------------#
#-Run each model-#
#----------------#

#Scenario 1: No DS
S1 <- jagsUI(data1, inits1, params, "PO.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S1, file="S1.Rdata")

#Scenario 2: 1 Transect
S2 <- jagsUI(data2, inits2, params, "ISDMalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S2, file="S2.Rdata")

#Scenario 3: Distance sampling robust
S3 <- jagsUI(data3, inits3, params, "ISDMalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S3, file="S3.Rdata")

#Scenario 4: Distance sampling robust
S4 <- jagsUI(data4, inits4, params, "ISDMalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S4, file="S4.Rdata")

#Scenario 5: ISDM robust DS & PO
S5 <- jagsUI(data5, inits5, params, "ISDMalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S5, file="S5.Rdata")

#Scenario 6: ISDM sparse DS & PO
S6 <- jagsUI(data6, inits6, params, "ISDMalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S6, file="S6.Rdata")

#Scenario 7: ISDM robust DS & sparse PO
S7 <- jagsUI(data7, inits7, params, "ISDMalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S7, file="S7.Rdata")

#Scenario 8: ISDM sparse DS & robust PO
S8 <- jagsUI(data8, inits8, params, "ISDMalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S8, file="S8.Rdata")

###

end <- Sys.time()


#######
#Scenario 2: 1 Transect
str(data9 <- list(x = x, G = (W*W), B = B[7], dst = yds[[7]]$dst,
                  y.ds = yds[[7]]$freq, coverage = yds[[7]]$p))

#Scenario 3: 5 Transects
str(data10 <- list(x = x, G = (W*W), B = B[6], dst = yds[[6]]$dst,
                  y.ds = yds[[6]]$freq, coverage = yds[[6]]$p))

#Scenario 4: 10 Transects
str(data11 <- list(x = x, G = (W*W), B = B[5], dst = yds[[5]]$dst,
                  y.ds = yds[[5]]$freq, coverage = yds[[5]]$p))

#Scenario 5: 15 Transects
str(data12 <- list(x = x, G = (W*W), B = B[4], dst = yds[[4]]$dst,
                  y.ds = yds[[4]]$freq, coverage = yds[[4]]$p))

#Scenario 6: 20 Transects
str(data13 <- list(x = x, G = (W*W), B = B[3], dst = yds[[3]]$dst,
                  y.ds = yds[[3]]$freq, coverage = yds[[3]]$p))

#Scenario 7: 25 Transects
str(data14 <- list(x = x, G = (W*W), B = B[2], dst = yds[[2]]$dst,
                  y.ds = yds[[2]]$freq, coverage = yds[[2]]$p))

#Scenario 8: 50 Transects
str(data15 <- list(x = x, G = (W*W), B = B[1], dst = yds[[1]]$dst,
                  y.ds = yds[[1]]$freq, coverage = yds[[1]]$p))

#----------------#
#-Initial values-#
#----------------#

#Scenario 2: 1 Transect
inits9 <- function(){list(N = yds[[7]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}

#Scenario 3: Distance sampling robust
inits10 <- function(){list(N = yds[[6]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}

#Scenario 4: Distance sampling robust
inits11 <- function(){list(N = yds[[5]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}

#Scenario 5: ISDM robust DS & PO
inits12 <- function(){list(N = yds[[4]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}

#Scenario 6: ISDM sparse DS & PO
inits13 <- function(){list(N = yds[[3]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}

#Scenario 7: ISDM robust DS & sparse PO
inits14 <- function(){list(N = yds[[2]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}

#Scenario 8: ISDM sparse DS & robust PO
inits15 <- function(){list(N = yds[[1]]$freq+1, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}



#Scenario 2: 1 Transect
S9 <- jagsUI(data9, inits9, params, "DSalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S2, file="S2.Rdata")

#Scenario 3: Distance sampling robust
S10 <- jagsUI(data10, inits10, params, "DSalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S3, file="S3.Rdata")

#Scenario 4: Distance sampling robust
S11 <- jagsUI(data11, inits11, params, "DSalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S4, file="S4.Rdata")

#Scenario 5: ISDM robust DS & PO
S12 <- jagsUI(data14, inits14, params, "DSalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S5, file="S5.Rdata")

#Scenario 6: ISDM sparse DS & PO
S13 <- jagsUI(data13, inits13, params, "DSalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S6, file="S6.Rdata")

#Scenario 7: ISDM robust DS & sparse PO
S14 <- jagsUI(data14, inits14, params, "DSalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S7, file="S7.Rdata")

#Scenario 8: ISDM sparse DS & robust PO
S15 <- jagsUI(data15, inits15, params, "DSalt2.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
# save(S8, file="S8.Rdata")













