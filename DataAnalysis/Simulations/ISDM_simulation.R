#NOTES:
#1 change PO sampling to sample at different rates check koshkina [check]
#2 check PO with baseline (alpha0) at 100% / idea for figure [check]
#3 reduce DS.S

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
# beta0 <- log(0.05) #Low Abundance (500 individuals; 0.01/pixel)
beta0 <- log(0.5) #Medium Abundance (5,000 individuals; 1/pixel)
# beta0 <- log(5) #High Abundance (50,000 individuals; 10/pixel)

#Effect parameter of enivornment on intensity
# beta1 <- 0 #No effect on intensity (i.e., homogenous)
# beta1 <- 0.5 #Weak effect on intensity X 1.14 (14% increase)
# beta1 <- 1 #Meadium effect on intensity X 1.66 (66% inrease)
beta1 <- 1.25 #Strong effect on intensity X 2.19 (119% increase)

#Intercept parameter for prensence only (PO) detection
# alpha0 <- logit(0.1) #Low detection (10%)
# alpha0 <- logit(0.5) #Medium detection (50%)
alpha0 <- logit(0.99) #High detection (100%)

#Effect parameter of environment on PO detection
# alpha1 <- 0 #No effect on PO detection
# alpha1 <- 2.25 #Weak effect on PO detection
# alpha1 <- 5 #Medium effect on PO detection
alpha1 <- 10 #Strong effect on PO detection

#Unobserved sampling error
error.R <- 0.5
error.S <- 0.3

#Scale parameter for DS detection
# sigma <- 1 #Low detection (10%)
sigma <- 2.5 #Medium detection (26%)
# sigma <- 5 #High detection (50%)

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
image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))

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

#----------------------------#
#-Simulate distance sampling-#
#----------------------------#

#Alternate transects
x.line.R <- rep(c(25,50,75), each = W) #Robust
y.line.R <- rep(seq(0.5, W-0.5, 1), 3) #Robust
x.line.S <- rep(c(25, 75), each = 60) #Sparse
y.line.S <- rep(c(seq(11, 40, 1), seq(61, 90, 1)), 2) #Sparse
#Visualize transect lines
points(x.line.R, y.line.R, col = 'grey', pch = 20)

image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))
points(u1, u2, col = 'black', pch = 20)
points(x.line.S, y.line.S, col = 'grey', pch = 20)

#Number of points in transect line
# J.R <- length(x.line.R)
# J.S <- length(x.line.S)
#Distance array for all points on transect line
# d.R <- array(NA, dim = c(W*W, J.R))
# d.S <- array(NA, dim = c(W*W, J.S))

#Distance to nearest transect
# dist.R <- NULL
# dist.S <- NULL
#Simulate above quantities
# for(g in 1:(W*W)){
#   for(jr in 1:J.R){
#     d.R[g,jr] <- sqrt((gr[g,1] - x.line.R[jr])^2 + (gr[g,2] - y.line.R[jr])^2)
#   }
#   for(js in 1:J.S){
#     d.S[g,js] <- sqrt((gr[g,1] - x.line.S[js])^2 + (gr[g,2] - y.line.S[js])^2)
#   }
#   dist.R[g] <- min(d.R[g,])
#   dist.S[g] <- min(d.S[g,])
# }
#Load distance to pixels for robust and sparse
load(file = "dist.R.Rdata")
load(file = "dist.S.Rdata")

#Detection probability for distance sampling
pds.R <- rep(NA, N)
pds.S <- rep(NA, N)
#Distance for each individual
dst.R <- dist.R[s]
dst.S <- dist.S[s]
#Individual presence/absence based on distance sampling
yds.R <- NULL
yds.S <- NULL
#Simulate distance sampling
pds.R <- exp(-dst.R * dst.R / (2 * sigma * sigma))
pds.S <- exp(-dst.S * dst.S / (2 * sigma * sigma))
yds.R <- rbinom(N, 1, pds.R)
yds.S <- rbinom(N, 1, pds.S)
yds.R[dst.R>12] <- 0
yds.S[dst.S>12] <- 0

#Coordinates of detected individuals
uxds.R <- u1[yds.R == 1]
uyds.R <- u2[yds.R == 1]
uxds.S <- u1[yds.S == 1]
uyds.S <- u2[yds.S == 1]
#Distance of detected individuals
dst.R <- dst.R[yds.R == 1]
dst.S <- dst.S[yds.S == 1]
#Pixel ID
pixds.R <- s[yds.R == 1]
pixds.S <- s[yds.S == 1]
#Distance class length
v <- 1
#Transect half-width
B <- 12
#Midpoint of distance class
mdpt <- seq(0.5, 13, 1)
#Number of distance class
nD <- length(mdpt)-1
#Pixels covered by distance sampling (transect half-width is 10 pixels)
coverage.R <- which(dist.R <= B)
coverage.S <- which(dist.S <= B)
#Distance class of detected individuals
dclass.R <- NULL
dclass.S <- NULL
for(i in 1:sum(yds.R)){
  for(k in 1:nD){
    if(mdpt[k] - 0.5 <= dst.R[i] && dst.R[i] < mdpt[k+1] - 0.5)
      dclass.R[i] <- k
  }
}
for(i in 1:sum(yds.S)){
  for(k in 1:nD){
    if(mdpt[k] - 0.5 <= dst.S[i] && dst.S[i] < mdpt[k+1] - 0.5)
      dclass.S[i] <- k
  }
}

regionC.R <- as.data.frame(coverage.R)
colnames(regionC.R) <- "pixel"
regionC.S <- as.data.frame(coverage.S)
colnames(regionC.S) <- "pixel"
C.R <- length(coverage.R)
C.S <- length(coverage.S)
dst.R <- dist.R[coverage.R] 
dst.S <- dist.S[coverage.S] 
yds.R <- as.data.frame(table(pixds.R))
yds.S <- as.data.frame(table(pixds.S))
yds.R$pixds.R <- as.numeric(as.character(yds.R$pixds.R))
colnames(yds.R) <- c("pixel", "freq")
tmp <- full_join(yds.R, regionC.R, by = "pixel")
tmp <- tmp%>%arrange(pixel)
yds.R <- tmp$freq
yds.R[is.na(yds.R)] <- 0
yds.S$pixds.S <- as.numeric(as.character(yds.S$pixds.S))
colnames(yds.S) <- c("pixel", "freq")
tmp <- full_join(yds.S, regionC.S, by = "pixel")
tmp <- tmp%>%arrange(pixel)
yds.S <- tmp$freq
yds.S[is.na(yds.S)] <- 0

#Visualize detected individuals from distance sampling
points(uxds.S, uyds.S, col = 'orange', pch = 20)

image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))
points(u1, u2, col = 'black', pch = 20)
points(x.line.R, y.line.R, col = 'grey', pch = 20)
points(uxds.R, uyds.R, col = 'orange', pch = 20)

#---------------------------------------#
#-Draw covariate value for PO detection-#
#---------------------------------------#
#Environmental covariate on PO detection (Multivariate normal)
#Mean of x-dim distribution
mu.x <- 50
#Mean of y-dim distribution
mu.y <- 50
#Variance of x-dim distribution
sigmax <- 0.5*abs(W) #Robust PO
#Variance of y-dim distribution
sigmay <- 0.5*abs(W)
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
ypa <- NULL
for(i in 1:N){
  ypa[i] <- rbinom(1, 1, ppo[s[i]])
}

#Pixel ID for true presence
pixpa <- s[ypa == 1]

#Coord of true presence
uxpa <- u1[ypa == 1]
uypa <- u2[ypa == 1]

#Unobserved sampling / PO pixel ID
pixpo.R <- sample(pixpa, length(pixpa)*error.R, replace = FALSE)
pixpo.S <- sample(pixpa, length(pixpa)*error.S, replace = FALSE)

#Coord of PO
uxpo.R <- u1[which(s%in%pixpo.R)]
uypo.R <- u2[which(s%in%pixpo.R)]
uxpo.S <- u1[which(s%in%pixpo.S)]
uypo.S <- u2[which(s%in%pixpo.S)]

#Number of PO detections per pixel
ypo.R <- as.data.frame(table(pixpo.R))
ypo.S <- as.data.frame(table(pixpo.S))
ypo.R$pixpo.R <- as.numeric(as.character(ypo.R$pixpo.R))
ypo.S$pixpo.S <- as.numeric(as.character(ypo.S$pixpo.S))
#Vector of pixels with PO detections
tmp <- rep(0, (W*W))
for(i in 1:length(ypo.R[,1])){
  tmp[ypo.R$pixpo[i]] <- ypo.R$Freq[i]
}
ypo.R <- tmp
tmp <- rep(0, (W*W))
for(i in 1:length(ypo.S[,1])){
  tmp[ypo.S$pixpo[i]] <- ypo.S$Freq[i]
}
ypo.S <- tmp
#Visualize covariate
image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))
points(u1, u2, col = 'black', pch = 20)
points(uxpa, uypa, col = "orange", pch = 20)
points(uxpo.R, uypo.R, col = 'red', pch = 20)

image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))
points(u1, u2, col = 'black', pch = 20)
points(uxpa, uypa, col = "orange", pch = 20)
points(uxpo.S, uypo.S, col = 'red', pch = 20)

#----------------------------------------------#
#-Compile BUGS data for each sampling scenario-#
#----------------------------------------------#

#Scenario 1: Presence only robust
str(data1 <- list(x = x, G = (W*W),
                   w = w.R,  y.po = ypo.R))

#Scenario 2: Presence only sparse
str(data2 <- list(x = x, G = (W*W),
                   w = w.S,  y.po = ypo.S))

#Scenario 3: Distance sampling robust
str(data3 <- list(x = x, G = (W*W), C = C.R, dst = dst.R,
                  nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass.R, nds = sum(yds.R),
                  y.ds = yds.R, coverage = coverage.R))

#Scenario 4: Distance sampling sparse
str(data4 <- list(x = x, G = (W*W), C = C.S, dst = dst.S,
                  nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass.S, nds = sum(yds.S),
                  y.ds = yds.S, coverage = coverage.S))

#Scenario 5: ISDM robust DS & PO
str(data5 <- list(x = x, G = (W*W), C = C.R, dst = dst.R,
                 nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass.R, nds = sum(yds.R), 
                 y.ds = yds.R, coverage = coverage.R, w = w.R, y.po = ypo.R))

#Scenario 6: ISDM sparse DS & PO
str(data6 <- list(x = x, G = (W*W), C = C.S, dst = dst.S,
                     nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass.S, nds = sum(yds.S), 
                     y.ds = yds.S, coverage = coverage.S, w = w.S, y.po = ypo.S))

#Scenario 7: ISDM robust DS & sparse PO
str(data7 <- list(x = x, G = (W*W), C = C.R, dst = dst.R,
                     nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass.R, nds = sum(yds.R), 
                     y.ds = yds.R, coverage = coverage.R, w = w.S, y.po = ypo.S))

#Scenario 8: ISDM sparse DS & robust PO
str(data8 <- list(x = x, G = (W*W), C = C.S, dst = dst.S,
                     nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass.S, nds = sum(yds.S), 
                     y.ds = yds.S, coverage = coverage.S, w = w.R, y.po = ypo.R))

#----------------#
#-Initial values-#
#----------------#

#Inital value for N for distance sampling
Nst.R <- yds.R + 1 #robust
Nst.S <- yds.S + 1 #sparse

#Scenario 1: Presence only robust
inits1 <- function(){list(beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25), 
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 2: Presence only sparse
inits2 <- function(){list(beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25), 
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 3: Distance sampling robust
inits3 <- function(){list(N = Nst.R, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}

#Scenario 4: Distance sampling robust
inits4 <- function(){list(N = Nst.S, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25))}

#Scenario 5: ISDM robust DS & PO
inits5 <- function(){list(N = Nst.R, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 6: ISDM sparse DS & PO
inits6 <- function(){list(N = Nst.S, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 7: ISDM robust DS & sparse PO
inits7 <- function(){list(N = Nst.R, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#Scenario 8: ISDM sparse DS & robust PO
inits8 <- function(){list(N = Nst.S, sigma = runif(1, sigma-0.25, sigma+0.25), 
                          beta1 = runif(1, beta1-0.25, beta1+0.25), beta0 = runif(1, beta0-0.25, beta0+0.25),
                          alpha1 = runif(1, alpha1-0.25, alpha1+0.25), alpha0 = runif(1, alpha0-0.25, alpha0+0.25))}

#------------#
#-Parameters-#
#------------#

params <- c("Ntot", "sigma", "alpha0", "alpha1", "beta0", "beta1")

#-------------#
#-MCMC values-#
#-------------#

nb <- 1000
ni <- 7000
nt <- 2
nc <- 3
na <- 1000

start <- Sys.time()
#----------------#
#-Run each model-#
#----------------#

#Scenario 1: Presence only robust
S1 <- jagsUI(data1, inits1, params, "PO.txt", n.thin=nt, 
              n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S1, file="S1.Rdata")

#Scenario 2: Presence only sparse
S2 <- jagsUI(data2, inits2, params, "PO.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S2, file="S2.Rdata")

#Scenario 3: Distance sampling robust
S3 <- jagsUI(data3, inits3, params, "DSalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S3, file="S3.Rdata")

#Scenario 4: Distance sampling robust
S4 <- jagsUI(data4, inits4, params, "DSalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S4, file="S4.Rdata")

#Scenario 5: ISDM robust DS & PO
S5 <- jagsUI(data5, inits5, params, "ISDMalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S5, file="S5.Rdata")

#Scenario 6: ISDM sparse DS & PO
S6 <- jagsUI(data6, inits6, params, "ISDMalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S6, file="S6.Rdata")

#Scenario 7: ISDM robust DS & sparse PO
S7 <- jagsUI(data7, inits7, params, "ISDMalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S7, file="S7.Rdata")

#Scenario 8: ISDM sparse DS & robust PO
S8 <- jagsUI(data8, inits8, params, "ISDMalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)
save(S8, file="S8.Rdata")

end <- Sys.time()

# image(rasterFromXYZ(cbind(gr, exp(beta0 + beta1*x))), col = topo.colors(20))
# image(rasterFromXYZ(cbind(gr, exp(DS$mean$beta0 + DS$mean$beta1*x))), col = topo.colors(20))
# image(rasterFromXYZ(cbind(gr, exp(PO$mean$beta0 + PO$mean$beta1*x))), col = topo.colors(20))
# image(rasterFromXYZ(cbind(gr, exp(ISDM$mean$beta0 + ISDM$mean$beta1*x))), col = topo.colors(20))
