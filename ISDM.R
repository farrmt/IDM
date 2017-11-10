#---------------------------------------------------#
#----Integrated Species Distribution Model----------#
#----Combines presence only (PO) and distance-------#
#----sampling data (DS). Buids off of simulation----#
#----in Dorazio (2014).-----------------------------#
#----Created by Matthew Farr------------------------#
#---------------------------------------------------#

#----------------#
#-Load libraries-#
#----------------#
library(mvtnorm)
library(rgeos)
library(sp)
library(raster)
library(jagsUI)

#----------------------------#
#-Define necessary functions-#
#----------------------------#

expit = function(eta) {1/(1+exp(-eta))}

#--------------------------#
#-Create sampling region B-#
#--------------------------#

#Width/length of area
W <- 100
#Pixel size
px <- 1
#Mid-point coordinates
grx <- seq(px/2, W - px/2, px)
#Create grid coordinates for sampling region
gr <- expand.grid(grx,grx)

#-------------------------#
#-Parameters' true values-#
#-------------------------#

#Intercept of intensity
beta0 <- - 2.5
#Effet of environment on intensity
beta1 <- 0.25
#Intercept of environment on PO detection
alpha0 <- - 5
#Effect of environment on PO detection
alpha1 <- 5
#Scale parameter of half-normal function for DS detection
sigma <- 2

#--------------------------#
#-Environmental covariates-#
#--------------------------#

#Simulate covariate on intensity (multivariate normal)
mu1.x = 0.75*(100)
mu1.y = 0.40*(100)
sigma1.x = 0.25*abs(100)
sigma1.y = 0.50*abs(100)
rho1.xy = 0.5
mu1 = c(mu1.x, mu1.y)
Sigma1 = matrix(c(sigma1.x^2, rep(rho1.xy*sigma1.x*sigma1.y, 2), sigma1.y^2), ncol=2)

mu2.x = 0.15*(100)
mu2.y = 0.80*(100)
sigma2.x = 0.50*abs(100)
sigma2.y = 0.25*abs(100)
rho2.xy = -0.4
mu2 = c(mu2.x, mu2.y)
Sigma2 = matrix(c(sigma2.x^2, rep(rho2.xy*sigma2.x*sigma2.y, 2), sigma2.y^2), ncol=2)

xcov = 0.4 * dmvnorm(gr, mean=mu1, sigma=Sigma1) + 0.6 * dmvnorm(gr, mean=mu2, sigma=Sigma2)
xcov = (xcov - mean(xcov))/sd(xcov)

#Visualize covariate on PO detection
image(rasterFromXYZ(cbind(gr,xcov)), col=topo.colors(10))

#Simulate covariate on PO detection (Multivariate normal)
mu.x <- 0.2*(100)
mu.y <- 0.1*(100)
sigma.x <- 0.3*abs(100)
sigma.y <- 0.15*abs(100)
rho.xy <- 0.1
mu <- c(mu.x, mu.y)
Sigmaxy <- matrix(c(sigma.x^2, rep(rho.xy*sigma.x*sigma.y, 2), sigma.y^2), ncol=2)

wcov <- dmvnorm(gr, mean=mu, sigma=Sigmaxy)
wcov <- (wcov - mean(wcov))/sd(wcov)

#Visualize covariate on PO detection
image(rasterFromXYZ(cbind(gr,wcov)), col=topo.colors(10))

#--------------------------------#
#-Simulate species' distribution-#
#--------------------------------#

#Probability of point in pixel
probs <- exp(beta1 * xcov + beta0)/sum(exp(beta1 * xcov + beta0))
#Simulate number of individuals in area
N <- rpois(1, sum(exp(beta1 * xcov + beta0)))
#Sample location of simulated individuals
pixel.id <- sample(1:(W*W), N, replace=TRUE, prob=probs)
#Assign X and Y cord
u1 <- gr[pixel.id,1]
u2 <- gr[pixel.id,2]

#-------------------#
#-Distance sampling-#
#-------------------#

#Transect lines for distance sampling
line <- structure(list(x = rep(c(25,50,75), each = W), y = rep(seq(0.5, W-0.5, 1), 3), .Names = c("x", "y")))
line <- cbind(line$x, line$y)
line <- as.matrix(line)
#X Coord
X <- line[,1]
#Y Coord
Y <- line[,2]
#Number of points in transect lines
J <- length(line[,1])
#Distance class length
v <- 1
#Transect half-width
B <- 12
#Midpoint of distance class
mdpt <- seq(0.5, 12, 1)
#Number of distance class
nD <- length(mdpt)
#Distance array for all points on transect lines
d <- array(NA, dim = c(W*W, J))
#Distance to nearest transect (3 transect)
dist <- NULL
#Distance to nearest transect (1 transect)
#dist1 <- NULL
for(g in 1:10000){
  for(j in 1:J){
    d[g,j] <- sqrt((gr[g,1] - X[j])^2 + (gr[g,2] - Y[j])^2)
  }
  dist[g] <- min(d[g,])
  #dist1[g] <- min(d[g,101:200])
}

#Detection probability for distance sampling
p.ds <- rep(NA, N)
#Distance for each individual
dst <- dist[pixel.id]
#Individual presence/absence
y.ds <- NULL
for(i in 1:N){
  p.ds[i] <- exp(-dst[i] * dst[i] / (2 * sigma * sigma))
  y.ds[i] <- rbinom(1, 1, p.ds[i])
}

#Coord of detected individuals
uxds <- u1[y.ds == 1]
uyds <- u2[y.ds == 1]
#Distance class of detected individuals
dst <- dst[y.ds == 1]
#Pixels covered by distance sampling (transect half-width is B:12 pixels)
coverage <- which(dist <= B)
#Distance class of detected individuals
dclass <- as.numeric(dst + 0.5)
#Number of detected individuals
y.ds <- rep(1, sum(y.ds))

#If using only 1 transect line instead of 3:
# uxds1 <- uxds[uxds > 40 & uxds < 60]
# uyds1 <- uyds[uxds > 40 & uxds < 60]
# dst1 <- dst[uxds > 40 & uxds < 60]
# coverage1 <- which(dist1 <= 12)
# dclass1 <- as.numeric(dst1 + 0.5)
# y.ds1 <- rep(1, length(dclass1))

#------------------------#
#-Presence only sampling-#
#------------------------#

#Detection probability of PO
p.po <- expit(alpha1*wcov + alpha0)

#Individuals detected in PO
y.po <- NULL
for(i in 1:N){
  y.po[i] <- rbinom(1, 1, p.po[pixel.id[i]])
}

#Coord of detected individuals
uxpo <- u1[y.po == 1]
uypo <- u2[y.po == 1]

#Pixel ID for PO
pixpo <- pixel.id[y.po == 1]
y.po <- as.data.frame(table(pixpo))
y.po$pixpo <- as.numeric(as.character(y.po$pixpo))

#Vector of pixels with PO detections
tmp <- rep(0, 10000)
for(i in 1:length(y.po[,1])){
  tmp[y.po$pixpo[i]] <- y.po$Freq[i]
}
y.po <- tmp

#------------------------#
#-Sampling visualization-#
#------------------------#

par(mar=c(3,3,3,6))
#Intensity of individuals in region B
image(rasterFromXYZ(cbind(gr,xcov)), col=topo.colors(10))
#Location of N individuals in region B
points(u1, u2, pch=20, col='black', cex = 0.8)
#Transect lines
lines(X, Y, col = 'grey')
#Detected individuals in DS for 3 transects
points(uxds, uyds, col = 'orange', pch = 20)
#Detected individuals in DS for 1 transects
#points(uxds1, uyds1, col = 'blue', pch = 20)
#Detected individuals in PO
points(uxpo, uypo, col = 'red', pch = 20)

#---------------------#
#-JAGS model for ISDM-#
#---------------------#

cat("
    model{ 
    
    #--------#
    #-Priors-#
    #--------#

    #Scale parameter DS
    sigma ~ dunif(0, 10)
    #Intercept intensity
    beta0 ~ dnorm(0, 0.01)
    #Effect intensity
    beta1 ~ dnorm(0, 0.01)
    #Intercept PO detection
    alpha0 ~ dnorm(0, 0.01)
    #Effect PO detection
    alpha1 ~ dnorm(0, 0.01)
    
    #------------#
    #-Likelihood-#
    #------------#

    for(g in 1:G){
    #Intensity of individuals
    intensity[g] <- exp(beta0 + beta1 * xcov[g])
    #Detection for PO
    logit(p.po[g]) <- alpha0 + wcov[g] * alpha1
    #Thinning of intensity
    v0[g] <- intensity[g] * p.po[g]
    #Observation process of PO
    y.po[g] ~ dpois(v0[g])
    }

    
    #Detection for DS
    for(k in 1:nD){
    p[k] <- exp(-mdpt[k]*mdpt[k]/(2*sigma*sigma))
    pi[k] <- v/B
    f[k] <- p[k] * pi[k]
    fc[k] <- f[k] / p.ds
    }
    
    p.ds <- sum(f[])

    for(i in 1:nds){
    dclass[i] ~ dcat(fc[1:nD])
    }
    
    #Observation process for DS
    y.ds ~ dbin(p.ds, N)
    N ~ dpois(lambda)
    lambda <- sum(intensity[coverage[]])

    #Number of individuals in region B
    Ntot <- sum(intensity[])

    }
    ",fill=TRUE,file="ISDM.txt")

#--------------#
#-Compile data-#
#--------------#

str(ISDMdata <- list(xcov = xcov, G = 10000,
                 nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass, nds = sum(y.ds),
                 y.ds = sum(y.ds), coverage = coverage,
                 wcov = wcov,  y.po = y.po))

#----------------#
#-Initial values-#
#----------------#

N.st <- sum(y.ds) + 1

ISDMinits <- function(){list(N = N.st, sigma=runif(1,1,3), beta1=runif(1,0,1), beta0=runif(1,-2,2), alpha0=runif(1,-10,0), alpha1 = runif(1,5,10))}

#--------------------#
#-Parameters to save-#
#--------------------#

ISDMparams <- c("beta0", "beta1", "alpha0", "alpha1", "Ntot", "lambda")

#---------------#
#-MCMC settings-#
#---------------#

nb <- 100000
ni <- 120000
nt <- 10
nc <- 3
na <- 1000

#----------------#
#-Run JAGS model-#
#----------------#

ISDM <- jagsUI(ISDMdata, ISDMinits, ISDMparams, "ISDM.txt", n.thin=nt, 
               n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#-------------------#
#-JAGS model for DS-#
#-------------------#

cat("
    model{ 
    
    #--------#
    #-Priors-#
    #--------#
    
    #Scale parameter DS
    sigma ~ dunif(0, 10)
    #Intercept intensity
    beta0 ~ dnorm(0, 0.01)
    #Effect intensity
    beta1 ~ dnorm(0, 0.01)
    
    #------------#
    #-Likelihood-#
    #------------#
    
    for(g in 1:G){
    #Intensity of individuals
    intensity[g] <- exp(beta0 + beta1 * xcov[g])
    }
    
    
    #Detection for DS
    for(k in 1:nD){
    p[k] <- exp(-mdpt[k]*mdpt[k]/(2*sigma*sigma))
    pi[k] <- v/B
    f[k] <- p[k] * pi[k]
    fc[k] <- f[k] / p.ds
    }
    
    p.ds <- sum(f[])
    
    for(i in 1:nds){
    dclass[i] ~ dcat(fc[1:nD])
    }
    
    #Observation process for DS
    y.ds ~ dbin(p.ds, N)
    N ~ dpois(lambda)
    lambda <- sum(intensity[coverage[]])
    
    #Number of individuals in region B
    Ntot <- sum(intensity[])
    
    }
    ",fill=TRUE,file="DS.txt")

#--------------#
#-Compile data-#
#--------------#

str(DSdata <- list(xcov = xcov, G = 10000,
                 nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass, nds = sum(y.ds),
                 y.ds = sum(y.ds), coverage = coverage))

#----------------#
#-Initial values-#
#----------------#

N.st <- sum(y.ds) + 1

DSinits <- function(){list(N = N.st, sigma=runif(1,1,3), beta1=runif(1,0,1), beta0=runif(1,-2,2))}

#--------------------#
#-Parameters to save-#
#--------------------#

DSparams <- c("beta0", "beta1", "Ntot", "lambda")

#---------------#
#-MCMC settings-#
#---------------#

#Same as above

#----------------#
#-Run JAGS model-#
#----------------#

DS <- jagsUI(DSdata, DSinits, DSparams, "DS.txt", n.thin=nt, 
               n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#-------------------#
#-JAGS model for PO-#
#-------------------#

cat("
    model{ 
    
    #--------#
    #-Priors-#
    #--------#
    
    #Intercept intensity
    beta0 ~ dnorm(0, 0.01)
    #Effect intensity
    beta1 ~ dnorm(0, 0.01)
    #Intercept PO detection
    alpha0 ~ dnorm(0, 0.01)
    #Effect PO detection
    alpha1 ~ dnorm(0, 0.01)
    
    #------------#
    #-Likelihood-#
    #------------#
    
    for(g in 1:G){
    #Intensity of individuals
    intensity[g] <- exp(beta0 + beta1 * xcov[g])
    #Detection for PO
    logit(p.po[g]) <- alpha0 + wcov[g] * alpha1
    #Thinning of intensity
    v0[g] <- intensity[g] * p.po[g]
    #Observation process of PO
    y.po[g] ~ dpois(v0[g])
    }
    
    
    #Detection for DS
    for(k in 1:nD){
    p[k] <- exp(-mdpt[k]*mdpt[k]/(2*sigma*sigma))
    pi[k] <- v/B
    f[k] <- p[k] * pi[k]
    fc[k] <- f[k] / p.ds
    }
    
    #Number of individuals in region B
    Ntot <- sum(intensity[])
    
    }
    ",fill=TRUE,file="PO.txt")

#--------------#
#-Compile data-#
#--------------#

str(POdata <- list(xcov = xcov, G = 10000, wcov = wcov,  y.po = y.po))

#----------------#
#-Initial values-#
#----------------#

POinits <- function(){list(beta1=runif(1,0,1), beta0=runif(1,-2,2), alpha0=runif(1,-10,0), alpha1 = runif(1,5,10))}

#--------------------#
#-Parameters to save-#
#--------------------#

POparams <- c("beta0", "beta1", "alpha0", "alpha1")

#---------------#
#-MCMC settings-#
#---------------#

#Same as above

#----------------#
#-Run JAGS model-#
#----------------#

PO <- jagsUI(POdata, POinits, POparams, "PO.txt", n.thin=nt, 
               n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)