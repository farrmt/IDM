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

expit <- function(eta) {1/(1+exp(-eta))}

#--------------------------#
#-Create sampling region B-#
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
beta0 <- - 1
#Effect parameter of enivornment on intensity
beta1 <- 1

#Intercept parameter for prensence only (PO) detection
alpha0 <- - 5
#Effect parameter of environment on PO detection
alpha1 <- 7.5

#Unobserved sampling error
error <- runif(1, 0.3, 0.4)

#Scale parameter of half-normal function for DS detection
sigma <- 2

#--------------------------#
#-Environmental covariates-#
#--------------------------#

#Simulate covariate on intensity (bivariate normal)
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

x = 0.4 * dmvnorm(gr, mean=mu1, sigma=Sigma1) + 0.6 * dmvnorm(gr, mean=mu2, sigma=Sigma2)
x = (x - mean(x))/sd(x)

#Visualize covariate on PO detection
image(rasterFromXYZ(cbind(gr,xcov)), col=topo.colors(10))

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





#-------------------#
#-Distance sampling-#
#-------------------#

#Transect lines for distance sampling
x.line <- rep(c(25,50,75), each = W)
y.line <- rep(seq(0.5, W-0.5, 1), 3)

#Visualize transect lines
points(x.line, y.line, col = 'grey', pch = 20)

#Number of points in transect line
J <- length(x.line)
#Distance array for all points on transect line
d <- array(NA, dim = c(W*W, J))

#Distance to nearest transect
dist <- NULL
#Simulate above quantities
for(g in 1:(W*W)){
  for(j in 1:J){
    d[g,j] <- sqrt((gr[g,1] - x.line[j])^2 + (gr[g,2] - y.line[j])^2)
  }
  dist[g] <- min(d[g,])
}

#Detection probability for distance sampling
p.ds <- rep(NA, N)
#Distance for each individual
dst <- dist[s]
#Individual presence/absence based on distance sampling
y.ds <- NULL
#Simulate distance sampling
for(i in 1:N){
  p.ds[i] <- exp(-dst[i] * dst[i] / (2 * sigma * sigma))
  y.ds[i] <- rbinom(1, 1, p.ds[i])
  if(dst[i] > 12)
    y.ds[i] == 0
}

# p.ds <- exp(-dst * dst / (2 * sigma * sigma))
# y.ds <- rbinom(N, 1, p.ds)
# y.ds[dst>12] <- 0

#Coordinates of detected individuals
uxds <- u1[y.ds == 1]
uyds <- u2[y.ds == 1]
#Distance of detected individuals
dst <- dst[y.ds == 1]
#Pixel ID
pixds <- s[y.ds == 1]
#Distance class length
v <- 1
#Transect half-width
B <- 12
#Midpoint of distance class
mdpt <- seq(0.5, 12, 1)
#Number of distance class
nD <- length(mdpt)
#Pixels covered by distance sampling (transect half-width is 10 pixels)
coverage <- which(dist <= B)
#Number of detected individuals
y.ds <- sum(y.ds)
#Distance class of detected individuals
dclass <- NULL
for(i in 1:y.ds){
  for(k in 1:nD){
    if(mdpt[k] - 0.5 <= dst[i] && dst[i] < mdpt[k+1] - 0.5)
      dclass[i] <- k
  }
}

#Visualize detected individuals from distance sampling
points(uxds, uyds, col = 'orange', pch = 20)

#---------------------------------------#
#-Draw covariate value for PO detection-#
#---------------------------------------#

#Simulate covariate on PO detection (Multivariate normal)
mu.x <- 0.3*(100)
mu.y <- 0.3*(100)
sigma.x <- 0.25*abs(100)
sigma.y <- 0.25*abs(100)
rho.xy <- 0.1
mu <- c(mu.x, mu.y)
Sigmaxy <- matrix(c(sigma.x^2, rep(rho.xy*sigma.x*sigma.y, 2), sigma.y^2), ncol=2)

w <- dmvnorm(gr, mean=mu, sigma=Sigmaxy)
w <- (w - mean(w))/sd(w)

#Visualize covariate on PO detection
image(rasterFromXYZ(cbind(gr,wcov)), col=topo.colors(10))

#------------------------#
#-Presence only sampling-#
#------------------------#

#Detection probability of PO
p.po <- expit(alpha1*wcov + alpha0) * error

#Individuals detected in PO
y.po <- NULL
for(i in 1:N){
  y.po[i] <- rbinom(1, 1, p.po[s[i]])
}

#Pixel ID for PO
pixpo <- s[y.po == 1]

#Coord of detected individuals
uxpo <- u1[y.po == 1]
uypo <- u2[y.po == 1]


#Number of PO detections per pixel
y.po <- as.data.frame(table(pixpo))
y.po$pixpo <- as.numeric(as.character(y.po$pixpo))
#Vector of pixels with PO detections
tmp <- rep(0, (W*W))
for(i in 1:length(y.po[,1])){
  tmp[y.po$pixpo[i]] <- y.po$Freq[i]
}
y.po <- tmp

#Visualize covariate
par(mar=c(3,3,3,6))
image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))
points(u1, u2, col = 'black', pch = 20)
points(uxpo, uypo, col = 'red', pch = 20)

#---------------#
#-Visualization-#
#---------------#
par(mar=c(3,3,3,6))
image(rasterFromXYZ(cbind(gr, x)), col = topo.colors(20))
points(u1, u2, col = 'black', pch = 20)
points(x.line, y.line, col = 'grey', pch = 20)
points(uxds, uyds, col = 'orange', pch = 20)
points(uxpo, uypo, col = 'red', pch = 20)

#----------------------------------#
#-Compile BUGS data for each model-#
#----------------------------------#

#Distance sampling model
str(dataDS <- list(x = x, G = (W*W),
                   nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass, nds = y.ds, y.ds = y.ds, coverage = coverage, pixds = pixds))

#Presence only model
str(dataPO <- list(x = x, G = (W*W),
                   w = w,  y.po = y.po))

#Integrated species distribution model
str(dataISDM <- list(x = x, G = (W*W),
                     nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass, nds = y.ds, y.ds = y.ds, coverage = coverage, pixds = pixds,
                     w = w,  y.po = y.po))

#----------------#
#-Initial values-#
#----------------#

#Inital value for N for distance sampling
N.st <- sum(y.ds) + 1

#Distance sampling model
initsDS <- function(){list(N = N.st, sigma=runif(1,1,3), beta1 = runif(1, 0.75, 1), beta0 = runif(1, -1, -0.75))}

#Presnece only model
initsPO <- function(){list(beta1 = runif(1, 0.75, 1), beta0 = runif(1, -1, -0.75), alpha1 = runif(1, 0, 2))}

#Integrated species distribution model
initsISDM <- function(){list(N = N.st, sigma=runif(1,1,3), beta1 = runif(1, 0.75, 1), beta0 = runif(1, -1, -0.75), alpha1 = runif(1, 0, 2))}

#------------#
#-Parameters-#
#------------#

#Distance sampling model
paramsDS <- c("Ntot", "beta0", "beta1")

#Presence only model
paramsPO <- c("Ntot", "alpha0", "alpha1", "beta0", "beta1")

#Integrated species distribution model
paramsISDM <- c("Ntot", "alpha0", "alpha1", "beta0", "beta1")

#-------------#
#-MCMC values-#
#-------------#
nb <- 1000
ni <- 7000
nt <- 2
nc <- 3
na <- 1000

#----------------#
#-Run each model-#
#----------------#

DS <- jagsUI(dataDS, initsDS, paramsDS, "DS.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

PO <- jagsUI(dataPO, initsPO, paramsPO, "PO.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

ISDM <- jagsUI(dataISDM, initsISDM, paramsISDM, "ISDM.txt", n.thin=nt, 
               n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

image(rasterFromXYZ(cbind(gr, exp(beta0 + beta1*x))), col = topo.colors(20))
image(rasterFromXYZ(cbind(gr, exp(DS$mean$beta0 + DS$mean$beta1*x))), col = topo.colors(20))
image(rasterFromXYZ(cbind(gr, exp(PO$mean$beta0 + PO$mean$beta1*x))), col = topo.colors(20))
image(rasterFromXYZ(cbind(gr, exp(ISDM$mean$beta0 + ISDM$mean$beta1*x))), col = topo.colors(20))
