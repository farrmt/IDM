#------------------------------#
#-Libraries used in simulation-#
#------------------------------#

library(mvtnorm)
library(jagsUI)

#------------------------------#
#-Functions used in simulation-#
#------------------------------#

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

#-------------------------------------#
#-Draw environmental covariate values-#
#-------------------------------------#

#Load x covariate instead
load(file = "xcov.Rdata")

#--------------------------------#
#-Parameter values for scenarios-#
#--------------------------------#

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
# alpha0 <- logit(0.25) #Medium detection (25%)
alpha0 <- logit(0.5) #High detection (50%)

#Effect parameter of environment on PO detection
# alpha1 <- 0 #No effect on PO detection
alpha1 <- 2.25 #Weak effect on PO detection X 0.86 (14% decrease)
# alpha1 <- 5 #Medium effect on PO detection X 0.66 (34% decrease)
# alpha1 <- 10 #Strong effect on PO detection X 0.4 (60% decrease)

#Unobserved sampling error
error <- 0.4

#Scale parameter for DS detection
# sigma <- 1 #Low detection (10%)
sigma <- 2.5 #Medium detection (26%)
# sigma <- 5 #High detection (50%)

#--------------------------------#
#-Values to save from simulation-#
#--------------------------------#

#Number of simulation iterations
iter <- 1

#True and estimated parameter values to save
Out <- array(NA, dim = c(iter, 7, 6), 
             dimnames = list(NULL, c("Truth", "DS", "PO", "ISDM", "DS.Rhat", "PO.Rhat", "ISDM.Rhat"),
                             c("N", "beta0", "beta1", "sigma", "alpha0", "alpha1")))

#------------------#
#-Begin Simulation-#
#------------------#

start.time <- Sys.time()
for(z in 1:iter){

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

#----------------------------#
#-Simulate distance sampling-#
#----------------------------#

#Transects
x.line <- rep(c(25,75), each = W)
y.line <- rep(seq(0.5, W-0.5, 1), 2)
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
p.ds <- exp(-dst * dst / (2 * sigma * sigma))
y.ds <- rbinom(N, 1, p.ds)
y.ds[dst>12] <- 0
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

#---------------------------------------#
#-Draw covariate value for PO detection-#
#---------------------------------------#
#Environmental covariate on PO detection (Multivariate normal)
#Mean of x-dim distribution
mu.x <- runif(1,20,80)
#Mean of y-dim distribution
mu.y <- runif(1,20,80)
#Variance of x-dim distribution
sigma.x <- 0.25*abs(100)
#Variance of y-dim distribution
sigma.y <- 0.25*abs(100)
#Covariance of x-dim and y-dim distributions
rho.xy <- 0.1
mu <- c(mu.x, mu.y)
#Covariance matrix
Sigmaxy <- matrix(c(sigma.x^2, rep(rho.xy*sigma.x*sigma.y, 2), sigma.y^2), ncol=2)

w <- dmvnorm(gr, mean=mu, sigma=Sigmaxy)
w <- (w - mean(w))/sd(w)

#----------------------------------#
#-Simulate opportunistic surveying-#
#----------------------------------#

#Detection probability of PO
p.po <- expit(alpha1*w + alpha0) * error
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

#----------------------------------#
#-Compile BUGS data for each model-#
#----------------------------------#

#Distance sampling model
str(dataDS <- list(x = x, G = (W*W), nD = nD, v = v, B = B, mdpt = mdpt, 
                   dclass = dclass, nds = y.ds, y.ds = y.ds, 
                   coverage = coverage, pixds = pixds))

#Presence only model
str(dataPO <- list(x = x, G = (W*W), w = w,  y.po = y.po))

#Integrated species distribution model
str(dataISDM <- list(x = x, G = (W*W), nD = nD, v = v, B = B, mdpt = mdpt,
                     dclass = dclass, nds = y.ds, y.ds = y.ds, 
                     coverage = coverage, pixds = pixds,
                     w = w,  y.po = y.po))

#----------------#
#-Initial values-#
#----------------#

#Inital value for N for distance sampling
N.st <- sum(y.ds) + 1

#Distance sampling model
initsDS <- function(){list(N = N.st, sigma=runif(1,1,3), 
                           beta1 = runif(1, 0.75, 1), beta0 = runif(1, -1, -0.75))}

#Presnece only model
initsPO <- function(){list(beta1 = runif(1, 0.75, 1), beta0 = runif(1, -1, -0.75), 
                           alpha0 = runif(1, -3, 0), alpha1 = runif(1, 0, 2))}

#Integrated species distribution model
initsISDM <- function(){list(N = N.st, sigma=runif(1,1,3),
                             beta1 = runif(1, 0.75, 1), beta0 = runif(1, -1, -0.75), 
                             alpha0 = runif(1, -3, 0), alpha1 = runif(1, 0, 2))}

#------------#
#-Parameters-#
#------------#

#Distance sampling model
paramsDS <- c("Ntot", "beta0", "beta1", "sigma")

#Presence only model
paramsPO <- c("Ntot", "alpha0", "alpha1", "beta0", "beta1")

#Integrated species distribution model
paramsISDM <- c("Ntot", "alpha0", "alpha1", "beta0", "beta1", "sigma")

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

#----------------------------------#
#-Save values from each simulation-#
#----------------------------------#

#True parameter values
Out[z,1,1] <- N
Out[z,1,2] <- beta0
Out[z,1,3] <- beta1
Out[z,1,4] <- sigma
Out[z,1,5] <- alpha0
Out[z,1,6] <- alpha1
  
#Estimated distance sampling parameter values
Out[z,2,1] <- DS$mean$Ntot
Out[z,2,2] <- DS$mean$beta0
Out[z,2,3] <- DS$mean$beta1
Out[z,2,4] <- DS$mean$sigma

Out[z,5,1] <- DS$Rhat$Ntot
Out[z,5,2] <- DS$Rhat$beta0
Out[z,5,3] <- DS$Rhat$beta1
Out[z,5,4] <- DS$Rhat$sigma

#Estimated presence only parameter values
Out[z,3,1] <- PO$mean$Ntot
Out[z,3,2] <- PO$mean$beta0
Out[z,3,3] <- PO$mean$beta1
Out[z,3,5] <- PO$mean$alpha0
Out[z,3,6] <- PO$mean$alpha1

Out[z,6,1] <- PO$Rhat$Ntot
Out[z,6,2] <- PO$Rhat$beta0
Out[z,6,3] <- PO$Rhat$beta1
Out[z,6,5] <- PO$Rhat$alpha0
Out[z,6,6] <- PO$Rhat$alpha1

#Estimated integrated species distribution modeling parameter values
Out[z,4,1] <- ISDM$mean$Ntot
Out[z,4,2] <- ISDM$mean$beta0
Out[z,4,3] <- ISDM$mean$beta1
Out[z,4,4] <- ISDM$mean$sigma
Out[z,4,5] <- ISDM$mean$alpha0
Out[z,4,6] <- ISDM$mean$alpha1

Out[z,7,1] <- ISDM$Rhat$Ntot
Out[z,7,2] <- ISDM$Rhat$beta0
Out[z,7,3] <- ISDM$Rhat$beta1
Out[z,7,4] <- ISDM$Rhat$sigma
Out[z,7,5] <- ISDM$Rhat$alpha0
Out[z,7,6] <- ISDM$Rhat$alpha1
}#End simulation
end.time <- Sys.time()
Time <- end.time - start.time

#-Save HPCC output-#
ID <- gsub(" ","_",Sys.time())
ID <- gsub(":", "-", ID)
output <- list(Out, Time)
heads <- c("Out", "Time")
output <- setNames(output, nm = heads)
save(output, file = paste("output", ID, ".R", sep=""))