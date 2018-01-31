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
#Load distance to pixels for robust and sparse
load(file = "dist.R.Rdata")
load(file = "dist.S.Rdata")

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
Out <- array(NA, dim = c(iter, 17, 6), 
             dimnames = list(NULL, NULL, c("N", "beta0", "beta1", "sigma", "alpha0", "alpha1")))

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

#Alternate transects
x.line.R <- rep(c(25,50,75), each = W) #Robust
y.line.R <- rep(seq(0.5, W-0.5, 1), 3) #Robust
x.line.S <- rep(c(25, 75), each = 60) #Sparse
y.line.S <- rep(c(seq(11, 40, 1), seq(61, 90, 1)), 2) #Sparse

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

#---------------------------------------#
#-Draw covariate value for PO detection-#
#---------------------------------------#
#Environmental covariate on PO detection (Multivariate normal)
#Mean of x-dim distribution
mu.x <- 50
#Mean of y-dim distribution
mu.y <- 50
#Variance of x-dim distribution
sigmax.R <- 0.5*abs(W) #Robust PO
sigmax.S <- 0.05*abs(W) #Sparse PO
#Variance of y-dim distribution
sigmay.R <- 0.5*abs(W)
sigmay.S <- 0.05*abs(W)
#Covariance of x-dim and y-dim distributions
rho.xy <- 0.1
mu <- c(mu.x, mu.y)
#Covariance matrix
Sigmaxy.R <- matrix(c(sigmax.R^2, rep(rho.xy*sigmax.R*sigmay.R, 2), sigmay.R^2), ncol=2)
Sigmaxy.S <- matrix(c(sigmax.S^2, rep(rho.xy*sigmax.S*sigmay.S, 2), sigmay.S^2), ncol=2)

w.R <- dmvnorm(gr, mean=mu, sigma=Sigmaxy.R)
w.R <- (w.R - mean(w.R))/sd(w.R)
w.S <- dmvnorm(gr, mean=mu, sigma=Sigmaxy.S)
w.S <- (w.S - mean(w.S))/sd(w.S)

#----------------------------------#
#-Simulate opportunistic surveying-#
#----------------------------------#

#Detection probability of PO
ppo.R <- expit(alpha1*w.R + alpha0) * error
ppo.S <- expit(alpha1*w.S + alpha0) * error
#Individuals detected in PO
ypo.R <- NULL
ypo.S <- NULL
for(i in 1:N){
  ypo.R[i] <- rbinom(1, 1, ppo.R[s[i]])
  ypo.S[i] <- rbinom(1, 1, ppo.S[s[i]])
}

#Pixel ID for PO
pixpo.R <- s[ypo.R == 1]
pixpo.S <- s[ypo.S == 1]

#Coord of detected individuals
uxpo.R <- u1[ypo.R == 1]
uypo.R <- u2[ypo.R == 1]
uxpo.S <- u1[ypo.S == 1]
uypo.S <- u2[ypo.S == 1]

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

#----------------#
#-Run each model-#
#----------------#

S <- list()

#Scenario 1: Presence only robust
S[[1]] <- jagsUI(data1, inits1, params, "PO.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#Scenario 2: Presence only sparse
S[[2]] <- jagsUI(data2, inits2, params, "PO.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#Scenario 3: Distance sampling robust
S[[3]] <- jagsUI(data3, inits3, params, "DSalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#Scenario 4: Distance sampling robust
S[[4]] <- jagsUI(data4, inits4, params, "DSalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#Scenario 5: ISDM robust DS & PO
S[[5]] <- jagsUI(data5, inits5, params, "ISDMalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#Scenario 6: ISDM sparse DS & PO
S[[6]] <- jagsUI(data6, inits6, params, "ISDMalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#Scenario 7: ISDM robust DS & sparse PO
S[[7]] <- jagsUI(data7, inits7, params, "ISDMalt.txt", n.thin=nt, 
             n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

#Scenario 8: ISDM sparse DS & robust PO
S[[8]] <- jagsUI(data8, inits8, params, "ISDMalt.txt", n.thin=nt, 
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

for(sc in 2:9){
  Out[z,sc,1] <- S[[sc-1]]$mean$Ntot
  Out[z,sc,2] <- S[[sc-1]]$mean$beta0
  Out[z,sc,3] <- S[[sc-1]]$mean$beta1
  Out[z,sc,4] <- S[[sc-1]]$mean$sigma
  Out[z,sc,5] <- S[[sc-1]]$mean$alpha0
  Out[z,sc,6] <- S[[sc-1]]$mean$alpha1
  
  Out[z,sc+8,1] <- S[[sc]]$Rhat$Ntot
  Out[z,sc+8,2] <- S[[sc]]$Rhat$beta0
  Out[z,sc+8,3] <- S[[sc]]$Rhat$beta1
  Out[z,sc+8,4] <- S[[sc]]$Rhat$sigma
  Out[z,sc+8,5] <- S[[sc]]$Rhat$alpha0
  Out[z,sc+8,6] <- S[[sc]]$Rhat$alpha1
}

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