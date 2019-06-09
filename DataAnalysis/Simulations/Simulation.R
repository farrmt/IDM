#----------------------------------------------------#
#----Simulation study for ISDM combining-------------#
#----distance sampling and opportunistic sampling----#
#----Created by Matthew Farr-------------------------#
#----------------------------------------------------#

#-----------#
#-Libraries-#
#-----------#

library(mvtnorm)
library(jagsUI)
library(dplyr)

#-------------------#
#-Working Directory-#
#-------------------#

setwd("./DataAnalysis/Simulations")

#-----------#
#-Functions-#
#-----------#

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
#-Create Transect Units-#
#-----------------------#

#Transect unit size
tu <- 10
#Mid-point coordinates of first dimension
tu.m <- seq(tu/2, W - tu/2, tu)
#Expand mid-point to entire grid
grt <- expand.grid(tu.m, tu.m)

#---------------------------#
#-Calculate pixel distances-#
#---------------------------#

#All possible transects
# x.line <- rep(tu.m, each = W)
# y.line <- rep(seq(0.5, W-0.5, 1), 10)
# line <- cbind(x.line, y.line)
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

#Load distance to pixels as above code is expensive
load(file = "dist.Rdata")

#-------------------------------------#
#-Draw environmental covariate values-#
#-------------------------------------#

#Code borrowed from Kéry & Royle 2016 pg. 534
#Variance-covariance matrix based on Euclidean distance
#V <- exp(-e2dist(gr, gr))
#Covariate values from a correlated multivariate normal (pg.534 AHM)
#w <- as.vector(t(chol(V))%*%rnorm(W^2))

#Load ecological covariate as above code is expensive
load(file = "wcov.Rdata")

#--------------------------------#
#-Values to save from simulation-#
#--------------------------------#

#Number of simulation iterations
iter <- 1

#True and estimated parameter values to save
Out <- array(NA, dim = c(iter, 25, 6), 
             dimnames = list(NULL, NULL, c("N", "lambda0", "beta1", "sigma", "p0", "alpha1")))

#------------------#
#-Begin Simulation-#
#------------------#

start.time <- Sys.time()
for(i in 1:iter){
  
#-----------------------#
#-Draw parameter values-#
#-----------------------#

#Intercept parameter for intensity function
lambda0 <- log(runif(1, 0.05, 1))
#Effect parameter of enivornment on intensity
beta1 <- runif(1, -1.25, 1.25)
#Intercept parameter for prensence only (PO) detection
p0.L <- logit(0.5)
p0.S <- logit(0.1)
#Effect parameter of environment on PO detection
alpha1 <- runif(1, 0, 2)
#Scale parameter for DS detection
sigma <- runif(1, 0.75, 1.25)

#--------------------------------#
#-Simulate number of individuals-#
#--------------------------------#

#Intensity function
intensity <- exp(lambda0 + beta1 * w)
#Probability of each pixel based on intensity function
probs <- intensity/sum(intensity)
#Simulate number of individuals in Region B
N <- rpois(1, sum(intensity))
#Sample location, s, of simulated individuals
s <- sample(1:(W*W), N, replace=TRUE, prob=probs)
#Assign X and Y coordinates
u1 <- gr[s,1]
u2 <- gr[s,2]

#-----------------------#
#-Sample Transect Units-#
#-----------------------#

#Number of transect scenarios
nts <- 6
#Samples per scenario
ns <- c(50,25,20,15,10,5)
#Sampled transects for each scenario
tsamp <- array(NA, dim = c(50,nts))
#First scenario
tsamp[,1] <- sample(1:100, ns[1], replace = FALSE)
#Remanding scenarios
for(j in 2:nts){
  tsamp[,j] <- c(sample(tsamp[1:ns[j-1],j-1], ns[j], replace = FALSE), rep(NA, ns[1] - ns[j]))
}

#----------------------------#
#-Simulate distance sampling-#
#----------------------------#

#Sampled pixels and corresponding distances
fill <- seq(1,5100,by=100)
coverage <- array(NA, dim = c(5000,nts))
dst <- array(NA, dim = c(5000,nts))
for(j in 1:nts){
  for(k in 1:ns[j]){
    coverage[fill[k]:(fill[k+1]-1),j] <- which((gr[,1] > (grt[tsamp[k,j],1] - 5)) & 
                                                 ((grt[tsamp[k,j],1] + 5) > gr[,1]) & 
                                                 (gr[,2] > (grt[tsamp[k,j],2] - 5)) & 
                                                 ((grt[tsamp[k,j],2] + 5) > gr[,2]))
    dst[fill[k]:(fill[k+1]-1),j] <- dist[coverage[fill[k]:(fill[k+1]-1),j]]
  }
}

#Individual (pixels) within region B
sds<- as.data.frame(s[s%in%coverage[,1]])
colnames(sds) <- "p"
#Detection probability for distance sampling
pids <- rep(NA, length(sds))
#Distance for each individual
dst50 <- as.data.frame(cbind(coverage[,1], dst[,1]))
colnames(dst50) <- c("p", "dst")
x50 <- inner_join(dst50, sds, by = "p")
dst.S <- x50$dst
#Individual presence/absence based on distance sampling
xds <- NULL
#Simulate distance sampling
pids <- exp(-dst.S * dst.S / (2 * sigma * sigma))
xds <- rbinom(length(sds[,1]), 1, pids)
x50$x <- xds
x50$uxds <- gr[x50[,1],1]
x50$uyds <- gr[x50[,1],2]
x50$dclass <- NA
x50 <- x50%>%arrange(p)

tmp2 <- as.data.frame(coverage[,2])
colnames(tmp2) <- c("p")
x25 <- inner_join(x50,tmp2,by="p")
x25 <- x25%>%arrange(p)

tmp2 <- as.data.frame(coverage[,3])
colnames(tmp2) <- c("p")
x20 <- inner_join(x50,tmp2,by="p")
x20 <- x20%>%arrange(p)

tmp2 <- as.data.frame(coverage[,4])
colnames(tmp2) <- c("p")
x15 <- inner_join(x50,tmp2,by="p")
x15 <- x15%>%arrange(p)

tmp2 <- as.data.frame(coverage[,5])
colnames(tmp2) <- c("p")
x10 <- inner_join(x50,tmp2,by="p")
x10 <- x10%>%arrange(p)

tmp2 <- as.data.frame(coverage[,6])
colnames(tmp2) <- c("p")
x5 <- inner_join(x50,tmp2,by="p")
x5 <- x5%>%arrange(p)

XD <- list(x50,x25,x20,x15,x10,x5)
#Size of subregion B
B <- NULL
#Number of obesrvation per pixel of subregion
x <- list()
#observed pixels
pixD <- list()
nobs <- list()
#Distance class
for(j in 1:nts){
  B[j] <- length(na.omit(coverage[,j]))
  x[[j]] <- as.data.frame(table(XD[[j]][XD[[j]][,3]==1,1]))
  x[[j]][,1] <- as.numeric(as.character(x[[j]][,1]))
  colnames(x[[j]]) <- c("p", "freq")
  coverB <- as.data.frame(coverage[1:B[j],j])
  colnames(coverB) <- "p"
  x[[j]] <- full_join(x[[j]], coverB, by = "p")
  x[[j]] <- x[[j]]%>%arrange(p)
  x[[j]]$dst <- dst[1:B[j],j]
  pixD[[j]] <- which(is.na(x[[j]][,2]))
  nobs[[j]] <- length(pixD[[j]])
}

x <- rapply(x, f=function(x) ifelse(is.na(x),0,x), how="replace" )


#---------------------------------------#
#-Draw covariate value for PO detection-#
#---------------------------------------#
#Code borrowed from Dorazio 2014
#Environmental covariate on PO detection (Multivariate normal)
#Mean of x-dim distribution
mu.x <- runif(1, 25, 75)
#Mean of y-dim distribution
mu.y <- runif(1, 25, 75)
#Variance of x-dim distribution
sigmax.L <- 0.75*abs(W)
sigmax.S <- 0.1*abs(W)
#Variance of y-dim distribution
sigmay.L <- 0.75*abs(W)
sigmay.S <- 0.1*abs(W)
#Covariance of x-dim and y-dim distributions
rho.xy <- 0.25
mu <- c(mu.x, mu.y)
#Covariance matrix
Sigmaxy.L <- matrix(c(sigmax.L^2, rep(rho.xy*sigmax.L*sigmay.L, 2), sigmay.L^2), ncol=2)
Sigmaxy.S <- matrix(c(sigmax.S^2, rep(rho.xy*sigmax.S*sigmay.S, 2), sigmay.S^2), ncol=2)

z.L <- dmvnorm(gr, mean=mu, sigma=Sigmaxy.L)
z.L <- (z.L - mean(z.L))/sd(z.L)
z.S <- dmvnorm(gr, mean=mu, sigma=Sigmaxy.S)
z.S <- (z.S - mean(z.S))/sd(z.S)

#----------------------------------#
#-Simulate opportunistic surveying-#
#----------------------------------#

#Detection probability of PO
p.L <- expit(p0.L + alpha1 * z.L)
p.S <- expit(p0.S + alpha1 * z.S)

#Individuals detected in PO
y.L <- y.S <- NULL
for(j in 1:N){
  y.L[j] <- rbinom(1, 1, p.L[s[i]])
  y.S[j] <- rbinom(1, 1, p.S[s[i]])
}

#Pixel ID for true presence
pixpo.L <- s[y.L == 1]
pixpo.S <- s[y.S == 1]

#Coord of true presence
uxpo.L <- u1[y.L == 1]
uypo.L <- u2[y.L == 1]
uxpo.S <- u1[y.S == 1]
uypo.S <- u2[y.S == 1]

#Unobserved sampling / PO pixel ID
error <- 0.3
pixpo.L <- sample(pixpo.L, length(pixpo.L)*error, replace = FALSE)
pixpo.S <- sample(pixpo.S, length(pixpo.S)*error, replace = FALSE)

#Number of PO detections per pixel
y.L <- as.data.frame(table(pixpo.L))
y.L$pixpo.L <- as.numeric(as.character(y.L$pixpo.L))
y.S <- as.data.frame(table(pixpo.S))
y.S$pixpo.S <- as.numeric(as.character(y.S$pixpo.S))

#Vector of pixels with PO detections
tmp <- rep(0, (W*W))
for(j in 1:length(y.L[,1])){
  tmp[y.L$pixpo.L[j]] <- y.L$Freq[j]
}
y.L <- tmp
tmp <- rep(0, (W*W))
for(j in 1:length(y.S[,1])){
  tmp[y.S$pixpo.S[j]] <- y.S$Freq[j]
}
y.S <- tmp

#----------------------------------------------#
#-Compile BUGS data for each sampling scenario-#
#----------------------------------------------#

#Scenario 1: PO Large DS 0%
str(data1 <- list(w = w, G = (W*W),
                  z = z.L,  y = y.L))

#Scenario 2: PO Large DS 5%
str(data2 <- list(w = w, G = (W*W), B = B[6], d = x[[6]]$dst,
                  x = x[[6]]$freq, coverage = x[[6]]$p,
                  z = z.L, y = y.L))

#Scenario 3: PO Large DS 10%
str(data3 <- list(w = w, G = (W*W), B = B[5], d = x[[5]]$dst,
                  x = x[[5]]$freq, coverage = x[[5]]$p,
                  z = z.L, y = y.L))

#Scenario 4: PO Large DS 15%
str(data4 <- list(w = w, G = (W*W), B = B[4], d = x[[4]]$dst,
                  x = x[[4]]$freq, coverage = x[[4]]$p,
                  z = z.L, y = y.L))

#Scenario 5: PO Large DS 20%
str(data5 <- list(w = w, G = (W*W), B = B[3], d = x[[3]]$dst,
                  x = x[[3]]$freq, coverage = x[[3]]$p,
                  z = z.L, y = y.L))

#Scenario 6: PO Large DS 25%
str(data6 <- list(w = w, G = (W*W), B = B[2], d = x[[2]]$dst,
                  x = x[[2]]$freq, coverage = x[[2]]$p,
                  z = z.L, y = y.L))

#Scenario 7: PO Small DS 0%
str(data7 <- list(w = w, G = (W*W),
                  z = z.S,  y = y.S))

#Scenario 8: PO Small DS 5%
str(data8 <- list(w = w, G = (W*W), B = B[6], d = x[[6]]$dst,
                  x = x[[6]]$freq, coverage = x[[6]]$p,
                  z = z.S, y = y.S))

#Scenario 9: PO Small DS 10%
str(data9 <- list(w = w, G = (W*W), B = B[5], d = x[[5]]$dst,
                  x = x[[5]]$freq, coverage = x[[5]]$p,
                  z = z.S, y = y.S))

#Scenario 10: PO Small DS 15%
str(data10 <- list(w = w, G = (W*W), B = B[4], d = x[[4]]$dst,
                   x = x[[4]]$freq, coverage = x[[4]]$p,
                   z = z.S, y = y.S))

#Scenario 11: PO Small DS 20%
str(data11 <- list(w = w, G = (W*W), B = B[3], d = x[[3]]$dst,
                   x = x[[3]]$freq, coverage = x[[3]]$p,
                   z = z.S, y = y.S))

#Scenario 12: PO Small DS 20%
str(data12 <- list(w = w, G = (W*W), B = B[2], d = x[[2]]$dst,
                   x = x[[2]]$freq, coverage = x[[2]]$p,
                   z = z.S, y = y.S))

#----------------#
#-Initial values-#
#----------------#

#Scenario 1: PO Large DS 0%
inits1 <- function(){list(N = rep(1, W*W), beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05), 
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}

#Scenario 2: PO Large DS 5%
inits2 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}

#Scenario 3: PO Large DS 10%
inits3 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}

#Scenario 4: PO Large DS 15%
inits4 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}

#Scenario 5: PO Large DS 20%
inits5 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}

#Scenario 6: PO Large DS 25%
inits6 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}

#Scenario 7: PO Small DS 0%
inits7 <- function(){list(N = rep(1, W*W), beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05), 
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}

#Scenario 8: PO Small DS 5%
inits8 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}

#Scenario 9: PO Small DS 10%
inits9 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}

#Scenario 10: PO Small DS 15%
inits10 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                           beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                           alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}

#Scenario 11: PO Small DS 20%
inits11 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                           beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                           alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}

#Scenario 12: PO Small DS 25%
inits12 <- function(){list(N = rep(1, W*W), sigma = runif(1, sigma-0.05, sigma+0.05), 
                           beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                           alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}

#------------#
#-Parameters-#
#------------#

params <- c("Ntot", "sigma", "p0", "alpha1", "lambda0", "beta1")

#-------------#
#-MCMC values-#
#-------------#

nb <- 10000
ni <- 16000
nt <- 2
nc <- 3
na <- 500

#----------------#
#-Run each model-#
#----------------#

S <- list()

#Scenario 1: PO Large DS 0%
S[[1]] <- jagsUI(data1, inits1, params, "PO.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)

#Scenario 2: PO Large DS 5%
S[[2]] <- jagsUI(data2, inits2, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)

#Scenario 3: PO Large DS 10%
S[[3]] <- jagsUI(data3, inits3, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 4: PO Large DS 15%
S[[4]] <- jagsUI(data4, inits4, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 5: PO Large DS 20%
S[[5]] <- jagsUI(data5, inits5, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 6: PO Large DS 25%
S[[6]] <- jagsUI(data6, inits6, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 7: PO Small DS 0%
S[[7]] <- jagsUI(data7, inits7, params, "PO.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 8: PO Small DS 5%
S[[8]] <- jagsUI(data8, inits8, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 9: PO Small DS 10%
S[[9]] <- jagsUI(data9, inits9, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 10: PO Small DS 15%
S[[10]] <- jagsUI(data10, inits10, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 11: PO Small DS 20%
S[[11]] <- jagsUI(data11, inits11, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 12: PO Small DS 25%
S[[12]] <- jagsUI(data12, inits12, params, "ISDM.txt", n.thin=nt, 
                  n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#----------------------------------#
#-Save values from each simulation-#
#----------------------------------#

#True parameter values
Out[i,1,1] <- N
Out[i,1,2] <- lambda0
Out[i,1,3] <- beta1
Out[i,1,4] <- sigma
Out[i,1,5] <- NA
Out[i,1,6] <- alpha1

for(sc in 2:13){
  tryCatch({Out[i,sc,1] <- S[[sc-1]]$mean$Ntot}, error = function(e){})
  tryCatch({Out[i,sc,2] <- S[[sc-1]]$mean$lambda0}, error = function(e){})
  tryCatch({Out[i,sc,3] <- S[[sc-1]]$mean$beta1}, error = function(e){})
  tryCatch({Out[i,sc,4] <- S[[sc-1]]$mean$sigma}, error = function(e){})
  tryCatch({Out[i,sc,5] <- S[[sc-1]]$mean$p0}, error = function(e){})
  tryCatch({Out[i,sc,6] <- S[[sc-1]]$mean$alpha1}, error = function(e){})
  
  tryCatch({Out[i,sc+12,1] <- S[[sc-1]]$Rhat$Ntot}, error = function(e){})
  tryCatch({Out[i,sc+12,2] <- S[[sc-1]]$Rhat$lambda0}, error = function(e){})
  tryCatch({Out[i,sc+12,3] <- S[[sc-1]]$Rhat$beta1}, error = function(e){})
  tryCatch({Out[i,sc+12,4] <- S[[sc-1]]$Rhat$sigma}, error = function(e){})
  tryCatch({Out[i,sc+12,5] <- S[[sc-1]]$Rhat$p0}, error = function(e){})
  tryCatch({Out[i,sc+12,6] <- S[[sc-1]]$Rhat$alpha1}, error = function(e){})
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

#------------#
#-References-#
#------------#

#Dorazio, R.M. (2014) Accounting for imperfect detection and survey bias in statistical analysis of presence-only data. Global Ecology and Biogeography, 23, 1472–1484.

#Kéry, M. & Royle, J.A. (2016) Applied hierarchical modeling in ecology: Analysis of distribution, abundance and species richness in R and BUGS (volume 1 – prelude and static models), Elsevier, Amsterdam.