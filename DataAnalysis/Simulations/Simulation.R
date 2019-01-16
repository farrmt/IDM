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

#Variance-covariance matrix based on Euclidean distance
#V <- exp(-e2dist(gr, gr))
#Covariate values from a correlated multivariate normal (pg.534 AHM)
#x <- as.vector(t(chol(V))%*%rnorm(W^2))

#Load x covariate as above code is expensive
load(file = "xcov.Rdata")

#--------------------------------#
#-Values to save from simulation-#
#--------------------------------#

#Number of simulation iterations
iter <- 1

#True and estimated parameter values to save
Out <- array(NA, dim = c(iter, 25, 6), 
             dimnames = list(NULL, NULL, c("N", "beta0", "beta1", "sigma", "alpha0", "alpha1")))

#------------------#
#-Begin Simulation-#
#------------------#

start.time <- Sys.time()
for(z in 1:iter){
  
#-----------------------#
#-Draw parameter values-#
#-----------------------#

#Intercept parameter for intensity function
beta0 <- log(runif(1, 0.05, 1))
#Effect parameter of enivornment on intensity
beta1 <- runif(1, -1.25, 1.25)
#Intercept parameter for prensence only (PO) detection
alpha0.L <- logit(0.5)
alpha0.S <- logit(0.1)
#Effect parameter of environment on PO detection
alpha1 <- runif(1, 0, 2)
#Scale parameter for DS detection
sigma <- runif(1, 0.75, 1.25)

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
for(i in 2:nts){
  tsamp[,i] <- c(sample(tsamp[1:ns[i-1],i-1], ns[i], replace = FALSE), rep(NA, ns[1] - ns[i]))
}

#----------------------------#
#-Simulate distance sampling-#
#----------------------------#

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

YD <- list(y50,y25,y20,y15,y10,y5)
#Size of subregion B
B <- NULL
#Number of obesrvation per pixel of subregion
yds <- list()
#observed pixels
pixD <- list()
nobs <- list()
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
  pixD[[i]] <- which(is.na(yds[[i]][,2]))
  nobs[[i]] <- length(pixD[[i]])
}

yds <- rapply(yds, f=function(x) ifelse(is.na(x),0,x), how="replace" )


#---------------------------------------#
#-Draw covariate value for PO detection-#
#---------------------------------------#
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

w.L <- dmvnorm(gr, mean=mu, sigma=Sigmaxy.L)
w.L <- (w.L - mean(w.L))/sd(w.L)
w.S <- dmvnorm(gr, mean=mu, sigma=Sigmaxy.S)
w.S <- (w.S - mean(w.S))/sd(w.S)

#----------------------------------#
#-Simulate opportunistic surveying-#
#----------------------------------#

#Detection probability of PO
ppo.L <- expit(alpha1*w.L + alpha0.L)
ppo.S <- expit(alpha1*w.S + alpha0.S)

#Individuals detected in PO
ypo.L <- ypo.S <- NULL
for(i in 1:N){
  ypo.L[i] <- rbinom(1, 1, ppo.L[s[i]])
  ypo.S[i] <- rbinom(1, 1, ppo.S[s[i]])
}

#Pixel ID for true presence
pixpo.L <- s[ypo.L == 1]
pixpo.S <- s[ypo.S == 1]

#Coord of true presence
uxpo.L <- u1[ypo.L == 1]
uypo.L <- u2[ypo.L == 1]
uxpo.S <- u1[ypo.S == 1]
uypo.S <- u2[ypo.S == 1]

#Unobserved sampling / PO pixel ID
error <- 0.3
pixpo.L <- sample(pixpo.L, length(pixpo.L)*error, replace = FALSE)
pixpo.S <- sample(pixpo.S, length(pixpo.S)*error, replace = FALSE)

#Number of PO detections per pixel
ypo.L <- as.data.frame(table(pixpo.L))
ypo.L$pixpo.L <- as.numeric(as.character(ypo.L$pixpo.L))
ypo.S <- as.data.frame(table(pixpo.S))
ypo.S$pixpo.S <- as.numeric(as.character(ypo.S$pixpo.S))

#Vector of pixels with PO detections
tmp <- rep(0, (W*W))
for(i in 1:length(ypo.L[,1])){
  tmp[ypo.L$pixpo.L[i]] <- ypo.L$Freq[i]
}
ypo.L <- tmp
tmp <- rep(0, (W*W))
for(i in 1:length(ypo.S[,1])){
  tmp[ypo.S$pixpo.S[i]] <- ypo.S$Freq[i]
}
ypo.S <- tmp

#----------------------------------------------#
#-Compile BUGS data for each sampling scenario-#
#----------------------------------------------#

#Scenario 1: PO Large DS 0%
str(data1 <- list(x = x, G = (W*W),
                  w = w.L,  y.po = ypo.L))

#Scenario 2: PO Large DS 5%
str(data2 <- list(x = x, G = (W*W), B = B[6], dst = yds[[6]]$dst,
                  y.ds = yds[[6]]$freq, coverage = yds[[6]]$p,
                  w = w.L, y.po = ypo.L))

#Scenario 3: PO Large DS 10%
str(data3 <- list(x = x, G = (W*W), B = B[5], dst = yds[[5]]$dst,
                  y.ds = yds[[5]]$freq, coverage = yds[[5]]$p,
                  w = w.L, y.po = ypo.L))

#Scenario 4: PO Large DS 15%
str(data4 <- list(x = x, G = (W*W), B = B[4], dst = yds[[4]]$dst,
                  y.ds = yds[[4]]$freq, coverage = yds[[4]]$p,
                  w = w.L, y.po = ypo.L))

#Scenario 5: PO Large DS 20%
str(data5 <- list(x = x, G = (W*W), B = B[3], dst = yds[[3]]$dst,
                  y.ds = yds[[3]]$freq, coverage = yds[[3]]$p,
                  w = w.L, y.po = ypo.L))

#Scenario 6: PO Large DS 25%
str(data6 <- list(x = x, G = (W*W), B = B[2], dst = yds[[2]]$dst,
                  y.ds = yds[[2]]$freq, coverage = yds[[2]]$p,
                  w = w.L, y.po = ypo.L))

#Scenario 7: PO Small DS 0%
str(data7 <- list(x = x, G = (W*W),
                  w = w.S,  y.po = ypo.S))

#Scenario 8: PO Small DS 5%
str(data8 <- list(x = x, G = (W*W), B = B[6], dst = yds[[6]]$dst,
                  y.ds = yds[[6]]$freq, coverage = yds[[6]]$p,
                  w = w.S, y.po = ypo.S))

#Scenario 9: PO Small DS 10%
str(data9 <- list(x = x, G = (W*W), B = B[5], dst = yds[[5]]$dst,
                  y.ds = yds[[5]]$freq, coverage = yds[[5]]$p,
                  w = w.S, y.po = ypo.S))

#Scenario 10: PO Small DS 15%
str(data10 <- list(x = x, G = (W*W), B = B[4], dst = yds[[4]]$dst,
                  y.ds = yds[[4]]$freq, coverage = yds[[4]]$p,
                  w = w.S, y.po = ypo.S))

#Scenario 11: PO Small DS 20%
str(data11 <- list(x = x, G = (W*W), B = B[3], dst = yds[[3]]$dst,
                  y.ds = yds[[3]]$freq, coverage = yds[[3]]$p,
                  w = w.S, y.po = ypo.S))

#Scenario 12: PO Small DS 20%
str(data12 <- list(x = x, G = (W*W), B = B[2], dst = yds[[2]]$dst,
                   y.ds = yds[[2]]$freq, coverage = yds[[2]]$p,
                   w = w.S, y.po = ypo.S))

#----------------#
#-Initial values-#
#----------------#

#Scenario 1: PO Large DS 0%
inits1 <- function(){list(beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05), 
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.L-0.05, alpha0.L+0.05))}

#Scenario 2: PO Large DS 5%
inits2 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.L-0.05, alpha0.L+0.05))}

#Scenario 3: PO Large DS 10%
inits3 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.L-0.05, alpha0.L+0.05))}

#Scenario 4: PO Large DS 15%
inits4 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                           beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                           alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.L-0.05, alpha0.L+0.05))}

#Scenario 5: PO Large DS 20%
inits5 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                           beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                           alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.L-0.05, alpha0.L+0.05))}

#Scenario 6: PO Large DS 25%
inits6 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.L-0.05, alpha0.L+0.05))}

#Scenario 7: PO Small DS 0%
inits7 <- function(){list(beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05), 
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.S-0.05, alpha0.S+0.05))}

#Scenario 8: PO Small DS 5%
inits8 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.S-0.05, alpha0.S+0.05))}

#Scenario 9: PO Small DS 10%
inits9 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                          beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                          alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.S-0.05, alpha0.S+0.05))}

#Scenario 10: PO Small DS 15%
inits10 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                           beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                           alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.S-0.05, alpha0.S+0.05))}

#Scenario 11: PO Small DS 20%
inits11 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                           beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                           alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.S-0.05, alpha0.S+0.05))}

#Scenario 12: PO Small DS 25%
inits12 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                           beta1 = runif(1, beta1-0.05, beta1+0.05), beta0 = runif(1, beta0-0.05, beta0+0.05),
                           alpha1 = runif(1, alpha1-0.05, alpha1+0.05), alpha0 = runif(1, alpha0.S-0.05, alpha0.S+0.05))}

#------------#
#-Parameters-#
#------------#

params <- c("Ntot", "sigma", "alpha0", "alpha1", "beta0", "beta1")

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
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

#Scenario 2: PO Large DS 5%
S[[2]] <- jagsUI(data2, inits2, params, "ISDM.txt", n.thin=nt, 
                 n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na)

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
Out[z,1,1] <- N
Out[z,1,2] <- beta0
Out[z,1,3] <- beta1
Out[z,1,4] <- sigma
Out[z,1,5] <- NA
Out[z,1,6] <- alpha1

for(sc in 2:13){
  tryCatch({Out[z,sc,1] <- S[[sc-1]]$mean$Ntot}, error = function(e){})
  tryCatch({Out[z,sc,2] <- S[[sc-1]]$mean$beta0}, error = function(e){})
  tryCatch({Out[z,sc,3] <- S[[sc-1]]$mean$beta1}, error = function(e){})
  tryCatch({Out[z,sc,4] <- S[[sc-1]]$mean$sigma}, error = function(e){})
  tryCatch({Out[z,sc,5] <- S[[sc-1]]$mean$alpha0}, error = function(e){})
  tryCatch({Out[z,sc,6] <- S[[sc-1]]$mean$alpha1}, error = function(e){})
  
  tryCatch({Out[z,sc+12,1] <- S[[sc-1]]$Rhat$Ntot}, error = function(e){})
  tryCatch({Out[z,sc+12,2] <- S[[sc-1]]$Rhat$beta0}, error = function(e){})
  tryCatch({Out[z,sc+12,3] <- S[[sc-1]]$Rhat$beta1}, error = function(e){})
  tryCatch({Out[z,sc+12,4] <- S[[sc-1]]$Rhat$sigma}, error = function(e){})
  tryCatch({Out[z,sc+12,5] <- S[[sc-1]]$Rhat$alpha0}, error = function(e){})
  tryCatch({Out[z,sc+12,6] <- S[[sc-1]]$Rhat$alpha1}, error = function(e){})
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
