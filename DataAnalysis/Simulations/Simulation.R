#----------------------------------------------------#
#----Simulation study for IDM combining--------------#
#----distance sampling and opportunistic sampling----#
#----Created by Matthew Farr-------------------------#
#----------------------------------------------------#

#-----------#
#-Libraries-#
#-----------#

library(RandomFields)
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

# w <- RFsimulate(model = RMgauss(scale = 0.5), x = gr[,1], y = gr[,2], grid = FALSE, seed = 1)$variable1
# w <- (w - mean(w))/sd(w)
load(file = "wcov.Rdata")

#--------------------------------#
#-Values to save from simulation-#
#--------------------------------#

#Number of simulation iterations
iter <- 1

#True and estimated parameter values to save
Out <- array(NA, dim = c(iter, 31, 6), 
             dimnames = list(NULL, NULL, c("N", "lambda0", "beta1", "sigma", "p0", "alpha1")))

#uncertainty values to save
Out2 <- array(NA, dim = c(iter, 9))

#------------------#
#-Begin Simulation-#
#------------------#

start.time <- Sys.time()
for(i in 1:iter){
  
  #-----------------------#
  #-Draw parameter values-#
  #-----------------------#
  
  #Intercept parameter for intensity function
  lambda0 <- runif(1, 0.05, 1)
  #Effect parameter of enivornment on intensity
  beta1 <- runif(1, -1.25, 1.25)
  #Intercept parameter for prensence only (PO) detection
  p0.L <- 0.5
  p0.S <- 0.1
  #Effect parameter of environment on PO detection
  alpha1 <- runif(1, 0, 2)
  #Scale parameter for DS detection
  sigma <- runif(1, 0.75, 1.25)
  
  #--------------------------------#
  #-Simulate number of individuals-#
  #--------------------------------#
  
  #Intensity function
  intensity <- exp(log(lambda0) + beta1 * w)
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
  nts <- 4
  #Samples per scenario
  ns <- c(20,15,10,5)
  #Sampled transects for each scenario
  tsamp <- array(NA, dim = c(20,nts))
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
  fill <- seq(1,2100,by=100)
  coverage <- array(NA, dim = c(2000,nts))
  dst <- array(NA, dim = c(2000,nts))
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
  dst2 <- as.data.frame(cbind(coverage[,1], dst[,1]))
  colnames(dst2) <- c("p", "dst")
  x20 <- inner_join(dst2, sds, by = "p")
  dst3 <- x20$dst
  #Individual presence/absence based on distance sampling
  xds <- NULL
  #Simulate distance sampling
  pids <- exp(-dst3 * dst3 / (2 * sigma * sigma))
  xds <- rbinom(length(sds[,1]), 1, pids)
  x20$x <- xds
  x20$uxds <- gr[x20[,1],1]
  x20$uyds <- gr[x20[,1],2]
  x20$dclass <- NA
  x20 <- x20%>%arrange(p)
  
  tmp2 <- as.data.frame(coverage[,2])
  colnames(tmp2) <- c("p")
  x15 <- inner_join(x20,tmp2,by="p")
  x15 <- x15%>%arrange(p)
  
  tmp2 <- as.data.frame(coverage[,3])
  colnames(tmp2) <- c("p")
  x10 <- inner_join(x20,tmp2,by="p")
  x10 <- x10%>%arrange(p)
  
  tmp2 <- as.data.frame(coverage[,4])
  colnames(tmp2) <- c("p")
  x5 <- inner_join(x20,tmp2,by="p")
  x5 <- x5%>%arrange(p)
  
  XD <- list(x20,x15,x10,x5)
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
  
  # z <- RFsimulate(model = RMgauss(scale = 2), x = gr[,1], y = gr[,2], grid = FALSE, seed = 89)$variable1
  # z <- (z - mean(z))/sd(z)
  load(file = "zcov.Rdata")
  
  #----------------------------------#
  #-Simulate opportunistic surveying-#
  #----------------------------------#
  
  #Detection probability of PO
  p.L <- expit(logit(p0.L) + alpha1 * z)
  p.S <- expit(logit(p0.S) + alpha1 * z)
  
  #Individuals detected in PO
  y.L <- y.S <- NULL
  for(j in 1:N){
    y.L[j] <- rbinom(1, 1, p.L[s[j]])
    y.S[j] <- rbinom(1, 1, p.S[s[j]])
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
                    z = z,  y = y.L))
  
  #Scenario 2: PO Large DS 5%
  str(data2 <- list(w = w, G = (W*W), B = B[4], d = x[[4]]$dst,
                    x = x[[4]]$freq, coverage = x[[4]]$p,
                    z = z, y = y.L))
  
  #Scenario 3: PO Large DS 10%
  str(data3 <- list(w = w, G = (W*W), B = B[3], d = x[[3]]$dst,
                    x = x[[3]]$freq, coverage = x[[3]]$p,
                    z = z, y = y.L))
  
  #Scenario 4: PO Large DS 15%
  str(data4 <- list(w = w, G = (W*W), B = B[2], d = x[[2]]$dst,
                    x = x[[2]]$freq, coverage = x[[2]]$p,
                    z = z, y = y.L))
  
  #Scenario 5: PO Large DS 20%
  str(data5 <- list(w = w, G = (W*W), B = B[1], d = x[[1]]$dst,
                    x = x[[1]]$freq, coverage = x[[1]]$p,
                    z = z, y = y.L))
  
  #Scenario 6: PO Small DS 0%
  str(data6 <- list(w = w, G = (W*W),
                    z = z,  y = y.S))
  
  #Scenario 7: PO Small DS 5%
  str(data7 <- list(w = w, G = (W*W), B = B[4], d = x[[4]]$dst,
                    x = x[[4]]$freq, coverage = x[[4]]$p,
                    z = z, y = y.S))
  
  #Scenario 8: PO Small DS 10%
  str(data8 <- list(w = w, G = (W*W), B = B[3], d = x[[3]]$dst,
                    x = x[[3]]$freq, coverage = x[[3]]$p,
                    z = z, y = y.S))
  
  #Scenario 9: PO Small DS 15%
  str(data9 <- list(w = w, G = (W*W), B = B[2], d = x[[2]]$dst,
                    x = x[[2]]$freq, coverage = x[[2]]$p,
                    z = z, y = y.S))
  
  #Scenario 10: PO Small DS 20%
  str(data10 <- list(w = w, G = (W*W), B = B[1], d = x[[1]]$dst,
                     x = x[[1]]$freq, coverage = x[[1]]$p,
                     z = z, y = y.S))
  
  #----------------#
  #-Initial values-#
  #----------------#
  
  #Scenario 1: PO Large DS 0%
  inits1 <- function(){list(beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05), 
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}
  
  #Scenario 2: PO Large DS 5%
  Nst2 <- function(){
    Nst2 <- data2$y + 1
    for(i in 1:data2$B){
      if(Nst2[data2$coverage[i]] < data2$x[i])
        Nst2[data2$coverage[i]] <- data2$x[i] + 1
    }
    return(Nst2)
  }
  
  inits2 <- function(){list(N = Nst2(), sigma = runif(1, sigma-0.05, sigma+0.05), 
                            beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}
  
  #Scenario 3: PO Large DS 10%
  Nst3 <- function(){
    Nst3 <- data3$y + 1
    for(i in 1:data3$B){
      if(Nst3[data3$coverage[i]] < data3$x[i])
        Nst3[data3$coverage[i]] <- data3$x[i] + 1
    }
    return(Nst3)
  }
  
  inits3 <- function(){list(N = Nst3(), sigma = runif(1, sigma-0.05, sigma+0.05), 
                            beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}
  
  #Scenario 4: PO Large DS 15%
  Nst4 <- function(){
    Nst4 <- data4$y + 1
    for(i in 1:data4$B){
      if(Nst4[data4$coverage[i]] < data4$x[i])
        Nst4[data4$coverage[i]] <- data4$x[i] + 1
    }
    return(Nst4)
  }
  
  inits4 <- function(){list(N = Nst4(), sigma = runif(1, sigma-0.05, sigma+0.05), 
                            beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}
  
  #Scenario 5: PO Large DS 20%
  Nst5 <- function(){
    Nst5 <- data5$y + 1
    for(i in 1:data5$B){
      if(Nst5[data5$coverage[i]] < data5$x[i])
        Nst5[data5$coverage[i]] <- data5$x[i] + 1
    }
    return(Nst5)
  }
  
  inits5 <- function(){list(N = Nst5(), sigma = runif(1, sigma-0.05, sigma+0.05), 
                            beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.L-0.05, p0.L+0.05))}
  
  #Scenario 6: PO Small DS 0%
  inits6 <- function(){list(beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05), 
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}
  
  #Scenario 7: PO Small DS 5%
  inits7 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                            beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}
  
  #Scenario 8: PO Small DS 10%
  inits8 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                            beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}
  
  #Scenario 9: PO Small DS 15%
  inits9 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                            beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                            alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}
  
  #Scenario 10: PO Small DS 20%
  inits10 <- function(){list(sigma = runif(1, sigma-0.05, sigma+0.05), 
                             beta1 = runif(1, beta1-0.05, beta1+0.05), lambda0 = runif(1, lambda0-0.05, lambda0+0.05),
                             alpha1 = runif(1, alpha1-0.05, alpha1+0.05), p0 = runif(1, p0.S-0.05, p0.S+0.05))}  
  #------------#
  #-Parameters-#
  #------------#
  
  params <- c("Ntot", "sigma", "p0", "alpha1", "lambda0", "beta1")
  
  #-------------#
  #-MCMC values-#
  #-------------#
  
  nb <- 1000
  ni <- 6000
  nt <- 1
  nc <- 3
  na <- 200
  
  #----------------#
  #-Run each model-#
  #----------------#
  
  S <- list()
  
  #Scenario 1: PO Large DS 0%
  print(1)
  S[[1]] <- jagsUI(data1, inits1, params, "PO.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 2.1: PO Large DS 5%; independent
  print(2.1)
  S[[2]] <- jagsUI(data2, inits2, params, "IDM1.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 2.2: PO Large DS 5%; dependent
  print(2.2)
  S[[3]] <- jagsUI(data2, inits2, params, "IDM2.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 3.1: PO Large DS 10%; independent
  print(3.1)
  S[[4]] <- jagsUI(data3, inits3, params, "IDM1.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 3.2: PO Large DS 10%; dependent
  print(3.2)
  S[[5]] <- jagsUI(data3, inits3, params, "IDM2.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 4.1: PO Large DS 15%; independent
  print(4.1)
  S[[6]] <- jagsUI(data4, inits4, params, "IDM1.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 4.2: PO Large DS 15%; dependent
  print(4.2)
  S[[7]] <- jagsUI(data4, inits4, params, "IDM2.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 5.1: PO Large DS 20%; independent
  print(5.1)
  S[[8]] <- jagsUI(data5, inits5, params, "IDM1.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 5.2: PO Large DS 20%; dependent
  print(5.2)
  S[[9]] <- jagsUI(data5, inits5, params, "IDM2.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 5.3: PO Large DS 20%; alternative structure dependent
  print(5.3)
  S[[10]] <- jagsUI(data5, inits5, params, "IDM3.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 6: PO Small DS 0%
  print(6)
  S[[11]] <- jagsUI(data6, inits6, params, "PO.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 7: PO Small DS 5%
  print(7)
  S[[12]] <- jagsUI(data7, inits7, params, "IDM1.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 8: PO Small DS 10%
  print(8)
  S[[13]] <- jagsUI(data8, inits8, params, "IDM1.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 9: PO Small DS 15%
  print(9)
  S[[14]] <- jagsUI(data9, inits9, params, "IDM1.txt", n.thin=nt, 
                   n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
  #Scenario 10: PO Small DS 20%
  print(10)
  S[[15]] <- jagsUI(data10, inits10, params, "IDM1.txt", n.thin=nt, 
                    n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na, parallel = TRUE)
  
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
  
  for(sc in 2:16){
    tryCatch({Out[i,sc,1] <- S[[sc-1]]$mean$Ntot}, error = function(e){})
    tryCatch({Out[i,sc,2] <- S[[sc-1]]$mean$lambda0}, error = function(e){})
    tryCatch({Out[i,sc,3] <- S[[sc-1]]$mean$beta1}, error = function(e){})
    tryCatch({Out[i,sc,4] <- S[[sc-1]]$mean$sigma}, error = function(e){})
    tryCatch({Out[i,sc,5] <- S[[sc-1]]$mean$p0}, error = function(e){})
    tryCatch({Out[i,sc,6] <- S[[sc-1]]$mean$alpha1}, error = function(e){})
    
    tryCatch({Out[i,sc+15,1] <- S[[sc-1]]$Rhat$Ntot}, error = function(e){})
    tryCatch({Out[i,sc+15,2] <- S[[sc-1]]$Rhat$lambda0}, error = function(e){})
    tryCatch({Out[i,sc+15,3] <- S[[sc-1]]$Rhat$beta1}, error = function(e){})
    tryCatch({Out[i,sc+15,4] <- S[[sc-1]]$Rhat$sigma}, error = function(e){})
    tryCatch({Out[i,sc+15,5] <- S[[sc-1]]$Rhat$p0}, error = function(e){})
    tryCatch({Out[i,sc+15,6] <- S[[sc-1]]$Rhat$alpha1}, error = function(e){})
  }
  
  for(kk in 1:9){
    Out2[i,kk] <- ifelse(S[[kk+1]]$q2.5$beta1 <= beta1 & S[[kk+1]]$q97.5$beta1 >= beta1, 1, 0)
  }
  
}#End simulation
end.time <- Sys.time()
Time <- end.time - start.time

#-Save HPCC output-#
ID <- gsub(" ","_",Sys.time())
ID <- gsub(":", "-", ID)
output <- list(Out, Time, Out2)
heads <- c("Out", "Time", "Out2")
output <- setNames(output, nm = heads)
save(output, file = paste("output", ID, ".R", sep=""))