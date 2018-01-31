#------------------------------#
#-Libraries used in simulation-#
#------------------------------#

library(mvtnorm)
library(jagsUI)
library(dplyr)

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
#-Values to save from simulation-#
#--------------------------------#

#Number of simulation iterations
iter <- 10

#True and estimated parameter values to save
Out <- array(NA, dim = c(iter, 5, 4), 
             dimnames = list(NULL, c("Truth", "DS", "DSalt", "DS.Rhat", "DSalt.Rhat"),
                             c("N", "beta0", "beta1", "sigma")))

#------------------#
#-Begin Simulation-#
#------------------#

start.time <- Sys.time()
for(z in 1:iter){
  
  #-----------------------#
  #-Draw parameter values-#
  #-----------------------#
  
  #Intercept parameter for intensity function
  beta0 <- runif(1, -1, -0.75)
  #Effect parameter of enivornment on intensity
  beta1 <- runif(1, 0.85, 1.15)
  #Scale parameter for DS detection
  sigma <- runif(1, 1.75, 2.25)
  
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
  y <- NULL
  #Simulate distance sampling
  p.ds <- exp(-dst * dst / (2 * sigma * sigma))
  y <- rbinom(N, 1, p.ds)
  y[dst>12] <- 0
  #Coordinates of detected individuals
  uxds <- u1[y == 1]
  uyds <- u2[y == 1]
  #Distance of detected individuals
  dst <- dst[y == 1]
  #Pixel ID
  pixds <- s[y == 1]
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
  y.ds <- sum(y)
  #Distance class of detected individuals
  dclass <- NULL
  for(i in 1:y.ds){
    for(k in 1:nD){
      if(mdpt[k] - 0.5 <= dst[i] && dst[i] < mdpt[k+1] - 0.5)
        dclass[i] <- k
    }
  }
  
  #-------------------------------#
  #-Alternative distance sampling-#
  #-------------------------------#
  
  regionC <- as.data.frame(coverage)
  colnames(regionC) <- "pixel"
  C <- length(coverage)
  pixel.ds <- s[s %in% coverage]
  N.ds <- length(pixel.ds)
  dstalt <- dist[coverage] 
  y.dsalt <- as.data.frame(table(pixds))
  y.dsalt$pixds <- as.numeric(as.character(y.dsalt$pixds))
  colnames(y.dsalt) <- c("pixel", "freq")
  tmp <- full_join(y.dsalt, regionC, by = "pixel")
  tmp <- tmp%>%arrange(pixel)
  y.dsalt <- tmp$freq
  y.dsalt[is.na(y.dsalt)] <- 0
  
  #----------------------------------#
  #-Compile BUGS data for each model-#
  #----------------------------------#
  
  #Distance sampling model
  str(dataDS <- list(x = x, G = (W*W), nD = nD, v = v, B = B, mdpt = mdpt, 
                     dclass = dclass, nds = y.ds, y.ds = y.ds, 
                     coverage = coverage, pixds = pixds))
  

  str(dataDSalt <- list(x = x, G = (W*W), C = C, dst = dstalt,
                     nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass, nds = y.ds,
                     y.ds = y.dsalt, coverage = coverage))
  
  #----------------#
  #-Initial values-#
  #----------------#
  
  #Inital value for N for distance sampling
  N.st <- y.ds + 1
  
  N.stalt <- y.dsalt + 1
  
  #Distance sampling model
  initsDS <- function(){list(N = N.st, sigma=runif(1,1,3), 
                             beta1 = runif(1, 0.75, 1), beta0 = runif(1, -1, -0.75))}
  
  initsDSalt <- function(){list(N = N.stalt, sigma=runif(1,1,3), 
                             beta1 = runif(1, 0.75, 1), beta0 = runif(1, -1, -0.75))}

  #------------#
  #-Parameters-#
  #------------#
  
  #Distance sampling model
  paramsDS <- c("Ntot", "beta0", "beta1", "sigma")
  
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

  DSalt <- jagsUI(dataDSalt, initsDSalt, paramsDS, "DSalt.txt", n.thin=nt,
               n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = FALSE)
  
  #----------------------------------#
  #-Save values from each simulation-#
  #----------------------------------#
  
  #True parameter values
  Out[z,1,1] <- N
  Out[z,1,2] <- beta0
  Out[z,1,3] <- beta1
  Out[z,1,4] <- sigma
  
  #Estimated distance sampling parameter values
  Out[z,2,1] <- DS$mean$Ntot
  Out[z,2,2] <- DS$mean$beta0
  Out[z,2,3] <- DS$mean$beta1
  Out[z,2,4] <- DS$mean$sigma
  
  Out[z,4,1] <- DS$Rhat$Ntot
  Out[z,4,2] <- DS$Rhat$beta0
  Out[z,4,3] <- DS$Rhat$beta1
  Out[z,4,4] <- DS$Rhat$sigma
  
  #Estimated presence only parameter values
  Out[z,3,1] <- DSalt$mean$Ntot
  Out[z,3,2] <- DSalt$mean$beta0
  Out[z,3,3] <- DSalt$mean$beta1
  Out[z,3,4] <- DSalt$mean$sigma

  Out[z,5,1] <- DS$Rhat$Ntot
  Out[z,5,2] <- DS$Rhat$beta0
  Out[z,5,3] <- DS$Rhat$beta1
  Out[z,5,4] <- DS$Rhat$sigma
  
}#End simulation
end.time <- Sys.time()
Time <- end.time - start.time

#-Save HPCC output-#
ID <- gsub(" ","_",Sys.time())
ID <- gsub(":", "-", ID)
output <- list(Out, Time)
heads <- c("Out", "Time")
output <- setNames(output, nm = heads)
save(output, file = paste("DS", ID, ".R", sep=""))