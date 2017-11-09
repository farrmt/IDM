library(mvtnorm)
library(rgeos)
library(sp)
library(raster)
library(jagsUI)

expit = function(eta) {1/(1+exp(-eta))}
logit = function(pp) { log(pp) - log(1-pp) }

# MF 10,000 pixels total

W <- 100
cell <- 1               # '2D bin width'
grx <- seq(cell/2, W - cell/2, cell) # mid-point coordinates
gr <- expand.grid(grx,grx)         # Create grid coordinates

beta1 <- 0.25
beta0 <- - 2.5
alpha1 <- 5
alpha0 <- - 5
sigma <- 2

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

#probs <- exp(beta1*xcov + beta0)/sum(exp(beta1*xcov + beta0)) # probability of point in pixel (sum = 1)
probs <- exp(beta1*xcov)/sum(exp(beta1*xcov))
#N <- rpois(1, sum(exp(beta1*xcov + beta0)))
N <- rpois(1, sum(exp(beta1*xcov)))
pixel.id <- sample(1:(W*W), N, replace=TRUE, prob=probs)
u1 <- gr[pixel.id,1]
u2 <- gr[pixel.id,2]

line <- structure(list(x = rep(c(25,50,75), each = W), y = rep(seq(0.5, W-0.5, 1), 3), .Names = c("x", "y")))
line <- cbind(line$x, line$y)
line <- as.matrix(line)

X <- line[,1]
Y <- line[,2]
J <- length(line[,1])

d <- array(NA, dim = c(W*W, J))
dist <- NULL
for(g in 1:10000){
  for(j in 1:J){
    d[g,j] <- sqrt((gr[g,1] - X[j])^2 + (gr[g,2] - Y[j])^2)
  }
  dist[g] <- min(d[g,])
}

p.ds <- rep(NA, N)
dst <- dist[pixel.id]
y.ds <- NULL
for(i in 1:N){
  p.ds[i] <- exp(-dst[i] * dst[i] / (2 * sigma * sigma))
  y.ds[i] <- rbinom(1, 1, p.ds[i])
}

uxds <- u1[y.ds == 1]
uyds <- u2[y.ds == 1]
dst <- dst[y.ds == 1]

mu.x <- 0.65*(100)
mu.y <- 0.65*(100)
sigma.x <- 0.15*abs(100)
sigma.y <- 0.30*abs(100)
rho.xy <- 0.1
mu <- c(mu.x, mu.y)
Sigmaxy <- matrix(c(sigma.x^2, rep(rho.xy*sigma.x*sigma.y, 2), sigma.y^2), ncol=2)

wcov <- dmvnorm(gr, mean=mu, sigma=Sigmaxy)
wcov <- (wcov - mean(wcov))/sd(wcov)
image(rasterFromXYZ(cbind(gr,wcov)), col=topo.colors(10))

#p.po <- expit(alpha1*wcov + alpha0)
p.po <- expit(alpha1*wcov)
image(rasterFromXYZ(cbind(gr,p.po)), col=topo.colors(20))
y.po <- NULL

for(i in 1:N){
  y.po[i] <- rbinom(1, 1, p.po[pixel.id[i]])
}

uxpo <- u1[y.po == 1]
uypo <- u2[y.po == 1]

par(mar=c(3,3,3,6))
image(rasterFromXYZ(cbind(u1,u2,p.po[pixel.id])), col=topo.colors(20))
par(mar=c(3,3,3,6))
image(rasterFromXYZ(cbind(uxpo,uypo,p.po[pixel.id[y.po == 1]])), col=topo.colors(20))

par(mar=c(3,3,3,6))
image(rasterFromXYZ(cbind(gr,xcov)), col=topo.colors(10))
points(u1, u2, pch=20, col='black', cex = 0.8)  # plot points
lines(X, Y, col = 'grey')
points(uxds, uyds, col = 'orange', pch = 20)
points(uxpo, uypo, col = 'red', pch = 20)

pixds <- pixel.id[y.ds == 1]
pixpo <- pixel.id[y.po == 1]

y.ds <- rep(1, sum(y.ds))
y.po <- rep(1, sum(y.po))

v <- 1
B <- 12
mdpt <- seq(0.5, 12, 1)
nD <- length(mdpt)
dclass <- as.numeric(dst + 0.5)

coverage <- which(dist <= 12)

str(data <- list(#y.ds = sum(y.ds), coverage = coverage,
  Habitat = xcov, G = 10000,
  wcov = wcov,  y.po = sum(y.po), npo = sum(y.po), pixpo = pixpo
  #nD = nD, v = v, B = B, mdpt = mdpt, dclass = dclass, nds = sum(y.ds)
))

cat("
    model{ 
    
    # Priors
    #sigma ~ dunif(0, 10)
    #beta0 ~ dnorm(0, 0.01)
    beta1 ~ dunif(0, 1)
    #alpha0 ~ dnorm(0, 0.01)
    alpha1 ~ dunif(0, 10)
    
    for(g in 1:G){   # g is the pixel index, there are G total pixels
    #intensity[g] <- exp(beta1*Habitat[g] + beta0)
    #logit(p.po[g]) <- alpha0 + wcov[g] * alpha1
    intensity[g] <- exp(beta1*Habitat[g])
    logit(p.po[g]) <- wcov[g] * alpha1
    v0[g] <- intensity[g] * p.po[g]
    }
    
    # for(i in 1:npo){
    # v0[i] <- intensity[pixpo[i]]*p.po[pixpo[i]]
    # }
    
    v1 <- sum(v0[])
    y.po ~ dpois(v1) 
    
    # for(k in 1:nD){
    # p[k] <- exp(-mdpt[k]*mdpt[k]/(2*sigma*sigma))
    # pi[k] <- v/B
    # f[k] <- p[k] * pi[k]
    # fc[k] <- f[k] / p.ds
    # }
    
    # p.ds <- sum(f[])
    # 
    # for(i in 1:nds){
    # dclass[i] ~ dcat(fc[1:nD])
    # }
    
    # y.ds ~ dbin(p.ds, N)
    # N ~ dpois(lambda)
    # lambda <- sum(intensity[coverage[]])
    Ntot <- sum(intensity[])
    }
    ",fill=TRUE,file="model1.txt")

N.st <- sum(y.ds) + 1

#inits <- function(){list(N = N.st, sigma=runif(1,1,3), beta1=runif(1,0,1), beta0=runif(1,-5,0), alpha0=runif(1,4,7), alpha1 = runif(1,0,3))}
inits <- function(){list(beta1=runif(1,0,1), alpha1 = runif(1,0,3))}
inits <- function(){list(N = N.st, sigma=runif(1,1,3), beta1=runif(1,0,1))}


params <- c("beta1", "sigma", "Ntot", "alpha1")

nb <- 2000
ni <- 12000
nt <- 1
nc <- 3
na <- 100

out2 <- jagsUI(data, inits, params, "model1.txt", n.thin=nt, 
               n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

