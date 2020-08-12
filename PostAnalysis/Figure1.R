#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("~/ISDM/PostAnalysis")

#-----------#
#-Libraries-#
#-----------#

library(dplyr)
library(RandomFields)
library(abind)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(extrafont)
loadfonts(quiet = TRUE)

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

rotate <- function(x)
{
  t(apply(x, 2, rev))
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
lambda0 <- log(0.5)
#Effect parameter of enivornment on intensity
beta1 <- 1
#Intercept parameter for prensence only (PO) detection
p0.L <- 0.5
p0.S <- 0.1
#Effect parameter of environment on PO detection
alpha1 <- 1
#Scale parameter for DS detection
sigma <- 0.65

#-------------------------------------#
#-Draw environmental covariate values-#
#-------------------------------------#

load(file = "~/ISDM/DataAnalysis/Simulations/wcov.Rdata")

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
x.line <- rep(tu.m, each = W)
y.line <- rep(seq(0.5, W-0.5, 1), 10)
line <- cbind(x.line, y.line)
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
load(file = "~/ISDM/DataAnalysis/Simulations/dist.Rdata")

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

load(file = "~/ISDM/DataAnalysis/Simulations/zcov.Rdata")

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


#--------------#
#-Figures 1A,B-#
#--------------#

#Locations of PO data
uxpo.L <- gr[pixpo.L,1]
uypo.L <- gr[pixpo.L,2]
uxpo.S <- gr[pixpo.S,1]
uypo.S <- gr[pixpo.S,2]

pxtype <- list(uxpo.L, uxpo.S)
pytype <- list(uypo.L, uypo.S)

#Create dataframe for Figures 1A,B
tst2 <- data.frame(gr, w)
colnames(tst2) <- c("X", "Y", "Covariate")
tst3 <- data.frame(pxtype[[1]], pytype[[1]])
colnames(tst3) <- c("X3","Y3")
tst4 <- data.frame(XD[[1]][XD[[1]][,3]==1,4:5])
colnames(tst4) <- c("X4","Y4")
tst5 <- list()
for(j in 1:ns[1]){
  tst5[[j]] <- as.data.frame(line[which((line[,1] > (grt[tsamp[j,1],1] - 5)) & 
                                          ((grt[tsamp[j,1],1] + 5) > line[,1]) & 
                                          (line[,2] > (grt[tsamp[j,1],2] - 5)) & 
                                          ((grt[tsamp[j,1],2] + 5) > line[,2])),])}

Fig1A <- ggplot(tst2, aes(x=X, y=Y)) + 
  geom_raster(aes(fill = Covariate)) +
  geom_point(data = tst3, aes(x=X3,y=Y3), col = "#006699", size = 1.125) +
  geom_point(data = tst4, aes(x=X4,y=Y4), col = "#990000", size = 1.125) +
  scale_fill_gradient2(low = "#FFFFFF", mid = "#CCCCCC", high = "#999999", 
                       labels = c(-2,0,2), breaks = c(-2,0,2),
                       limits=c(-2.5,2.5), oob = scales::squish,
                       guide = guide_colorbar(barheight = 0.5, frame.colour = "black", ticks.colour = "black")) +
  scale_x_continuous(limits = c(0,100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  ggtitle("High Presence-only") +
  theme(text = element_text(family = "Times New Roman", size = 16),
        plot.title = element_text(hjust = 0.5, family = "Times New Roman", size = 16),
        axis.title=element_blank(),
        axis.line=element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_blank(),
        legend.key.width=unit(0.5, "in"),
        legend.position="bottom",
        legend.justification="center",
        legend.margin = margin(-10,0,10,0),
        plot.margin = unit(c(-0.25,0,0,0), "in"))

for(i in 1:20){
  Fig1A <- Fig1A + geom_line(data = tst5[[i]], aes(x=x.line,y=y.line), col = "#990000", size = 1)
}

Fig1A <- ggplotGrob(Fig1A)

tst3 <- data.frame(pxtype[[2]], pytype[[2]])
colnames(tst3) <- c("X3","Y3")

Fig1B <- ggplot(tst2, aes(x=X, y=Y)) + 
  geom_raster(aes(fill = Covariate)) +
  geom_point(data = tst3, aes(x=X3,y=Y3), col = "#006699", size = 1.125) +
  geom_point(data = tst4, aes(x=X4,y=Y4), col = "#990000", size = 1.125) +
  scale_fill_gradient2(low = "#FFFFFF", mid = "#CCCCCC", high = "#999999", 
                       labels = c(-2,0,2), breaks = c(-2,0,2),
                       limits=c(-2.5,2.5), oob = scales::squish,
                       guide = guide_colorbar(barheight = 0.5, frame.colour = "black", ticks.colour = "black")) +
  scale_x_continuous(limits = c(0,100), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0,100), expand = c(0, 0)) +
  ggtitle("Low Presence-only") +
  theme(text = element_text(family = "Times New Roman", size = 16),
        plot.title = element_text(hjust = 0.5, family = "Times New Roman", size = 16),
        axis.title=element_blank(),
        axis.line=element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        legend.title=element_blank(),
        legend.key.width=unit(0.5, "in"),
        legend.position="bottom",
        legend.justification="center",
        legend.margin = margin(-10,0,10,0),
        plot.margin = unit(c(-0.25,0,0,0), "in"))

for(i in 1:20){
  Fig1B <- Fig1B + geom_line(data = tst5[[i]], aes(x=x.line,y=y.line), col = "#990000", size = 1)
}

Fig1B <- ggplotGrob(Fig1B)

#-------------#
#-Figure 1C-F-#
#-------------#

#List of all filenames
filenames <- list.files(path = "~/ISDM/DataAnalysis/Simulations/SimulationOutput", pattern = "output", full.names = TRUE)

#Load first file
load(filenames[1])

#Initialize vector for all output
#####
#Out[,1,] are true values
#Out[,2,] estimates for PO large DS 0%
#Out[,3,] estimates for PO large DS 5%; independent
#Out[,4,] estimates for PO large DS 5%; dependent
#Out[,5,] estimates for PO large DS 10%; independent
#Out[,6,] estimates for PO large DS 10%; dependent
#Out[,7,] estimates for PO large DS 15%; independent
#Out[,8,] estimates for PO large DS 15%; dependent                                                                                   
#Out[,9,] estimates for PO large DS 20%; independent
#Out[,10,] estimates for PO large DS 20%; dependent
#Out[,11,] estimates for PO large DS 20%; dependent alt
#Out[,12,] estimates for PO small DS 0%
#Out[,13,] estimates for PO small DS 5%
#Out[,14,] estimates for PO small DS 10%
#Out[,15,] estimates for PO small DS 15%
#Out[,16,] estimates for PO small DS 20%
#Out[,17,] Rhat for PO large DS 0%
#Out[,18,] Rhat for PO large DS 5%; independent
#Out[,19,] Rhat for PO large DS 5%; dependent
#Out[,20,] Rhat for PO large DS 10%; independent
#Out[,21,] Rhat for PO large DS 10%; dependent
#Out[,22,] Rhat for PO large DS 15%; independent
#Out[,23,] Rhat for PO large DS 15%; dependent
#Out[,24,] Rhat for PO large DS 20%; independent
#Out[,25,] Rhat for PO large DS 20%; dependent
#Out[,26,] Rhat for PO large DS 20%; dependent alt
#Out[,27,] Rhat for PO small DS 0%
#Out[,28,] Rhat for PO small DS 5%
#Out[,29,] Rhat for PO small DS 10%
#Out[,30,] Rhat for PO small DS 15%
#Out[,31,] Rhat for PO small DS 20%
#####
Out <- output$Out
Out2 <- output$Out2

#Time vector
Time <- output$Time

#Harvest parameters from files and remove model runs with Rhat > 1.1
for(i in 2:length(filenames)){
  load(filenames[i])
  for(j in 1:length(output$Out[,1,1])){
    if(max(output$Out[j,c(18,20,22,24,28:31),2:6], na.rm = TRUE) < 1.1){
    #if(max(output$Out[j,c(18:25,27:31),2:6], na.rm = TRUE) < 1.1){
      Out <- abind(Out, output$Out[j,,], along = 1)
      Out2 <- rbind(Out2, output$Out2[j,])
    }
  }
  Time <- c(Time, output$Time)
}

#Remove first sample if Rhat > 1.1
if(max(Out[1,c(18,20,22,24,28:31),2:6], na.rm = TRUE) < 1.1){
  Out <- Out[-1,,]
  Out2 <- Out2[-1,]
  }

#Sample 1000 iterations
set.seed(123)
iter <- sort(sample(dim(Out)[1], 1000, replace = FALSE))
Out <- Out[iter,,]

name <- c("0", "5", "10", "15", "20")
name <- as.character(name)
name <- factor(name, levels=unique(name))

truth <- Out[,rep(1,10),]
truth[,1:5,5] <- 0.5
truth[,6:10,5] <- 0.1

y75 <- apply(Out[,c(2,3,5,7,9,12:16),] - truth, MARGIN = c(2,3), FUN = quantile, probs = 0.75, na.rm = TRUE)
y50 <- apply(Out[,c(2,3,5,7,9,12:16),] - truth, MARGIN = c(2,3), FUN = quantile, probs = 0.5, na.rm = TRUE)
y25 <- apply(Out[,c(2,3,5,7,9,12:16),] - truth, MARGIN = c(2,3), FUN = quantile, probs = 0.25, na.rm = TRUE)

ymax <-  ((y75 - y25) * 1.5) + y75
ymin <- y25 - ((y75 - y25) * 1.5)

lambda0 <- data.frame(ymax[,2], y75[,2], y50[,2], y25[,2], ymin[,2])
colnames(lambda0) <- c("ymax", "y75", "y50", "y25", "ymin")

beta1 <- data.frame(ymax[,3], y75[,3], y50[,3], y25[,3], ymin[,3])
colnames(beta1) <- c("ymax", "y75", "y50", "y25", "ymin")

p0 <- data.frame(ymax[,5], y75[,5], y50[,5], y25[,5], ymin[,5])
colnames(p0) <- c("ymax", "y75", "y50", "y25", "ymin")

alpha1 <- data.frame(ymax[,6], y75[,6], y50[,6], y25[,6], ymin[,6])
colnames(alpha1) <- c("ymax", "y75", "y50", "y25", "ymin")

Fig1C <- ggplotGrob(ggplot(lambda0[c(1:5),], aes(name)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", fill = "#A6A6A6", size = 0.75) +
                      geom_hline(yintercept = 0, col = "black", size = 1) +
                      coord_cartesian(ylim = c(-1.25, 1)) +
                      theme_few() +
                      theme(text = element_text(family = "Times New Roman", size = 16),
                            axis.title.y = element_text(margin = margin(0,-3,0,0)),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            plot.margin = unit(c(-0.25,0.0625,-0.1875,0), "in")) +
                      labs(y = "Intercept", x = ""))

Fig1D <- ggplotGrob(ggplot(lambda0[c(6:10),], aes(name)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", fill = "#A6A6A6", size = 0.75) +
                      geom_hline(yintercept = 0, col = "black", size = 1) +
                      coord_cartesian(ylim = c(-1.25, 1)) +
                      theme_few() +
                      theme(text = element_text(family = "Times New Roman", size = 16),
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank(),
                            plot.margin = unit(c(-0.25,0.0625,-0.1875,-0.0625), "in")) +
                      labs(y = "", x = ""))

Fig1E <- ggplotGrob(ggplot(beta1[c(1:5),], aes(name)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", fill = "#A6A6A6", size = 0.75) +
                      geom_hline(yintercept = 0, col = "black", size = 1) +
                      coord_cartesian(ylim = c(-0.175,0.175)) +
                      theme_few() +
                      theme(text = element_text(family = "Times New Roman", size = 16),
                            axis.title.y = element_text(margin = margin(0,-3,0,0)),
                            plot.margin = unit(c(-0.25,0.0625,-0.1875,0), "in")) +
                      labs(y = "Effect", x = ""))

Fig1F <- ggplotGrob(ggplot(beta1[c(6:10),], aes(name)) + 
                      geom_boxplot(aes(ymin = ymin, lower = y25, middle = y50, upper = y75, ymax = ymax), 
                                   stat = "identity", fill = "#A6A6A6", size = 0.75) +
                      geom_hline(yintercept = 0, col = "black", size = 1) +
                      coord_cartesian(ylim = c(-0.175,0.175)) +
                      theme_few() +
                      theme(text = element_text(family = "Times New Roman", size = 16),
                            plot.margin = unit(c(-0.25,0.0625,-0.1875,-0.0625), "in")) +
                      labs(y = "", x = ""))

#Set widths for Figures 1F,C
Fig1F$widths <- Fig1D$widths
Fig1C$widths <- Fig1E$widths

#Add letters to figures
Figure1A <- arrangeGrob(Fig1A, top = grid::textGrob("A", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = -0.5625, hjust = -0.5,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))
Figure1B <- arrangeGrob(Fig1B, top = grid::textGrob("B", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = -0.5625, hjust = -0.5,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))
Figure1C <- arrangeGrob(Fig1C, top = grid::textGrob("C", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = -3.125,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))
Figure1D <- arrangeGrob(Fig1D, top = grid::textGrob("D", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = -3.25,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))
Figure1E <- arrangeGrob(Fig1E, top = grid::textGrob("E", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = -3.375,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))
Figure1F <- arrangeGrob(Fig1F, top = grid::textGrob("F", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 0, hjust = -3.6875,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))

#X axis label
axisy <- grid::textGrob("Distance sampling coverage \n (percent of study area)", gp=grid::gpar(fontsize = 16, fontfamily = "Times New Roman"))
axisy$vjust <- 0.375

#Save Figure 1
tiff(file = "~/ISDM/PostAnalysis/Figure1.tiff", res = 600, width = 6.5, height = 8, units = "in")
grid.arrange(arrangeGrob(Figure1A, Figure1B, ncol = 2, nrow = 1), 
             arrangeGrob(Figure1C, Figure1D, Figure1E, Figure1F, ncol = 2, nrow = 2),
             bottom = axisy,
             heights = c(1,1.25))
dev.off()
