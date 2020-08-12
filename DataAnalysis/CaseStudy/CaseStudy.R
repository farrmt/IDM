#---------------------------------------------#
#----Integrated Species Distribution Model----#
#----Case study of black-backed jackal in-----#
#----the Masai Mara National Reserve.---------#
#----Created by Matthew Farr------------------#
#---------------------------------------------#

#-----------------------#
#-Set Working Directory-#
#-----------------------#

setwd("~/IDM/DataAnalysis/CaseStudy")

#-----------#
#-Libraries-#
#-----------#

library(jagsUI)

#-----------#
#-Load Data-#
#-----------#

load(file = "~/IDM/DataFormatting/IDM.Rdata")

#-------------#
#-Attach Data-#
#-------------#

attach(IDMdata)

#--------------#
#-Compile Data-#
#--------------#

data <- list(y = y.po, x = y.ds, dst = dst, w = Bias,
             D = D, Dend = Dend, Bend = Bend, G = G, TT = TT,
             c1 = c1, c2 = c2, c3 = c3, scale = scale,
             Gstart = Gstart, Gend = Gend, cover = cover, border = border[1:555,1], region = region,
             water = water[1:555,1], Lion = Lion[1:555,1:21], NDVI = NDVI[1:555,1:21])

#--------------------#
#-Parameters to Save-#
#--------------------#

params <- c("gamma0", "gamma1", "p0", "alpha1",
            "beta6", "beta5", "beta4", "beta3", "beta2", "beta1", "log_lambda0",
            "tau", "tauz", "Density")

#---------------#
#-Inital Values-#
#---------------#

inits <- function(){list(log_lambda0 = exp(runif(1, 0.0001, 0.001)),
                         beta1 = runif(1, 0.8, 1), beta2 = runif(1, 1.5, 2),
                         beta3 = runif(1, -1.4, -1), beta4 = runif(1, 0.1, 0.4),
                         beta5 = runif(1, -0.1, 0), beta6 = runif(1, -0.1, 0),
                         gamma1 = runif(1, -1, 0),gamma0 = runif(1, 5.3, 5.4), 
                         alpha1 = runif(1, 12, 14), p0 = runif(1, 0.25, 0.5),
                         sig = runif(1, 1, 3), sigz = runif(1, 1, 3))}

#---------------#
#-MCMC settings-#
#---------------#

nb <- 5000
ni <- 20000
nt <- 5
nc <- 3
na <- 1000

#-----------#
#-Run Model-#
#-----------#

#WARNING: Model takes ~65 hours to converge on parrallel computing. 
#This model was run on CentOS7 Linux using a high performance computer.

IDM <- jagsUI(data, inits, params, "IDM.txt", n.thin=nt, 
               n.chains=nc, n.burnin=nb, n.iter=ni, n.adapt=na,  parallel = TRUE)

save(IDM, file = "CaseStudy.Rdata")
