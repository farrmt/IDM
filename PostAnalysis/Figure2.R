#-------------------------------------------------------#
#----Script to create all subcomponents of Figure 2.----#
#----Created by Matthew Farr----------------------------#
#-------------------------------------------------------#

#-----------------------#
#-Set working directory-#
#-----------------------#

setwd("~/ISDM/PostAnalysis")

#-----------#
#-Libraries-#
#-----------#

library(cowplot)
library(raster)
library(ggplot2)
library(ggthemes)
library(gridExtra)
library(extrafont)
loadfonts(quiet = TRUE)

#------------#
#-Load Files-#
#------------#

#Model ouput from integrated model
load("~/ISDM/DataAnalysis/CaseStudy/CaseStudy.Rdata") #FIX FILE 
#Formatted data used in analysis
load("~/ISDM/DataFormatting/ISDM.Rdata")

#----------#
#-Figure2A-#
#----------#

#Load Figure 2A. Figure subcomponent was created in ArcGis.
Fig2A <- ggdraw() + 
  draw_image("~/ISDM/PostAnalysis/Figure2A.tif", scale = 1)

#-----------#
#-Figure 2B-#
#-----------#

#Load in empty raster to store density values
predD <- estD <- stack("../DataFormatting/ArcGIS/BaseData/tst2.tif")

#Mean across iterations 
Density <- apply(ISDM$sims.list$Density, MARGIN = 2, mean)

#Initialize predicted density
predDen <- NULL

for(i in 1:555){
  values(estD)[which(values(estD)==i)] <- Density[i]
  values(predD)[which(values(predD)==i)] <- predDen[i] <- exp(ISDM$mean$log_lambda0 +
                                                  ISDM$mean$beta1*ISDMdata$border[i] + 
                                                  ISDM$mean$beta2 * ISDMdata$region[i] +
                                                  ISDM$mean$beta3 * ISDMdata$region[i] * ISDMdata$border[i] +
                                                  ISDM$mean$beta4 * ISDMdata$water[i] +
                                                  ISDM$mean$beta5 * mean(ISDMdata$Lion[i,]) +
                                                  ISDM$mean$beta6 * mean(ISDMdata$NDVI[i,])
  ) * 400
}

#Estimated density
plot(estD)
writeRaster(estD, "EstDensity_Revised2.tif", overwrite=TRUE) #FIX file path

#Predicted density
plot(predD)
writeRaster(predD, "PredDensity_Revised2.tif", overwrite=TRUE) #FIX file path

#This file is now edited in ArcMap

#Load edited map
Fig2B <- ggdraw() + 
  draw_image("~/ISDM/PostAnalysis/Figure2B_Revised.tif", scale = 1)

#-----------#
#-Figure 2C-#
#-----------#

#Add mean and CIs to dataframe
beta.val <- cbind(c(ISDM$q2.5$beta1, ISDM$q2.5$beta2, ISDM$q2.5$beta3, ISDM$q2.5$beta4, ISDM$q2.5$beta6, ISDM$q2.5$beta5),
                  c(ISDM$q25$beta1, ISDM$q25$beta2, ISDM$q25$beta3, ISDM$q25$beta4, ISDM$q25$beta6, ISDM$q25$beta5),
                  c(ISDM$mean$beta1, ISDM$mean$beta2, ISDM$mean$beta3, ISDM$mean$beta4, ISDM$mean$beta6, ISDM$mean$beta5),
                  c(ISDM$q75$beta1, ISDM$q75$beta2, ISDM$q75$beta3, ISDM$q75$beta4, ISDM$q75$beta6, ISDM$q75$beta5),
                  c(ISDM$q97.5$beta1, ISDM$q97.5$beta2, ISDM$q97.5$beta3, ISDM$q97.5$beta4, ISDM$q97.5$beta6, ISDM$q97.5$beta5))

covnames <- c("Border", "Region", "Border x Region", "Water", "NDVI", "Lion")

values <- data.frame(covnames, beta.val)
colnames(values) <- c("covariate", "lower.alpha", "l25.alpha", "mean.alpha", "u75.alpha", "upper.alpha")

values$covariate <- factor(as.character(values$covariate), levels = unique(values$covariate))


Fig2C <- ggplotGrob(ggplot(values) + 
  geom_errorbar(aes(x = covariate, ymin = mean.alpha, ymax = mean.alpha),
                width = 0.375) +
  geom_errorbar(aes(x = covariate, ymin = lower.alpha, ymax = upper.alpha),
                width = 0, size = 1.25) +
  geom_errorbar(aes(x = covariate, ymin = l25.alpha, ymax = u75.alpha),
                width = 0, size = 3) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
  geom_hline(yintercept = 0, alpha = 0.75) +
  theme_few() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "in"),
        text = element_text(family = "Times New Roman", size = 16),
        panel.background = element_rect(fill = "transparent", color = NA),
        axis.text.x = element_text(angle = 50, hjust = 0.5, vjust = 0.5),
        legend.position = "none") +
  scale_x_discrete(labels = c("Border" = "Border", "Region" = "Regime", "Border x Region" = "Border x\nRegime",
                              "Water" = "Water", "NDVI" = "NDVI", "Lion" = "Lion")) +
  labs(y ="Covariate effect (log scale)", x = expression()))

#-----------#
#-Figure 2D-#
#-----------#

#Create distances between 0 and 1000 meters
dist <- seq(0, 1000, 1)
dist <- scale(dist)

#Initialize vectors for each region
disturbed <- matrix(NA, nrow = 1001, ncol = 3000)
undisturbed <- matrix(NA, nrow = 1001, ncol = 3000)

#Generate predicted values
for(i in 1:1001){
  for(k in 1:3000){
    disturbed[i,k] <- exp(ISDM$sims.list$log_lambda0[k] + ISDM$sims.list$beta2[k] + ISDM$sims.list$beta1[k] * dist[i] + ISDM$sims.list$beta3[k] * dist[i]) * 400
    undisturbed[i,k] <- exp(ISDM$sims.list$log_lambda0[k] + ISDM$sims.list$beta1[k] * dist[i]) * 400
  }
}

#Summarize over MCMC iterations
values2B <- data.frame(apply(disturbed, MARGIN = 1, mean), t(apply(disturbed, MARGIN = 1, quantile, probs = c(0.025, 0.975))),
                   apply(undisturbed, MARGIN = 1, mean), t(apply(undisturbed, MARGIN = 1, quantile, probs = c(0.025, 0.975))))
colnames(values2B) <- c("mean-1", "q2.5-1", "q97.5-1", "mean-2", "q2.5-2", "q97.5-2")
dist <- seq(0, 1000, 1)

Fig2D <- ggplotGrob(ggplot(values2B) +
  geom_ribbon(aes(x = dist, ymin = values2B$'q2.5-1', ymax = values2B$'q97.5-1'), alpha = 0.5, fill = "#CCCCCC") +
  geom_line(aes(x = dist, y = values2B$'mean-1'), color = "black", size = 1.25) +
  geom_ribbon(aes(x = dist, ymin = values2B$'q2.5-2', ymax = values2B$'q97.5-2'), alpha = 0.75, fill = "#CCCCCC") +
  geom_line(aes(x = dist, y = values2B$'mean-2'), linetype = 2, color = "black", size = 1.25) +
  labs(y = expression(Density~km^2), x = "Distance from border") +
  scale_x_continuous(limits = c(0,1000), expand = c(0, 0)) +
  theme_few() +
  theme(plot.margin = unit(c(0, 0, 0, 0.125), "in"),
        panel.background = element_rect(fill = "transparent", color = NA),
        text = element_text(family = "Times New Roman", size = 16),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = "none"))

#Set heights of Figure 2D
Fig2D$heights <- Fig2C$heights
#Move x-axis title
Fig2D$layout[12,1] <- 8

#----------#
#-Figure 2-#
#----------#

#Add letters to figures
Figure2A <- arrangeGrob(Fig2A, top = grid::textGrob("A", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 1.125, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))
Figure2B <- arrangeGrob(Fig2B, top = grid::textGrob("B", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 1.125, hjust = 0,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))
Figure2C <- arrangeGrob(Fig2C, top = grid::textGrob("C", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 1.75, hjust = -3.0625,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))
Figure2D <- arrangeGrob(Fig2D, top = grid::textGrob("D", x = unit(0, "in"), 
                                                    y = unit(0, "in"), just=c("left","top"), vjust = 1.75, hjust = -4.75,
                                                    gp=grid::gpar(fontsize=18, fontfamily = "Times New Roman", fontface = 2)))

#Layout for panels of Figure 2
lay = rbind(rep(1,26), rep(1,26), rep(2,26), rep(2,26), rep(2,26), rep(2,26), c(NA, rep(3,14), rep(4,10), NA), c(NA, rep(3,14), rep(4,10), NA), c(NA, rep(3,14), rep(4,10), NA))

#Save Figure 2
tiff(file = "~/ISDM/PostAnalysis/Figure2_revised2.tiff", res = 600, width = 6.5, height = 9, units = "in")
grid.arrange(Figure2A, Figure2B, Figure2C, Figure2D, layout_matrix = lay)
dev.off()
