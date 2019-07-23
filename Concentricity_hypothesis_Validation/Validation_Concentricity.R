# Script to analize the distance between the center of the two circunferences that adjust the proteins.
# Author: Yasel Garces

# Load Libraries
library(ggplot2)
library(cowplot)
library(plyr)
library(gdata)

# Change Directory
setwd('/home/yasel/Dropbox/Paper Viroplasmas/Programs/Validation_concentricity')
# Load data
distance<-read.csv("DistanceCenterProtein.csv")
distance<-subset(distance,(Combination!="PDI") & (Combination!="dsRNA"))
distance$Combination<-drop.levels(distance$Combination)

# Boxplot considering the factor variable 
distance$Combination<-factor(distance$Combination,levels = levels(distance$Combination)[c(1,2,3,4,6,5,7,8)])
p<-ggplot(data = distance, aes(x = Combination, y = Distance, fill=Combination))+ geom_boxplot(notch = TRUE)+
  theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1, size=13),
        axis.text.y = element_text(hjust = 1, size=13))+
          geom_hline(yintercept = c(0.14142,0.052),color="red")+
  scale_y_continuous(breaks = sort(c(seq(min(distance$Distance), 
                                         max(distance$Distance), length.out=5), 0.14142,0.052)),
                     labels = function (x) round(x,3),
                     name=expression(paste("Displacement Distance"," ",(paste(mu,m)))))
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/ValidationDisplacementDistance.pdf",width = 5.5,height = 4)
p
dev.off()

# Successful number of concentric adjusts
ConcentricAdjust<-ddply(distance, c("Combination"), function(x) sum(x$Distance<0.14142))
ConcentricTotal<-ddply(distance, c("Combination"), function(x) length(x$Distance))

by(distance,
   distance[,"Combination"],
   function(x) {
     Median=median(x$Distance)
     qnt <- quantile(x$Distance, probs=c(.25, .75), na.rm = TRUE)
     H <- 1.57 * IQR(x$Distance)/sqrt(length(x$Distance))
     list(Median,CI=c(Median-H,Median+H))
   })
