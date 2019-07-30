#=======================================================================================
#' Study to validate the viral elements concentricity hypothesis based on 
#' the distance between the centers of the adjusted circumferences of both proteins.
#' @param DistanceCenterProtein.csv: File with the next structure:
#'             Column 1 (Combination): Name of the protein combined with NSP2
#'             Column 2 (Distance): Distance between the center of the adjusted 
#'                         circumferences to NSP2 and the accompanying protein.
#' @return Boxplot with the distance between the center of the ciircumferences by 
#' protein combination.
#' @return Table with the performance of the algorithm (see details below) 
#' @return Differences in population medians for all the combinations of proteins
# (difference in location) through a Mann-Whitney test.                        
#' @author Yasel Garces (88yasel@gmail.com)
#=======================================================================================
# Load Libraries
library(ggplot2)
library(cowplot)
library(plyr)
library(gdata)
theme_set(theme_cowplot()) # Change theme

# Set Directory
setwd('/home/yasel/TRABAJO/IBt/Viroplasms/GitHub Codes Paper  (No Mover)/Nanoscale_organization_of_rotavirus_replication_machineries/Concentricity_hypothesis_Validation')
# Load data
distance<-read.csv("DistanceCenterProtein.csv")

# Change the levels organization of the "Combination" factor variable.
distance$Combination<-factor(distance$Combination,levels = levels(distance$Combination)[c(1,2,3,4,6,5,7,8)])
# Boxplot with the distance between the center of the ciircumferences by protein combination. 
ggplot(data = distance, aes(x = Combination, y = Distance, fill=Combination))+ geom_boxplot(notch = TRUE)+
  theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1, size=13),
        axis.text.y = element_text(hjust = 1, size=13))+
          geom_hline(yintercept = c(0.14142,0.052),color="red")+
  scale_y_continuous(breaks = sort(c(seq(min(distance$Distance), 
                                         max(distance$Distance), length.out=5), 0.14142,0.052)),
                     labels = function (x) round(x,3),
                     name=expression(paste("Displacement Distance"," ",(paste(mu,m)))))

# When the distance between the centers in under 0.14142 microns, we assume that these 
# distributions are concentrics. See the paper for details. Based on this, it is possible 
# to compute the TP rate. 
ConcentricAdjust<-ddply(distance, c("Combination"), function(x) sum(x$Distance<0.14142))
ConcentricTotal<-ddply(distance, c("Combination"), function(x) length(x$Distance))
performance<-data.frame('Combination'=paste('NSP2',ConcentricAdjust$Combination,sep = '-'),
           'Success'=ConcentricAdjust$V1,'Total'=ConcentricTotal$V1)
performance

# Compute the differences in population medians for all the combinations of proteins
# (difference in location) through a Mann-Whitney test. 
by(distance,
   distance[,"Combination"],
   function(x) {
     Median=median(x$Distance)
     qnt <- quantile(x$Distance, probs=c(.25, .75), na.rm = TRUE)
     H <- 1.57 * IQR(x$Distance)/sqrt(length(x$Distance))
     list(Median,CI=c(Median-H,Median+H))
   })
