#' This script analyses the distribution of the viral proteins into 
#' the viroplasm taking NSP5 as the reference protein. For details about 
#' the full research consult the below article:
#' Garcés et al. Nanoscale organization of rotavirus replication machineries. 
#' eLife 2019;8:e42906. https://elifesciences.org/articles/42906, 
#' doi: 10.7554/eLife.42906.
#'
#' @param This script uses the data collected for each protein combination that
#' are stored in the csv files:
#' 1- NSP5vsVP4.csv
#' 2- NSP5vsVP6.csv
#' This files have the following structure:
#' Column 1: Distance between the distribution of both proteins.
#' Column 2: Radius of the circumference that adjust the central protein.
#' Column 3: Radius of the circumference that adjust the other protein.
#' @return A set of graphics and statistics. See below for details.
#' **Note** This script is practically the same than "AnalysisNSP2.R" but 
#' as a difference it include the regression 
#' model analysis.
#' @author Yasel Garces (88yasel@gmail.com)  
#===================================================================================
#===================================================================================
## FUNCTIONS 
# Load the functions saved in Functions.R
source('/home/yasel/TRABAJO/IBt/Viroplasms/GitHub Codes Paper  (No Mover)/Nanoscale_organization_of_rotavirus_replication_machineries/R Codes/Functions.R')
#===================================================================================
#===================================================================================
# # Load Libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(ggsignif)
theme_set(theme_cowplot()) # Change theme
# Set work directory
setwd('/home/yasel/TRABAJO/IBt/Viroplasms/GitHub Codes Paper  (No Mover)/Nanoscale_organization_of_rotavirus_replication_machineries/R Codes/ResultsCSV/NSP5/')
# Load data
VP4<-read.csv('NSP5vsVP4.csv')
VP6<-read.csv('NSP5vsVP6.csv')
## Data manipulation
# Add a new variable with the name of the NSP5 accompanying protein
VP4<-mutate(VP4, Protein='VP4')
VP6<-mutate(VP6, Protein='VP6')

# Merge the data and convert the "Protein" column to a factor variable
viroData<-rbind(VP6,VP4)
viroData$Protein<-as.factor(viroData$Protein)
# Convert from pixels to microns (this is based on our experimetal design, for 
# other experiments you need to take care about how to do this conversion).
viroData$Distance=viroData$Distance/100
viroData$ratioNSP5=viroData$ratioNSP5/100
viroData$ratioOther=viroData$ratioOther/100
#===================================================================================
### Exploratory analysis of the results obtained by the algorithm VPs-DLSFC ###
# Plot a histogram with the number of samples per condition.
p<-ggplot(viroData,aes(Protein))+geom_histogram(stat = "count", fill=I("red"), 
                                                col=I("red"), 
                                                alpha=I(.2))+
  theme(legend.position = "none",axis.text.x = element_text(vjust = 0.5))+
  geom_text(stat = "count", aes(label = ..count.., y = ..count..+3))
# Save the graphic
pdf(file = "HistogramNSP5.pdf",width = 2,height = 2)
p
dev.off()
#===================================================================================
# VP6 and VP4 spatial distribution taking NSP5 as reference protein.
# Create a convinient data frame that allow to graphic with a boxplot the 
# radii's distribution of all the combinations off proteins.
allProteins<-mutate(viroData,Comparation="Other Protein")
forNSP5<-select(allProteins,one_of(c("ratioNSP5","Protein")))
forNSP5<-mutate(forNSP5,Comparation="NSP5")
allProteins<-data.frame(Comparation=c(forNSP5$Comparation,allProteins$Comparation),
                        Protein=c(as.character(forNSP5$Protein),as.character(allProteins$Protein)),
                        Value=c(forNSP5$ratioNSP5,allProteins$ratioOther))
# Orders the factor levels
allProteins$Protein<-factor(allProteins$Protein,levels = levels(allProteins$Protein)[c(2,1,3)])
# Boxplot for the radius of the fitting circumferences. In each experimental condition
# (NSP5-VP4, NSP5-VP6) we plot two boxes, the red box is for the radius of NSP5 (reference protein), 
# and the blue boxes represent the radius of the accompanying VP components (names)
p<-ggplot(allProteins, aes(x = Protein, y = Value,fill=Comparation)) +
  geom_boxplot(notch=TRUE)+ylab(expression(
    paste("Radius of the adjusted circumference"," ", (paste(mu,m)))))+
  theme(legend.position = "bottom")+
  scale_y_continuous(breaks = seq(from = 0.2, to = 1.5, by =0.1),limits = c(0.2,1.2))
p
# Save the graphic
pdf(file = "Radius_With_NSP5.pdf",width = 6,
    height = 5)
p
dev.off()
# Filters the "allProtein" data frame to only consider the radii of NSP5
NSP5compare<-subset(allProteins,Comparation=="NSP5")
# Create a boxplot with the NSP5's radii and save the graphic.
pdf(file = "Radius_Only_NSP5.pdf",width = 3,
    height = 3)
ggplot(NSP5compare,aes(x = Protein, y = Value))+geom_boxplot(notch = TRUE)+
  geom_signif(test="wilcox.test", comparisons = combn(levels(NSP5compare$Protein),2, simplify = F)[-4],
              step_increase = 0.2,map_signif_level=TRUE)+
  theme(legend.position = "none")
dev.off()
#===================================================================================
# Boxplot for the distance between NSP5 and {VP6,VP4}.
# Orders the factor levels
viroData$Protein<-factor(viroData$Protein,levels = levels(viroData$Protein)[c(2,1,3)])
p<-ggplot(viroData, aes(x = Protein, y = Distance,fill=Protein)) +
  geom_boxplot(notch=TRUE,show.legend=FALSE)+
  ylab(expression(paste("Distance to NSP5"," ",(paste(mu,m)))))+
  xlab("Protein")+scale_y_continuous(breaks = seq(from = -.5, to = .5, by =0.1),limits = c(0,.5))+
  geom_signif(test="wilcox.test", comparisons = combn(levels(viroData$Protein),2, simplify = F)[-5],
              step_increase = 0.05,map_signif_level=TRUE)
p
# Save figure
pdf(file = "Distance_With_NSP5.pdf",width = 5,
    height = 4)
p
dev.off()
#===================================================================================
# Linear regression fitting taking the radius of NSP5 as independent variable
# and the radius of VP6 and VP4 as dependent variable. The gray shadow represents the 
# confidence interval at a level of 95%.
p<-ggplot(viroData, aes(x = ratioNSP5, y = ratioOther,color=Protein)) + geom_point(show.legend=FALSE) + 
  facet_grid(~ Protein)+ xlab("Ratio NSP2")+
  geom_smooth(method='lm',formula=y~x-1,show.legend=FALSE)+
  scale_y_continuous(name=expression(paste("Protein Radius"," ", 
                                           (paste(mu,m)))), 
                     breaks=seq(0, 0.9,0.1),limits = c(0.3,0.9))+
  scale_x_continuous(name=expression(paste("Radius NSP5"," ", 
                                           (paste(mu,m)))))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5))
p
pdf(file = "LR_With_NSP5.pdf",width = 5,
    height = 4)
p
dev.off()
#===================================================================================
# Residuals errors for each linear regression model. The gray line represent the
# regression model, the points over the line are the predicted values,
# and the dots filled with a color gradient are the real values. The errors between 
# the predicted and real values, are represented as a gradient of colors as follows 
# (from lowest/coldest to highest/warmest).
# NSP2 vs VP6
data<-data.frame(x=VP6$ratioNSP5/100,y=VP6$ratioOther/100)
LM_VP6<-LMbyProtein(data,"VP6","NSP5")
# NSP5 vs VP4
data<-data.frame(x=VP4$ratioNSP5/100,y=VP4$ratioOther/100)
LM_VP4<-LMbyProtein(data,"VP4","NSP5")

## Main results of the linear regression model
coefficients<-rbind(LM_VP6$coef,LM_VP4$coef)
colnames(coefficients)<-c("Estimate","Std.Error","t.value","p.value")
coefficients<-data.frame(Protein =c("VP6","VP4"),coefficients)
coefficients$Protein<-factor(coefficients$Protein,levels = levels(coefficients$Protein)[c(2,1)])
coefficients
# --------------------------------------------
# Save the graphs
## LM_VP6
pdf(file = "LM_VP6_NSP5.pdf",width = 4.5,height = 4)
LM_VP6$Res_plot
dev.off()
## LM_VP4
pdf(file = "LM_VP4_NSP5.pdf",width = 4.5,height = 4)
LM_VP4$Res_plot
dev.off()
