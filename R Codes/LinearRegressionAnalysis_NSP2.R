#' This script contains the linear regression analysis considering the radius 
#' of NSP2 as the independent variable.
#' Garcés et al. Nanoscale organization of rotavirus replication machineries. 
#' eLife 2019;8:e42906. https://elifesciences.org/articles/42906, 
#' doi: 10.7554/eLife.42906.
#' 
#' @param This script uses the data collected for each protein combination that
#' are stored in the csv files:
#' 1- NSP2rojo-VP4verde.csv
#' 2- NSP2rojo-VP6verde.csv
#' 3- NSP2rojo-VP760verde.csv
#' 4- NSP2rojo-VP7159verde.csv
#' 5- NSP2rojo-VP1verde.csv
#' 6- NSP2rojo-VP2verde.csv
#' 7- NSP4rojo-NSP2verde.csv
#' 8- NSP5rojo-NSP2verde.csv
#' All these files have the following structure:
#' Column 1: Distance between the distribution of both proteins.
#' Column 2: Radius of the circumference that adjust the central protein.
#' Column 3: Radius of the circumference that adjust the other protein.
#' @return A set of graphics and statistics. See below for details.
#' @author Yasel Garces (88yasel@gmail.com)
#===================================================================================
#===================================================================================
## FUNCTIONS 
# Load the functions saved in Functions.R
source('/home/yasel/TRABAJO/IBt/Viroplasms/GitHub Codes Paper  (No Mover)/Nanoscale_organization_of_rotavirus_replication_machineries/R Codes/Functions.R')
#===================================================================================
# Load Libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(ggsignif)
theme_set(theme_cowplot()) # Change theme

# Set work directory
setwd('/home/yasel/TRABAJO/IBt/Viroplasms/GitHub Codes Paper  (No Mover)/Nanoscale_organization_of_rotavirus_replication_machineries/R Codes/ResultsCSV/NSP2')
# Load data
VP4<-read.csv('NSP2rojo-VP4verde.csv')
VP6<-read.csv('NSP2rojo-VP6verde.csv')
VP760<-read.csv('NSP2rojo-VP760verde.csv')
VP7159<-read.csv('NSP2rojo-VP7159verde.csv')
NSP4<-read.csv('NSP4rojo-NSP2verde.csv')
NSP5<-read.csv('NSP5rojo-NSP2verde.csv')
VP1<-read.csv('NSP2rojo-VP1verde.csv')
VP2<-read.csv('NSP2rojo-VP2verde.csv')
#===================================================================================
## Data manipulation
# Add a new variable with the name of the NSP2 accompanying protein
VP4<-mutate(VP4, Protein='VP4')
VP6<-mutate(VP6, Protein='VP6')
VP760<-mutate(VP760, Protein='VP760')
VP7159<-mutate(VP7159, Protein='VP7159')
NSP4<-mutate(NSP4, Protein='NSP4')
NSP5<-mutate(NSP5, Protein='NSP5')
VP1<-mutate(VP1, Protein='VP1')
VP2<-mutate(VP2, Protein='VP2')

# Remove the outliers with the objective to obtain a more accurate model.
VP4<-remove_outliers(VP4)
VP6<-remove_outliers(VP6)
VP760<-remove_outliers(VP760)
VP7159<-remove_outliers(VP7159)
NSP4<-remove_outliers(NSP4)
NSP5<-remove_outliers(NSP5)
VP1<-remove_outliers(VP1)
VP2<-remove_outliers(VP2)

# Merge the data and convert the "Protein" column to a factor variable
viroData<-rbind(VP6,NSP4,NSP5,VP4,VP7159,VP760,VP1,VP2)
viroData$Protein<-as.factor(viroData$Protein)
# Convert from pixels to microns (this is based on our experimetal design, for 
# other experiments you need to take care about how to do this conversion).
viroData$Distance=viroData$Distance/100
viroData$ratioNSP2=viroData$ratioNSP2/100
viroData$ratioOther=viroData$ratioOther/100

# Orders the factor levels
viroData$Protein<-factor(viroData$Protein,levels = levels(viroData$Protein)[c(2,1,3,4,6,5,8,7)])

# Plot a histogram with the number of samples per condition.
pdf(file = "Histogram.pdf",width = 6,height = 6)
ggplot(viroData,aes(Protein,fill=Protein))+geom_histogram(stat = "count")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 0.5))+
  geom_text(stat = "count", aes(label = ..count.., y = ..count..+2))
dev.off()
#===================================================================================
# Linear regression fitting taking the radius of NSP2 as independent variable
# and the radius of all the other proteins as dependent variable. The gray shadow represents the 
# confidence interval at a level of 95%.
p<-ggplot(viroData, aes(x = ratioNSP2, y = ratioOther,color=Protein)) + geom_point(show.legend=FALSE) + 
  facet_grid(~ Protein)+ xlab("Ratio NSP2")+
  geom_smooth(method='lm',formula=y~x-1,show.legend=FALSE)+
  scale_y_continuous(name=expression(paste("Protein Radius"," ", 
                                           (paste(mu,m)))), 
                     breaks=seq(0, max(viroData$ratioOther)+0.4,0.1))+
  scale_x_continuous(name=expression(paste("Radius NSP2"," ", 
                                           (paste(mu,m)))))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5))
p
# Save the graphic 
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/LinearRegre_Ratio.pdf",width = 7.5,height = 5.5)
p
dev.off()
#===================================================================================
# Hypothesis test (Test of Equal or Given Proportions)
# The ‘Mann-Whitney’ test was carry out to evaluate the differences between
# the radii of the dependent and independent variables.
KW_NSP5<-wilcox.test(NSP5$ratioOther/100,NSP5$ratioNSP2/100,conf.int = TRUE)
KW_NSP4<-wilcox.test(NSP4$ratioOther/100,NSP4$ratioNSP2/100,conf.int = TRUE)
KW_VP6<-wilcox.test(VP6$ratioOther/100,VP6$ratioNSP2/100,conf.int = TRUE)
KW_VP4<-wilcox.test(VP4$ratioOther/100,VP4$ratioNSP2/100,conf.int = TRUE)
KW_VP760<-wilcox.test(VP760$ratioOther/100,VP760$ratioNSP2/100,conf.int = TRUE)
KW_VP7159<-wilcox.test(VP7159$ratioOther/100,VP7159$ratioNSP2/100,conf.int = TRUE)
KW_VP1<-wilcox.test(VP1$ratioOther/100,VP1$ratioNSP2/100,conf.int = TRUE)
KW_VP2<-wilcox.test(VP2$ratioOther/100,VP2$ratioNSP2/100,conf.int = TRUE)

# p-value of the Mann-Whitney test
p_value<-c(KW_NSP5$p.value, KW_NSP4$p.value, KW_VP1$p.value,KW_VP2$p.value, 
           KW_VP6$p.value, KW_VP4$p.value, KW_VP760$p.value, KW_VP7159$p.value)
# W (the value of the test statistic) 
W<-c(KW_NSP5$statistic, KW_NSP4$statistic, KW_VP1$statistic,KW_VP2$statistic,
     KW_VP6$statistic, KW_VP4$statistic, KW_VP760$statistic, KW_VP7159$statistic)
# Confidence interval
lessCoef<-c(KW_NSP5$conf.int[1], KW_NSP4$conf.int[1], KW_VP1$conf.int[1],KW_VP2$conf.int[1],
            KW_VP6$conf.int[1], KW_VP4$conf.int[1], KW_VP760$conf.int[1], KW_VP7159$conf.int[1])
greatCoef<-c(KW_NSP5$conf.int[2], KW_NSP4$conf.int[2], KW_VP1$conf.int[2],KW_VP2$conf.int[2],
             KW_VP6$conf.int[2], KW_VP4$conf.int[2], KW_VP760$conf.int[2], KW_VP7159$conf.int[2])
# difference in location
Diff_Location<-c(KW_NSP5$estimate, KW_NSP4$estimate, KW_VP1$estimate,KW_VP2$estimate,
                KW_VP6$estimate, KW_VP4$estimate, KW_VP760$estimate, KW_VP7159$estimate)
# Create a data frame with all the previous information
wilcoxTest<-data.frame(Protein = c("NSP5","NSP4","VP1","VP2","VP6","VP4","VP760",
                                   "VP7159"),W,Diff_Location,lessCoef,greatCoef,p_value)
wilcoxTest$Protein<-factor(wilcoxTest$Protein,levels = levels(wilcoxTest$Protein)[c(2,1,3,4,6,5,8,7)])

# Plot the results of the Mann-Whitney test.
pdf(file = "DifferenceLocation.pdf",width = 4.5,height = 5)
ggplot(wilcoxTest,aes(x = Protein,y = Diff_Location,color=Protein))+
  geom_errorbar(ymin = lessCoef, ymax = greatCoef,width = 0.3 , size=1)+
  geom_point(size=3)+  
  scale_y_continuous(breaks = seq(-0.1,0.7,0.05), labels = seq(-0.1,0.7,0.05),limits = c(-0.1,0.7))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Difference Location")+geom_text(aes(y = greatCoef+0.12,
                                            label = format(p_value, scientific = TRUE)),angle=90)
dev.off()
#===================================================================================
# Linear regression analysis for all the protein combinations.
# NSP2 vs NSP5
data<-data.frame(x=NSP5$ratioNSP2/100,y=NSP5$ratioOther/100)
LM_NSP5<-LMbyProtein(data,"NSP5","NSP2")
# NSP2 vs NSP4
data<-data.frame(x=NSP4$ratioNSP2/100,y=NSP4$ratioOther/100)
LM_NSP4<-LMbyProtein(data,"NSP4","NSP2")
# NSP2 vs VP6
data<-data.frame(x=VP6$ratioNSP2/100,y=VP6$ratioOther/100)
LM_VP6<-LMbyProtein(data,"VP6","NSP2")
# NSP2 vs VP4
data<-data.frame(x=VP4$ratioNSP2/100,y=VP4$ratioOther/100)
LM_VP4<-LMbyProtein(data,"VP4","NSP2")
# NSP2 vs VP760
data<-data.frame(x=VP760$ratioNSP2/100,y=VP760$ratioOther/100)
LM_VP760<-LMbyProtein(data,"VP760","NSP2")
# NSP2 vs VP7159
data<-data.frame(x=VP7159$ratioNSP2/100,y=VP7159$ratioOther/100)
LM_VP7159<-LMbyProtein(data,"VP7159","NSP2")
# NSP2 vs VP1
data<-data.frame(x=VP1$ratioNSP2/100,y=VP1$ratioOther/100)
LM_VP1<-LMbyProtein(data,"VP1")
# NSP2 vs VP2
data<-data.frame(x=VP2$ratioNSP2/100,y=VP2$ratioOther/100)
LM_VP2<-LMbyProtein(data,"VP2","NSP2")
#===================================================================================
## Coefficients of the linear regressions
coefficients<-rbind(LM_NSP5$coef,LM_NSP4$coef, LM_VP1$coef,LM_VP2$coef,
                    LM_VP6$coef,LM_VP4$coef,LM_VP760$coef,LM_VP7159$coef)
colnames(coefficients)<-c("Estimate","Std.Error","t.value","p.value")
coefficients<-data.frame(Protein=c("NSP5","NSP4","VP1","VP2",
                                   "VP6","VP4","VP760","VP7159"),coefficients,
                         row.names = c("NSP5","NSP4","VP1","VP2",
                                      "VP6","VP4","VP760","VP7159"))
# Plot the slope of the regression model for each protein.
coefficients$Protein<-factor(coefficients$Protein,levels = levels(coefficients$Protein)[c(2,1,3,4,6,5,8,7)])
# Plot the coefficients of the linear regresions
pdf(file = "SlopeEstimationCI.pdf",width = 4.5,height = 5)
ggplot(coefficients,aes(x = Protein,y = Estimate,color=Protein))+
  geom_errorbar(ymin = coefficients$Estimate - coefficients$Std.Error, 
                ymax = coefficients$Estimate + coefficients$Std.Error,
                width = 0.4 , size=0.4)+geom_point(size=1.5)+  
  scale_y_continuous(breaks = seq(0.45,2,0.1), labels = seq(0.45,2,0.1),limits = c(0.45,2))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Slope Estimation & CI")+geom_text(aes(y = Estimate - coefficients$Std.Error- 0.25,
                                            label = format(Estimate, scientific = FALSE)),angle=90)
dev.off()
# Plot the p.value of the regression models for each protein.
pdf(file = "PValueRM.pdf",width = 4,height = 4)
ggplot(coefficients, aes(y=Protein,x=p.value, label=round(p.value,4),color=Protein)) +
  geom_point(stat='identity', size=8)+
  geom_text(color="white", size=2)+theme(panel.grid.major.y =
                                           element_line(colour="gray", size=1))+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  )+theme(legend.position = "none")+
  xlab("P-Value")+ylab("Proteins")
dev.off()
#===================================================================================
# Save the plots
## LM_NSP5
pdf(file = "LM_NSP5.pdf",width = 4.5,height = 4)
LM_NSP5$Res_plot
dev.off()
## LM_NSP4
pdf(file = "LM_NSP4.pdf",width = 4.5,height = 4)
LM_NSP4$Res_plot
dev.off()
## LM_VP6
pdf(file = "LM_VP6.pdf",width = 4.5,height = 4)
LM_VP6$Res_plot
dev.off()
## LM_VP4
pdf(file = "LM_VP4.pdf",width = 4.5,height = 4)
LM_VP4$Res_plot
dev.off()
## LM_VP760
pdf(file = "LM_VP760.pdf",width = 4.5,height = 4)
LM_VP760$Res_plot
dev.off()
## LM_VP7159
pdf(file = "LM_VP7159.pdf",width = 4.5,height = 4)
LM_VP7159$Res_plot
dev.off()
## LM_VP1
pdf(file = "LM_VP1.pdf",width = 4.5,height = 4)
LM_VP1$Res_plot
dev.off()
## LM_VP2
pdf(file = "LM_VP2.pdf",width = 4.5,height = 4)
LM_VP2$Res_plot
dev.off()
######################################################################
# Model Representation
######################################################################
library(geomnet)
library(ggforce)
# Data frame with all the proteins and the coefficient of the linear regressions.
mydata<-data.frame(
  Protein = c("NSP5","NSP2","NSP4","VP2","VP1","VP6","VP4","VP760","VP7159"),
  Radius = c(0.6,
             1, 
             1.1,
             1.6,
             2.1,
             2.2,
             3.2,
             3.7,
             3.8
  ))
mydata<-mutate(mydata,Start=rep(0,9),End=rep(2*pi,9))
mydata$Protein<-factor(mydata$Protein,levels = 
                         levels(mydata$Protein)[c(3,1,2,4,5,7,6,9,8)])
mydata$Radius<-as.numeric(mydata$Radius)

p<-ggplot()+geom_point()+
  geom_arc(aes(x0=0, y0=0, r=Radius, start=Start,end=End, 
               color=Protein), data=mydata,size = 5)+
  geom_text(data=mydata, mapping=aes(x=runif(9,-1,1), 
                                     y=Radius-runif(9,0,0.1), 
                                     label=Protein), size=4)+
  scale_x_discrete()+scale_y_discrete()+ylab("")+xlab("")
# Save the plot
pdf(file = "Scheme.pdf",width = 6,height = 5)
p
dev.off()
