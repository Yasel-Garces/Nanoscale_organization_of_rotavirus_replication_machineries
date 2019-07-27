#' This script analyses the distribution of the viral proteins into 
#' the viroplasm taking NSP4 as the reference protein. For details about 
#' the full research consult the paper
#''Garcés et al. Nanoscale organization of rotavirus replication machineries. 
#'eLife 2019;8:e42906. https://elifesciences.org/articles/42906, doi: 10.7554/eLife.42906.
#'
#'  @param This script uses the data collected for each protein combination that
#' are stored in the csv files:
#' 1- NSP4vsVP6.csv
#' This files have the following structure:
#' Column 1: Distance between the distribution of both proteins.
#' Column 2: Radius of the circumference that adjust the central protein.
#' Column 3: Radius of the circumference that adjust the other protein.
#' @return A set of graphics and statistics. See below for details.
#' **Note** This script is practically the same than "AnalysisNSP2.R", 
#' but as a difference it include the regression model analysis.
#' @author Yasel Garces (88yasel@gmail.com)  
#===================================================================================
#===================================================================================
## FUNCTIONS 
LMbyProtein <- function(data,Protein){
  # This function adjusts a least-square linear model to a set of points (data) and
  # return a residual plot, a table with the coeficients of the linear regression and 
  # the residuals values.
  # INPUT: 
  #     data: data frame with columns (x,y).
  #     Protein: Name of the viral protein that is going to be adjusted as dependet 
  #              variable
  # OUTPUT: A list with the next elements:
  #     res_plot: Residual of the linear regression model.
  #     coef: Coefficients of the linear regression model.
  #     residual: residual valuess.
  # AUTHOR: Yasel Garces (88yasel@gmail.com)
  #-------------------------------------------------------------
  # Adjust the linear regression model to the data
  fit<-lm(y ~ x - 1,data)
  # Summary the resut
  t<-summary(fit)
  
  data$predicted <- predict(fit)   # Save the predicted values
  data$residuals <- residuals(fit) # Save the residual values
  
  # Residuals graphic
  # Use the residuals values to make an aesthetic adjustment 
  # (e.g. red colour when residual in very high) to highlight points 
  # which are poorly predicted by the model.
  Residuals<-abs(data$residuals)
  res_plot<-ggplot(data, aes(x = x, y = y)) +
    geom_smooth(method = "lm",formula = y~x-1, se = FALSE, color = "lightgrey") +
    geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
    geom_point(aes(color = Residuals)) + # size also mapped
    geom_point(aes(y = predicted),shape = 1)+
    scale_color_continuous(low = "black", high = "red") +
    guides(color = guide_colorbar())+
    scale_x_continuous(breaks = seq(0,max(data$x)+0.1,by = 0.1))+
    scale_y_continuous(breaks = seq(0,max(data$predicted)+0.1,by = 0.1))+
    geom_label(x = min(data$x), hjust =0, y = max(data$predicted)-0.05,
               label = paste("RSE=", abbreviate(as.character(t$sigma),5),"\n",
                             "R-squared=", abbreviate(as.character(t$r.squared),5)))+
    ylab(Protein)+xlab("NSP4")+scale_fill_continuous(guide = guide_legend()) +
    theme(legend.key.width = unit(2.6, 'lines'), legend.position="bottom",
          axis.text.y =element_text(size=15),
          axis.text.x =element_text(size=15))
  
  # Coefficients of the linear regression model
  coef<-t$coefficients
  # Residual values
  residual<-t$residuals
  # List to return
  list(Res_plot=res_plot,coef=coef,res=residual)
} 
#===================================================================================
#===================================================================================
# Load Libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
theme_set(theme_cowplot()) # Change theme
# Set work directory
setwd('/home/yasel/TRABAJO/IBt/Viroplasms/GitHub Codes Paper  (No Mover)/Nanoscale_organization_of_rotavirus_replication_machineries/R Codes/ResultsCSV/NSP4/')
# Load data
viroData<-read.csv('NSP4vsVP6.csv')

## Data manipulation
# Add a new variable with the name of the NSP4 accompanying protein
# and convert to a factor variable
viroData<-mutate(viroData, Protein='VP6')
viroData$Protein<-as.factor(viroData$Protein)
# Convert from pixels to microns (this is based on our experimetal design, for 
# other experiments you need to take care about how to do this conversion).
viroData$Distance=viroData$Distance/100
viroData$ratioNSP4=viroData$ratioNSP4/100
viroData$ratioOther=viroData$ratioOther/100
#===================================================================================
# Exploratory analysis of the results obtained by the algorithm VPs-DLSFC.
# VP6 spatial distribution taking NSP4 as reference protein.
# Create a convinient data frame that allow to graphic with a boxplot 
# the radii's distribution of all the combinations off proteins.
allProteins<-mutate(viroData,Comparison="Other Protein")
forNSP4<-select(allProteins,one_of(c("ratioNSP4","Protein")))
forNSP4<-mutate(forNSP4,Comparison="NSP4")
allProteins<-data.frame(Comparison=c(forNSP4$Comparison,allProteins$Comparison),
                        Protein=c(as.character(forNSP4$Protein),as.character(allProteins$Protein)),
                        Value=c(forNSP4$ratioNSP4,allProteins$ratioOther))

# Boxplot for the radii of the fitting circumference of NSP4 and VP6.
p<-ggplot(allProteins, aes(x = Protein, y = Value,fill=Comparison)) +
  geom_boxplot(notch=TRUE)+ylab(expression(
    paste("Radius of the adjusted circumference"," ", (paste(mu,m)))))+
  theme(legend.position = "bottom")+
  scale_y_continuous(breaks = seq(from = 0.3, to = 1, by =0.1),limits = c(0.3,1))
p
# Two-sample Mann-Whitney hypothesis test, considering as variables the radii of
# NSP4 in contrast with the radii of VP6.
W_test<-wilcox.test(viroData$ratioOther,viroData$ratioNSP4,conf.int = TRUE)
# Save the plot
pdf(file = "Radius_With_NSP4.pdf",width = 6.8,
    height = 5)
p
dev.off()
#===================================================================================
# Boxplot for the distance between NSP4 and VP6.
p<-ggplot(viroData, aes(x = Protein, y = Distance,fill=Protein)) +
  geom_boxplot(notch=TRUE,show.legend=FALSE)+
  ylab(expression(paste("Distance to NSP4"," ",(paste(mu,m)))))+
  xlab("Protein")+scale_y_continuous(breaks = seq(from = 0, to = .3, by =0.05),limits = c(0,.3))
p
# Save figure
pdf(file = "Distance_With_NSP4.pdf",width = 6.8,
    height = 5)
p
dev.off()
#===================================================================================
# Linear regression fitting taking the radius of NSP4 as independent variable
# and the radius of VP6 as dependent variable. The gray shadow represents the 
# confidence interval at a level of 95%.
p<-ggplot(viroData, aes(x = ratioNSP4, y = ratioOther,color=Protein)) + geom_point(show.legend=FALSE) + 
  facet_grid(~ Protein)+ xlab("Ratio NSP2")+
  geom_smooth(method='lm',formula=y~x-1,show.legend=FALSE)+
  scale_y_continuous(name=expression(paste("Protein Radius"," ", 
                                           (paste(mu,m)))), 
                     breaks=seq(0, 1,0.1),limits = c(0.3,1))+
  scale_x_continuous(name=expression(paste("Radius NSP4"," ", 
                                           (paste(mu,m)))))+
  theme(legend.position = "none",axis.text.x = element_text(vjust = 0.5))
p
# Save figure
pdf(file = "LR_With_NSP4.pdf",width = 6.8,
    height = 5)
p
dev.off()
#===================================================================================
# Residuals errors for each linear regression model. The gray line represent the
# regression model, the points over the line are the predicted values,
# and the dots filled with a color gradient are the real values. The errors between 
# the predicted and real values, are represented as a gradient of colors as follows 
# (from lowest/coldest to highest/warmest).

# Linear regression analysis with NSP4 as independent variable and VP6 as dependent.
data<-data.frame(x=viroData$ratioNSP4,y=viroData$ratioOther)
LM_VP6<-LMbyProtein(data,"VP6")
# Save the graph
pdf(file = "Residuals_NSP4.pdf",width = 5.5,height = 5)
LM_VP6$Res_plot
dev.off()
#################################################################t