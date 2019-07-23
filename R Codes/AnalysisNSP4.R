# This is a very simple aproach to compute a very very tiny analysis for the viroplasmas.
#####################################################################
##-----------------------------------------------------##
LMbyProtein<- function(data,Protein){
  fit<-lm(y ~ x - 1,data)
  t<-summary(fit)
  
  data$predicted <- predict(fit)   # Save the predicted values
  data$residuals <- residuals(fit) # Save the residual values
  
  # Residuals graphic
  # Use the residuals to make an aesthetic adjustment 
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
  
  # Coeff
  coef<-t$coefficients
  # Residuals
  residual<-t$residuals
  
  list(Res_plot=res_plot,coef=coef,res=residual)
} 
#####################################################################

# Libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)

# Load data
setwd('/home/yasel/Dropbox/Paper Viroplasmas/Programs/R/ResultsCSV/NSP4/')

VP6<-read.csv('NSP4vsVP6.csv')

## Data manipulation
# Create a factor variable with the protein name
VP6<-mutate(VP6, Protein='VP6')

# Merge the data in a same data frame
viroData<-VP6
viroData$Protein<-as.factor(viroData$Protein)
viroData$Distance=viroData$Distance/100
viroData$ratioNSP4=viroData$ratioNSP4/100
viroData$ratioOther=viroData$ratioOther/100
##################################################################################
##################################################################################
allProteins<-mutate(viroData,Comparation="Other Protein")
forNSP4<-select(allProteins,one_of(c("ratioNSP4","Protein")))
forNSP4<-mutate(forNSP4,Comparation="NSP4")
allProteins<-data.frame(Comparation=c(forNSP4$Comparation,allProteins$Comparation),
                        Protein=c(as.character(forNSP4$Protein),as.character(allProteins$Protein)),
                        Value=c(forNSP4$ratioNSP4,allProteins$ratioOther))


p<-ggplot(allProteins, aes(x = Protein, y = Value,fill=Comparation)) +
  geom_boxplot(notch=TRUE)+ylab(expression(
    paste("Radius of the adjusted circumference"," ", (paste(mu,m)))))+
  theme(legend.position = "bottom")+
  scale_y_continuous(breaks = seq(from = 0.3, to = 1, by =0.1),limits = c(0.3,1))
p

W_test<-wilcox.test(VP6$ratioOther/100,VP6$ratioNSP4/100,conf.int = TRUE)
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Radius_With_NSP4.pdf",width = 6.8,
    height = 5)
p
dev.off()
# ==============================================================

# ==============================================================
# Boxplot (distance to NSP4) ===================================
p<-ggplot(viroData, aes(x = Protein, y = Distance,fill=Protein)) +
  geom_boxplot(notch=TRUE,show.legend=FALSE)+
  ylab(expression(paste("Distance to NSP4"," ",(paste(mu,m)))))+
  xlab("Protein")+scale_y_continuous(breaks = seq(from = 0, to = .3, by =0.05),limits = c(0,.3))
p
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Distance_With_NSP4.pdf",width = 6.8,
    height = 5)
p
dev.off()
#----------------------
# Lineal dependence between NSP4 and the others proteins. That is for the ratio of the circle
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
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/LR_With_NSP4.pdf",width = 6.8,
    height = 5)
p
dev.off()

# NSP2 vs VP6
data<-data.frame(x=VP6$ratioNSP4/100,y=VP6$ratioOther/100)
LM_VP6<-LMbyProtein(data,"VP6")

pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Residuals_NSP4.pdf",width = 5.5,height = 5)
LM_VP6$Res_plot
dev.off()
#################################################################t