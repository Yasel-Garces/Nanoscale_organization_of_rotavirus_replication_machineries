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
    ylab(Protein)+xlab("NSP5")+scale_fill_continuous(guide = guide_legend()) +
    theme(legend.key.width = unit(2.6, 'lines'), legend.position="bottom",
          axis.text.y =element_text(size=15),
          axis.text.x =element_text(size=15))
  
  # Coeff
  coef<-t$coefficients
  # Residuals
  residual<-t$residuals
  
  list(Res_plot=res_plot,coef=coef,res=residual)
} 

# Remove outliers in data frame
remove_outliers <- function(data, na.rm = TRUE, ...) {
  x<-data$ratioOther
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  data$ratioOther<-y
  na.omit(data,cols="ratioOther")
}
#####################################################################

# Libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(ggsignif)

# Load data
setwd('/home/yasel/Dropbox/Paper Viroplasmas/Programs/R/ResultsCSV/NSP5/')

VP4<-read.csv('NSP5vsVP4.csv')
VP6<-read.csv('NSP5vsVP6.csv')

## Data manipulation
# Create a factor variable with the protein name
VP4<-mutate(VP4, Protein='VP4')
VP6<-mutate(VP6, Protein='VP6')

# Merge the data in a same data frame
viroData<-rbind(VP6,VP4)
viroData$Protein<-as.factor(viroData$Protein)
viroData$Distance=viroData$Distance/100
viroData$ratioNSP5=viroData$ratioNSP5/100
viroData$ratioOther=viroData$ratioOther/100
##################################################################################
####################################################
# Histogram
p<-ggplot(viroData,aes(Protein))+geom_histogram(stat = "count", fill=I("red"), 
                                                col=I("red"), 
                                                alpha=I(.2))+
  theme(legend.position = "none",axis.text.x = element_text(vjust = 0.5))+
  geom_text(stat = "count", aes(label = ..count.., y = ..count..+3))
pdf(file = "HistogramNSP5.pdf",width = 2,height = 2)
p
dev.off()
#---------------------------------------------------
allProteins<-mutate(viroData,Comparation="Other Protein")
forNSP5<-select(allProteins,one_of(c("ratioNSP5","Protein")))
forNSP5<-mutate(forNSP5,Comparation="NSP5")
allProteins<-data.frame(Comparation=c(forNSP5$Comparation,allProteins$Comparation),
                        Protein=c(as.character(forNSP5$Protein),as.character(allProteins$Protein)),
                        Value=c(forNSP5$ratioNSP5,allProteins$ratioOther))

allProteins$Protein<-factor(allProteins$Protein,levels = levels(allProteins$Protein)[c(2,1,3)])
p<-ggplot(allProteins, aes(x = Protein, y = Value,fill=Comparation)) +
  geom_boxplot(notch=TRUE)+ylab(expression(
    paste("Radius of the adjusted circumference"," ", (paste(mu,m)))))+
  theme(legend.position = "bottom")+
  scale_y_continuous(breaks = seq(from = 0.2, to = 1.5, by =0.1),limits = c(0.2,1.2))
p
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Radius_With_NSP5.pdf",width = 6,
    height = 5)
p
dev.off()

NSP5compare<-subset(allProteins,Comparation=="NSP5")
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Radius_Only_NSP5.pdf",width = 3,
    height = 3)
ggplot(NSP5compare,aes(x = Protein, y = Value))+geom_boxplot(notch = TRUE)+
  geom_signif(test="wilcox.test", comparisons = combn(levels(NSP5compare$Protein),2, simplify = F)[-4],
              step_increase = 0.2,map_signif_level=TRUE)+
  theme(legend.position = "none")
dev.off()
# ==============================================================

# ==============================================================
# Boxplot (distance to NSP5) ===================================
viroData$Protein<-factor(viroData$Protein,levels = levels(viroData$Protein)[c(2,1,3)])
p<-ggplot(viroData, aes(x = Protein, y = Distance,fill=Protein)) +
  geom_boxplot(notch=TRUE,show.legend=FALSE)+
  ylab(expression(paste("Distance to NSP5"," ",(paste(mu,m)))))+
  xlab("Protein")+scale_y_continuous(breaks = seq(from = -.5, to = .5, by =0.1),limits = c(0,.5))+
  geom_signif(test="wilcox.test", comparisons = combn(levels(viroData$Protein),2, simplify = F)[-5],
              step_increase = 0.05,map_signif_level=TRUE)
p
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Distance_With_NSP5.pdf",width = 5,
    height = 4)
p
dev.off()
#----------------------
# Lineal dependence between NSP2 and the others proteins. That is for the ratio of the circle
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
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/LR_With_NSP5.pdf",width = 5,
    height = 4)
p
dev.off()
#################################################################
# Inferential analysis (test between the radius of the proteins.)
# Resume the data
# Test
KW_VP6<-wilcox.test(VP6$ratioOther/100,VP6$ratioNSP5/100,conf.int = TRUE)
KW_VP4<-wilcox.test(VP4$ratioOther/100,VP4$ratioNSP5/100,conf.int = TRUE)

# p-value
p_value=c(KW_VP6$p.value,KW_VP4$p.value)
# W
W=c(KW_VP6$statistic, KW_VP4$statistic)
# DF
Diff_Location=c(KW_VP6$estimate, KW_VP4$estimate)
# CI
lessCoef<-c(KW_VP6$conf.int[1], KW_VP4$conf.int[1])
greatCoef<-c(KW_VP6$conf.int[2], KW_VP4$conf.int[2])
# Create data frame
wilcoxTest<-data.frame(Protein = c("VP6","VP4"),W,Diff_Location,lessCoef,greatCoef,p_value)
wilcoxTest$Protein<-factor(wilcoxTest$Protein,levels = levels(wilcoxTest$Protein)[c(2,1)])
# Plot the results of Wilcox test.
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/DifferenceLocationNSP5.pdf",width = 4,height = 3)
ggplot(wilcoxTest,aes(x = Protein,y = Diff_Location,color=Protein))+
  geom_errorbar(ymin = lessCoef, ymax = greatCoef,width = 0.3 , size=1)+
  geom_point(size=3)+  
  scale_y_continuous(breaks = seq(0.05,0.3,0.05), labels = seq(0.05,0.3,0.05),limits = c(0.05,0.3))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Difference Location")+geom_text(aes(y = Diff_Location,x=c(0.7,1.7),
                                            label = format(p_value, scientific = TRUE)))
dev.off()
#################################################################
# Study of the linear regression.
# NSP2 vs VP6
data<-data.frame(x=VP6$ratioNSP5/100,y=VP6$ratioOther/100)
LM_VP6<-LMbyProtein(data,"VP6")

# NSP5 vs VP4
data<-data.frame(x=VP4$ratioNSP5/100,y=VP4$ratioOther/100)
LM_VP4<-LMbyProtein(data,"VP4")
##-----------------------------------------------------##
## Coefficients of the linear regression
coefficients<-rbind(LM_VP6$coef,LM_VP4$coef)
colnames(coefficients)<-c("Estimate","Std.Error","t.value","p.value")
coefficients<-data.frame(Protein =c("VP6","VP4"),coefficients)
coefficients$Protein<-factor(coefficients$Protein,levels = levels(coefficients$Protein)[c(2,1)])

pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/SlopeEstimationCINSP5.pdf",width = 4,height = 3)
ggplot(coefficients,aes(x = Protein,y = Estimate,color=Protein))+
  geom_errorbar(ymin = coefficients$Estimate - coefficients$Std.Error, 
                ymax = coefficients$Estimate + coefficients$Std.Error,
                width = 0.4 , size=0.4)+geom_point(size=1.5)+  
  scale_y_continuous(breaks = seq(1.25,1.7,0.1), labels = seq(1.25,1.7,0.1),limits = c(1.25,1.7))+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+
  ylab("Slope Estimation & CI")+geom_text(aes(y = Estimate,x=c(0.7,1.7),
                                              label = format(Estimate, scientific = FALSE)))
dev.off()
# --------------------------------------------
# Save the plots
## LM_VP6
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/LM_VP6_NSP5.pdf",width = 4.5,height = 4)
LM_VP6$Res_plot
dev.off()
## LM_VP4
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/LM_VP4_NSP5.pdf",width = 4.5,height = 4)
LM_VP4$Res_plot
dev.off()
