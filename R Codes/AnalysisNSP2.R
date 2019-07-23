# Main program for the statistic analysis of the viroplasm's rotavirus.

# Load Libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(ggsignif)
library(RColorBrewer)
library(gplots)
library(gdata)
# --------------

# Change work directory and Load data
setwd('/home/yasel/Dropbox/Paper Viroplasmas/Programs/R/ResultsCSV/NSP2/')

VP4<-read.csv('NSP2rojo-VP4verde.csv')
VP6<-read.csv('NSP2rojo-VP6verde.csv')
VP7_Mon<-read.csv('NSP2rojo-VP760verde.csv')
VP7_Tri<-read.csv('NSP2rojo-VP7159verde.csv')
NSP4<-read.csv('NSP4rojo-NSP2verde.csv')
NSP5<-read.csv('NSP5rojo-NSP2verde.csv')
VP1<-read.csv('NSP2rojo-VP1verde.csv')
VP2<-read.csv('NSP2rojo-VP2verde.csv')
##-----------------------------------------------
## Data manipulation
# Create a factor variable with the protein name
#dsRNA<-mutate(dsRNA, Protein='dsRNA')
#PDI<-mutate(PDI, Protein='PDI')
VP4<-mutate(VP4, Protein='VP4')
VP6<-mutate(VP6, Protein='VP6')
VP7_Mon<-mutate(VP7_Mon, Protein='VP7_Mon')
VP7_Tri<-mutate(VP7_Tri, Protein='VP7_Tri')
NSP4<-mutate(NSP4, Protein='NSP4')
NSP5<-mutate(NSP5, Protein='NSP5')
VP1<-mutate(VP1, Protein='VP1')
VP2<-mutate(VP2, Protein='VP2')

###################################################
# Because exist differents numbers of samples for 
# each protein combinations, and we want to compare
# the spatial distribution of differents proteins, 
# we select properly 40 samples of each combination. 
# The selections criterion was that in the final set
# there not exist statistical differences between 
# the protein NSP2 in all combinations. 

# Selection of the samples
NSP5<-tail(arrange(NSP5,desc(NSP5$ratioNSP2)),35)
NSP4<-tail(arrange(NSP4,desc(NSP4$ratioNSP2)),40)
# VP4<-tail(arrange(VP4,desc(VP4$ratioNSP2)),32)
# VP7159<-tail(arrange(VP7159,desc(ratioNSP2)),32)
VP7_Mon<-tail(arrange(VP7_Mon,desc(VP7_Mon$ratioNSP2)),60)
VP1<-head(arrange(VP1,desc(VP1$ratioNSP2)),40)
VP2<-head(arrange(VP2,desc(VP2$ratioNSP2)),40)
VP6<-head(arrange(VP6,desc(VP6$ratioNSP2)),60)

##-----------------------------------------------
# Merge the data in a same data frame
viroData<-rbind(VP6,NSP4,NSP5,VP4,VP7_Tri,VP7_Mon,VP1,VP2)
viroData$Protein<-as.factor(viroData$Protein)
viroData$Distance=viroData$Distance/100
viroData$ratioNSP2=viroData$ratioNSP2/100
viroData$ratioOther=viroData$ratioOther/100
## END DATA PREPARATION ##

##################################################################################
####################################################
allProteins<-mutate(viroData,Comparation="Other Protein")
forNSP2<-select(allProteins,one_of(c("ratioNSP2","Protein")))
forNSP2<-mutate(forNSP2,Comparation="NSP2")
allProteins<-data.frame(Comparation=c(forNSP2$Comparation,allProteins$Comparation),
                        Protein=c(as.character(forNSP2$Protein),as.character(allProteins$Protein)),
                        Value=c(forNSP2$ratioNSP2,allProteins$ratioOther))

allProteins$Protein<-factor(allProteins$Protein,levels = levels(allProteins$Protein)[c(2,1,3,4,6,5,7,8)])

p<-ggplot(allProteins, aes(x = Protein, y = Value,fill=Comparation)) +
  geom_boxplot(notch=TRUE)+ylab(expression(
    paste("Circumference's Radius"," ", (paste(mu,m)))))+
  theme(legend.position = "bottom")+
  scale_y_continuous(breaks = seq(from = 0.2, to = 1.5, by =0.1),limits = c(0.2,1.2))
p
# Save
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Boxplot_NSP2vsOthers.pdf",width = 6.8,height = 5)
p+theme(text = element_text(size=20),axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()

# Just NSP2 in all combinations
compareNSP2<-subset(allProteins,Comparation=="NSP2")
p<-ggplot(compareNSP2,aes(x = Protein,y = Value,fill="red"))+geom_boxplot(notch = TRUE)+
  geom_signif(test="wilcox.test", comparisons = combn(levels(compareNSP2$Protein),2, simplify = F)[-4],
              step_increase = 0.2,map_signif_level=TRUE)+theme(legend.position = "none")+ylab("Value")+xlab("NSP2")+
  scale_y_continuous(breaks = seq(0,1,0.15),labels = seq(0,1,0.15))
# Save
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/WilcoxTest_NSP2_all.pdf",width = 8,height = 8)
p+theme(text = element_text(size=20),axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()
# ==============================================================
# ==============================================================
# Boxplot (distance to NSP2) ===================================
VP1andVP2<-subset(viroData,Protein %in% c("VP1","VP2"))
viroData<-subset(viroData,Protein!="VP1" & Protein!="VP2")
viroData$Protein<-drop.levels(viroData$Protein)
viroData$Protein<-factor(viroData$Protein,levels = levels(viroData$Protein)[c(2,1,4,3,5,6)])

# Just in viroData ----------------------------------------------
p<-ggplot(viroData, aes(x = Protein, y = Distance,fill=Protein)) +
  geom_boxplot(notch=TRUE,show.legend=FALSE)+
  ylab(expression(paste("Distance to NSP2"," ",(paste(mu,m)))))+
  xlab("Protein")+scale_y_continuous(breaks = seq(from = -.5, to = .8, by =0.1))
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Boxplot_Distance_NSP2.pdf",width = 4.8,height = 3)
p+theme(text = element_text(size=20),axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18))
dev.off()

pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Boxplot_Wilcox_Proteins.pdf",width = 6.8,height = 7)
p+geom_signif(test="wilcox.test", comparisons = combn(levels(viroData$Protein),2, simplify = F)[-5],
              step_increase = 0.01,map_signif_level=TRUE)+
  theme(text = element_text(size=20),axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()
# Just in VP1andVP2 ----------------------------------------------
VP1andVP2$Protein<-drop.levels(VP1andVP2$Protein)
VP1andVP2$Protein<-factor(VP1andVP2$Protein,levels = levels(VP1andVP2$Protein)[c(2,1)])
p<-ggplot(VP1andVP2, aes(x = Protein, y = Distance,fill=Protein)) +
  geom_boxplot(notch=TRUE,show.legend=FALSE)+
  ylab(expression(paste("Distance to NSP2"," ",(paste(mu,m)))))+
  xlab("Protein")+scale_y_continuous(breaks = seq(from = -.5, to = .6, by =0.05))
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/Boxplot_Distance_VP1_VP2.pdf",width = 4.8,height = 3)
p+geom_signif(test="wilcox.test", comparisons = combn(levels(VP1andVP2$Protein),2, simplify = F)[-5],
               step_increase = 0.05,map_signif_level=TRUE)+
  theme(text = element_text(size=15),axis.text.x = element_text(size = 15),
                                                                  axis.text.y = element_text(size = 15))
dev.off()

# Cluster ======================================0
# Merge data
ClusterViro<- summarise(group_by(viroData,Protein), Distance = mean(Distance), SdDistance=sd(ratioNSP2), 
                        ratioNSP2=mean(ratioNSP2),Ratio=mean(ratioOther))
colnames(ClusterViro)<-c("Protein","Mean(Dist)","Sd(Dist)","Mean(NSP2)","Mean(Other)")

NewData<-as.data.frame(ClusterViro[,-1]) 
row.names(NewData)<-as.character(ClusterViro$Protein)

hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")

d <- distfunc(NewData)
fit <- hclustfunc(d)
plot(fit)
groups <- cutree(fit, k=5) 

cols <- brewer.pal(max(groups), "Set1")
# Save data
pdf(file = "/home/yasel/Dropbox/Paper Viroplasmas/HeatMapNSP2.pdf",width = 5,height = 5)
heatmap.2(as.matrix(NewData),dendrogram="row",trace="none", margin=c(8,9), 
          hclust=hclustfunc, distfun=distfunc, RowSideColors=cols[groups],srtCol=45)
dev.off()

