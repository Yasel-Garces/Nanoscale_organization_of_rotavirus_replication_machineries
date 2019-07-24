#' This script analyses the distribution of the viral proteins into 
#' the viroplasm taking NSP2 as the reference protein. For details about 
#' the full research consult the paper
#' "Nanoscale organization of rotavirus replication machineries", eLife.
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
#' 
#' @author Yasel Garces (88yasel@gmail.com)
#===================================================================================
# Load Libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(plotly)
library(ggsignif)
library(RColorBrewer)
library(gplots)
library(gdata)
theme_set(theme_cowplot()) # Change theme
# ----------------------------------------------------------------------------------
# Change work directory and Load data
setwd('/home/yasel/TRABAJO/IBt/Viroplasms/GitHub Codes Paper  (No Mover)/Nanoscale_organization_of_rotavirus_replication_machineries/R Codes/ResultsCSV/NSP2')
# Load data
VP4<-read.csv('NSP2rojo-VP4verde.csv')
VP6<-read.csv('NSP2rojo-VP6verde.csv')
VP7_Mon<-read.csv('NSP2rojo-VP760verde.csv')
VP7_Tri<-read.csv('NSP2rojo-VP7159verde.csv')
NSP4<-read.csv('NSP4rojo-NSP2verde.csv')
NSP5<-read.csv('NSP5rojo-NSP2verde.csv')
VP1<-read.csv('NSP2rojo-VP1verde.csv')
VP2<-read.csv('NSP2rojo-VP2verde.csv')
# ----------------------------------------------------------------------------------
## Data manipulation
# Create a factor variable with the protein's name
VP4<-mutate(VP4, Protein='VP4')
VP6<-mutate(VP6, Protein='VP6')
VP7_Mon<-mutate(VP7_Mon, Protein='VP7_Mon')
VP7_Tri<-mutate(VP7_Tri, Protein='VP7_Tri')
NSP4<-mutate(NSP4, Protein='NSP4')
NSP5<-mutate(NSP5, Protein='NSP5')
VP1<-mutate(VP1, Protein='VP1')
VP2<-mutate(VP2, Protein='VP2')
# ----------------------------------------------------------------------------------
# Because exist differents numbers of samples for each protein combinations, 
# and we want to compare the spatial distribution of differents proteins, 
# we select properly 40 samples of each combination. The selections criterion 
# was that in the final set there not exist statistical differences between 
# the protein NSP2 in all combinations. 
# **Note** The parameters for this selection can change for other images dataset
# Selection of the samples
NSP5<-tail(arrange(NSP5,desc(NSP5$ratioNSP2)),35)
NSP4<-tail(arrange(NSP4,desc(NSP4$ratioNSP2)),40)
VP7_Mon<-tail(arrange(VP7_Mon,desc(VP7_Mon$ratioNSP2)),60)
VP1<-head(arrange(VP1,desc(VP1$ratioNSP2)),40)
VP2<-head(arrange(VP2,desc(VP2$ratioNSP2)),40)
VP6<-head(arrange(VP6,desc(VP6$ratioNSP2)),60)
# ----------------------------------------------------------------------------------
# Merge the data in a same data frame
viroData<-rbind(VP6,NSP4,NSP5,VP4,VP7_Tri,VP7_Mon,VP1,VP2)
# Convert the protein's name to a factor variable
viroData$Protein<-as.factor(viroData$Protein)
# Convert from pixels to microns (this is based on our experimetal design, for 
# other experiments you need to take care about how to do this conversion).
viroData$Distance=viroData$Distance/100
viroData$ratioNSP2=viroData$ratioNSP2/100
viroData$ratioOther=viroData$ratioOther/100
#===================================================================================
# Exploratory analysis of the results obtained by the algorithm VPs-DLSFC
# Create a convinient data frame that allow to boxplot the radii distribution of 
# all the combinations off proteins.
allProteins<-mutate(viroData,Comparison="Other Protein")
forNSP2<-select(allProteins,one_of(c("ratioNSP2","Protein")))
forNSP2<-mutate(forNSP2,Comparison="NSP2")
allProteins<-data.frame(Comparison=c(forNSP2$Comparison,allProteins$Comparison),
                        Protein=c(as.character(forNSP2$Protein),as.character(allProteins$Protein)),
                        Value=c(forNSP2$ratioNSP2,allProteins$ratioOther))
# Order the factor levels
allProteins$Protein<-factor(allProteins$Protein,
                            levels = levels(allProteins$Protein)[c(2,1,3,4,6,5,7,8)])
# Boxplot for the radius of the fitting circumferences. In each experimental condition 
# we plot two boxes, the red box is for the radius of NSP2 (reference protein), and 
# the blue box represents the radius of the accompanying VP components (names                                                                                                  in x-axis).
p<-ggplot(allProteins, aes(x = Protein, y = Value,fill=Comparison)) +
  geom_boxplot(notch=TRUE)+ylab(expression(
    paste("Circumference's Radius"," ", (paste(mu,m)))))+
  theme(legend.position = "bottom")+
  scale_y_continuous(breaks = seq(from = 0.2, to = 1.5, by =0.1),limits = c(0.2,1.2))
p
# Save the figure
pdf(file = "Boxplot_NSP2vsOthers.pdf",width = 6.8,height = 5)
p+theme(text = element_text(size=20),axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()
# ----------------------------------------------------------------------------------
# Graphical representation of the Mann-Whitney hypothesis test in two-by-two comparison
# between the radii distributions of NSP2 in each experiment (Mann and Whitney, 1947). To avoid
# confusions, for each box, under the NSP2 name in x-axis, we included the name of the accompanying
# protein that specifies in which experiment we obtain that distribution of NSP2. Each combination is
# linked by a line, and the result of the test is up of the line.

# Filter the information of NSP2
compareNSP2<-subset(allProteins,Comparison=="NSP2")
# Plot
p<-ggplot(compareNSP2,aes(x = Protein,y = Value,fill="red"))+geom_boxplot(notch = TRUE)+
  geom_signif(test="wilcox.test", comparisons = combn(levels(compareNSP2$Protein),2, simplify = F)[-4],
              step_increase = 0.2,map_signif_level=TRUE)+theme(legend.position = "none")+
  ylab("Value")+xlab("NSP2")+
  scale_y_continuous(breaks = seq(0,1,0.15),labels = seq(0,1,0.15))
# Save
pdf(file = "WilcoxTest_NSP2_all.pdf",width = 8,height = 8)
p+theme(text = element_text(size=20),axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()
# ----------------------------------------------------------------------------------
# Boxplot and results of the Mann-Whitney hypothesis test for the distance between 
# each viral element and NSP2. Each combination of the Mann-Whitney test is linked 
# by a line, and the result of the test it is above the line. Note that this test 
# reports significant differences between the distribution of the distance to 
# NSP2 of two different VP components.
# Filter the data to take into account all the proteins except for VP1 and VP2
viroData<-subset(viroData,Protein!="VP1" & Protein!="VP2")
viroData$Protein<-drop.levels(viroData$Protein)
# Order the factor levels
viroData$Protein<-factor(viroData$Protein,levels = levels(viroData$Protein)[c(2,1,4,3,5,6)])
# Plot whitout Mann-Whitney test
p<-ggplot(viroData, aes(x = Protein, y = Distance,fill=Protein)) +
  geom_boxplot(notch=TRUE,show.legend=FALSE)+
  ylab(expression(paste("Distance to NSP2"," ",(paste(mu,m)))))+
  xlab("Protein")+scale_y_continuous(breaks = seq(from = -.5, to = .8, by =0.1))
pdf(file = "Boxplot_Distance_NSP2.pdf",width = 4.8,height = 3)
# Save the figure
p+theme(text = element_text(size=20),axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18))
dev.off()
# Plot with a Mann-Whitney test.
pdf(file = "Boxplot_Wilcox_Proteins.pdf",width = 6.8,height = 7)
p+geom_signif(test="wilcox.test", comparisons = combn(levels(viroData$Protein),2, simplify = F)[-5],
              step_increase = 0.01,map_signif_level=TRUE)+
  theme(text = element_text(size=20),axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))
dev.off()
# ----------------------------------------------------------------------------------
# Distance of VP1 and VP2 to NSP2 and result of the Mann-Whitney test. Because the 
# distributions of NSP2 in combination with VP1 and VP2 are statistically different 
# to the other NSP2 distributions (see the publication), we show these two cases 
# independently.
# Filter the data to only consider the information relative to VP1 and VP2.
VP1andVP2<-subset(viroData,Protein %in% c("VP1","VP2"))
VP1andVP2$Protein<-drop.levels(VP1andVP2$Protein)
VP1andVP2$Protein<-factor(VP1andVP2$Protein,levels = levels(VP1andVP2$Protein)[c(2,1)])
# Plot
p<-ggplot(VP1andVP2, aes(x = Protein, y = Distance,fill=Protein)) +
  geom_boxplot(notch=TRUE,show.legend=FALSE)+
  ylab(expression(paste("Distance to NSP2"," ",(paste(mu,m)))))+
  xlab("Protein")+scale_y_continuous(breaks = seq(from = -.5, to = .6, by =0.05))
pdf(file = "Boxplot_Distance_VP1_VP2.pdf",width = 4.8,height = 3)
p+geom_signif(test="wilcox.test", comparisons = combn(levels(VP1andVP2$Protein),2, simplify = F)[-5],
               step_increase = 0.05,map_signif_level=TRUE)+
  theme(text = element_text(size=15),axis.text.x = element_text(size = 15),
                                                                  axis.text.y = element_text(size = 15))
dev.off()
# ----------------------------------------------------------------------------------
# Hierarchically clustered heatmap for the standard deviation of the distance to 
# NSP2, the mean distance to NSP2, the mean radius of NSP2, and the mean radius 
# of the accompanying protein layers, NSP5, NSP4, VP6, VP4, and VP7. 
# Summarise the data considering the factor variable "protein"
ClusterViro<- summarise(group_by(viroData,Protein), Distance = mean(Distance), SdDistance=sd(ratioNSP2), 
                        ratioNSP2=mean(ratioNSP2),Ratio=mean(ratioOther))
# Colnames
colnames(ClusterViro)<-c("Protein","Mean(Dist)","Sd(Dist)","Mean(NSP2)","Mean(Other)")
# Convert to data frame
NewData<-as.data.frame(ClusterViro[,-1]) 
row.names(NewData)<-as.character(ClusterViro$Protein)
# Function for the hierarchical cluster analysis
hclustfunc <- function(x) hclust(x, method="complete")
# Function to compute the euclidean distance
distfunc <- function(x) dist(x,method="euclidean")
d <- distfunc(NewData) # Compute the distance
fit <- hclustfunc(d) # Cluster
plot(fit) # Plot only he cluster
groups <- cutree(fit, k=5)  # Cut the cluster into k groups 
cols <- brewer.pal(max(groups), "Set1") # Set a color palette
# Save the heatmap as pdf
pdf(file = "HeatMapNSP2.pdf",width = 5,height = 5)
heatmap.2(as.matrix(NewData),dendrogram="row",trace="none", margin=c(8,9), 
          hclust=hclustfunc, distfun=distfunc, RowSideColors=cols[groups],srtCol=45)
dev.off()