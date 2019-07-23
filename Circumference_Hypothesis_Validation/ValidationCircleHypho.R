#' Comparison between the major and minor semi axis of the ellipse with the radius 
#' of the adjusted circumference. Under the assumption that the protein P can be 
#' approximate by a circumference, we expect that does not exist significative 
#' statistical differences between the radius of the circumference (r) and each 
#' one of the ellipse semi-axis (a, b).
#' @param "Result.csv": Dataset with the next structure:
#' Column 1 (Combination): Name of the protein combined with NSP2
#' Column 2 (Protein): Adjusted proteins through the algorithms DLSFC and DLFE
#' Column 3 (MajorSemiAxis): Major semi-axis of the ellipse
#' Column 4 (MinorSemiAxis): Minor semi-axis of the ellipse
#' Column 5 (RadiusCirc): Circumference radius
#' @return Shapiro-Wilk test of normality
#' @return Results of the regression linear model between the radii of the 
#' circumferences and the ellipse semi axis.
#' @return Results of the ‘Mann-Whitney’ test
#' @author Yasel garces (88yasel@gmail.com)
#========================================================================
# Libraries
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Change directory
setwd("/home/yasel/Dropbox/Paper Viroplasmas/Programs/Validation_Circumference_Hypothesis")

# Load data
semiaxis<-read.csv("Result.csv")
#========================================================================
# For each protein combination, we accomplish a Shapiro–Wilk hypothesis test to study if
# any of the three distributions functions (circumference radius (r), semi-major axis (a) and
# semi-minor axis (b)) came from a normally distributed population
# Major semi axis
by(semiaxis,
   semiaxis[,c("Combination","Protein")],
   function(x) shapiro.test(x$MajorSemiAxis))
# Minor semi axis
by(semiaxis,
   semiaxis[,c("Combination","Protein")],
   function(x) shapiro.test(x$MinorSemiAxis))
# Circumference radius
by(semiaxis,
   semiaxis[,c("Combination","Protein")],
   function(x) shapiro.test(x$RadiusCirc))
#========================================================================
# Adjust a linear regressio model for each protein. Independet variable = Radius
# of the adjusted circumference; Dependent variable= Major Semi Axis/Minor Semi axis
notNSP2<-subset(semiaxis,Protein != "NSP2")
# MajorSemiAxis ~ RadiusCirc
ggplot(notNSP2,aes(x = RadiusCirc,y = MajorSemiAxis,color=Protein))+geom_point()+
  facet_grid(.~Combination)+geom_smooth(method='lm',formula=y~x)+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))
by(semiaxis,
   semiaxis[,c("Combination","Protein")],
   function(x) summary(lm(x$MajorSemiAxis ~ x$RadiusCirc -1)))
# MinorSemiAxis ~ RadiusCirc
ggplot(notNSP2,aes(x = RadiusCirc,y = MinorSemiAxis,color=Protein))+geom_point()+
  facet_grid(.~Combination)+geom_smooth(method='lm',formula=y~x)+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))
by(semiaxis,
   semiaxis[,c("Combination","Protein")],
   function(x) summary(lm(x$MinorSemiAxis ~ x$RadiusCirc -1)))
#========================================================================
# Hypothesis test (Test of Equal or Given Proportions)
# The ‘Mann-Whitney’ test was carry out to evaluate the differences between
# the radii of the circumferences and each of the ellipses semi-axis.
# Minor semi axis
by(semiaxis,
   semiaxis[,c("Combination","Protein")],
   function(x) wilcox.test(x$MinorSemiAxis,x$RadiusCirc,conf.int = TRUE))
# Major semi axis
by(semiaxis,
   semiaxis[,c("Combination","Protein")],
   function(x) wilcox.test(x$MajorSemiAxis,x$RadiusCirc,conf.int = TRUE))
#========================================================================