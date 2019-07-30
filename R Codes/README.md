# Source R Codes
## General Information
This set of scripts allow to analyze statistically the spatial distribution of 7 viral proteins related to the viroplasms. The observed viral components are spatially organized as 5 concentric layers, in which NSP5 localizes at the center of the VPs, surrounded by a layer of NSP2 and NSP4 proteins, followed by an intermediate zone comprised of the VP1, VP2, VP6. In the outermost zone, we observed a ring of VP4 and finally a layer of VP7.    
To be more clear and to provide an easies way to reproduce our results, the scripts were split according to the reference protein, for example, "AnalysisNSP2.R" carry up the statistical analysis of the results obtained from the images with NSP2 as reference protein.    
For more details consult:    
Garcés et al. Nanoscale organization of rotavirus replication machineries. eLife 2019;8:e42906. https://elifesciences.org/articles/42906, doi: 10.7554/eLife.42906.

## Files available in this repository
1. **AnalysisNSP2.R**: This script analyses the distribution of the viral proteins into the viroplasm taking NSP2 as the reference protein, and return a set of graphics and statistics (see comments in the code for details).     
2. **AnalysisNSP4.R**: This script analyses the distribution of the viral proteins into the viroplasm taking NSP4 as the reference protein, and return a set of graphics and statistics (see comments in the code for details).     
3. **AnalysisNSP5.R**: This script analyses the distribution of the viral proteins into the viroplasm taking NSP5 as the reference protein, and return a set of graphics and statistics (see comments in the code for details).     
4. **Functions.R**: Set of functions that have been used in the above scripts.      
5. **LinearRegressionAnalysis_NSP2.R**: This script contains the linear regression analysis considering the radius of NSP2 as the independent variable.    
6. **ResultsCSV**: This directory contains the results of the VP-DLSFC segmentation algorithm by reference protein (NSP2, NSP4, and NSP5). All the csv files have the same structure:     
	* Column 1: Distance between the distribution of both proteins.    
	* Column 2: Radius of the circumference that adjust the central protein.     
	* Column 3: Radius of the circumference that adjust the other protein.     

In all R scripts press ```Ctrl+Enter``` to send an individual line of code from the editor to the console or to send a block of highlighted code to the console. ```Ctrl+Shift+Enter``` to send the entire script to the console.