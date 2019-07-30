# Concentricity Hypothesis Validation
## General Information
In this section we carried up a study to validate the viral elements concentricity hypothesis based on the distance between the centers of the adjusted circumferences of both proteins.     

In each experiment, we have a couple of proteins (normally NSP2 as a central protein). For each combination (for example NSP2-NSP4 and NSP2-VP6) we adjust the least square circumference through the algorithm DLSFC to each protein independently and compute the distance between their centers.     

For more details consult [our article](https://elifesciences.org/articles/42906):     
Garcés et al. Nanoscale organization of rotavirus replication machineries. eLife 2019;8:e42906. https://elifesciences.org/articles/42906, doi: 10.7554/eLife.42906.     

## Files available in this repository
1. **Run_Concentricity_Validation.m**: Main program to validate the hypothesis that the central and the accompainy protein are concentric.     
```Run_Concentricity_Validation()```: It is only necessary to specify the directories where the images for the validation are. Also, it is possible to control other variables, for details see the comments inside the script.    
2. **fit_circumference_LSFC.m**: Adjust a circumference to a set of points in R^2 through the algorithm DLSFC (see the article for details).     
```fit_circumference_LSFC(x,y)``` where (x,y) are points in R^2.     
3. **Concentricity_Validation.R**: Statistical analysis to validate the viral elements concentricity hypothesis based on the distance between the centers of the adjusted circumferences of both proteins.     
Example: ```Ctrl+Enter``` to send an individual line of code from the editor to the console or to send a block of highlighted code to the console. ```Ctrl+Shift+Enter``` to send the entire script to the console.    
4. **DistanceCenterProtein.csv**: Data file with the structure:
	* Column 1 (Combination): Name of the protein combined with NSP2.     	
	* Column 2 (Distance): Distance between the center of the adjusted circumferences to NSP2 and the accompanying protein.