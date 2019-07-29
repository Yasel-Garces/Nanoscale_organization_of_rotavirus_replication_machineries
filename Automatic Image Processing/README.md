# Automatic Image Processing through the algorithm VP-DLSFC

## General Information
The use of primitive models for the segmentation of the SRM images has many benefits, like computational efficiency, simple programmation, and understandable information about the objects that were segmented. In this regard, we developed a simple method for fitting a circumference to scattered data, which we called “Direct Least Squares Fitting Circumference” (DLSFC). This approach considers basic mathematical analysis tools for the computation of the extreme value of a continuous function with N variables.

The proposed segmentation algorithm is composed by 3 steps ([Appendix 1 Figure 1](https://elifesciences.org/articles/42906)). The first is a manual pre-segmentation that guarantees the existence of one and only one viroplasm per image ([Appendix 1 Figure 1A](https://elifesciences.org/articles/42906)). In the second step, we adjust a circumference to the reference protein (these are paired experiments, where NSP2 is present in all the combinations as reference protein) through the algorithm DLSFC ([Appendix 1 Figure 1B](https://elifesciences.org/articles/42906)), and finally the radius of the accompanying protein is adjusted using the [equation (12)](https://elifesciences.org/articles/42906). Note that, the center of the accompanying protein is the same that the center of the reference protein (concentric model). Out of the manual pre segmentation step, this algorithm is automatic, deterministic, not iterative and with a linear computational complexity.

This repository collects all the source codes used to automatically analyze and quantify the relative distribution of seven viral proteins within and around the VPs. This research was published in eLife [link to the article](https://elifesciences.org/articles/42906) and can be cited as:

Garcés et al. Nanoscale organization of rotavirus replication machineries. eLife 2019;8:e42906. https://elifesciences.org/articles/42906, doi: 10.7554/eLife.42906.

## Functions available in this repository
This repository contains 4 matlab scripts, these are:
1. **Run_Image_Analysis.m**: Main script to read and process all the images for all the hours post-infection.    
``Run_Image_Analysis()``: Only it is necessary to specify inside the function the path where the images are. For more details see the code.
2. **radius_accompanying_protein.m**: Runs the last step of the algorithm VP-DLSFC. The center of the reference protein is taken as the center of the accompanying protein, and then the radius of the fit circumference for this second protein is computed.    
```radius_accompanying_protein([26, -24],x,y)``` where (x,y) are points in R^2.
3. **fit_circumference_LSFC.m**: Adjust a circumference to a set of points in R^2 through the algorithm DLSFC (see [the article](https://elifesciences.org/articles/42906) for details).    
```fit_circumference_LSFC(x,y)``` where (x,y) are points in R^2.
4. **csvwriteh.m**: This function write matrix to a csv file with header. This code was obtained from [this blog](https://www.nesono.com/node/415).    
```csvwriteh( filename, data, header )``` where "FILENAME" is filename for csv output file, "DATA" is a matrix with data for csv file and "HEADER" is a cell array with names per column.
