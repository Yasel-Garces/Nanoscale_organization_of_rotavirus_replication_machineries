# Circumference Hypothesis Validation
## General Information
To proof the hypothesis that the viral elements of the VPs can be approximated through circumferences, we carried up a series of experiments based on the comparison between the circumference obtained by the algorithm DLSFC, and a least squared ellipse resulting of the "Direct Least Square Fitting Ellipse" (DLSFE) algorithm ([Fitzgibbon et al., 1999](https://ieeexplore.ieee.org/document/765658)). Under the assumption that the protein P can be approximate by a circumference, we expect that does not exist significative statistical differences between the radius of the circumference and each one of the ellipse semi-axis. Note that this approach is independent of which protein we choose to test the circularity hypothesis.

For each protein combination, we accomplish a Shapiro–Wilk hypothesis test to study if any of the three distributions functions (circumference radius (r), semi-major axis (a) and semi-minor axis (b)) came from a normally distributed population [Shapiro and Wilk, 1965](https://www.jstor.org/stable/2333709?seq=1#page_scan_tab_contents). The results revealed that, our data does not seem to follow a Gaussian distribution, and as consequence, it is impossible to use a parametric test (for example t-student). Based on our data conditions, we considered the Mann-Whitney test [Mann and Whitney, 1947](https://projecteuclid.org/euclid.aoms/1177730491); [Hollander et al., 2013](https://www.wiley.com/en-us/Nonparametric+Statistical+Methods%2C+3rd+Edition-p-9780470387375)) to evaluate the differences between the radii of the circumferences and each of the ellipses semi-axis. This is a nonparametric test that can be applied when the observations are independent (our variables are independent because the radius and the semi-axes were obtained by different algorithms), the variables are ordinal (also our variables meet this requirement because they are numeric), and finally, we want to test if the radius of the circumferences and each one of the ellipses semi-axis, are equal (null hypothesis) or are not (alternative hypothesis).     

For more details consult [our article](https://elifesciences.org/articles/42906):     
Garcés et al. Nanoscale organization of rotavirus replication machineries. eLife 2019;8:e42906. https://elifesciences.org/articles/42906, doi: 10.7554/eLife.42906.     

## Files available in this repository
1. **Run_Validation.m**: Main script to run the circumference hypothesis validation.    
```Run_Validation()```: It is only necessary to specify the directories where the images for the validation are. Also, it is possible to control other variables, for details see the comments inside the script.    
2. **fit_circumference_LSFC.m**: Adjust a circumference to a set of points in R^2 through the algorithm DLSFC (see the article for details).     
```fit_circumference_LSFC(x,y)``` where (x,y) are points in R^2.     
3. **fit_ellipse_LSFE**: It adjust a ellipse to points in R^2 using the algorithm Direct Least Square Fitting Ellipse ([DLSFE](http://autotrace.sourceforge.net/WSCG98.pdf)): R. Halir and J. Flusser. Numerically stable direct least-squares fitting of ellipses. In Sixth International Conference of Computer Graphics and Visualization, 1998.     
```fit_ellipse_LSFE(x,y)``` where (x,y) are the data points.     
4. **Circumference Hypothesis Validation.R**: Statistical analysis of the comparison between the major and minor semi-axis of the ellipse with the radius of the adjusted circumference.     
Example: ```Ctrl+Enter``` to send an individual line of code from the editor to the console or to send a block of highlighted code to the console. ```Ctrl+Shift+Enter``` to send the entire script to the console.    
5. **Result.csv**: File with the result of the validation:      
    * Column 1 (Combination): Name of the protein combined with NSP2.     
    * Column 2 (Protein): Adjusted proteins through the algorithms DLSFC and DLFE.     
    * Column 3 (MajorSemiAxis): Major semi-axis of the ellipse.      
    * Column 4 (MinorSemiAxis): Minor semi-axis of the ellipse.      
    * Column 5 (RadiusCirc): Circumference radius.     
