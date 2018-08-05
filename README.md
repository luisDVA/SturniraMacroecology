# SturniraMacroecology
## Co-occurrence and character convergence in two Neotropical bats
Code, data and data definition

The code is set up in a workflow that keeps data, code, and outputs in separate directories and uses the 'here' package to create relative paths. To reproduce the analyses, set up the following folders in your respective R working directory.


### Working directory

 * analysis code
   * store all R scripts here
 * data
   * store all csv files here
 * lmmEffects
   * model outputs will be written here 
 * figures (optional)
   * figures will be exported here
 * filtered
   * store scripts with _no_sing suffix here    
 * clustered
   * store scripts with _clustered suffix here
 * variability_tests
   * write outputs from tests of equality of coefficients of variation here

### File descriptions 
*Data*

cranBiocl.csv -- Cranial measurements  
flwbiocl.csv -- Forearm length (FL) and body mass (W) data  
spOccPatterns.csv -- Occurrence pattern for each locality (sympatry or allopatry)  

*R scripts*

01_LMMforearm.R -- Linear mixed model and residual randomization tests for forearm length  
02_LMMbodymass.R -- Linear mixed model and residual randomization tests for body mass  
03_LMMisosize.R -- Linear mixed model and residual randomization tests for skull isosize  
04_LMMheadlength.R -- Linear mixed model and residual randomization tests for head length (condylobasal length)  
05_effectsTable.R -- Collating and reporting the linear mixed model estimates  
06_anovaTable.R -- Reporting the ANOVAs from the linear mixed models  
07_spatialFig1.R -- Plotting the specimen localities  
08_measPlotsClean.R -- Descriptive plots of the specimen measurements  
09_effectsPlots.R -- Effects plot for species trait means  
10_cv_equality.R -- Test for equality in coefficients of variation  
Baur_isosize_fn.R -- Function to calculate skull isosizes, from Baur & Leuenberger (2011); sourced by 03_LMMisosize.R

*filtered*

LMM scripts with an added step to exclude localities where n=1

*clustered*

LMM scripts using spatial clusters instead of localities to define allopatry/sympatry.

---

SuppTables.Rmd will run the necessary analyses and create a binary output file with Supplementary Data S1-S8.
