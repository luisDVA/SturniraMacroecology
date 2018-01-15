# SturniraMacroecology
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

### File descriptions 
*Data*

cranBiocl.csv -- Cranial measurements
flwbiocl.csv -- Forearm length (FL) and body mass (W) data
fullDataF.csv -- Specimens with both cranial and external (FL,W) data
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
Baur_isosize_fn.R -- Function to calculate skull isosizes, from Baur & Leuenberger (2011), sourced by 03_LMMisosize.R
