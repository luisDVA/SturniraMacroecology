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

# File descriptions 

cranBiocl.csv - Cranial measurements

flwbiocl.csv - Forearm length (FL) and body mass (W) data

fullDataF.csv - Specimens with both cranial and external (FL,W) data

spOccPatterns.csv - Occurrence pattern for each locality (sympatry or allopatry)
