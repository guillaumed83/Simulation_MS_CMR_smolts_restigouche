# Restigouche smolt monitoring program

This repositery contains the scripts necessary to conduct the simulation analysis presented in the MS 
"Estimating Atlantic salmon smolt abundance in a large Canadian catchment using multiple rotary screw traps" 

## Background
A set of models are tested to evaluate their suitability to the Restigouche smolt CMR dataset.

## Files available
there are 2 folders: 
1. **Models_simulation**: Contains the various model files needed for the analysisdescribed in the MS - models are labelled similarly M1-M5.  
2. **Scripts**: Contains various scripts to generate model inputs and running the analysis

+ * *Final_data_simulation_replicates_logN.R* * documents how the different datasets and initial values used in the simulation/jags are generated.

+ all files starting with * *Final_run_simulations_xxx* * are needed to run the inference in jags

+ * *Final_analyse_simul_LogN.R* * documents how summary statistics and figures are generated



