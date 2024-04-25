# Computational modeling of YAP phosphorylation within adhesions

This repository contains the MATLAB codes to run a particle-based simulation for YAP phosphorylation within adhesions in a 3D simulation domain. The model involves the diffusion of YAP, YAP-adhesion binding, YAP phosphorylation within the adhesion, and (p)YAP release from the adhesion. The simulation also takes into account adhesion turnover and the formation of new adhesions.

The main function is written in "stochastic_model_yap/stochastic_model.m". Script "run.m" sets the input parameters, calls the main function to run the simulation, and plots the pYAP ratio vs time (averaged over three replicates). Each input parameter has been defined in the scripts and functions. To view the 3d simulation domain at specific intervals, you can uncomment the line where "PlotPosition" function is called in "stochastic_model_yap/stochastic_model.m". The figures will be printed in the folder “results”.

The parameter "lifetime" in "stochastic_model_yap/stochastic_model.m" can be adjusted to run the simulations for different values of lifetime (e.g. 60 s, 135 s). To simulate with no adhesion turnover, the value of lifetime can be set to a very large number (larger than the expected simulation time).  

## The output folder contains: 
-figures for time evolution of pYAP ratio
-A mat file with the output data for all steps




## Simulation results
The results used to generate the figures in the manuscript are also copied to the main folder. Script "plot_results" can be used to make all the figures in the manuscript.




