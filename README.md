# YAP-phosphyrylation-within-adhesion-stochastic-model

This repository contains the MATLAB codes to run a particle-based simulation for YAP phosphorylation within adhesions in a 3D simulation domain. The model involves the diffusion of YAP, YAP-adhesion binding, YAP phosphorylation within the adhesion, and (p)YAP release from the adhesion. The simulation also takes into account adhesion turnover and the formation of new adhesions.

The main function is written in "stochastic_model_yap/stochastic_model.m". Script "run.m" sets the input parameters, calls the main function to run the simulation, and plots the pYAP ratio vs time (averaged over three replicates). Each input parameter has been defined in the scripts and functions. To view the 3d simulation domain at specific intervals, you can uncomment the line where "PlotPosition" function is called in "stochastic_model_yap/stochastic_model.m". The figures will be printed in the folder “results”.

## The output folder contains: 
-figures for time evolution of pYAP ratio (similar to figure S6,S7 left columns in the SI) for 10 cases; 5 black and 5 blue data points for the corresponding lifetime in figure 4a. What you see in figure 4a is the pYAP ratio averaged over three iterations (iteration = 3 in the code now but you can set it to 1 for your tests).
-A mat file with the output data for all steps

Parameter "lifetime" in "stochastic_model_yap/stochastic_model.m" can be adjusted to run the simulations for different values of lifetime (e.g. 60 s, 135 s). To simulate with no adhesion turnover, the value of lifetime can be set to a very large number (larger than the expected simulation time).  


## Simulation results
The results used to generate the figiures in the manuscript is also 




