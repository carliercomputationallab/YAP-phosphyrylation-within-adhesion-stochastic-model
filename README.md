# Computational modeling of YAP phosphorylation within adhesions

This repository contains MATLAB code to execute a particle-based simulation focused on YAP phosphorylation within adhesions in a 3D simulation domain. The model includes YAP diffusion, YAP-adhesion binding, YAP phosphorylation within the adhesion, and (p)YAP release from the adhesion. Additionally, the simulation takes into account adhesion turnover and the generation of new adhesions.

The main function is implemented in "stochastic_model_yap/stochastic_model.m". The script "run.m" configures the input parameters, calls the main function to execute the simulation, and plots the pYAP ratio versus time (averaged over three replicates). Each input parameter is clearly defined within the scripts and functions. To visualize the 3D simulation domain at specific intervals, you can uncomment the line calling the "PlotPosition" function in "stochastic_model_yap/stochastic_model.m". The resulting figures will be saved in the "results" folder.

The parameter "lifetime" in "stochastic_model_yap/stochastic_model.m" can be adjusted to conduct simulations for different lifetime values (e.g., 60 s, 135 s). To simulate with no adhesion turnover, set the value of lifetime to a sufficiently large number (greater than the expected simulation time). 

## Contents of the output folder
- Figures illustrating the time evolution of pYAP ratio
- A .mat file containing the output data for all steps

## Plotting figures in the manuscript
The simulation results used to generate the figures in the manuscript can be found in the main folder. The script "plot_results" can be used to plot all figures in the manuscript and supplemental information. 



