# Neuroscience-Project
This project attempts to evaluate the ability of a mathematical model to predict the location of electrodes based on the potential decay between them\\\
The following is an explanation of the scripts and datasets found in the repository\\\
Electrode_Positons.mat - a 3d vector which describes the location of each of the electrodes in the array\\\
MCRackRecording.m - An auxiliary function which extracts the data from data extracted through MCRack software\\\
Near_Sol.mat - A struct which contains the optimization data when the optimization is initialized as close as possible to the true solution as possible\\\
Optimization_Sol.mat - A struct which contains the optimization data when the optimization is initialized at random points within the boundary conditions\\\
Peak_Data.mat - A matrix containing the minimum of the spikes which characterize the recorded potential from each of the electrodes\\\
Potential_Matrix.mat - A matrix containing a potential used to characterize the potential measured between each to electrodes after averaging\\\
Regression_Fits.mat - A struct which contains the data for the regression analysis described in the paper.\\\
analysis.m - The main script\\\
layout_100_12x12_120Channels.mat - a dataset accessed by MCRackRecording.m decribing the positions of the electrodes\\\
nearest_upper_square.m and randinterval.m are auxiliary functions accessed in the main script\\\
