# Error-Growth-Borrus-2021
Files relate to the work Marshall Borrus (msborrus@gmail.com) did for the paper Midlatitude Error Growth in Atmosphere GCMs (https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2021GL096126).

# Files

## AM4_Master.m
Should calculate EGR, saturation time, mean values, dry and moist N^2 values, and creates the AM4_Data variable which is used for plotting.
## AM4_Master_2.m
Does the same as AM4_Master for the updated AM4 runs - but does additional plotting and analysis. 
## Dycore_Master.m
Does most of what AM4_Master does, but for the dycore data, which is shared differently.
## Master_Plotter.m
Generates the figures in the paper given the correct input data (stored on box and OAK). 
## SCRATCH_TcompleteCode.m
Allows you to create a surface plot of latitude, temperature, and saturation time for RMS of Temperature.  
## SCRATCH_UandTcompleteCode.m
Allows you to create a surface plot of latitude, temperature, and saturation time for RMS of horizontal winds U 
## Saturation_Time_Shortcut.m
Calculates saturation time without all the other plots and calculations generated from AM4_Master#.m
## dycore_interp.m
This file interpolates the raw output from dycore runs, and is meant to be run in parallel. 
## SST Play
Allows you to change SST profiles of AM4 runs
