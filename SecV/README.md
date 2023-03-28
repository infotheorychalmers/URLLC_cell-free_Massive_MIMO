
This folder contains the files used to produce the plots in Section V. In particular, it contains:

1. .m files to perform numerical routines and obtain CDF points of the network availability; 
2. .m files to post-process the obtained numerical results and plot them; 
<!--3. .txt to run simulations in Chalmers' cluster; and folders with all the saved data from the run simulations.-->

1.1 *Functions to generate CDF points of the network availability*:

**Main function**: get_CDF_points.m 

NOTE: Set variable DEBUG = 1 to perform sanity checks in your computer. Otherwise, the function is ready to run in Chalmers' cluster. 

This function computes the average error probability of a randomly chosen UE that has been deployed, as the rest of UEs in the network, randomly over the scenario.

The scenario is generated using the function generateSetup.m

The channel estimates are obtained using the function functionChannelEstimates.m (or functionChannelEstimatesDL.m in the case of DL pilots)

The combiners/precoders are obtained using the function functionCombinerPrecoder.m

Other auxiliary functions:

golden_search.m (To optimize over s)
saddlepoint_approximation.m (computes the saddlepoint approximation of the average error probability) 
functionRlocalscattering.m (used inside generateSetup.m to generate the channel correlation matrix) 

get_CDF_points.m can generate CDF points of the network availability for three different network settings:

- CENTRALIZED: cell-free setting where the signal processing is performed at a CPU.
- DISTRIBUTED: cell-free setting where the channel estimation and spatial processing is performed at the AP, and the CPU performs extra signal processing only based on the statistics of the channel and the results obtained at each AP after performing the spatial processing.
- CELLULAR: Classical Massive MIMO cellular network. It also includes the case of small cellular networks.

1.2. *Functions to study channel hardening in a single-cell scenario*:

hardening_study_cell_free.m
hardening_study_cellular.m

2. *Auxiliary functions plot_*.m used to plot and generate .csv containing the data points used to generate the plots in the paper*

<!--3. Scripts ready to run in Chalmer's cluster. (job_*) -->



