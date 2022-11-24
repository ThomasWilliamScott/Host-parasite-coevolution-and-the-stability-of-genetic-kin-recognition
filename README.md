# Host-parasite-coevolution-and-the-stability-of-genetic-kin-recognition
This repository contains all supplementary code for "Host-parasite coevolution and the stability of genetic kin recognition" (Scott, Grafen &amp; West, unpublished)

This repository contains three types of file:

  1) Scripts for generating data (specifically, these scripts numerically implement our mamthematical model). 
  2) Data files (these are the ".mat" files, which are saved outputs of the data generating scripts).
  3) Scripts for generating figures (specifically, these scripts load up the data files and use them to generate the figures).
  
There are three data generating scripts: "Script_for_generating_parameter_sweep_data.m", "Script_for_generating_single_trial_data.m", "Script_for_generating_initial_genotype_frequencies.m". Running the "Script_for_generating_parameter_sweep_data.m" script will generate long-term (equilibrium) results for a range of d (parasite virulence) and lag (parasite evolutionary lag) values, and save these results in matrices. Running the "Script_for_generating_single_trial_data.m" script will generate over-time (dynamical) results for a specific set of parameter values. Both of these scripts call the "Script_for_generating_initial_genotype_frequencies.m" in order to generate the initial genotype frequencies that start off each run.

The "parameter sweep" and "single trial" data files are named with slightly different formats. An example "parameter sweep" data file is: "alpha=1_scenario=2_trial=3_T=200000_tag=10_theta=0.25_b=0.3_c=0.1_dmin=0_lagmin=0_dmax=1_lagmax=100_dint=0.05_lagint=5_muC=0_mu=0.001majNeutral=0.9_majResist=0.9_majHelp=0.1_majChoice=1.mat". An example "single trial" data file is: "Single_trial_alpha=0_scenario=1_trial=1_T=75000_tag=10_theta=0.25_b=0.3_c=0.1_d=0.6_lag=20_muC=0_mu=0.001majNeutral=0.9_majResist=0.9_majHelp=0.1_majChoice=0.mat". The data files are named with all of the parameter values in the title. Multiple trials have been saved for each parameter combination. When plotting equilibrium results, we take averages over multiple trials, as this averages over the slight stochasticity involved in generating initial genotype frequencies. 

There are 10 data generating figures, each titled with the format: "Generate_Figure_...". Code is given for all figures that appear in the main text and supplementary information, except for conceptual figures that were generated on powerpoint.
