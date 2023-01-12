# anc_pmcmc
Application of mcstate pmcmc to ANC data from various settings

## sim
-sim_main.R: Generates simulated data, fits pMCMC, and summarises results. 

-data_gen.R: Generates simulated data

-run_pmcmc.R: Fits the mcstate pMCMC

-Folder sim_datasets: 3 example simulated datasets


## nnp
-nnp_submit.R: Submits the pmcmc to the cluster

-nnp_summarize.R: Creates plots to summarize results

-nnp_readdata.R: Processes raw data for use with the pmcmc

-run_pmcmc.R: Fits the mcstate pMCMC


## wkenya
-wkenya_submit.R: Submits pmcmc to cluster

-ANC_cMIS_summary.R: Produce plots

-run_pmcmc.R: Fits the mcstate pMCMC


## shared: Support files used across projects
-odinmodelmatchedstoch.R: Odin model with stochastic update of EIR

-odin_model_stripped_matched.R: Odin model without the stochastic update, used to generate simulated data.

-plot_particle_filter.R: Plots the pmcmc trajectory against the observed prevalence.

-addCIs.R: Add binomial confidence intervals to observed prevalence data