library("odin")
library("odin.dust")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library(dplyr)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
library(lubridate)


##Test move of seasonality parameter look up to mpl
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')
source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')


test_mpl <- model_param_list_create(country = 'Burkina Faso',
                                    admin_unit = 'Cascades')

data_sim_comptest3 <- readRDS('sim/sim_datasets/data_sim3.RDS')
##Test run_pmcmc function##
test_run <- run_pmcmc(data = data_sim_comptest3,
                      n_particles = 10,
                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                      max_EIR=1000,
                      max_steps = 1e7,
                      atol = 1e-5,
                      rtol = 1e-6,
                      n_steps = 10,
                      n_threads = 2,
                      lag_rates = 10,
                      country = 'Burkina Faso',
                      admin_unit = 'Cascades')
plot_particle_filter(test_run$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)

data_raw_bf_pg_banfora <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_bf_pg_banfora.RDS')
source('shared/run_pmcmc.R')
test_run <- run_pmcmc(data = data_raw_bf_pg_banfora,
                      n_particles = 10,
                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                      max_EIR=1000,
                      max_steps = 1e7,
                      atol = 1e-5,
                      rtol = 1e-6,
                      n_steps = 10,
                      n_threads = 2,
                      lag_rates = 10,
                      country = 'Burkina Faso',
                      admin_unit = 'Cascades',
                      seasonality_on = 1)
test_run <- run_pmcmc(data = data_raw_bf_pg_banfora,
                      n_particles = 10,
                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                      max_EIR=1000,
                      max_steps = 1e7,
                      atol = 1e-5,
                      rtol = 1e-6,
                      n_steps = 10,
                      n_threads = 2,
                      lag_rates = 10,
                      country = 'Burkina Faso',
                      admin_unit = 'Cascades',
                      seasonality_on = 0,
                      state_check = 1)
windows(10,8)
plot_particle_filter(test_run$history,true_history=data_raw_bf_pg_banfora,times=data_raw_bf_pg_banfora$t)
1 - coda::rejectionRate(as.mcmc(test_run$mcmc)) ##Acceptance rate
coda::effectiveSize(as.mcmc(test_run$mcmc)) ##ESS
cov(result_32_200$pars) ##Covariance
summary(as.mcmc(test_run$mcmc)) ##Summarize mcmc run
plot(as.mcmc(test_run$mcmc)) ##Plot traces and distributions
