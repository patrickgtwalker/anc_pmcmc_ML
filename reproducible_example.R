library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)

library(zoo)
library(lubridate)
library(dplyr)


##Run locally
data_sim_comptest3 <- readRDS('sim/sim_datasets/data_sim3.RDS')

source('shared/run_pmcmc.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')

sim_local_test<- run_pmcmc(data = data_sim_comptest3,
                      n_particles = 10,
                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                      max_EIR=1000,
                      max_steps = 1e7,
                      atol = 1e-5,
                      rtol = 1e-6,
                      n_steps = 10,
                      n_threads = 2,
                      lag_rates = 10,
                      seasonality_on = 0)


