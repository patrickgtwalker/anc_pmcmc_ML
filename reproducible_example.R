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


##Run on cluster
root <- "T:/jth/contexts.new"
sources <- c("shared/run_pmcmc.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts.new", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "8Core",cores = 2, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx,config = config_1)

sim_cluster_test <- obj$enqueue(run_pmcmc(data = data_sim_comptest3,
                           n_particles = 10,
                           proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                           max_EIR=1000,
                           max_steps = 1e7,
                           atol = 1e-5,
                           rtol = 1e-6,
                           n_steps = 10,
                           n_threads = 2,
                           lag_rates = 10,
                           seasonality_on = 0))
sim_cluster_test$status()
sim_cluster_test$log()
