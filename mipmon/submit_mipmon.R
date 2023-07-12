library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library(tidyverse)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(zoo)
library(lubridate)

#Required functions
source('shared/run_pmcmc.R')
source('nnp/in_development/run_pmcmc_mg.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')


test_run_mg_mipmon <- run_pmcmc_mg(data_raw_pg = mipmon_data_pg[[1]],
                                 data_raw_mg = mipmon_data_mg[[1]],
                                 n_particles = 10,
                                 proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                 max_EIR=1000,
                                 max_steps = 1e7,
                                 atol = 1e-5,
                                 rtol = 1e-6,
                                 n_steps = 20,
                                 n_threads = 2,
                                 lag_rates = 10,
                                 country = 'Mozambique',
                                 admin_unit = 'Maputo',
                                 seasonality_on = 0,
                                 state_check = 0)
windows(10,7)
plot_particle_filter(test_run_mg_mipmon$history,true_history=mipmon_data_pg[[1]],times=mipmon_data_pg[[1]]$t)


##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "nnp/in_development/run_pmcmc_pg.R",
             "nnp/in_development/run_pmcmc_mg.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts_may23", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "20Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj$login()
obj$cluster_load(TRUE)

test_run <- obj$enqueue_bulk(district, function(country){
  lapply(modeling_smc(district=i,country='Burkina Faso',region='Region'))
},region,country)


mipmon_mg_bulk_seas <- obj$enqueue_bulk(1:4, function(i,data_pg,data_mg){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 1,
               country = 'Mozambique',
               admin_unit = 'Maputo',
               state_check = 0)
},data_pg=mipmon_data_pg,data_mg=mipmon_data_mg)
mipmon_mg_bulk_seas$status() #'springy_antlion' 
mipmon_mg_bulk_seas$tasks$e05a353cd298e3b8e92627a6a454267d$log()
mipmon_mg_bulk_seas <- obj$task_bundle_get('springy_antlion')
obj$task_bundle_info()
mipmon_mg_bulk_seas_results <- lapply(1:4, function(id){
  mipmon_mg_bulk_seas$tasks[[id]]$result()
})
mipmon_mg_bulk_seas_results[[1]]$history['prev_05',,]
