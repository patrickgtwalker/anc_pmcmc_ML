remotes::install_github("mrc-ide/mcstate", upgrade = TRUE)
remotes::install_github("mrc-ide/dust", upgrade = TRUE)
remotes::install_github("mrc-ide/mode", upgrade = TRUE)
remotes::install_github("mrc-ide/odin.dust@8aef08d", upgrade = TRUE)
remotes::install_github("mrc-ide/didehpc", upgrade = TRUE)
library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library(dplyr)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)
library(zoo)
library(lubridate)

#Required functions
source('shared/run_pmcmc.R')
source('nnp/in_development/run_pmcmc_pg.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')

##Import NNP data (incomplete data sets)
##Primigrav
data_raw_ng_pg_asa <- readRDS('nnp/data/data_raw_ng_pg_asa.RDS')
data_raw_ng_pg_ifenorth <- readRDS('nnp/data/data_raw_ng_pg_ifenorth.RDS')
data_raw_ng_pg_ejigbo <- readRDS('nnp/data/data_raw_ng_pg_ejigbo.RDS')
data_raw_ng_pg_moro <- readRDS('nnp/data/data_raw_ng_pg_moro.RDS')

data_raw_bf_pg_banfora <- readRDS('nnp/data/data_raw_bf_pg_banfora.RDS')
data_raw_bf_pg_orodara <- readRDS('nnp/data/data_raw_bf_pg_orodara.RDS')
data_raw_bf_pg_gaoua <- readRDS('nnp/data/data_raw_bf_pg_gaoua.RDS')

data_raw_mz_pg_guro <- readRDS('nnp/data/data_raw_mz_pg_guro.RDS')
data_raw_mz_pg_chemba <- readRDS('nnp/data/data_raw_mz_pg_chemba.RDS')
data_raw_mz_pg_changara <- readRDS('nnp/data/data_raw_mz_pg_changara.RDS')

nnp_pg_list <- list(data_raw_bf_pg_banfora,data_raw_bf_pg_gaoua,data_raw_bf_pg_orodara,
                    data_raw_mz_pg_changara,data_raw_mz_pg_chemba,data_raw_mz_pg_guro,
                    data_raw_ng_pg_asa,data_raw_ng_pg_ejigbo,data_raw_ng_pg_ifenorth,data_raw_ng_pg_moro)


country <- c('Burkina Faso','Burkina Faso','Burkina Faso',
             'Mozambique','Mozambique','Mozambique',
             'Nigeria','Nigeria','Nigeria','Nigeria')
admin <- c('Cascades','Sud-Ouest','Haut-Bassins',
           'Tete','Sofala','Manica',
           'Kwara','Osun','Osun','Kwara')

##Test run_pmcmc function##
source('nnp/in_development/run_pmcmc_pg.R')

test_run_pg_corr <- run_pmcmc_pg(data = data_raw_bf_pg_banfora,
                          n_particles = 10,
                          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                          max_EIR=1000,
                          max_steps = 1e7,
                          atol = 1e-5,
                          rtol = 1e-6,
                          n_steps = 3,
                          n_threads = 2,
                          lag_rates = 10,
                          country = 'Burkina Faso',
                          admin_unit = 'Cascades',
                          seasonality_on = 0,
                          state_check = 0)
windows(10,7)
plot_particle_filter(test_run_pg_corr$history,true_history=data_raw_bf_pg_banfora,times=data_raw_bf_pg_banfora$t)

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "nnp/in_development/run_pmcmc_pg.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts.temp", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust@8aef08d",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "32Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
config_dide <- didehpc::didehpc_config(template = "8Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj$login()
obj$cluster_load(TRUE)

##Submit bulk non-seasonality runs
nnp_pgcorr_bulk_std_test <- obj$enqueue_bulk(1:10, function(i,data_site){
  run_pmcmc_pg(data = data_site[[i]],
            n_particles = 10,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 10,
            n_threads = 8,
            lag_rates = 10,
            seasonality_on = 0,
            state_check = 0)
},data_site=nnp_pg_list)
nnp_pgcorr_bulk_std_test$status() #'draughty_cardinal'
# obj$unsubmit(nnp_pgcorr_bulk_std_test$ids)

nnp_pgcorr_bulk_std <- obj$enqueue_bulk(1:10, function(i,data_site){
  run_pmcmc_pg(data = data_site[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 4,
               lag_rates = 10,
               seasonality_on = 0,
               state_check = 0)
},data_site=nnp_pg_list)
nnp_pgcorr_bulk_std$status() #'falconiform_amurstarfish'
nnp_pgcorr_bulk_std <- obj$task_bundle_get('falconiform_amurstarfish')
nnp_pgcorr_bulk_std$tasks[[1]]$log()
obj$unsubmit(nnp_pgcorr_bulk_std$ids)
nnp_pgcorr_bulk_std_short <- obj$task_bundle_get('patriarchical_africanparadiseflycatcher')
nnp_pgcorr_bulk_std_short$status()
plot_particle_filter(nnp_pgcorr_bulk_std_short$tasks[[1]]$result()$history,true_history=nnp_pg_list[[1]],times=nnp_pg_list[[1]]$t)
nnp_pgcorr_bulk_std_results <- lapply(1:10, function(id){
  nnp_pgcorr_bulk_std$tasks[[id]]$result()
})

nnp_ng_pgcorr_bulk_std <- obj$enqueue_bulk(7:10, function(i,data_site){
  run_pmcmc_pg(data = data_site[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 1,
               lag_rates = 10,
               seasonality_on = 0,
               state_check = 0)
},data_site=nnp_pg_list)
nnp_ng_pgcorr_bulk_std$status()#'pseudoregal_katydid' submitted 2 Mar 10:50am
nnp_ng_pgcorr_bulk_std_results <- lapply(1:4, function(id){
  nnp_ng_pgcorr_bulk_std$tasks[[id]]$result()
})
nnp_pgcorr_bulk_std_results_update <- append(nnp_pgcorr_bulk_std_results[1:6],nnp_ng_pgcorr_bulk_std_results)

nnp_ng_pgcorr_bulk_std$times()
obj$login() 
obj$cluster_load(TRUE)

nnp_pgcorr_bulk_seas_test <- obj$enqueue_bulk(1, function(i,data_site,country,admin){
  run_pmcmc_pg(data = data_site[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 8,
               lag_rates = 10,
               seasonality_on = 1,
               country = country[i],
               admin_unit = admin[i],
               state_check = 0)
},data_site=nnp_pg_list,country=country,admin=admin)
nnp_pgcorr_bulk_seas_test$status() #'slimline_siamesecat'

nnp_pgcorr_bulk_seas_test$tasks[[1]]$log()
obj$unsubmit(nnp_pgcorr_bulk_seas_test$ids)

nnp_pgcorr_bulk_seas <- obj$enqueue_bulk(1:10, function(i,data_site,country,admin){
  run_pmcmc_pg(data = data_site[[i]],
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
               country = country[i],
               admin_unit = admin[i],
               state_check = 0)
},data_site=nnp_pg_list,country=country,admin=admin)
obj$login()
nnp_pgcorr_bulk_seas <- obj$task_bundle_get('nonhistoric_dwarfmongoose')
nnp_pgcorr_bulk_seas$status() #'nonhistoric_dwarfmongoose'
nnp_pgcorr_bulk_seas$times()
nnp_pgcorr_bulk_seas_results <- lapply(1:10, function(id){
  nnp_pgcorr_bulk_seas$tasks[[id]]$result()
})
nnp_ng_pgcorr_bulk_seas <- obj$enqueue_bulk(7:10, function(i,data_site,country,admin){
  run_pmcmc_pg(data = data_site[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 1,
               lag_rates = 10,
               seasonality_on = 1,
               country = country[i],
               admin_unit = admin[i],
               state_check = 0)
},data_site=nnp_pg_list,country=country,admin=admin)
nnp_ng_pgcorr_bulk_seas$status() #'jesting_zanzibardaygecko' submitted 2 Mar 10:48am
nnp_ng_pgcorr_bulk_seas_results <- lapply(1:4, function(id){
  nnp_ng_pgcorr_bulk_seas$tasks[[id]]$result()
})
obj$unsubmit(nnp_ng_pgcorr_bulk_seas$ids)
nnp_pgcorr_bulk_seas_results_update <- append(nnp_pgcorr_bulk_seas_results[1:6],nnp_ng_pgcorr_bulk_seas_results)
View(nnp_pgcorr_bulk_seas_results_update)
