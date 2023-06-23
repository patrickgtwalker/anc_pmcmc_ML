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
source('nnp/in_development/run_pmcmc_mg.R')
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

##Multigrav
data_raw_ng_mg_asa <- readRDS('nnp/data/data_raw_ng_mg_asa.RDS')
data_raw_ng_mg_ifenorth <- readRDS('nnp/data/data_raw_ng_mg_ifenorth.RDS')
data_raw_ng_mg_ejigbo <- readRDS('nnp/data/data_raw_ng_mg_ejigbo.RDS')
data_raw_ng_mg_moro <- readRDS('nnp/data/data_raw_ng_mg_moro.RDS')

data_raw_bf_mg_banfora <- readRDS('nnp/data/data_raw_bf_mg_banfora.RDS')
data_raw_bf_mg_orodara <- readRDS('nnp/data/data_raw_bf_mg_orodara.RDS')
data_raw_bf_mg_gaoua <- readRDS('nnp/data/data_raw_bf_mg_gaoua.RDS')

data_raw_mz_mg_guro <- readRDS('nnp/data/data_raw_mz_mg_guro.RDS')
data_raw_mz_mg_chemba <- readRDS('nnp/data/data_raw_mz_mg_chemba.RDS')
data_raw_mz_mg_changara <- readRDS('nnp/data/data_raw_mz_mg_changara.RDS')

nnp_mg_list <- list(data_raw_bf_mg_banfora,data_raw_bf_mg_gaoua,data_raw_bf_mg_orodara,
                    data_raw_mz_mg_changara,data_raw_mz_mg_chemba,data_raw_mz_mg_guro,
                    data_raw_ng_mg_asa,data_raw_ng_mg_ejigbo,data_raw_ng_mg_ifenorth,data_raw_ng_mg_moro)

country <- c('Burkina Faso','Burkina Faso','Burkina Faso',
             'Mozambique','Mozambique','Mozambique',
             'Nigeria','Nigeria','Nigeria','Nigeria')
admin <- c('Cascades','Sud-Ouest','Haut-Bassins',
           'Tete','Sofala','Manica',
           'Kwara','Osun','Osun','Kwara')
##Test run_pmcmc function##
source('nnp/in_development/run_pmcmc_mg.R')

test_run_mg_corr <- run_pmcmc_mg(data_raw_pg = data_raw_bf_pg_banfora,
                                 data_raw_mg = data_raw_bf_mg_banfora,
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
plot_particle_filter(test_run_mg_corr$history,true_history=data_raw_bf_pg_banfora,times=data_raw_bf_pg_banfora$t)
dim(t(test_run_mg_corr$history[1,,-1]))

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "nnp/in_development/run_pmcmc_pg.R",
             "nnp/in_development/run_pmcmc_mg.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts.temp", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust@8aef08d",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "32Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
config_single <- didehpc::didehpc_config(template = "GeneralNodes",cores =1, parallel = FALSE,wholenode = FALSE, cluster = 'fi--didemrchnb')
config_dide <- didehpc::didehpc_config(template = "8Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj <- didehpc::queue_didehpc(ctx,config = config_single)
obj$login()
obj$cluster_load(TRUE)

##Submit bulk non-seasonality runs
nnp_mgcorr_bulk_std_test <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
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
},data_pg=nnp_pg_list,data_mg=nnp_mg_list)
nnp_mgcorr_bulk_std_test$status() #'banal_eagle'
nnp_mgcorr_bulk_std_test$tasks[[1]]$log()
obj$unsubmit(nnp_pgcorr_bulk_std_test$ids)

nnp_mgcorr_bulk_std <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
               n_particles = 200,
               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
               max_EIR=1000,
               max_steps = 1e7,
               atol = 1e-5,
               rtol = 1e-6,
               n_steps = 1000,
               n_threads = 8,
               lag_rates = 10,
               seasonality_on = 0,
               state_check = 0)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list)
nnp_mgcorr_bulk_std$status() #'trifling_bongo'
nnp_mgcorr_bulk_std <- obj$task_bundle_get('trifling_bongo')
nnp_mgcorr_bulk_std_results <- lapply(1:10, function(id){
  nnp_mgcorr_bulk_std$tasks[[id]]$result()
})

nnp_ng_mgcorr_bulk_std <- obj$enqueue_bulk(7:10, function(i,data_pg,data_mg){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
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
},data_pg=nnp_pg_list,data_mg=nnp_mg_list)
nnp_ng_mgcorr_bulk_std$status() #'flannel_easternglasslizard' submitted 2 Mar 10:48am

nnp_mgcorr_bulk_seas <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
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
               country = country[i],
               admin_unit = admin[i],
               state_check = 0)
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_mgcorr_bulk_seas$status() #ittybitty_atlanticblackgoby - submitted 28 Feb 9:27am
#'amateurish_elver' - submitted 27 Feb 10:25am - [1] completed
nnp_mgcorr_bulk_seas$tasks[[1]]$log()
nnp_mgcorr_bulk_seas_270223 <- obj$task_bundle_get('amateurish_elver')
nnp_mgcorr_bulk_seas_270223$times()
nnp_mgcorr_bulk_bf_seas$status() #'chemophobic_macaque'
obj$unsubmit(nnp_mgcorr_bulk_seas$ids)
obj$login()
obj$cluster_load(TRUE)
nnp_mgcorr_bulk_bf_seas$tasks[[1]]$log()
obj$config

nnp_mgcorr_bulk_seas_single <- obj$enqueue_bulk(1:10, function(i,data_pg,data_mg,country,admin){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
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
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_mgcorr_bulk_seas_single$status() #'monotonous_betafish' - submitted 1 Mar 9:35am
nnp_mgcorr_bulk_seas_single$times()
nnp_mgcorr_bulk_seas_single <- obj$task_bundle_get('monotonous_betafish')
nnp_mgcorr_bulk_seas_results <- lapply(1:6, function(id){
  nnp_mgcorr_bulk_seas_single$tasks[[id]]$result()
})
obj$unsubmit(nnp_mgcorr_bulk_seas_single$ids[7:10])

nnp_ng_mgcorr_bulk_seas_single <- obj$enqueue_bulk(7:10, function(i,data_pg,data_mg,country,admin){
  run_pmcmc_mg(data_raw_pg = data_pg[[i]],
               data_raw_mg = data_mg[[i]],
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
},data_pg=nnp_pg_list,data_mg=nnp_mg_list,country=country,admin=admin)
nnp_ng_mgcorr_bulk_seas_single$status() #'antimonarchal_flee' submitted 2 Mar 10:40am
nnp_ng_mgcorr_bulk_seas_single <- obj$task_bundle_get('antimonarchal_flee')
nnp_ng_mgcorr_bulk_seas_results <- lapply(1:4, function(id){
  nnp_ng_mgcorr_bulk_seas_single$tasks[[id]]$result()
})
obj$task_bundle_list()
nnp_mgcorr_bulk_seas_results_update <- append(nnp_mgcorr_bulk_seas_results,nnp_ng_mgcorr_bulk_seas_results)
