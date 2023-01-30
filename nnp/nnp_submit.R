remotes::install_github("mrc-ide/mcstate", upgrade = TRUE)
remotes::install_github("mrc-ide/dust", upgrade = TRUE)
remotes::install_github("mrc-ide/mode", upgrade = TRUE)
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
test_run_std <- run_pmcmc(data = data_raw_bf_pg_banfora,
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
                      seasonality_on = 1,
                      state_check = 0)
plot_particle_filter(test_run_std$history,true_history=data_raw_bf_pg_banfora,times=data_raw_bf_pg_banfora$t)


##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
sources_old <- c("test_scripts/oldversions/run_pmcmc.R",
             "test_scripts/oldversions/model_parameters.R",
             "test_scripts/oldversions/equilibrium-init-create-stripped.R",
             'test_scripts/oldversions/utils.R')

ctx <- context::context_save("T:/jth/contexts.new", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))

ctx <- context::context_save("T:/jth/contexts.temp", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust@8aef08d",'mrc-ide/mcstate')))

ctx_old <- context::context_save("T:/jth/contexts", sources = sources_old,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','odin'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

config_16 <- didehpc::didehpc_config(cores = 16, parallel = TRUE)
obj_16 <- didehpc::queue_didehpc(ctx,config = config_16)
obj_16$login()
obj_16$cluster_load(TRUE)

config_16 <- didehpc::didehpc_config(cores = 16, parallel = TRUE)
obj_16 <- didehpc::queue_didehpc(ctx,config = config_16)
obj_16$login()
obj_16$cluster_load(TRUE)

config_1 <- didehpc::didehpc_config(template = "8Core",cores =8, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj$enqueue(sessionInfo())
check <- obj$task_get("6efd57d099127ba7761ce60902008e19")
check2 <- obj$enqueue(sessionInfo())
check$status()
check2$status()
check2$result()
obj$login()
obj$config
obj$cluster_load(TRUE)

config_12 <- didehpc::didehpc_config(cores = 12, parallel = TRUE)
obj_12 <- didehpc::queue_didehpc(ctx,config = config_12)
obj_12$login()
obj_12$config
obj_12$cluster_load(TRUE)

##Fit to primigrav data to <5 yo
nnp_pg_bulk_test <- obj_12$enqueue_bulk(1:10, function(i,data_site){
  run_pmcmc(data = data_site[[i]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 10,
            n_threads = 1,
            lag_rates = 10,
            country = 'Burkina Faso',
            admin_unit = 'Cascades',
            seasonality_on = 0,
            state_check = 0)
},data_site=nnp_pg_list)

nnp_pg_bulk_test$status()
obj_1$unsubmit(nnp_pg_bulk_test$ids)
nnp_pg_bulk_test
nnp_pg_bulk_test$times()

nnp_pg_bulk_seas_test <- obj_1$enqueue_bulk(1:10, function(i,data_site){
  run_pmcmc(data = data_site[[i]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 10,
            n_threads = 1,
            lag_rates = 10,
            country = 'Burkina Faso',
            admin_unit = 'Cascades',
            seasonality_on = 1,
            state_check = 0)
},data_site=nnp_pg_list)

nnp_pg_bulk_seas_minitest <- obj$enqueue(run_pmcmc(data = nnp_pg_list[[1]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 10,
            n_threads = 1,
            lag_rates = 10,
            country = 'Burkina Faso',
            admin_unit = 'Cascades',
            seasonality_on = 1,
            state_check = 0))
obj$login()
nnp_pg_bulk_seas_minitest$status()
nnp_pg_bulk_seas_minitest$log()
nnp_pg_bulk_seas_minitest$times()

test_old <- obj$enqueue(run_pmcmc(data = nnp_pg_list[[1]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 16))
test_old$status()
test_old$log()
check_vers <- function(x){return(list(x,packageVersion(x)))}
packages <- c('odin.dust','mode','dust')
lapply(packages,check_vers)
check_versions <- obj$enqueue(lapply(packages,check_vers))
check_versions$status()
check_versions$log()

check_versions$result()
obj$login()
nnp_pg_bulk_seas_test$status()
nnp_pg_bulk_seas_test$tasks$`443bf11594b11993ff1944593cc0aa7b`$times()
nnp_pg_bulk_seas_test$tasks$`443bf11594b11993ff1944593cc0aa7b`$log()
obj_1$unsubmit(nnp_pg_bulk_seas_test$ids)

nnp_pg_result_list <- lapply(1:10, function(id){
  nnp_pg_bulk$tasks[[id]]$result()
})

##Submit bulk non-seasonality runs
nnp_pg_bulk_std_mzbf <- obj$enqueue_bulk(1:6, function(i,data_site){
  run_pmcmc(data = data_site[[i]],
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
},data_site=nnp_pg_list)
##'gothic_halibut'
nnp_pg_bulk_std_mzbf$status()
nnp_pg_bulk_std_mzbf$times()
std_mzbf_result_list <- lapply(1:6, function(id){
  nnp_pg_bulk_std_mzbf$tasks[[id]]$result()
})

##Submit MZ and BF seasonal model
obj$login()
obj$cluster_load(TRUE)
nnp_pg_bulk_seas_mzbf <- obj$enqueue_bulk(1:6, function(i,data_site,country,admin){
  run_pmcmc(data = data_site[[i]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 8,
            lag_rates = 10,
            country = country[6],
            admin_unit = admin[6],
            seasonality_on = 1,
            state_check = 0)
},data_site=nnp_pg_list,country=country,admin=admin)
#'highhanded_wildcat'
nnp_pg_bulk_seas_mzbf$status()
nnp_pg_bulk_seas_mzbf$times()
seas_mzbf_result_list <- lapply(1:6, function(id){
  nnp_pg_bulk_seas_mzbf$tasks[[id]]$result()
})
