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
library(reshape2)

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts.temp", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust@8aef08d",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "32Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj$login()
obj$cluster_load(TRUE)
names()

test_tanz_submit <- obj$enqueue_bulk(1:2, function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 10,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 100,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_data_list_14to17)

test_tanz_submit$status()
test_tanz_submit$tasks$c788477447acd71ee8be7aba3a6ba47c$log()

tanz_lt20_2015to2017_submit <- obj$enqueue_bulk(1:26, function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_data_list_15to17)
tanz_lt20_2015to2017_submit$status()#'nonspiritous_lobo'
tanz_lt20_2015to2017_submit <- obj$task_bundle_get('nonspiritous_lobo')
test_result <- tanz_lt20_2015to2017_submit$tasks$`0519fdb81fdca79d98dcc1fa2b19cd3e`$result()

tanz_lt20_2015to2017_results <- lapply(1:26, function(id){
  tanz_lt20_2015to2017_submit$tasks[[id]]$result()
})

names(tanz_lt20_2015to2017_results) <- gsub(' Region','',names(tanz_data_list_15to17$lt20))


##Submit District level Lake Malawi
obj$login()
obj$cluster_load(TRUE)
obj$config

tanz_lakemal_lt20_2015to2017_submit <- obj$enqueue_bulk(1:14, function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_lakemal_list_15to17)

tanz_lakemal_lt20_2015to2017_submit$status() #'neutral_asianelephant'
tanz_lakemal_lt20_2015to2017_submit <- obj$task_bundle_get('neutral_asianelephant')
tanz_lakemal_lt20_2015to2017_submit$tasks$`754d1d793898e799ba97f07b42a06f73`$log()
tanz_lakemal_lt20_2015to2017_submit$times()
as.vector(tanz_lakemal_lt20_2015to2017_submit$status())

##Re-run tasks 8 and 13
obj$login()
obj$cluster_load(TRUE)
tanz_lakemal_lt20_2015to2017_submit_2 <- obj$enqueue_bulk(c(8,11), function(x,obs_list) {
  print(names(obs_list$lt20[x]))
  run_pmcmc(data = obs_list$lt20[[x]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 4,
            lag_rates = 10,
            preyears = 5,
            seasonality_on = 0,
            state_check = 0,
            seasonality_check = 0)},obs_list=tanz_lakemal_list_15to17)
obj$task_bundle_list()
bundle_info <- obj$task_bundle_info()
tanz_lakemal_lt20_2015to2017_submit_2$status()
tanz_lakemal_lt20_2015to2017_submit_2 <- obj$task_bundle_get('spiritless_germanshepherd')
tanz_lakemal_lt20_2015to2017_submit$tasks[[8]] <- tanz_lakemal_lt20_2015to2017_submit_2$tasks[[1]]
tanz_lakemal_lt20_2015to2017_submit$tasks[[8]]$log()
tanz_lakemal_lt20_2015to2017_submit$tasks[[11]] <- tanz_lakemal_lt20_2015to2017_submit_2$tasks[[2]]
tanz_lakemal_lt20_2015to2017_submit$status()
tanz_lakemal_lt20_2015to2017_results <- lapply(1:14, function(id){
  tanz_lakemal_lt20_2015to2017_submit$tasks[[id]]$result()
})
names(tanz_lakemal_lt20_2015to2017_results) <- names(tanz_lakemal_list_15to17$lt20)

tanz_lakemal_lt20_2015to2017_results$`Ludewa District Council`$mcmc$log_prior
