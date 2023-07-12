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
source('nnp/in_development/run_pmcmc_pg.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')

##Get saved processed data
WA_pg_data_list <- readRDS('trial/Data/WA_pg_data_list.rds')
WA_sg_data_list <- readRDS('trial/Data/WA_sg_data_list.rds')
WA_all_data_list <- readRDS('trial/Data/WA_all_data_list.rds')

names(WA_pg_data_list)
admins <- c(`Burkina Faso` = 'Plateau-Central',
            Gambia = 'Upper River',
            Ghana = 'Upper East',
            Mali = 'Koulikoro') ##There were 3 Mali sites in the trial. Chose the 
                                ##admin area in the middle
##test
test_run_pg_wa <- run_pmcmc_pg(data_raw = WA_pg_data_list[[1]],
                               n_particles = 10,
                               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                               max_EIR=1000,
                               max_steps = 1e7,
                               atol = 1e-5,
                               rtol = 1e-6,
                               n_steps = 20,
                               n_threads = 2,
                               lag_rates = 10,
                               country = names(WA_pg_data_list[1]),
                               admin_unit = admins[[names(WA_pg_data_list[1])]],
                               seasonality_on = 1,
                               state_check = 0)
windows(10,7)
plot_particle_filter(test_run_pg_wa$history,true_history=WA_pg_data_list[[1]],times=c(1:nrow(WA_pg_data_list[[1]])))


##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "nnp/in_development/run_pmcmc_pg.R",
             "nnp/in_development/run_pmcmc_mg.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
packages <- list(loaded = 'plyr',attached = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'))
ctx <- context::context_save("T:/jth/contexts_may23", sources = sources,
                             packages = packages,
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "32Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj$login()
obj$cluster_load(TRUE)

wa_pg_bulk_seas <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
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
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_seas$status()#'didactic_rainbowfish'
wa_pg_bulk_seas$tasks$`4c0c66f66ddac58dbeb57dc8a754cff5`$log()
obj$unsubmit(wa_pg_bulk_seas$ids)
wa_pg_bulk_seas$tasks$f7ca68b6142eb5a7b89d8effea87e586$log()
wa_pg_bulk_seas$tasks$`821caff3e86e94b8eb66e2c3a1348a41`$log()
wa_pg_seas_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_seas$tasks[[id]]$result()
})

wa_pg_bulk_seas_2 <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
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
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_seas_2 <- obj$task_bundle_get('elder_oryx')
wa_pg_bulk_seas_2$status()#'elder_oryx'
wa_pg_bulk_seas_2$tasks[[3]]$log()

wa_pg_seas_2_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_seas_2$tasks[[id]]$result()
})

admins[[names(WA_pg_data_list[1])]]

wa_pg_bulk_std <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
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
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_std <- obj$task_bundle_get('defunct_wryneck')
wa_pg_bulk_std$status()#''defunct_wryneck''
wa_pg_bulk_std$tasks[[1]]$log()

wa_pg_std_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_std$tasks[[id]]$result()
})


wa_pg_bulk_seas_padded <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
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
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L,
               start_pf_time = 30*4)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_seas_padded<-obj$task_bundle_get('unhumane_hammerheadshark')
wa_pg_bulk_seas_padded$status()#'unhumane_hammerheadshark'
#obj$unsubmit(wa_pg_bulk_seas_padded$ids)
wa_pg_seas_padded_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_seas_padded$tasks[[id]]$result()
})

wa_pg_bulk_seas_padded$tasks[[4]]$log()

data_gamb_all <- rbind(WA_pg_data_list$Gambia,WA_sg_data_list$Gambia)%>%
  group_by(month)%>%
  summarise(tested = sum(tested),
            positive = sum(positive))

wa_pg_gamb_seas <- obj$enqueue(run_pmcmc_pg(data_raw = data_gamb_all,
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
               country = 'Gambia',
               admin_unit = 'Upper River',
               state_check = 0,
               seed = 2L,
               start_pf_time = 30))
wa_pg_gamb_seas$id #"4cd1aad5848796687a13e1234594c54f"
wa_pg_gamb_seas <- obj$task_get('4cd1aad5848796687a13e1234594c54f')
wa_pg_gamb_seas$status()
wa_pg_gamb_seas$log()

data_bf_all <- rbind(WA_pg_data_list$`Burkina Faso`,WA_sg_data_list$`Burkina Faso`)%>%
  group_by(month)%>%
  summarise(tested = sum(tested),
            positive = sum(positive))

wa_pg_bfall_seas <- obj$enqueue(run_pmcmc_pg(data_raw = data_bf_all,
                                            n_particles = 200,
                                            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                            max_EIR=1000,
                                            max_steps = 1e7,
                                            atol = 1e-5,
                                            rtol = 1e-6,
                                            n_steps = 250,
                                            n_threads = 4,
                                            lag_rates = 10,
                                            seasonality_on = 1,
                                            seasonality_check = 1,
                                            country = 'Burkina Faso',
                                            admin_unit = 'Plateau-Central',
                                            state_check = 0,
                                            seed = 2L,
                                            start_pf_time = 30*4))
wa_pg_bfall_seas$status()
wa_pg_bfall_seas$log()

wa_pg_bf_seas <- obj$enqueue(run_pmcmc_pg(data_raw = WA_pg_data_list[[1]],
                                             n_particles = 200,
                                             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                             max_EIR=1000,
                                             max_steps = 1e7,
                                             atol = 1e-5,
                                             rtol = 1e-6,
                                             n_steps = 500,
                                             n_threads = 4,
                                             lag_rates = 10,
                                             seasonality_on = 1,
                                             seasonality_check = 1,
                                             country = 'Burkina Faso',
                                             admin_unit = 'Plateau-Central',
                                             state_check = 0,
                                             seed = 2L,
                                             start_pf_time = 30*4))
wa_pg_bf_seas$status()
wa_pg_bf_seas$times(unit_elapsed = 'mins')
wa_pg_bf_seas$log()

wa_pg_bfall_seas_16 <- obj$enqueue(run_pmcmc_pg(data_raw = data_bf_all,
                                             n_particles = 200,
                                             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                             max_EIR=1000,
                                             max_steps = 1e7,
                                             atol = 1e-5,
                                             rtol = 1e-6,
                                             n_steps = 250,
                                             n_threads = 16,
                                             lag_rates = 10,
                                             seasonality_on = 1,
                                             seasonality_check = 1,
                                             country = 'Burkina Faso',
                                             admin_unit = 'Plateau-Central',
                                             state_check = 0,
                                             seed = 2L,
                                             start_pf_time = 30*4))
wa_pg_bfall_seas_16$status()
wa_pg_bfall_seas_16$times(unit_elapsed = 'mins')
wa_pg_bfall_seas_16$log()

wa_pg_bf_seas_16 <- obj$enqueue(run_pmcmc_pg(data_raw = WA_pg_data_list[[1]],
                                          n_particles = 200,
                                          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                          max_EIR=1000,
                                          max_steps = 1e7,
                                          atol = 1e-5,
                                          rtol = 1e-6,
                                          n_steps = 500,
                                          n_threads = 16,
                                          lag_rates = 10,
                                          seasonality_on = 1,
                                          seasonality_check = 1,
                                          country = 'Burkina Faso',
                                          admin_unit = 'Plateau-Central',
                                          state_check = 0,
                                          seed = 2L,
                                          start_pf_time = 30*4))
wa_pg_bf_seas$status()
wa_pg_bf_seas$log()

wa_pg_bfall_seas_pc_2 <- run_pmcmc_pg(data_raw = data_bf_all,
                                             n_particles = 50,
                                             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                             max_EIR=1000,
                                             max_steps = 1e7,
                                             atol = 1e-5,
                                             rtol = 1e-6,
                                             n_steps = 20,
                                             n_threads = 2,
                                             lag_rates = 10,
                                             seasonality_on = 1,
                                             seasonality_check = 1,
                                             country = 'Burkina Faso',
                                             admin_unit = 'Plateau-Central',
                                             state_check = 0,
                                             seed = 2L,
                                             start_pf_time = 30*4)
wa_pg_bfall_seas_pc_2$history[1,,]
wa_pg_bfall_seas_pc$seas_history
wa_pg_bfall_seas_pc_2$history[1,,]

wa_pg_bf_seas_16_u5prev <- obj$enqueue(run_pmcmc(data_raw = WA_pg_data_list[[1]],
                                             n_particles = 200,
                                             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                             max_EIR=1000,
                                             max_steps = 1e7,
                                             atol = 1e-5,
                                             rtol = 1e-6,
                                             n_steps = 500,
                                             n_threads = 16,
                                             lag_rates = 10,
                                             seasonality_on = 1,
                                             seasonality_check = 0,
                                             country = 'Burkina Faso',
                                             admin_unit = 'Plateau-Central',
                                             state_check = 0,
                                             seed = 2L))
wa_pg_bf_seas_16_u5prev$status()
wa_pg_bf_seas_16_u5prev$log()

wa_pg_bfall_seas_16_u5prev <- obj$enqueue(run_pmcmc(data_raw = data_bf_all,
                                                 n_particles = 200,
                                                 proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                 max_EIR=1000,
                                                 max_steps = 1e7,
                                                 atol = 1e-5,
                                                 rtol = 1e-6,
                                                 n_steps = 500,
                                                 n_threads = 16,
                                                 lag_rates = 10,
                                                 seasonality_on = 1,
                                                 seasonality_check = 0,
                                                 country = 'Burkina Faso',
                                                 admin_unit = 'Plateau-Central',
                                                 state_check = 0,
                                                 seed = 2L))
wa_pg_bfall_seas_16_u5prev$status()
wa_pg_bfall_seas_16_u5prev$log()

data_bf_all_standardized <- data_bf_all%>%
  rename(positive_old = positive,
         tested_old = tested)%>%
  mutate(prev=positive_old/tested_old,
         tested = 200,
         positive = round(prev*tested))
data_bf_all_standardized_cis <- addCIs(df=data_bf_all_standardized,Ys=data_bf_all_standardized$positive,data_bf_all_standardized$tested)
wa_pg_bfall_seas_16_sampsz <- obj$enqueue(run_pmcmc_pg(data_raw = data_bf_all_standardized,
                                                       n_particles = 200,
                                                       proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                       max_EIR=1000,
                                                       max_steps = 1e7,
                                                       atol = 1e-5,
                                                       rtol = 1e-6,
                                                       n_steps = 250,
                                                       n_threads = 16,
                                                       lag_rates = 10,
                                                       seasonality_on = 1,
                                                       seasonality_check = 0,
                                                       country = 'Burkina Faso',
                                                       admin_unit = 'Plateau-Central',
                                                       state_check = 0,
                                                       seed = 2L,
                                                       start_pf_time = 30*4))
wa_pg_bfall_seas_16_sampsz$status()
wa_pg_bfall_seas_16_sampsz$log()
#obj$unsubmit(wa_pg_bfall_seas_16_sampsz$id)
obj$login()
##Check seasonality equilibrium
check_seas_pg_wa <- run_pmcmc_pg(data_raw = WA_pg_data_list[[1]],
                               n_particles = 10,
                               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                               max_EIR=1000,
                               max_steps = 1e7,
                               atol = 1e-5,
                               rtol = 1e-6,
                               n_steps = 20,
                               n_threads = 2,
                               lag_rates = 10,
                               country = names(WA_pg_data_list[1]),
                               admin_unit = admins[[names(WA_pg_data_list[1])]],
                               seasonality_on = 1,
                               seasonality_check = 1)
check_seas_pg_wa$seas_history


####
wa_pg_bulk_seas_090623 <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
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
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L,
               start_pf_time = 30*4)
},data_pg=WA_pg_data_list,admins=admins)
wa_pg_bulk_seas_090623$status() #''disastrous_andeancondor''
wa_pg_bulk_seas_090623 <- obj$task_bundle_get('disastrous_andeancondor')

wa_pg_bulk_seas_090623$tasks[[1]]$log()

wa_pg_bulk_seas_090623_result_list <- lapply(1:4, function(id){
  wa_pg_bulk_seas_090623$tasks[[id]]$result()
})
wa_pg_bulk_seas_090623_result_list[[1]]$history[1,1,1]
wa_all_bulk_seas_090623 <- obj$enqueue_bulk(1:4, function(i,data_pg,admins){
  run_pmcmc_pg(data_raw = data_pg[[i]],
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
               country = names(data_pg[i]),
               admin_unit = admins[[names(data_pg[i])]],
               state_check = 0,
               seed = 2L,
               start_pf_time = 30*4)
},data_pg=WA_all_data_list,admins=admins)
wa_all_bulk_seas_090623 <- obj$task_bundle_get('citric_snowmonkey')
wa_all_bulk_seas_090623$status() #'citric_snowmonkey'

wa_all_bulk_seas_090623_result_list <- lapply(1:4, function(id){
  wa_all_bulk_seas_090623$tasks[[id]]$result()
})
