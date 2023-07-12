library(foreign)
library(sas7bdat)
library(tidyverse)
library(readxl)
library(devtools)
devtools::install_github("mrc-ide/umbrella")
library(umbrella)
library(ggplot2)
library(lubridate)
library(zoo)
library(terra)
library(stringr)
library(geofacet)
library(reshape2)

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))
#See tanz_read_data.R for how to pull and organize Tanga data
test_tanga <- sifter::run_pmcmc(data_raw = tanga_data_list_15to22[[1]], 
                                n_particles = 10,
                                proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                max_EIR=1000,
                                max_steps = 1e7,
                                atol = 1e-5,
                                rtol = 1e-6,
                                n_steps = 3,
                                n_threads = 2,
                                lag_rates = 10,
                                country = 'Tanzania',
                                admin_unit = 'Tanga',
                                seasonality_on = 1,
                                state_check = 0,
                                seasonality_check = 0,
                                preyears = 5,
                                start_pf_time = 30*4
)

##Configure cluster settings
ctx_sifter <- context::context_save("T:/jth/contexts.sift",
                                    packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde','RecordLinkage'),
                                    package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust",'mrc-ide/mcstate','jt-hicks/sifter@issue-5')))
config_1 <- didehpc::didehpc_config(template = "24Core",cores =6, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
# config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj_sift <- didehpc::queue_didehpc(ctx_sifter,config = config_1)
obj_sift$cluster_load(TRUE)
obj_sift$login()

##Test Longer mosquito emergence run on cluster
test_run_mozemerge <- obj_sift$enqueue(sifter::run_pmcmc(data_raw = sifter::data_sim, #I've added data_sim to the package for an easy test
                                     n_particles = 200,
                                     proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                     max_EIR=1000,
                                     max_steps = 1e7,
                                     atol = 1e-6,
                                     rtol = 1e-6,
                                     n_steps = 100,
                                     n_threads = 2,
                                     lag_rates = 10,
                                     country = 'Burkina Faso',
                                     admin_unit = 'Cascades',
                                     seasonality_on = 1,
                                     state_check = 0,
                                     seasonality_check = 0,
                                     stoch_param = 'betaa'))
test_run_mozemerge$id #"032ce779172b0b58f4d523ba8fc7e55c"
test_run_mozemerge$status()
test_run_mozemerge_result <- test_run_mozemerge$result()
plot(as.mcmc(test_run_mozemerge_result$mcmc))

test_run_mozemerge_longer <- obj_sift$enqueue(sifter::run_pmcmc(data_raw = sifter::data_sim, #I've added data_sim to the package for an easy test
                                                         n_particles = 200,
                                                         proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                         max_EIR=1000,
                                                         max_steps = 1e7,
                                                         atol = 1e-6,
                                                         rtol = 1e-6,
                                                         n_steps = 1000,
                                                         n_threads = 6,
                                                         lag_rates = 10,
                                                         # country = 'Burkina Faso',
                                                         # admin_unit = 'Cascades',
                                                         seasonality_on = 0,
                                                         state_check = 0,
                                                         seasonality_check = 0,
                                                         stoch_param = 'betaa'))
test_run_mozemerge_longer$id #"28bac243ae79ce6a8913605ab3762678"
test_run_mozemerge_longer$status()
