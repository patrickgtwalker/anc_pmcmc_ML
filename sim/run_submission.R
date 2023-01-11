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

##Generate simulated data##
source('data_gen.R')
source('model_parameters.R')
source('equilibrium-init-create-stripped.R')

data_sim_comptest <- data_gen(EIR_volatility = 0.8, init_EIR = 20)
keep <- data_sim_comptest
keep2 <- data_sim_comptest
keep3 <- data_sim_comptest
plot(keep$EIR_true)
plot(keep2$EIR_true)
plot(keep3$EIR_true)

saveRDS(keep,'sim_datasets/data_sim1.RDS')
saveRDS(keep2,'sim_datasets/data_sim2.RDS')
saveRDS(keep3,'sim_datasets/data_sim3.RDS')

data_sim_comptest3 <- readRDS('sim_datasets/data_sim3.RDS')

##Test run_pmcmc function##
source('run_pmcmc.R')
test_run <- run_pmcmc(data = data_sim_comptest3,
                      n_particles = 10,
                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                      max_EIR=1000,
                      max_steps = 1e7,
                      atol = 1e-5,
                      rtol = 1e-6,
                      n_steps = 10,
                      n_threads = 2)
plot_particle_filter(test_run$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("run_pmcmc.R",
             "model_parameters.R","equilibrium-init-create-stripped.R")


ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))
##Set up and run scenarios on 64 node core
config_64 <- didehpc::didehpc_config(cores = 64, parallel = TRUE)
obj_64 <- didehpc::queue_didehpc(ctx,config = config_64)

obj_64$cluster_load(TRUE)
obj_64$config
obj_64$login()

run_64_10_test <- obj_64$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 10,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 10,
                                      n_threads = 64))
run_64_10_test$status()
run_64_10_test$log()

run_64_10 <- obj_64$enqueue(run_pmcmc(data = data_sim_comptest3,
                                           n_particles = 10,
                                           proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                           max_EIR=1000,
                                           max_steps = 1e7,
                                           atol = 1e-5,
                                           rtol = 1e-6,
                                           n_steps = 1000,
                                           n_threads = 64))
run_64_10$status()
run_64_10$log()

run_64_50 <- obj_64$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 50,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 64))
run_64_50$status()
run_64_50$log()
# obj_64$unsubmit(c(run_64_50$id,run_64_200$id))

run_64_200 <- obj_64$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 200,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 64))
run_64_200$status()
run_64_200$log()

##Set up and run on 32 node core##
config_32 <- didehpc::didehpc_config(cores = 32, parallel = TRUE)
obj_32 <- didehpc::queue_didehpc(ctx,config = config_32)
obj_32$unsubmit(c(run_32_10$id,run_32_50$id,run_32_200$id))

obj_32$cluster_load(TRUE)
obj_32$config
obj_32$login()

run_32_10 <- obj_32$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 10,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 32))
run_32_10$status()
run_32_10$log()

run_32_50 <- obj_32$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 50,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 32))
run_32_50$status()
run_32_50$log()

run_32_200 <- obj_32$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 200,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 32))
run_32_200$status()
run_32_200$log()

##Set up and run on 4 node core##
config_4 <- didehpc::didehpc_config(cores = 4, parallel = TRUE)
obj_4 <- didehpc::queue_didehpc(ctx,config = config_4)
obj_4$unsubmit(c(run_4_10$id,run_4_50$id,run_4_200$id))

obj_4$cluster_load(TRUE)
obj_4$config
obj_4$login()

run_4_10 <- obj_4$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 10,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 4))
run_4_10$status()
run_4_10$log()

run_4_50 <- obj_4$enqueue(run_pmcmc(data = data_sim_comptest3,
                                    n_particles = 50,
                                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                    max_EIR=1000,
                                    max_steps = 1e7,
                                    atol = 1e-5,
                                    rtol = 1e-6,
                                    n_steps = 1000,
                                    n_threads = 4))
run_4_50$status()
run_4_50$log()

run_4_200 <- obj_4$enqueue(run_pmcmc(data = data_sim_comptest3,
                                    n_particles = 200,
                                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                    max_EIR=1000,
                                    max_steps = 1e7,
                                    atol = 1e-5,
                                    rtol = 1e-6,
                                    n_steps = 1000,
                                    n_threads = 4))
run_4_200$status()
run_4_200$log()

##Set up and run on 1 node ##
config_1 <- didehpc::didehpc_config(cores = 1, parallel = TRUE)
obj_1 <- didehpc::queue_didehpc(ctx,config = config_1)
obj_1$unsubmit(c(run_1_10$id,run_1_50$id,run_1_200$id))

obj_1$cluster_load(TRUE)
obj_1$config
obj_1$login()

run_1_10 <- obj_1$enqueue(run_pmcmc(data = data_sim_comptest3,
                                    n_particles = 10,
                                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                    max_EIR=1000,
                                    max_steps = 1e7,
                                    atol = 1e-5,
                                    rtol = 1e-6,
                                    n_steps = 1000,
                                    n_threads = 1))
run_1_10$status()
run_1_10$log()

run_1_50 <- obj_1$enqueue(run_pmcmc(data = data_sim_comptest3,
                                    n_particles = 50,
                                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                    max_EIR=1000,
                                    max_steps = 1e7,
                                    atol = 1e-5,
                                    rtol = 1e-6,
                                    n_steps = 1000,
                                    n_threads = 1))
run_1_50$status()
run_1_50$log()

run_1_200 <- obj_1$enqueue(run_pmcmc(data = data_sim_comptest3,
                                    n_particles = 200,
                                    proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                    max_EIR=1000,
                                    max_steps = 1e7,
                                    atol = 1e-5,
                                    rtol = 1e-6,
                                    n_steps = 1000,
                                    n_threads = 1))
run_1_200$status()
run_1_200$log()

config_8 <- didehpc::didehpc_config(cores = 8, parallel = TRUE)
obj_8 <- didehpc::queue_didehpc(ctx,config = config_8)

obj_8$cluster_load(TRUE)
obj_8$config
obj_8$login()
maybe <- obj_8$task_get('25b2d62e22e868d6f9793464e412efcb')
maybe$expr()
maybe$log()
data_sim_comptest3_t36 <- data_sim_comptest3[data_sim_comptest3$t<=36*30,]
run_8_comp_test <- obj_8$enqueue(run_pmcmc(data = data_sim_comptest3_t36,
                                           n_particles = 128,
                                           proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                           max_EIR=1000,
                                           max_steps = 1e7,
                                           atol = 1e-5,
                                           rtol = 1e-6,
                                           n_steps = 10000,
                                           n_threads = 8))
run_8_comp_test$status()
run_8_comp_test$log()

##NNP submission
