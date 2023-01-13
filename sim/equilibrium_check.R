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
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)

##Set default theme
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA)))

#Required functions
source('sim/data_gen.R')
source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')

##Generate simulated data##
windows(10,8)
data_sim_comptest <- data_gen(EIR_volatility = 0.8, init_EIR = 20)
plot(data_sim_comptest$EIR_true)

##Three previously run simulated data sets are saved in the folder
## 'anc_pmcmc/sim/sim_datasets
data_sim_comptest3 <- readRDS('sim/sim_datasets/data_sim3.RDS')

##Test run_pmcmc function##
test_run <- run_pmcmc(data = data_sim_comptest3,
                      n_particles = 10,
                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                      max_EIR=1000,
                      max_steps = 1e7,
                      atol = 1e-5,
                      rtol = 1e-6,
                      n_steps = 10,
                      n_threads = 2,
                      state_check = 0)
plot_particle_filter(test_run$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)
matplot(data_sim_comptest3$t, t(test_run$history['Dout', , -1]), type = "l",
        xlab = "Time", ylab = "Dout")
