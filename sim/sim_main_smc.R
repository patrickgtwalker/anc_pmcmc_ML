library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library('mcstate')
#library(didehpc)
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



source('shared/utils.R')
source('sim/data_gen.R')
source('shared/data_gen_moz.R')
source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')

##### STANDARD INPUTS #####################
genRandWalk <- function(x,vol,randWalk) {
  if (x == 0)    return (randWalk)
  else return(genRandWalk(x-1,vol,c(randWalk,min(exp(log(randWalk[length(randWalk)])+rnorm(1)*vol),max_EIR))))
}

init_EIR=100
max_EIR=1000
prop_treated = 0.4
init_age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
time= 5*365
EIR_step=30
out_step=0.1
EIR_volatility=0.1
het_brackets<-5
EIR_times=seq(0,time,by=EIR_step)
EIR_vals=genRandWalk(length(EIR_times)-1,EIR_volatility,init_EIR)

### RUN CURRENT EIR STRIPPED MODEL #########
model_file<-"shared/odin_model_stripped_matched.R"
mpl <- model_param_list_create(init_EIR = init_EIR,
                               init_ft = prop_treated,
                               EIR_times=EIR_times,
                               EIR_vals=EIR_vals
)

pars <- equilibrium_init_create_stripped(age_vector = init_age,
                                         init_EIR = init_EIR,
                                         ft = prop_treated,
                                         model_param_list = mpl,
                                         het_brackets = het_brackets)

##The malaria model but only on human side (i.e. no mosquitoes to worry about)
generator <- odin(model_file)
state_use <- pars[names(pars) %in% coef(generator)$name]

# create model with initial values
mod <- generator(user = state_use, use_dde = TRUE)

tt <- seq(0, time, out_step)

# run the simulation to base the data
start.time <- Sys.time()
mod_run <- mod$run(tt, step_max_n = 1e7,
                   atol = 1e-5,
                   rtol = 1e-5)


# shape output
out <- mod$transform_variables(mod_run)

### RUN MODEL WITH SMC ####
weibull <- function(time, alpha = 3.4, beta = 39.34) {
  pow <- -(time/beta)^alpha
  y <- exp(pow)
  return(y)
}
get_smc_profile<-function(start_sim,end_sim,SMC_start,nround,gap){
  prop_prof<-c(rep(1-weibull(1:gap),nround-1),1-weibull(1:100))
  return(list(
    SMC_times=c(start_sim,seq(SMC_start,SMC_start+length(prop_prof)-1,by=1),end_sim),
    SMC_vals=c(1,prop_prof,1)
  ))
}
gap=30
nround=4
plot(1-c(rep(1-weibull(1:gap),nround-1),1-weibull(1:100)))
lines(1-c(rep(1-weibull(1:gap),nround-1),1-weibull(1:100)))

plot(smc_prof$SMC_times,smc_prof$SMC_vals)
smc_prof<-get_smc_profile(0,1800,900,4,30)
smc_cov<-0.9
smc_model_file<-"shared/odin_model_stripped_matched_smc.R"
generator_smc <- odin(smc_model_file)
state_use_smc<-state_use
state_use_smc$SMC_times<-smc_prof$SMC_times
state_use_smc$SMC_vals<-smc_prof$SMC_vals
state_use_smc$smc_cov<-smc_cov
mod_smc <- generator_smc(user = state_use_smc, use_dde = TRUE)
mod_run_smc <- mod_smc$run(tt, step_max_n = 1e7,
                           atol = 1e-5,
                           rtol = 1e-5)
out_smc <- mod_smc$transform_variables(mod_run_smc)
out_smc$de
# plot data and generate data
plot(out$t,out$inc,col="white",ylim=c(0,0.004))
lines(out$t,out$inc,col="blue",lwd=4)
lines(out_smc$t,out_smc$inc,col="red",lwd=4)

plot(out$t,out$incunder5,col="white",ylim=c(0,0.01))
lines(out$t,out$incunder5,col="blue",lwd=4)
lines(out_smc$t,out_smc$incunder5,col="red",lwd=4)
out$
plot(out$t,out$incunder5,col="white",ylim=c(0,0.01),xlim=c(800,1100))
lines(out$t,out$incunder5,col="blue",lwd=4)
abline(v=seq(900,1020,by=30),lty=2)
lines(out_smc$t,out_smc$incunder5,col="red",lwd=4)


### RUN BETA VARYING MODEL
model_file_beta<-"shared/MiP_odin_model_nodelay.R"
mpl_initial <- sifter::model_param_list_create(init_EIR = init_EIR,
                                               init_ft = prop_treated
)

pars_initial <- sifter::equilibrium_init_create_stripped(age_vector = init_age,
                                                         init_EIR = init_EIR,
                                                         ft = prop_treated,
                                                         model_param_list = mpl_initial,
                                                         het_brackets = het_brackets)
init_betaa <- pars_initial$betaa_eq


betaa_times<-seq(0,4*365,by=30)
volatility<-0.1

### just a random walk on logscale
log_betaa_vals <- genRandWalk(length(betaa_times)-1,volatility,log(init_betaa))

betaa_vals <- exp(log_betaa_vals)
##set up the simulation for the simualted data 
time<- 4*365
out_step=1

mpl <- sifter::model_param_list_create(init_EIR = init_EIR,
                                       init_ft = prop_treated,
                                       betaa_times=betaa_times,
                                       betaa_vals=betaa_vals,
                                       lag_rates = 10
)

pars <- sifter::equilibrium_init_create_stripped(age_vector = init_age,
                                                 init_EIR = init_EIR,
                                                 ft = prop_treated,
                                                 model_param_list = mpl,
                                                 het_brackets = het_brackets)
generator_beta <- odin(model_file_beta)
state_use_beta <- pars[names(pars) %in% coef(generator_beta)$name]
# create model with initial values
mod_beta <- generator_beta(user = state_use_beta, use_dde = TRUE)
tt <- seq(0, time, out_step)

# run the simulation to base the data
mod_run_beta <- mod_beta$run(tt, step_max_n = 1e7,
                   atol = 1e-5,
                   rtol = 1e-5)

# shape output
out_beta <- mod_beta$transform_variables(mod_run_beta)
smc_prof
model_file_beta_smc<-"shared/MiP_odin_model_nodelay_smc.R"
generator_beta_smc <- odin(model_file_beta_smc)
state_use_beta_smc<-state_use_beta
state_use_beta_smc$SMC_times<-smc_prof$SMC_times
state_use_beta_smc$SMC_vals<-smc_prof$SMC_vals
state_use_beta_smc$smc_cov<-1
mod_beta_smc <- generator_beta_smc(user = state_use_beta_smc, use_dde = TRUE)
mod_run_beta_smc <- mod_beta_smc$run(tt, step_max_n = 1e7,
                             atol = 1e-5,
                             rtol = 1e-5)
out_beta_smc <- mod_beta_smc$transform_variables(mod_run_beta_smc)
out_beta_smc$clin_inc[1,,,]
# plot data and generate data
 plot(out_beta$t,out_beta_smc$inc/out_beta$inc,type='l',col="red",lwd=4,xlim=c(900,1200))
 lines(out_beta_smc$t,out_beta_smc$inc,col="blue",lwd=4,xlim=c(900,1200),ylim=c(0,0.004))
 
 plot(out_beta$t,out_beta$inc1,type='l',col="red",lwd=4,xlim=c(900,1200),ylim=c(0,0.001))
 lines(out_beta_smc$t,out_beta_smc$inc1,col="blue",lwd=4)
 plot(out_beta$t,out_beta_smc$inc_smc,type='l',col="red",lwd=4,xlim=c(900,1200),ylim=c(0,0.004))
 
 out_beta_smc$cov_vect[1,,]
 