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

source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')


season_model <- odin::odin("shared/odin_model_stripped_seasonal.R")
#Provide age categories, proportion treated, and number of heterogeneity brackets
init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
prop_treated <- 0.4
het_brackets <- 5
country <- 'Burkina Faso'
admin_unit <- 'Cascades'
init_EIR <- 50
EIR_SD <- 1

max_EIR <- 1000
state_check <- 0
lag_rates <- 10
start_stoch <- 1
time_origin <- 1
seasonality_on <- 1

#Create model parameter list. Also loads seasonality profile data file to match to desired admin_unit and country
mpl <- model_param_list_create(init_age = init_age,
                                  pro_treated = prop_treated,
                                  het_brackets = het_brackets,
                                  country = country,
                                  admin_unit = admin_unit,
                               max_EIR = max_EIR,
                               state_check = state_check,
                               lag_rates = lag_rates,
                               start_stoch = start_stoch,
                               time_origin = time_origin,
                               seasonality_on = seasonality_on)

mpl <- append(mpl,list(EIR_SD = EIR_SD))#Only need this because it's expected as output
## Run equilibrium function
state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                          init_EIR = init_EIR,
                                          ft = prop_treated,
                                          model_param_list = mpl,
                                          het_brackets = het_brackets)
##run seasonality model
state_use <- state[names(state) %in% coef(season_model)$name]

# create model with initial values
mod <- season_model$new(user = state_use, use_dde = TRUE)

# tt <- c(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")))
tt <- c(0:(10*365))

# run seasonality model
mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)

# shape output
out <- mod$transform_variables(mod_run)
# windows(10,8)

plot(out$t,out$prev,type='l')
plot(out$prev ~ out$t, t = "l", ylim = c(0, 0.8), xlab = "Year", ylab = "Prevalence")

true_val <- data.frame(t=out$t,
                       prev=out$prev,
                       inc05=out$inc05)%>%
  mutate(date = as.Date('2010-01-01')+t,
         month = as.yearmon(date))
saveRDS(true_val,'sim/true_val.rds')
avg_true_prev <- true_val %>%
  group_by(month)%>%
  summarise(prev=mean(prev))%>%
  mutate(date = as.Date(as.yearmon(month), frac = 0.5))
saveRDS(avg_true_prev,'sim/avg_true_prev.rds')
plot(avg_true_prev$prev ~ avg_true_prev$date, t = "l", ylim = c(0, 0.8), xlab = "Year", ylab = "Prevalence")

##Simulate 5 years of data from this true run
obs_times <- data.frame(month=as.yearmon(seq(as.Date('2015-01-01'),by='month',length.out = 12*5)))

sim_obs <- merge(avg_true_prev[,c('month','date','prev')],obs_times,by='month')

sim_obs$tested <- round(rnorm(nrow(sim_obs),50,10))
sim_obs$positive <- rbinom(nrow(sim_obs),sim_obs$tested,sim_obs$prev)
saveRDS(sim_obs,'sim/sim_obs.rds')

plot(true_val$prev ~ true_val$date, t = "l", ylim = c(0, 0.8), xlab = "Year", ylab = "Prevalence")
points(sim_obs$date, sim_obs$positive/sim_obs$tested, pch = 19, col = "darkred")

sim_obs_peak <- sim_obs[sim_obs$month>=as.yearmon('Oct 2015'),]
saveRDS(sim_obs_peak,'sim/sim_obs_peak.rds')
plot(true_val$prev ~ true_val$date, t = "l", ylim = c(0, 0.8), xlab = "Year", ylab = "Prevalence")
points(sim_obs_peak$date, sim_obs_peak$positive/sim_obs_peak$tested, pch = 19, col = "darkred")

sim_obs_trough <- sim_obs[sim_obs$month>=as.yearmon('Apr 2015'),]
saveRDS(sim_obs_trough,'sim/sim_obs_trough.rds')
plot(true_val$prev ~ true_val$date, t = "l", ylim = c(0, 0.8), xlab = "Year", ylab = "Prevalence")
points(sim_obs_trough$date, sim_obs_trough$positive/sim_obs_trough$tested, pch = 19, col = "darkred")

test_run_peak_std <- run_pmcmc(data = sim_obs_peak,
                          n_particles = 10,
                          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                          max_EIR=1000,
                          max_steps = 1e7,
                          atol = 1e-5,
                          rtol = 1e-6,
                          n_steps = 5,
                          n_threads = 2,
                          lag_rates = 10,
                          country = 'Burkina Faso',
                          admin_unit = 'Cascades',
                          seasonality_on = 0,
                          state_check = 0)
plot_particle_filter(test_run_peak_std$history,true_history=sim_obs_peak,times=sim_obs_peak$month)

test_run_peak_seas <- run_pmcmc(data = sim_obs_peak,
                               n_particles = 10,
                               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                               max_EIR=1000,
                               max_steps = 1e7,
                               atol = 1e-5,
                               rtol = 1e-6,
                               n_steps = 10,
                               n_threads = 2,
                               lag_rates = 10,
                               country = 'Burkina Faso',
                               admin_unit = 'Cascades',
                               preyears = 5,
                               seasonality_on = 1,
                               state_check = 0)
plot_particle_filter(test_run_peak_seas$history,true_history=sim_obs_peak,times=sim_obs_peak$month)

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts.temp", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust@8aef08d",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "24Core",cores =8, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj$cluster_load(TRUE)
obj$login()
obj$config
test <- obj$enqueue(run_pmcmc(data = sim_obs_peak,
                               n_particles = 10,
                               proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                               max_EIR=1000,
                               max_steps = 1e7,
                               atol = 1e-5,
                               rtol = 1e-6,
                               n_steps = 10,
                               n_threads = 8,
                               lag_rates = 10,
                               country = 'Burkina Faso',
                               admin_unit = 'Cascades',
                               seasonality_on = 1,
                               state_check = 0))
test$status()

sim_peak_std <- obj$enqueue(run_pmcmc(data = sim_obs_peak,
                              n_particles = 200,
                              proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                              max_EIR=1000,
                              max_steps = 1e7,
                              atol = 1e-5,
                              rtol = 1e-6,
                              n_steps = 1000,
                              n_threads = 8,
                              lag_rates = 10,
                              country = 'Burkina Faso',
                              admin_unit = 'Cascades',
                              seasonality_on = 0,
                              state_check = 0))
sim_peak_std$status()
sim_peak_std$id #"59287feb9176d06c8320be2948156166"
sim_peak_std$log()
sim_trough_std <- obj$enqueue(run_pmcmc(data = sim_obs_trough,
                                      n_particles = 200,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 8,
                                      lag_rates = 10,
                                      country = 'Burkina Faso',
                                      admin_unit = 'Cascades',
                                      seasonality_on = 0,
                                      state_check = 0))
sim_trough_std$status()
sim_trough_std$id #"cbad64d6b0bfa082596da5004accef01"
sim_peak_seas <- obj$enqueue(run_pmcmc(data = sim_obs_peak,
                                      n_particles = 200,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 8,
                                      lag_rates = 10,
                                      country = 'Burkina Faso',
                                      admin_unit = 'Cascades',
                                      seasonality_on = 1,
                                      state_check = 0))
sim_peak_seas$status()
sim_peak_seas$id #"1d84e3d602faeae2c4866fcf97b9eaf7"
sim_trough_seas <- obj$enqueue(run_pmcmc(data = sim_obs_trough,
                                        n_particles = 200,
                                        proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                        max_EIR=1000,
                                        max_steps = 1e7,
                                        atol = 1e-5,
                                        rtol = 1e-6,
                                        n_steps = 1000,
                                        n_threads = 8,
                                        lag_rates = 10,
                                        country = 'Burkina Faso',
                                        admin_unit = 'Cascades',
                                        seasonality_on = 1,
                                        state_check = 0))
sim_trough_seas$status()
sim_trough_seas$id #"9b958d300ccf3fe95701769fd33bb5dc"
