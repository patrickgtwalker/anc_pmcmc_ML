library("odin")
library("odin.dust")
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
library(lubridate)

source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')
source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')

init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
prop_treated <- 0.4
het_brackets <- 5

max_EIR <- 1000
state_check <- 1
lag_rates <- 10
country <- 'Burkina Faso'
admin_unit <- 'Cascades'
seasonality_on <- 1

mpl_pf <- model_param_list_create(init_age = init_age,
                                  pro_treated = prop_treated,
                                  het_brackets = het_brackets,
                                  max_EIR = max_EIR,
                                  state_check = state_check,
                                  lag_rates = lag_rates,
                                  country = country,
                                  admin_unit = admin_unit,
                                  seasonality_on = seasonality_on)
# print(mpl_pf$state_check)
# print(mpl_pf$ssa0)
##Model with no seasonality
model <- odin("shared/odin_model_stripped_matched.R")
##Model with seasonality component, but set to 1
model_seas <- odin("shared/odin_model_stripped_seasonal.R")


init_EIR <- 10
EIR_vol <- 0.2
mpl <- append(mpl_pf,list(EIR_SD = EIR_vol,
                          init_EIR = init_EIR,
                          EIR_times = c(0,0.5*365),
                          EIR_vals = c(10,10)))

state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                          EIR = init_EIR,
                                          ft = prop_treated,
                                          model_param_list = mpl,
                                          het_brackets = het_brackets,
                                          state_check = mpl$state_check)
##run seasonality model first if seasonality_on == 1
state_use <- state[names(state) %in% coef(model)$name]
state_use_seas <- state[names(state) %in% coef(model_seas)$name]
# print(state_use)
# create model with initial values
mod <- model$new(user = state_use, use_dde = TRUE)
mod_seas <- model_seas$new(user = state_use_seas, use_dde = TRUE)
# print('generated seasonal model')
# tt <- c(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")))
tt <- seq(0, 0.5*365,length.out=500)

# run model
mod_run <- mod$run(tt)
mod_seas_run <- mod_seas$run(tt)
# print('ran seasonal model')
# View(mod_run)
# shape output
out <- mod$transform_variables(mod_run)
# windows(10,8)
plot(out$t,out$prev_all,type = 'l')

out_seas <- mod_seas$transform_variables(mod_seas_run)
# windows(10,8)
plot(out_seas$t,out_seas$prev_all,type = 'l')

# View(out)
print(out$mv_init[1])
print(out$mv0_init[1])
init4pmcmc <- transform_init(out)
