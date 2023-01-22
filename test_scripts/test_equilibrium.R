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
lag_rates <- 10
country <- 'Burkina Faso'
admin_unit <- 'Cascades'
seasonality_on <- 1 ## Run a deterministic seasonal model before running the stochastic model to get more realistic immunity levels
                    ## If seasonality_on = 0, runs the stochastic model based on the standard equilibrium solution
state_check <- 1 ## Run equilibrium checks
                  # If state_check = 1, returns expected deriv values which should equal 0
                  # If state_check = 1 and seasonality_on = 1, then the deterministic seasonal model is still run, but theta2 is forced to 1, forcing a constant seasonality profile
                  # If state_check = 0, no values are printed

mpl_pf <- model_param_list_create(init_age = init_age,
                                  pro_treated = prop_treated,
                                  het_brackets = het_brackets,
                                  max_EIR = max_EIR,
                                  state_check = state_check,
                                  lag_rates = lag_rates,
                                  country = country,
                                  admin_unit = admin_unit,
                                  seasonality_on = seasonality_on)

##Model with no seasonality
# model <- odin("shared/odin_model_stripped_matched.R")
##Model with seasonality component, but set to 1
model_seas <- odin("shared/odin_model_stripped_seasonal.R")


init_EIR <- 10
EIR_vol <- 0.2
mpl <- append(mpl_pf,list(EIR_SD = EIR_vol,
                          init_EIR = init_EIR,
                          EIR_times = c(0,0.5*365),
                          EIR_vals = c(10,10)))

state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                          init_EIR = init_EIR,
                                          ft = prop_treated,
                                          model_param_list = mpl,
                                          het_brackets = het_brackets,
                                          state_check = mpl$state_check)

state_use_seas <- state[names(state) %in% coef(model_seas)$name] #Only selects parameters needed by the odin model

# create model with initial values
mod_seas <- model_seas$new(user = state_use_seas, use_dde = TRUE)

# tt <- c(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days"))) #matched run_pmcmc.R
tt <- seq(0, 0.5*365,length.out=500)

# run model
mod_seas_run <- mod_seas$run(tt, verbose=FALSE,step_size_max=9)

out_seas <- mod_seas$transform_variables(mod_seas_run)
# windows(10,8)
plot(out_seas$t,out_seas$prev_all,type = 'l')
plot(out_seas$t,out_seas$inc,type = 'l')
plot(out_seas$t,out_seas$Sv,type = 'l')
plot(out_seas$t,out_seas$Iv,type = 'l')
plot(out_seas$t,out_seas$Ev,type = 'l')

##Equilibrium checks: All below should = 0
out_seas$S[1,,]-state$init_S[,]
out_seas$D[1,,]-state$init_D[,]
out_seas$T[1,,]-state$init_T[,]
out_seas$A[1,,]-state$init_A[,]
out_seas$U[1,,]-state$init_U[,]
out_seas$P[1,,]-state$init_P[,]
out_seas$Sv[1]-(state$init_Sv*state$mv0)
out_seas$Iv[1]-(state$init_Iv*state$mv0)
out_seas$Ev[1]-(state$init_Ev*state$mv0)
out_seas$IB_init[1,,]-state$init_IB
out_seas$ID_init[1,,]-state$init_ID
out_seas$ICA_init[1,,]-state$init_ICA
out_seas$FOI_init[1,,]-state$FOI_eq
out_seas$EIR_init[1,,]-state$EIR_eq
out_seas$PL[1]-state$init_PL
out_seas$LL[1]-state$init_LL
out_seas$EL[1]-state$init_EL
out_seas$mv[1]-state$mv0
out_seas$rel_foi_init[1,]-state$rel_foi
out_seas$foi_age_init[1,]-state$foi_age
out_seas$FOIv_init[1]-state$FOIv_eq
out_seas$FOIvij_init[1,,]-state$FOIvij_eq
out_seas$cA_init[1,,]-state$cA
# View(out)
print(out$mv_init[1])
print(out$mv0_init[1])
init4pmcmc <- transform_init(out)

##Test original equilibrium function
source('test_scripts/equilibrium-init-create.R')
mpl_pf <- model_param_list_create(init_age = init_age,
                                  pro_treated = prop_treated,
                                  het_brackets = het_brackets,
                                  max_EIR = max_EIR,
                                  state_check = state_check,
                                  lag_rates = lag_rates,
                                  country = country,
                                  admin_unit = admin_unit,
                                  seasonality_on = seasonality_on,
                                  num_int = 1)
state <- equilibrium_init_create(age_vector = mpl_pf$init_age,
                                          EIR = init_EIR,
                                          ft = prop_treated,
                                          model_param_list = mpl_pf,
                                          het_brackets = het_brackets)


##Original version
source('test_scripts/1_original_version_fns.R')
original_matched <- run_model_matched()
plot(original_matched$t,original_matched$prev,type = 'l')
matched_equilibrium <- run_model_matched()
original_matched$S[2,,,]-matched_equilibrium$init_S[,,1]
original_matched$D[2,,,]-matched_equilibrium$init_D[,,1]
original_matched$T[2,,,]-matched_equilibrium$init_T[,,1]
original_matched$A[2,,,]-matched_equilibrium$init_A[,,1]
original_matched$U[2,,,]-matched_equilibrium$init_U[,,1]
original_matched$P[2,,,]-matched_equilibrium$init_P[,,1]
original_matched$FOI[2,,,]-matched_equilibrium$init_P[,,1]
