library(odin)
library(ggplot2)
library(sifter)
library(zoo)
library(dplyr)
library(torch)
library(luz)
library(tidyverse)


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
source('ML/ML_functions.R')


net <- nn_module(
  initialize = function(d_in, d_hidden, d_out) {
    self$net <- nn_sequential(
      nn_linear(d_in, d_hidden),
      nn_relu(),
      nn_linear(d_hidden, d_hidden),
      nn_relu(),
      nn_linear(d_hidden, d_out),
    )
  },
  forward = function(x) {
    self$net(x)
  }
)


n_sims<-20
volatility<-0.8
init_EIR<-20
duration=20*365
out_step=30

sims_compendium_train<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_train<-sims_compendium_train%>%
  mutate(log_EIR=log(EIR_true),
         inc_true_1000=inc_true*1000,
         prev_true_zeroed=replace(prev_true,prev_true<0,1e-12),
         log_odds_prev=log(prev_true_zeroed/(1-prev_true_zeroed)))

sims_compendium_test<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_test<-sims_compendium_test%>%
  mutate(log_EIR=log(EIR_true),
         inc_true_1000=inc_true*1000,
         prev_true_zeroed=replace(prev_true,prev_true<0,1e-12),
         log_odds_prev=log(prev_true_zeroed/(1-prev_true_zeroed)))

sims_compendium_valid<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_valid<-sims_compendium_valid%>%
  mutate(log_EIR=log(EIR_true),
         inc_true_1000=inc_true*1000,
         prev_true_zeroed=replace(prev_true,prev_true<0,1e-12),
         log_odds_prev=log(prev_true_zeroed/(1-prev_true_zeroed)))

sims_compendium_test$step<-sims_compendium_test$t/30
sims_compendium_train$step<-sims_compendium_train$t/30
sims_compendium_valid$step<-sims_compendium_valid$t/30
#sims_compendium_train$step_back<-(max(sims_compendium_train$t)-sims_compendium_train$t)/30+1
#sims_compendium_test$step_back<-(max(sims_compendium_test$t)-sims_compendium_test$t)/30+1
#sims_compendium_valid$step_back<-(max(sims_compendium_valid$t)-sims_compendium_valid$t)/30+1
#sims_compendium_train$step_back<-(max(sims_compendium_train$t)-sims_compendium_train$t)/30+1

dir.create("./sim_compendia")
#saveRDS(sims_compendium_train,"./train_first_relu.RDS")
#saveRDS(sims_compendium_test,"./test_first_relu.RDS")
#saveRDS(sims_compendium_valid,"./valid_first_relu.RDS")


predict_prev_inf_lags<-generate_preds_valid_lag(model=net,fit_var="prev_true",pred_var = "log_EIR",t_var="step",middle=T,
                                                epochs=50, sims_compendium_train=sims_compendium_train,sims_compendium_test = sims_compendium_test,sims_compendium_valid = sims_compendium_valid
)
dir.create("./fitted_models")
#saveRDS(predict_prev_inf_lags,"./fitted_first_relu.RDS")
predict_prev_inf_lags$compare_plot

