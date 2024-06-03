# Load required packages
if (!require("pacman", character.only = TRUE)) {
  install.packages("pacman")
}
pacman::p_load(odin,ggplot2,zoo,dplyr,torch,luz,tidyverse)



##Set default theme
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA)))

#Required functions - these are all specific to the mechanistic models
source('sim/data_gen.R')
source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')

### my attempt to get the model to speak to torch...
source('ML/ML_functions.R')

### my current model 
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

### n_sims is number of simulations, with a random walk on log-scale starting at init_EIR with volatility

##NB - init EIR is something we will eventually want to vary...
n_sims<-200
init_EIR<-20
volatility<-0.8

## 20 years of monthly data 
duration=20*365
out_step=30

### training data 
sims_compendium_train<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_train<-sims_compendium_train%>%
  mutate(log_EIR=log(EIR_true))

### testing data 
sims_compendium_test<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_test<-sims_compendium_test%>%
  mutate(log_EIR=log(EIR_true))

### validation data 
sims_compendium_valid<-generate_sim_compendium(n_sims,volatility,init_EIR,duration,out_step)
sims_compendium_valid<-sims_compendium_valid%>%
  mutate(log_EIR=log(EIR_true))

### add a step of integer increments rather than t currently in days (so 30)
sims_compendium_test$step<-sims_compendium_test$t/30
sims_compendium_train$step<-sims_compendium_train$t/30
sims_compendium_valid$step<-sims_compendium_valid$t/30

dir.create("./sim_compendia")
saveRDS(sims_compendium_train,"./train_first_relu_200.RDS")
saveRDS(sims_compendium_test,"./test_first_relu_200.RDS")
saveRDS(sims_compendium_valid,"./valid_first_relu_200.RDS")

### fit model - default is 100 hidden states...
predict_prev_inf_lags<-generate_preds_valid_lag(model=net,fit_var="prev_true",pred_var = "log_EIR",t_var="step",middle=T,
                                                epochs=50, sims_compendium_train=sims_compendium_train,sims_compendium_test = sims_compendium_test,sims_compendium_valid = sims_compendium_valid
)
dir.create("./fitted_models")
predict_prev_inf_lags$compare_plot
luz_save(predict_prev_inf_lags$fit_model,"./200_model.RDS")



