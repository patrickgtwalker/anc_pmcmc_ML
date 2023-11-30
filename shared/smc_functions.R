weibull <- function(time, alpha = 3.4, beta = 39.34) {
  pow <- -(time/beta)^alpha
  y <- exp(pow)
  return(y)
}

### for running SMC ###
get_smc_profile<-function(start_sim,end_sim,SMC_start,nround,gap,years, alpha , beta){
  single_year<-c(rep(1-weibull(1:gap, alpha , beta),nround-1),1-weibull(1:100, alpha , beta))
  prop_prof<-rep(single_year,years)
  prop_times<-as.vector(sapply(seq(SMC_start,SMC_start+years*365-1,by=365),function(time){
    time:(time+length(single_year)-1)
  })
  )
  
  return(list(
    SMC_times=c(start_sim,prop_times,end_sim),
    SMC_vals=c(1,prop_prof,1)
  ))
}


get_smc_sims<-function(init_EIR,beta_draws,start_smc,nrounds,gap,years, alpha = 3.4, beta = 39.34){
  prop_treated = 0.4
  init_age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  time= 5*365
  EIR_step=30
  out_step=0.1
  EIR_volatility=0.1
  het_brackets<-5
  model_file_beta<-"shared/MiP_odin_model_nodelay.R"
  mpl <- sifter::model_param_list_create(init_EIR = init_EIR,
                                         init_ft = prop_treated,
                                         betaa_times=beta_draws$t,
                                         lag_rates = 10
  )
  pars <- sifter::equilibrium_init_create_stripped(age_vector = init_age,
                                                   init_EIR = init_EIR,
                                                   ft = prop_treated,
                                                   model_param_list = mpl,
                                                   het_brackets = het_brackets)
  
  
  generator_beta <- odin(model_file_beta)
  state_use_beta <- pars[names(pars) %in% coef(generator_beta)$name]
  model_file_beta_smc<-"shared/MiP_odin_model_nodelay_smc.R"
  generator_beta_smc <- odin(model_file_beta_smc)
  state_use_beta_smc<-state_use_beta
  smc_prof<-get_smc_profile(min(beta_draws$t),max(beta_draws$t),start_smc,nrounds,gap,years, alpha , beta)
  state_use_beta_smc$SMC_times<-smc_prof$SMC_times
  state_use_beta_smc$SMC_vals<-smc_prof$SMC_vals
  state_use_beta_smc$smc_cov<-1
  
  
  plot_df<-data.frame(draw=numeric(),t=numeric(),inc=numeric(),inc05=numeric(),inc_target=numeric(),prev_child=numeric(),
                      inc_smc=numeric(),inc05_smc=numeric(),inc_target_smc=numeric(),prev_child_smc=numeric(),SMC_prot=numeric(),
                      Yval=numeric(),Yval_smc=numeric())
  tt=min(beta_draws$t):max(beta_draws$t)
  draws<-seq_along(beta_draws)[-c(1,2)]
  for(i in draws){
    print(i)
    state_use_beta$betaa_vals<-beta_draws[,i]
    state_use_beta$betaa_vals
    # create model with initial values
    mod_beta <- generator_beta$new(user = state_use_beta, use_dde = TRUE)
    # run the simulation to base the data
    mod_run_beta <- mod_beta$run(tt, step_max_n = 1e7,
                                 atol = 1e-5,
                                 rtol = 1e-5)
    state_use_beta_smc<-state_use_beta
    state_use_beta_smc$SMC_times<-smc_prof$SMC_times
    state_use_beta_smc$SMC_vals<-smc_prof$SMC_vals
    state_use_beta_smc$smc_cov<-1
    mod_smc <- generator_beta_smc$new(user = state_use_beta_smc, use_dde = TRUE)
    mod_run_smc <- mod_smc$run(tt, step_max_n = 1e7,
                               atol = 1e-5,
                               rtol = 1e-5)
    out_beta_smc <- mod_smc$transform_variables(mod_run_smc)
    # shape output
    out_beta <- mod_beta$transform_variables(mod_run_beta)
    plot_df<-plot_df%>%add_row(draw=i-2,t=out_beta$t,inc=out_beta$inc,inc05=out_beta$inc05,inc_target=out_beta$inc_smc,prev_child=out_beta$prev,
                               inc_smc=out_beta_smc$inc,inc05_smc=out_beta_smc$inc05,inc_target_smc=out_beta_smc$inc_smc,prev_child_smc=out_beta_smc$prev,
                               SMC_prot=out_beta_smc$SMC_prot)
  }
  return(plot_df)
}