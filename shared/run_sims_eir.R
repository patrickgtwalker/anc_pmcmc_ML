library(ggplot2)
library(bayesplot)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
library(scales)
library("ggsci")

theme_set(theme_minimal())

source('shared/utils.R')
source('sim/data_gen.R')
source('shared/data_gen_moz.R')
source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')

### RUN MODEL WITH SMC ####
Lindi<-readRDS("./smc_stuff/lindi_eirs_4pw.rds")
Mtwara<-readRDS("./smc_stuff/Mtwara_eirs_4pw.rds")
Ruvuma<-readRDS("./smc_stuff/ruvuma_eirs_4pw.rds")

### proph_profile of SMC ###
weibull <- function(time, alpha = 3.4, beta = 39.34) {
  pow <- -(time/beta)^alpha
  y <- exp(pow)
  return(y)
}
### for running SMC ###
get_smc_profile<-function(start_sim,end_sim,SMC_start,nround,gap,years){
  single_year<-c(rep(1-weibull(1:gap),nround-1),1-weibull(1:100))
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
check_prof<-get_smc_profile(1,1000,100,3,30,2)
plot(check_prof$SMC_times,check_prof$SMC_vals)

prop_prof=1:300
check<-seq(1,1000,by=365)


Lindi$EIR_history
draws<-seq_along(Lindi$EIR_history)[-c(1,length(Lindi$EIR_history))]

Lindi$EIR_history$date[37]
Lindi$EIR_history$t[37]
get_smc_sims<-function(EIR_draws,smc_cov,start_smc,nrounds,gap,years){
  prop_treated = 0.4
  init_age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  time= 5*365
  EIR_step=30
  out_step=0.1
  EIR_volatility=0.1
  het_brackets<-5
  model_file<-"shared/odin_model_stripped_matched.R"
  mpl <- model_param_list_create(init_EIR = 10,
                                 init_ft = prop_treated,
                                 EIR_times=EIR_draws$EIR_history$t
  )
  
  
  
  generator <- odin(model_file)
  
  model_file_smc<-"shared/odin_model_stripped_matched_smc_check.R"
  generator_smc <- odin(model_file_smc)

  
  plot_df<-data.frame(draw=numeric(),t=numeric(),inc=numeric(),inc05=numeric(),inc_target=numeric(),prev_child=numeric(),
                      inc_smc=numeric(),inc05_smc=numeric(),inc_target_smc=numeric(),prev_child_smc=numeric(),SMC_prot=numeric()
                      )
  tt=min(EIR_draws$EIR_history$t):max(EIR_draws$EIR_history$t)
  draws<-seq_along(Lindi$EIR_history)[-c(1,length(Lindi$EIR_history))]
  for(i in draws){
    print(i)
    pars <- equilibrium_init_create_stripped(age_vector = init_age,
                                             init_EIR = EIR_draws$init_EIR$init_EIR[i-1],
                                             ft = prop_treated,
                                             model_param_list = mpl,
                                             het_brackets = het_brackets)
    state_use<- pars[names(pars) %in% coef(generator)$name]
    state_use$EIR_vals<-EIR_draws$EIR_history[,i]
    state_use$EIR_times<-EIR_draws$EIR_history$t
    smc_prof<-get_smc_profile(min(EIR_draws$EIR_history$t),max(EIR_draws$EIR_history$t),start_smc,nrounds,gap,years)
    state_use_smc<-state_use
    state_use_smc$SMC_times<-smc_prof$SMC_times
    state_use_smc$SMC_vals<-smc_prof$SMC_vals
    state_use_smc$smc_cov<-smc_cov
    
    # create model with initial values
    mod <- generator$new(user = state_use, use_dde = TRUE)
    # run the simulation to base the data
    mod_run <- mod$run(tt, step_max_n = 1e7,
                                 atol = 1e-5,
                                 rtol = 1e-5)
    state_use_smc<-state_use
    state_use_smc$SMC_times<-smc_prof$SMC_times
    state_use_smc$SMC_vals<-smc_prof$SMC_vals
    state_use_smc$smc_cov<-1
    mod_smc <- generator_smc$new(user = state_use_smc, use_dde = TRUE)
    mod_run_smc <- mod_smc$run(tt, step_max_n = 1e7,
                               atol = 1e-5,
                               rtol = 1e-5)
    out_smc <- mod_smc$transform_variables(mod_run_smc)
    # shape output
    out<- mod$transform_variables(mod_run)
    plot_df<-plot_df%>%add_row(draw=i-2,t=out$t,inc=out$inc,inc05=out$inc05,inc_target=out$inc_smc,prev_child=out$prev,
                               inc_smc=out_smc$inc,inc05_smc=out_smc$inc05,inc_target_smc=out_smc$inc_smc,prev_child_smc=out_smc$prev,
                               SMC_prot=out_smc$SMC_prot)
  }
  return(plot_df)
}

SMC_start<-500
nrounds<-4
gap<-30
years<-5


lindi_sims<-get_smc_sims(Lindi,1,1476,4,30,5)
lindi_5_sims<-get_smc_sims(Lindi,1,1476,5,30,5)

lindi_summary<-lindi_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

lindi_summary_5<-lindi_5_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

mtwara_sims<-get_smc_sims(Mtwara,1,1476,4,30,5)
mtwara_5_sims<-get_smc_sims(Mtwara,1,1476,5,30,5)

ruvuma_sims<-get_smc_sims(Ruvuma,1,1476,4,30,5)
ruvuma_5_sims<-get_smc_sims(Ruvuma,1,1476,5,30,5)

mtwara_summary<-mtwara_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

mtwara_summary_5<-mtwara_5_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

ruvuma_summary<-ruvuma_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

ruvuma_summary_5<-ruvuma_5_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))


windows(height=10,width=20)
lindi_plot<-ggplot(lindi_summary,aes(x=as.Date("2015/1/1")-min(Lindi$EIR_history$t)+t,y=incu5_med*1000*30))+geom_line(col="#FF410DFF",lwd=1.5)+
  #geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
   geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=1.5)+
  geom_ribbon(aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30),fill="#6EE2FFFF",alpha=0.5)+
  geom_ribbon(aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30),fill="#FF410DFF",alpha=0.5)+
  xlab("")+ylab("Cases per 1000 children per month")+
  geom_ribbon(data=lindi_summary_5,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30),fill="green",alpha=0.5)+
  geom_line(data=lindi_summary_5,aes(y=incu5_smc_med*1000*30),col="green",lwd=1.5)+
  xlab("")+ylab("Cases per 1000 children per month")

mtwara_plot<-ggplot(mtwara_summary,aes(x=as.Date("2015/1/1")-min(Mtwara$EIR_history$t)+t,y=incu5_med*1000*30))+geom_line(col="#FF410DFF",lwd=1.5)+
  #geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=1.5)+
  geom_ribbon(aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30),fill="#FF410DFF",alpha=0.5)+
  geom_ribbon(aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30),fill="#6EE2FFFF",alpha=0.5)+
  xlab("")+ylab("Cases per 1000 children per month")+
  geom_ribbon(data=mtwara_summary_5,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30),fill="green",alpha=0.5)+
  geom_line(data=mtwara_summary_5,aes(y=incu5_smc_med*1000*30),col="green",lwd=1.5)+
  xlab("")+ylab("Cases per 1000 children per month")

ruvuma_plot<-ggplot(ruvuma_summary,aes(x=as.Date("2015/1/1")-min(Ruvuma$EIR_history$t)+t,y=incu5_med*1000*30))+geom_line(col="#FF410DFF",lwd=1.5)+
  #geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=1.5)+
  geom_ribbon(aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30),fill="#FF410DFF",alpha=0.5)+
  geom_ribbon(aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30),fill="#6EE2FFFF",alpha=0.5)+
  xlab("")+ylab("Cases per 1000 children per month")+
  geom_ribbon(data=ruvuma_summary_5,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30),fill="green",alpha=0.5)+
  geom_line(data=ruvuma_summary_5,aes(y=incu5_smc_med*1000*30),col="green",lwd=1.5)+
  xlab("")+ylab("Cases per 1000 children per month")

gridExtra::grid.arrange(lindi_plot,mtwara_plot,ruvuma_plot)

 # geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+ylab("Cases per 1000 children under 5 per year")
geom_vline(xintercept=seq(SMC_start,SMC_start+nrounds*gap-1,by=gap),lty=2)+
  

