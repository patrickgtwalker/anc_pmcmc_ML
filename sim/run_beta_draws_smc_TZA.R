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

### proph_profile of SMC ###
weibull <- function(time, alpha = 3.4, beta = 39.34) {
  pow <- -(time/beta)^alpha
  y <- exp(pow)
  return(y)
}
BF_draws<-as.data.frame(readRDS("./smc_stuff/bf_betaa_flat3.rds"))
BF_times<-readRDS("./smc_stuff/bf_smc_times.rds")
gambia_draws<-as.data.frame(readRDS("./smc_stuff/gambia_betaa_flat3.rds"))
gambia_times<-readRDS("./smc_stuff/gambia_smc_times.rds")
mali_draws<-as.data.frame(readRDS("./smc_stuff/mali_betaa_flat3.rds"))
mali_times<-readRDS("./smc_stuff/mali_smc_times.rds")
ghana_draws<-as.data.frame(readRDS("./smc_stuff/ghana_betaa_flat3.rds"))
ghana_times<-readRDS("./smc_stuff/ghana_smc_times.rds")
plot(weibull(1:100))
lines(weibull(1:100,alpha=7,beta=40))
gambia_sims_step<-get_smc_sims(1.712536,gambia_draws,gambia_times[1]+30,3,30,1,7,40)
BF_sims_step<-get_smc_sims(49.81896,BF_draws,BF_times[1]+30,3,30,1,7,40)
mali_sims_step<-get_smc_sims(7.419249,mali_draws,mali_times[1]+30,3,30,1,7,40)
ghana_sims_step_check<-get_smc_sims(111.0336,ghana_draws,ghana_times[1]+30,4,30,1,7,40)

ghana_sims_step_2yr<-get_smc_sims(111.0336,ghana_draws,ghana_times[1]+30,4,30,2,7,40)
BF_sims_step_2yr<-get_smc_sims(49.81896,BF_draws,BF_times[1]+30,4,30,2,7,40)
mali_sims_step_2yr<-get_smc_sims(7.419249,mali_draws,mali_times[1]+30,4,30,2,7,40)
gambia_sims_step_2yr<-get_smc_sims(1.712536,gambia_draws,gambia_times[1]+30,4,30,2,7,40)

ghana_summary_step_2yr<-ghana_sims_step_2yr%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))
BF_summary_step_2yr<-BF_sims_step_2yr%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))
mali_summary_step_2yr<-mali_sims_step_2yr%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

gambia_summary_step_2yr<-gambia_sims_step_2yr%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))


step_df_long_2yr<-data.frame(time=date_zero+time_vect-30,
                         inc= c(BF_summary_step_2yr$incu5_med,gambia_summary_step_2yr$incu5_med,mali_summary_step_2yr$incu5_med,ghana_summary_step_2yr$incu5_med),
                         inc_smc=c(BF_summary_step_2yr$incu5_smc_med,gambia_summary_step_2yr$incu5_smc_med,mali_summary_step_2yr$t,ghana_summary_step_2yr$incu5_smc_med),
                         site=c(rep("Burkina Faso",length(BF_summary_step_2yr$t)),rep("Gambia",length(gambia_summary_step_2yr$t)),rep("Mali",length(mali_summary_step_2yr$t)),rep("Ghana",length(ghana_summary_step_2yr$t)))
)
step_df_long_2yr$site<-factor(step_df_long_2yr$site,levels=c("Ghana","Burkina Faso","Mali","Gambia"))
ggplot(step_df_long_2yr,aes(x=time,y=inc*1000*30))+
  geom_line()+facet_wrap(~site)

ggplot(ghana_summary_step_2yr,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  #geom_ribbon(fill="#FF410DFF",lwd=2,aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30))+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+mali_times,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_ribbon(fill="#6EE2FFFF",lwd=2,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30))+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per month")

gambia_summary_step<-gambia_sims_step%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

gambia_step_effiacy<-gambia_sims_step%>%filter(t>(gambia_times[1]+30)&t<(gambia_times[1]+30+365))%>%
  group_by(draw)%>%
  summarise(efficacy=1-sum(inc_target_smc)/sum(inc_target))
summary(gambia_step_effiacy)
BF_summary_step<-BF_sims_step%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))
BF_step_effiacy<-BF_sims_step%>%filter(t>(BF_times[1]+30)&t<(BF_times[1]+30+365))%>%
  group_by(draw)%>%
  summarise(efficacy=1-sum(inc_target_smc)/sum(inc_target))

BF_step_effiacy
summary(BF_step_effiacy)
mali_summary_step<-mali_sims_step%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))
mali_step_effiacy<-mali_sims_step%>%filter(t>(mali_times[1]+30)&t<(mali_times[1]+30+365))%>%
  group_by(draw)%>%
  summarise(efficacy=1-sum(inc_target_smc)/sum(inc_target))
summary(mali_step_effiacy)
ghana_summary_step<-ghana_sims_step%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))
ghana_step_effiacy<-ghana_sims_step%>%filter(t>ghana_times[1]+30)%>%
  group_by(draw)%>%
  summarise(efficacy=1-sum(inc_target_smc)/sum(inc_target))

ggplot(mali_summary_step,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  #geom_ribbon(fill="#FF410DFF",lwd=2,aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30))+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+mali_times,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_ribbon(fill="#6EE2FFFF",lwd=2,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30))+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per month")


ggplot(BF_summary_step,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  #geom_ribbon(fill="#FF410DFF",lwd=2,aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30))+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+BF_times,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_ribbon(fill="#6EE2FFFF",lwd=2,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30))+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per month")

summary(gambia_step_effiacy)

date_zero<-min(gambia_draws$date)-min(gambia_draws$t)
ggplot(gambia_summary_step,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+gambia_times,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per year")

ggplot(ghana_summary_step,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+ghana_times,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per year")

time_vect=c(BF_summary_step$t,gambia_summary_step$t,mali_summary_step$t,ghana_summary_step$t)
step_df_long<-data.frame(time=date_zero+time_vect-30,
                         inc= c(BF_summary_step$incu5_med,gambia_summary_step$incu5_med,mali_summary_step$incu5_med,ghana_summary_step$incu5_med),
                         inc_smc=c(BF_summary_step$incu5_smc_med,gambia_summary_step$incu5_smc_med,mali_summary_step$t,ghana_summary_step$incu5_smc_med),
                         site=c(rep("Burkina Faso",length(BF_summary_step$t)),rep("Gambia",length(gambia_summary_step$t)),rep("Mali",length(mali_summary_step$t)),rep("Ghana",length(ghana_summary_step$t)))
                         )
 


step_df_long$site<-factor(step_df_long$site,levels=c("Ghana","Burkina Faso","Mali","Gambia"))
 ggplot(step_df_long,aes(x=time,y=inc*1000*30))+
  geom_line()+facet_wrap(~site)

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

  

BF_draws<-as.data.frame(readRDS("./smc_stuff/bf_betaa_flat3.rds"))
BF_times<-readRDS("./smc_stuff/bf_smc_times.rds")
gambia_draws<-as.data.frame(readRDS("./smc_stuff/gambia_betaa_flat3.rds"))
gambia_times<-readRDS("./smc_stuff/gambia_smc_times.rds")
mali_draws<-as.data.frame(readRDS("./smc_stuff/mali_betaa_flat3.rds"))
mali_times<-readRDS("./smc_stuff/mali_smc_times.rds")
ghana_draws<-as.data.frame(readRDS("./smc_stuff/ghana_betaa_flat3.rds"))
ghana_times<-readRDS("./smc_stuff/ghana_smc_times.rds")

smc_times<-as.data.frame(readRDS("./smc_stuff/ghana_smc_times.rds"))
beta_draws
start_smc<-592
beta_draws

date_zero<-min(BF_draws$date)-min(BF_draws$t)
BF_sims<-get_smc_sims(49.81896,BF_draws,BF_times[1]+30,3,30,1)
BF_summary<-BF_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))


BF_plot<-ggplot(BF_summary,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  #geom_ribbon(fill="#FF410DFF",lwd=2,aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30))+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+BF_times,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_ribbon(fill="#6EE2FFFF",lwd=2,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30))+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per month")

gambia_sims<-get_smc_sims(1.712536,gambia_draws,gambia_times[1]+30,3,30,1)
gambia_summary<-gambia_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

gambia_effiacy<-gambia_sims%>%filter(t<min(t)+365)%>%
  group_by(draw)%>%
  summarise(efficacy=1-sum(inc_target_smc)/sum(inc_target))

summary(gambia_effiacy$efficacy)
sum(gambia_summary$incu5_smc_med[])/
sum(gambia_summary$incu5_med)

1-sum(mali_summary$incu5_smc_med)/
  sum(mali_summary$incu5_med)

date_zero<-min(gambia_draws$date)-min(gambia_draws$t)
gambia_plot<-ggplot(gambia_summary,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+gambia_times,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per year")

mali_sims<-get_smc_sims(7.419249,mali_draws,mali_times[1]+30,3,30,1)
mali_summary<-mali_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

date_zero<-min(mali_draws$date)-min(mali_draws$t)
mali_plot<-ggplot(mali_summary,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+mali_times,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per year")

ghana_sims<-get_smc_sims(111.0336,ghana_draws,ghana_times[1]+30,4,30,1)
ghana_summary<-ghana_sims%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))

date_zero<-min(ghana_draws$date)-min(ghana_draws$t)
ghana_plot<-ggplot(ghana_summary,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero+ghana_times[1:4],lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  #geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per year")

windows(height=20,width=30)
ggpubr::ggarrange(ghana_plot,BF_plot,mali_plot,gambia_plot,ncol=2,nrow=2,
                  labels=c("Ghana","BF","Mali","Gambia"))
  plot(gambia_draws$date,gambia_draws$X1)
gambia_draws$d
  gap<-30
nrounds<-4
plot(weibull(1:100))
abline(v=30)
weibull(30)
plot(get_smc_profile(0,500,50,nrounds,gap)$SMC_vals)
check_SMC<-get_smc_profile(min(beta_draws$t),max(beta_draws$t),592,nrounds,gap)
plot(check_SMC$SMC_vals)
check_SMC$SMC_vals


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

init_EIR = 300
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
smc_prof<-get_smc_profile(min(beta_draws$t),max(beta_draws$t),start_smc,nrounds,gap)
state_use_beta_smc$SMC_times<-smc_prof$SMC_times
state_use_beta_smc$SMC_vals<-smc_prof$SMC_vals
state_use_beta_smc$smc_cov<-1
state_use_beta$betaa_vals<-beta_draws[,3]
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

dim(out_beta_smc$Y)
out_beta_smc$t
out_beta_smc$Y[100,4,,2]
(out_beta_smc$Y[100,4,3,2])/sum(sum(out_beta$Y[100,4,3]))

out_beta_smc$phi[100,4,3,2]
out_beta$phi[100,4,3]

out_beta$
redo<-get_smc_sims(1,beta_draws[,1:10],SMC_start,nrounds,gap)
ggplot(check,aes(x=date_zero+t,y=))+geom_line(col="#FF410DFF",aes(by=as.factor(draw)),alpha=0.2)
redo$
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
time<- 4*365
out_step=1
  
  median_inc_05<-check%>%
  group_by(t)%>%
  summarise(med_inc05=mean(inc_target),med_inc05_smc=mean(inc_target_smc))

  median_inc_05<-check%>%
    group_by(t)%>%
    summarise(med_inc05=mean(inc05),med_inc05_smc=mean(inc05_smc))
  
  head(check)
  check$prev_child
  redo%>%
  filter(t<530+365)%>%
  group_by(draw)%>%
  summarise(inc_diff=1-sum(inc05_smc)/sum(inc05))
 
  ggplot(redo,aes(x=date_zero+t,y=inc_target_smc/inc_target))+geom_line(col="#FF410DFF",aes(by=as.factor(draw)),alpha=0.2)
    
   
ggplot(redo,aes(x=date_zero+t,y=inc_target*1000*30))+geom_line(col="#FF410DFF",aes(by=as.factor(draw)),alpha=0.2)+
  geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=seq(SMC_start,SMC_start+nrounds*gap-1,by=gap),lty=2)+xlim(date_zero+check$t[1],date_zero+check$t[1]+365)+
  geom_line(aes(y=inc_target_smc*1000*30,by=as.factor(draw)),col="#6EE2FFFF",alpha=0.2)+
  geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+ylab("Cases per 1000 children under 5 per year")

ggplot(plot_df,aes(x=t,y=inc05))+geom_line(col="#FF410DFF",aes(by=as.factor(draw)),alpha=0.2)+
  geom_line(data=median_inc_05,aes(y=med_inc05),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=seq(SMC_start,SMC_start+nrounds*gap-1,by=gap),lty=2)+
geom_line(aes(y=inc05_smc,by=as.factor(draw)),col="#6EE2FFFF",alpha=0.2)+
geom_line(data=median_inc_05,aes(y=med_inc05_smc),lwd=3,col="#6EE2FFFF")+
  xlim()
