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
source('sim/run_pmcmc.R')
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
                      n_threads = 2)
plot_particle_filter(test_run$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)

##Set up cluster##
root <- "T:/jth/contexts" ##Edit to your contexts path
sources <- c("sim/run_pmcmc.R",
             "shared/model_parameters.R","shared/equilibrium-init-create-stripped.R")


ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

##Set up and run on 32 node core##
config_32 <- didehpc::didehpc_config(cores = 32, parallel = TRUE)
obj_32 <- didehpc::queue_didehpc(ctx,config = config_32)

obj_32$cluster_load(TRUE)
obj_32$config
obj_32$login()

run_32_200 <- obj_32$enqueue(run_pmcmc(data = data_sim_comptest3,
                                      n_particles = 200,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 1000,
                                      n_threads = 32))
run_32_200$status()
run_32_200$log()

##Produce some diagnostics
1 - coda::rejectionRate(as.mcmc(result_32_200$mcmc)) ##Acceptance rate
coda::effectiveSize(as.mcmc(result_32_200$mcmc)) ##ESS
cov(result_32_200$pars) ##Covariance
summary(as.mcmc(result_32_200$mcmc)) ##Summarize mcmc run
plot(as.mcmc(result_32_200$mcmc)) ##Plot traces and distributions

##Plot particle filter results vs simulated data
plot_particle_filter(result_32_200$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)

##Below creates various figures to show prevalence, EIR, and incidence trajectories
##Used for ASTMH presentation
prev_cis <- addCIs(data_sim_comptest3,data_sim_comptest3$positive,data_sim_comptest3$tested)
##Remove 5% burn-in and reformat prevalence output
prev_history <- data.frame(t(result_32_200$history['prev', 51:1000, -1]))
long_prev_history <- prev_history%>%
  mutate(time=c(1:nrow(prev_history))*30)%>%
  melt(id='time')

windows(7,5)
prev_plot <- ggplot(long_prev_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=prev_cis,aes(x=t,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=prev_cis,aes(x=t,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_plot_true <- ggplot(long_prev_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0)+
  geom_line(data=data_sim_comptest3,aes(x=t,y=prev_true),size=1,col = "#1F78B4",alpha=1)+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_plot_true
prev_data_only <- ggplot(prev_cis)+
  geom_point(aes(x=t,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(aes(x=t,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_data_andtrue <- ggplot(prev_cis)+
  geom_point(aes(x=t,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(aes(x=t,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  geom_line(data=data_sim_comptest3,aes(x=t,y=prev_true),size=1,col = "#1F78B4",alpha=1)+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_plot_blank <- ggplot(long_prev_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0)+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Time (days)',y='Prevalence')
prev_plot_blank
eir_history <- data.frame(t(result_32_200$history['EIR', 51:1000, -1]))
long_eir_history <- eir_history%>%
  mutate(time=c(1:nrow(eir_history))*30)%>%
  melt(id='time')
windows(7,5)
true_eir <- data.frame(EIR_true=rep(data_sim_comptest3$EIR_true,rep(30,61)))
true_eir$time <- 31:1860
eir_plot <- ggplot(long_eir_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=true_eir[1:1800,],aes(x=time,y=EIR_true),size=1,
            col = "#1F78B4")+
  scale_y_continuous(limits=c(0,max(c(long_eir_history$value))))+
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')
eir_plot
windows(7,5)
true_eir_plot <- ggplot(true_eir[1:1800,])+
  geom_line(aes(x=time,y=EIR_true), size=1, col = "#1F78B4")+
  scale_y_log10(limits=c(min(long_eir_history$value),max(long_eir_history$value)),breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')
true_eir_plot
eir_plot_blank <- ggplot(long_eir_history)+
  geom_line(aes(x=time,y=value,group=variable),col = "#A6CEE3",alpha=0)+
  scale_y_continuous(limits=c(0,max(c(long_eir_history$value))))+
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')
eir_plot_blank

true_eir_plot <- ggplot(true_eir[1:1800,])+
  geom_line(aes(x=time,y=EIR_true),size=1,
            col = "#666666")+
  scale_y_log10(limits=c(min(long_eir_history$value),max(long_eir_history$value)),breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')
true_eir_plot

inc_history <- data.frame(t(result_32_200$history['inc', 51:1000, -1]))
long_inc_history <- inc_history%>%
  mutate(time=c(1:nrow(eir_history))*30)%>%
  melt(id='time')
inv_med <- long_inc_history%>%
  group_by(time)%>%
  summarise(med_inc=median(value))
windows(7,5)
inc_plot <- ggplot(long_inc_history)+
  geom_line(aes(x=time,y=value*10000,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_line(data=data_sim_comptest3,aes(x=t+30,y=inc_true*10000),size=1,col = "#1F78B4",alpha=1)+
  scale_y_continuous(limits=c(0,max(long_inc_history$value)*10000))+
  labs(x='Time (days)',y='Incidence (<5 year old)\nper 10,000 persons')
inc_plot
inc_plot_blank <- ggplot(long_inc_history)+
  geom_line(aes(x=time,y=value*10000,group=variable),col = "#A6CEE3",alpha=0)+
  geom_line(data=inv_med,aes(x=time,y=med_inc*10000),col = "#1F78B4",alpha=0)+
  scale_y_continuous(limits=c(0,max(long_inc_history$value)*10000))+
  labs(x='Time (days)',y='Incidence (<5 year old)\nper 10,000 persons')
inc_plot_blank

inc_plot_true <- ggplot(long_inc_history)+
  geom_line(aes(x=time,y=value*10000,group=variable),col = "#A6CEE3",alpha=0)+
  geom_line(data=data_sim_comptest3,aes(x=t+30,y=inc_true*10000),size=1,col = "#1F78B4",alpha=1)+
  scale_y_continuous(limits=c(0,max(long_inc_history$value)*10000))+
  labs(x='Time (days)',y='Incidence (<5 year old)\nper 10,000 persons')
inc_plot_true

windows(15,5)
eir_plot_blank + inc_plot_blank + prev_plot_blank
windows(15,5)
true_eir_plot + inc_plot_blank + prev_plot_blank
true_eir_plot + inc_plot_true + prev_plot_blank
true_eir_plot + inc_plot_true + prev_plot_true
true_eir_plot + inc_plot_true + prev_data_andtrue
true_eir_plot + inc_plot_true + prev_data_only
true_eir_plot + inc_plot_true + prev_plot
true_eir_plot + inc_plot + prev_plot
eir_plot + inc_plot_true+ prev_plot
eir_plot + inc_plot + prev_plot
