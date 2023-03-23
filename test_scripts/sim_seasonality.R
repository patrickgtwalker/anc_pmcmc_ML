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
library(reshape2)

source('shared/run_pmcmc.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')
source('shared/utils.R')

gen_seasonal_sim <- function(model_path = "shared/odin_model_stripped_seasonal.R",
                             init_age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80),
                             prop_treated = 0.4,
                             het_brackets = 5,
                             country = 'Burkina Faso',
                             admin_unit = 'Cascades',
                             init_EIR = 50,
                             EIR_SD = 1,
                             max_EIR = 1000,
                             state_check = 0,
                             lag_rates = 10,
                             start_stoch = 1,
                             time_origin = 1,
                             seasonality_on = 1,
                             time_length = 10,
                             start_date = '2010-01-01',
                             length_out = 5){
  season_model <- odin::odin(model_path)
  #Provide age categories, proportion treated, and number of heterogeneity brackets
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
  tt <- c(0:(time_length*365))
  
  # run seasonality model
  mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)
  
  # shape output
  out <- mod$transform_variables(mod_run)
  true_val <- data.frame(t=out$t,
             prev=out$prev,
             prev_all = out$prev_all,
             inc05=out$inc05,
             inc=out$inc)%>%
    mutate(date = as.character(as.Date(start_date)+t),
           month = as.yearmon(date))
  ##Simulate 5 years of data from this true run
  
  out_start <- as.Date(paste0(year(max(true_val$date))-length_out+1,'-01-01'))
  # print(out_start)
  obs_times <- data.frame(date=as.character(as.Date(as.yearmon(seq(out_start,by='month',length.out = 12*length_out)),frac=0.5)))
  # print('obs_times')
  # print(class(obs_times$date))
  # print(obs_times)
  # print('true_val')
  # print(class(true_val$date))
  # print(true_val)
  # sim_obs <- merge(true_val,obs_times,by='date')
  sim_obs <- left_join(obs_times, true_val, by="date")%>%
    mutate(date = as.Date(date))
  # print(sim_obs)
  sim_obs$tested <- round(rnorm(nrow(sim_obs),1000,10))
  sim_obs$positive <- rbinom(nrow(sim_obs),sim_obs$tested,sim_obs$prev)

  # Find the first peak and trough during the first year of interest
  first_peak <- true_val[true_val$inc==max(true_val[(365*(time_length-length_out)):(365*(time_length-length_out+1)),]$inc),]$month
  first_trough <- true_val[true_val$inc==min(true_val[(365*(time_length-length_out)):(365*(time_length-length_out+1)),]$inc),]$month
  sim_obs_peak <- sim_obs[sim_obs$month>=as.yearmon(first_peak),]
  sim_obs_trough <- sim_obs[sim_obs$month>=as.yearmon(first_trough),]
  sim_obs_peakplus3 <- sim_obs[sim_obs$month>=as.yearmon(first_peak)+3/12,]
  sim_obs_troughplus3 <- sim_obs[sim_obs$month>=as.yearmon(first_trough)+3/12,]
  return(list(true_val = true_val,
              sim_obs = sim_obs,
              sim_obs_peak = sim_obs_peak,
              sim_obs_peakplus3 = sim_obs_peakplus3,
              sim_obs_trough = sim_obs_trough,
              sim_obs_troughplus3 = sim_obs_troughplus3))
}
plot(out$t,out$prev,type='l')
plot(out$prev ~ out$t, t = "l", ylim = c(0, 0.8), xlab = "Year", ylab = "Prevalence")

sim_casc_med <- gen_seasonal_sim(init_EIR = 10)
sim_casc_hi <- gen_seasonal_sim(init_EIR = 100)
sim_casc_low <- gen_seasonal_sim(init_EIR = 1)
sim_nord_med <- gen_seasonal_sim(admin_unit = 'Nord',init_EIR = 10)
sim_nord_hi <- gen_seasonal_sim(admin_unit = 'Nord',init_EIR = 100)
sim_nord_low <- gen_seasonal_sim(admin_unit = 'Nord',init_EIR = 1)
addCIs()
test_peak <- addCIs(sim_casc_med$sim_obs_peak,sim_casc_med$sim_obs_peak$positive,sim_casc_med$sim_obs_peak$tested)
ggplot(test_peak)+
  geom_point(aes(x=t,y=mean))+
  geom_errorbar(aes(x=t,ymin=lower,ymax=upper),width=0)+
  scale_y_continuous(limits=c(0,1))
saveRDS(sim_casc_med,'sim/sim_datasets/sim_casc_med.rds')
saveRDS(sim_casc_hi,'sim/sim_datasets/sim_casc_hi.rds')
saveRDS(sim_casc_low,'sim/sim_datasets/sim_casc_low.rds')
saveRDS(sim_nord_med,'sim/sim_datasets/sim_nord_med.rds')
saveRDS(sim_nord_hi,'sim/sim_datasets/sim_nord_hi.rds')
saveRDS(sim_nord_low,'sim/sim_datasets/sim_nord_low.rds')

sim_casc_med <- readRDS('sim/sim_datasets/sim_casc_med.rds')
sim_casc_hi <- readRDS('sim/sim_datasets/sim_casc_hi.rds')
sim_casc_low <- readRDS('sim/sim_datasets/sim_casc_low.rds')
sim_nord_med <- readRDS('sim/sim_datasets/sim_nord_med.rds')
sim_nord_hi <- readRDS('sim/sim_datasets/sim_nord_hi.rds')
sim_nord_low <- readRDS('sim/sim_datasets/sim_nord_low.rds')

sim_seas_list <- list(sim_casc_low,sim_casc_med,sim_casc_hi,
                      sim_nord_low,sim_nord_med,sim_nord_hi)

source('shared/run_pmcmc.R')
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

test_run_seas_check <- run_pmcmc(data = sim_obs_peak,
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
                                state_check = 0,
                                seasonality_check = 1)
test_run_seas_check$seas_history[[1]]$t
windows(15,5)
plot(test_run_seas_check$seas_history[[1]]$t,test_run_seas_check$seas_history[[1]]$prev)
plot(test_run_seas_check$history['prev',1,])
test_seas_traj <- sapply(1:10,function(x) test_run_seas_check$seas_history[[x]]$prev)
matplot(x=test_run_seas_check$seas_history[[1]]$t,y=test_seas_traj,type='l',col='black')
##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("shared/run_pmcmc.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R",
             'shared/utils.R')
ctx <- context::context_save("T:/jth/contexts.temp", sources = sources,
                             packages = c('dplyr','statmod','coda','zoo','lubridate','stringi','dde'),
                             package_sources = conan::conan_sources(c('mrc-ide/mode','mrc-ide/dust',"mrc-ide/odin.dust@8aef08d",'mrc-ide/mcstate')))
config_1 <- didehpc::didehpc_config(template = "32Core",cores =4, parallel = TRUE,wholenode = FALSE, cluster = 'fi--didemrchnb')
config_dide <- didehpc::didehpc_config(template = "8Core",cores =1, parallel = TRUE,wholenode = FALSE, cluster = 'fi--dideclusthn')
obj <- didehpc::queue_didehpc(ctx,config = config_1)
obj <- didehpc::queue_didehpc(ctx,config = config_dide)
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
                              preyears = 5,
                              seasonality_on = 0,
                              state_check = 0))
sim_peak_std$status()
sim_peak_std$id #"59287feb9176d06c8320be2948156166"
                #"64f6b6cebd6bf0047d96148881833c9f" preyears = 5
sim_peak_std$log()
plot_particle_filter(sim_peak_std$result()$history,true_history=sim_obs_peak,times=sim_obs_peak$date)
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
                                      preyears = 5,
                                      seasonality_on = 0,
                                      state_check = 0))
sim_trough_std$status()
sim_trough_std$id #"cbad64d6b0bfa082596da5004accef01"
                  #"c0393cfdf2a362fef5b64f05fc686cb4" <- preyears = 5
plot_particle_filter(sim_trough_std$result()$history,true_history=sim_obs_trough,times=sim_obs_trough$date)
obj$login()
obj$cluster_load(TRUE)
obj$config
sim_peak_seas <- obj$enqueue(run_pmcmc(data = sim_obs_peak,
                                      n_particles = 200,
                                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                      max_EIR=1000,
                                      max_steps = 1e7,
                                      atol = 1e-5,
                                      rtol = 1e-6,
                                      n_steps = 100,
                                      n_threads = 8,
                                      lag_rates = 10,
                                      country = 'Burkina Faso',
                                      admin_unit = 'Cascades',
                                      preyears = 5,
                                      seasonality_on = 1,
                                      state_check = 0,
                                      seasonality_check = 1))

sim_peak_seas$status()
sim_peak_seas$id #"1d84e3d602faeae2c4866fcf97b9eaf7"
                #"ddf96e97ee1a452c4d7c8dcb50a0f848" <- run with seasonality check on
                #"11e51e93986f647985f68a9acb443e75" <- preyears=5
obj$unsubmit(c("c15892374daea956b1d237134b24d0fb","ba84a35f68cad4121dc570a2cfb41464",sim_peak_seas$id))
sim_peak_seas_results <- sim_peak_seas$result()
sim_trough_seas <- obj$enqueue(run_pmcmc(data = sim_obs_trough,
                                        n_particles = 200,
                                        proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                        max_EIR=1000,
                                        max_steps = 1e7,
                                        atol = 1e-5,
                                        rtol = 1e-6,
                                        n_steps = 100,
                                        n_threads = 8,
                                        lag_rates = 10,
                                        country = 'Burkina Faso',
                                        admin_unit = 'Cascades',
                                        preyears = 5,
                                        seasonality_on = 1,
                                        state_check = 0,
                                        seasonality_check = 1))
sim_trough_seas$status()
sim_trough_seas$id #"9b958d300ccf3fe95701769fd33bb5dc"
                  #"af0e0ac225c866176c9fb4f51ad19364" <-Run with seasonality check on
                  #"e3ff6834fa8f7a7db1c4bf59c8fe7ddc" <- preyears=5
windows(15,5)
plot(sim_trough_seas$result()$seas_history[[1]]$t,sim_trough_seas$result()$seas_history[[1]]$prev)
plot(sim_trough_seas$result()$history['prev',1,])
obj$task_get("59287feb9176d06c8320be2948156166")$result()
results<-list(sim_peak_std$result(),sim_trough_std$result(),sim_peak_seas$result(),sim_trough_seas$result())
results<-list(obj$task_get("59287feb9176d06c8320be2948156166")$result(),obj$task_get("cbad64d6b0bfa082596da5004accef01")$result(),obj$task_get("ddf96e97ee1a452c4d7c8dcb50a0f848")$result(),obj$task_get("af0e0ac225c866176c9fb4f51ad19364")$result())
create_diag_figs <- function(result,model,start){
  print('acceptance rate')
  print(1 - coda::rejectionRate(as.mcmc(result$mcmc)))
  print('effective size')
  print(coda::effectiveSize(as.mcmc(result$mcmc)))
  
  title <- paste0('Diagnostic plots for ',model,' model at ',start)
  
  diag <- ((bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'log_prior')+mcmc_dens(result$mcmc[51:1000,],pars = 'log_prior'))/
             (bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'log_likelihood')+mcmc_dens(result$mcmc[51:1000,],pars = 'log_likelihood'))/
             (bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'log_posterior')+mcmc_dens(result$mcmc[51:1000,],pars = 'log_posterior'))/
             (bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'EIR_SD')+mcmc_dens(result$mcmc[51:1000,],pars = 'EIR_SD'))/
             (bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'log_init_EIR')+mcmc_dens(result$mcmc[51:1000,],pars = 'log_init_EIR'))) + 
    plot_layout(guides = "collect") + plot_annotation(title = title)
  
  
  return(diag)
}

windows(10,7)
create_diag_figs(results[[1]],'Standard','Peak')
windows(10,7)
create_diag_figs(results[[2]],'Standard','Trough')
windows(10,7)
create_diag_figs(results[[3]],'Seasonal','Peak')
windows(10,7)
create_diag_figs(results[[4]],'Seasonal','Trough')
i=16
test_seas_peak_prev <- bind_rows(lapply(1:100,function(x){
  as.data.frame(sim_peak_seas$result()$seas_history[[x]][c('t','prev')])
  }
  ),.id = "column_label")
test_seas_peak_prev$date <- as.Date(min(as.Date(sim_obs_peak$month))-(max(test_seas_peak_prev$t)-test_seas_peak_prev$t))
test_seas_peak_inc <- bind_rows(lapply(1:100,function(x){
  as.data.frame(obj$task_get("ddf96e97ee1a452c4d7c8dcb50a0f848")$result()$seas_history[[x]][c('t','inc')])
}
),.id = "column_label")
test_seas_peak_inc$date <- as.Date(min(as.Date(sim_obs_peak$month))-(max(test_seas_peak_inc$t)-test_seas_peak_inc$t))

ggplot(test_seas_peak_prev)+
  geom_line(aes(x=t,y=prev,group=column_label))
matplot(x=test_run_seas_check$seas_history[[1]]$t,y=test_seas_traj,type='l',col='black')


df <- data.frame(time = integer(),
                 median = numeric(),
                 mean = numeric(),
                 upper = numeric(),
                 lower = numeric(),
                 model = character(),
                 start = character(),
                 date = Date())
df_sample <- data.frame(time = integer(),
                        value = numeric(),
                        model = character(),
                        start = character(),
                        date = Date())
models <- list('Standard','Standard','Seasonal','Seasonal')
starts <- list('Peak','Trough','Peak','Trough')
dates_list <- list(sim_obs_peak$date,sim_obs_trough$date,sim_obs_peak$date,sim_obs_trough$date)
for(i in 1:length(results)){
  # prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -length(dates_list[[i]])]))
  prev_history <- data.frame(t(results[[i]]$history['prev', , -length(dates_list[[i]])]))
  
  
  long_prev_sum <- prev_history%>%
    mutate(t=c(1:nrow(prev_history)))%>%
    melt(id='t')%>%
    rename(time=t)%>%
    group_by(time)%>%
    summarise(median=median(value),
              mean=mean(value),
              upper=quantile(value,probs=0.975),
              lower=quantile(value,probs=0.025))%>%
    mutate(model = models[[i]],
           start = starts[[i]],
           date = dates_list[[i]])
  df <- rbind(df,long_prev_sum)
  
  prev_sample <- prev_history[, sample(ncol(prev_history), 100)] %>%
    mutate(t=c(1:nrow(prev_history)))%>%
    melt(id='t')%>%
    rename(time=t)%>%
    mutate(model = models[[i]],
           start = starts[[i]],
           date = rep(dates_list[[i]],100))
  df_sample <- rbind(df_sample,prev_sample)
  
  
}


###Summarise and make figures
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))
true_prev_plot <- ggplot(true_val)+
  geom_line(aes(x=date,y=prev),linewidth=1)+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Prevalence')+
  theme(
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"),
        axis.title.x = element_blank()
  )
true_prev_plot
ggsave(paste0('Q:/anc_pmcmc/sim/figures/true_prev_060223.pdf'),plot=true_prev_plot, width = 7, height = 5)
true_inc_plot <- ggplot(true_val)+
  geom_line(aes(x=date,y=inc*1000),linewidth=1)+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,10),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Incidence\nper 1,000 people')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank()
  )
true_inc_plot
ggsave(paste0('Q:/anc_pmcmc/sim/figures/true_inc_060223.pdf'),plot=true_inc_plot, width = 7, height = 5)
windows(7,5)
true_inc_prev <- true_prev_plot+true_inc_plot + plot_layout(nrow=2)
ggsave(paste0('Q:/anc_pmcmc/sim/figures/true_previnc_060223.pdf'),plot=true_inc_prev, width = 7, height = 5)

sim_obs <- addCIs(sim_obs,sim_obs$positive,sim_obs$tested)
obs_prev_plot <- ggplot(true_val)+
  geom_line(aes(x=date,y=prev),linewidth=1)+
  geom_point(data=sim_obs,aes(x=date,y=mean),pch = 19,color="#999999")+
  geom_errorbar(data=sim_obs,aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Prevalence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank()
  )
obs_prev_plot
ggsave(paste0('Q:/anc_pmcmc/sim/figures/obs_prev_060223.pdf'),plot=obs_prev_plot, width = 7, height = 5)
obs_inc_prev <- obs_prev_plot+true_inc_plot + plot_layout(nrow=2)
ggsave(paste0('Q:/anc_pmcmc/sim/figures/obs_prev_060223-2.pdf'),plot=obs_inc_prev, width = 7, height = 5)

est_prev_plot_peak_std <- ggplot(true_val)+
  geom_line(aes(x=date,y=prev),linewidth=1)+
  geom_point(data=sim_obs,aes(x=date,y=mean),pch = 19,color="#999999")+
  geom_errorbar(data=sim_obs,aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
  geom_line(data=df_sample[df_sample$model=='Standard'&df_sample$start=='Peak',],aes(x=date,y=value,group=variable),alpha=0.1,color="#377EB8")+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df[df$model=='Standard'&df$start=='Peak',],aes(x=date,y=median),linewidth=1,color="#377EB8")+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Prevalence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank()
  )
est_prev_plot_peak_std
est_prev_plot_trough_std <- ggplot(true_val)+
  geom_line(aes(x=date,y=prev),linewidth=1)+
  geom_point(data=sim_obs,aes(x=date,y=mean),pch = 19,color="#999999")+
  geom_errorbar(data=sim_obs,aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
  geom_line(data=df_sample[df_sample$model=='Standard'&df_sample$start=='Trough',],aes(x=date,y=value,group=variable),alpha=0.1,color="#1B9E77")+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df[df$model=='Standard'&df$start=='Trough',],aes(x=date,y=median),linewidth=1,color="#1B9E77")+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Prevalence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank()
  )
est_prev_plot_trough_std
test_seas_peak_prev <- bind_rows(lapply(1:100,function(x){
  as.data.frame(sim_peak_seas$result()$seas_history[[x]][c('t','prev')])
}
),.id = "column_label")

est_prev_plot_peak_seas <- ggplot(true_val)+
  geom_line(aes(x=date,y=prev),linewidth=1)+
  geom_point(data=sim_obs,aes(x=date,y=mean),pch = 19,color="#999999")+
  geom_errorbar(data=sim_obs,aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
  geom_line(data=df_sample[df_sample$model=='Seasonal'&df_sample$start=='Peak',],aes(x=date,y=value,group=variable),alpha=0.1,color="#377EB8")+
  geom_line(data=test_seas_peak_prev,aes(x=date,y=prev,group=column_label))+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df[df$model=='Seasonal'&df$start=='Peak',],aes(x=date,y=median),linewidth=1,color="#377EB8")+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Prevalence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank()
  )
est_prev_plot_peak_seas
est_prev_plot_trough_seas <- ggplot(true_val)+
  geom_line(aes(x=date,y=prev),linewidth=1)+
  geom_point(data=sim_obs,aes(x=date,y=mean),pch = 19,color="#999999")+
  geom_errorbar(data=sim_obs,aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
  geom_line(data=df_sample[df_sample$model=='Seasonal'&df_sample$start=='Trough',],aes(x=date,y=value,group=variable),alpha=0.1,color="#1B9E77")+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df[df$model=='Seasonal'&df$start=='Trough',],aes(x=date,y=median),linewidth=1,color="#1B9E77")+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Prevalence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank()
  )
est_prev_plot_trough_seas

df_inc <- data.frame(time = integer(),
                 median = numeric(),
                 mean = numeric(),
                 upper = numeric(),
                 lower = numeric(),
                 model = character(),
                 start = character(),
                 date = Date())
df_inc_sample <- data.frame(time = integer(),
                        value = numeric(),
                        model = character(),
                        start = character(),
                        date = Date())
models <- list('Standard','Standard','Seasonal','Seasonal')
starts <- list('Peak','Trough','Peak','Trough')
dates_list <- list(sim_obs_peak$date,sim_obs_trough$date,sim_obs_peak$date,sim_obs_trough$date)
for(i in 1:length(results)){
  inc_history <- data.frame(t(results[[i]]$history['inc', , -length(dates_list[[i]])]))
  # inc_history <- data.frame(t(results[[i]]$history['inc', 51:1000, -length(dates_list[[i]])]))
  
  
  long_inc_sum <- inc_history%>%
    mutate(t=c(1:nrow(inc_history)))%>%
    melt(id='t')%>%
    rename(time=t)%>%
    group_by(time)%>%
    summarise(median=median(value),
              mean=mean(value),
              upper=quantile(value,probs=0.975),
              lower=quantile(value,probs=0.025))%>%
    mutate(model = models[[i]],
           start = starts[[i]],
           date = dates_list[[i]])
  df_inc <- rbind(df_inc,long_inc_sum)
  
  inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
    mutate(t=c(1:nrow(inc_history)))%>%
    melt(id='t')%>%
    rename(time=t)%>%
    mutate(model = models[[i]],
           start = starts[[i]],
           date = rep(dates_list[[i]],100))
  df_inc_sample <- rbind(df_inc_sample,inc_sample)
  
  
}

est_inc_plot_peak_std <- ggplot(true_val)+
  geom_line(aes(x=date,y=inc*1000),linewidth=1)+
  geom_line(data=df_inc_sample[df_inc_sample$model=='Standard'&df_inc_sample$start=='Peak',],aes(x=date,y=value*1000,group=variable),alpha=0.1,color="#377EB8")+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df_inc[df_inc$model=='Standard'&df_inc$start=='Peak',],aes(x=date,y=median*1000),linewidth=1,color="#377EB8")+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,10),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Incidence\nper 1,000 people')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank())
est_inc_plot_peak_std
windows(7,5)
est_prev_plot_peak_std + est_inc_plot_peak_std + plot_layout(nrow=2)
est_inc_plot_trough_std <- ggplot(true_val)+
  geom_line(aes(x=date,y=inc*1000),linewidth=1)+
  geom_line(data=df_inc_sample[df_inc_sample$model=='Standard'&df_inc_sample$start=='Trough',],aes(x=date,y=value*1000,group=variable),alpha=0.1,color="#1B9E77")+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df_inc[df_inc$model=='Standard'&df_inc$start=='Trough',],aes(x=date,y=median*1000),linewidth=1,color="#1B9E77")+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,10),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Incidence\nper 1,000 people')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank())
est_inc_plot_trough_std
windows(10,7)
est_inc_plot_peak_seas <- ggplot(true_val)+
  geom_line(aes(x=date,y=inc*1000),linewidth=1)+
  geom_line(data=df_inc_sample[df_inc_sample$model=='Seasonal'&df_inc_sample$start=='Peak',],aes(x=date,y=value*1000,group=variable),alpha=0.1,color="#377EB8")+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df_inc[df_inc$model=='Seasonal'&df_inc$start=='Peak',],aes(x=date,y=median*1000),linewidth=1,color="#377EB8")+
  geom_line(data=test_seas_peak_inc,aes(x=date,y=inc*1000,group=column_label))+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,10),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Incidence\nper 1,000 people')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank())
est_inc_plot_peak_seas


est_peak_prev <- est_prev_plot_trough_std + est_prev_plot_peak_std + true_inc_plot + true_inc_plot + plot_layout(ncol=2, nrow=2)
windows(14,5)
est_peak_prev
ggsave(paste0('Q:/anc_pmcmc/sim/figures/est_prev_060223.pdf'),plot=est_peak_prev, width = 14, height = 5)

est_previnc <- est_prev_plot_trough_std + est_prev_plot_peak_std + est_inc_plot_trough_std + est_inc_plot_peak_std + plot_layout(ncol=2, nrow=2)
windows(14,5)
est_previnc
ggsave(paste0('Q:/anc_pmcmc/sim/figures/est_previnc_060223.pdf'),plot=est_previnc, width = 14, height = 5)

est_inc_plot_peak_seas <- ggplot(true_val)+
  geom_line(aes(x=date,y=inc*1000),linewidth=1)+
  geom_line(data=df_inc_sample[df_inc_sample$model=='Seasonal'&df_inc_sample$start=='Peak',],aes(x=date,y=value,group=variable),alpha=0.1,color="#377EB8")+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df_inc[df_inc$model=='Seasonal'&df_inc$start=='Peak',],aes(x=date,y=median),linewidth=1,color="#377EB8")+
  scale_x_date(date_labels = "%b %Y")+
  coord_cartesian(ylim=c(0, 0.015),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Incidence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
  )
est_inc_plot_peak_seas
est_inc_plot_trough_seas <- ggplot(true_val)+
  geom_line(aes(x=date,y=inc))+
  geom_line(data=df_inc_sample[df_inc_sample$model=='Seasonal'&df_inc_sample$start=='Trough',],aes(x=date,y=value,group=variable),alpha=0.1,color="#1B9E77")+
  # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  geom_line(data=df_inc[df_inc$model=='Seasonal'&df_inc$start=='Trough',],aes(x=date,y=median),linewidth=1,color="#1B9E77")+
  scale_x_date(date_labels = "%b %Y")+
  coord_cartesian(ylim=c(0, 0.015),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Incidence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
  )
est_inc_plot_trough_seas

##Bulk submission
sim_seas_list
obj$login()
obj$cluster_load(TRUE)
obj$config
admin_list <- c(rep('Cascades',3),rep('Nord',3))
sim_seas_peak_bulk <- obj$enqueue_bulk(1:6, function(x,obs_list,admin_list) {run_pmcmc(data = obs_list[[x]]$sim_obs_peak,
                                         n_particles = 200,
                                         proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                         max_EIR=1000,
                                         max_steps = 1e7,
                                         atol = 1e-5,
                                         rtol = 1e-6,
                                         n_steps = 1000,
                                         n_threads = 4,
                                         lag_rates = 10,
                                         country = 'Burkina Faso',
                                         admin_unit = admin_list[x],
                                         preyears = 5,
                                         seasonality_on = 1,
                                         state_check = 0,
                                         seasonality_check = 1)},obs_list=sim_seas_list,admin_list=admin_list)
sim_seas_peak_bulk$status()#'adamantium_slothbear'
sim_seas_peak_bulk <- obj$task_bundle_get('adamantium_slothbear')
obj$unsubmit(sim_seas_peak_bulk$ids)
sim_seas_peak_bulk$tasks$`468ee673981b4ddba8c66bf293f845ad`$log()
sim_seas_peak_result_list <- lapply(1:6, function(id){
  sim_seas_peak_bulk$tasks[[id]]$result()
})
windows(10,7)
create_diag_figs(sim_seas_peak_result_list[[3]],'Seasonal','Peak')
ggsave('')
plot(sim_seas_list[[1]]$sim_obs_trough$month,sim_seas_list[[1]]$sim_obs_trough$prev)
plot(sim_seas_list[[2]]$sim_obs_trough$month,sim_seas_list[[2]]$sim_obs_trough$prev)
plot(sim_seas_list[[3]]$sim_obs_trough$month,sim_seas_list[[3]]$sim_obs_trough$prev)
sim_seas_trough_bulk <- obj$enqueue_bulk(1:6, 
                                       function(x,obs_list,admin_list) {
                                         run_pmcmc(data = obs_list[[x]]$sim_obs_trough,
                                                   n_particles = 200,
                                                   proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                   max_EIR=1000,
                                                   max_steps = 1e7,
                                                   atol = 1e-5,
                                                   rtol = 1e-6,
                                                   n_steps = 1000,
                                                   n_threads = 4,
                                                   lag_rates = 10,
                                                   country = 'Burkina Faso',
                                                   admin_unit = admin_list[x],
                                                   preyears = 5,
                                                   seasonality_on = 1,
                                                   state_check = 0,
                                                   seasonality_check = 1)},obs_list=sim_seas_list,admin_list=admin_list)
sim_seas_trough_bulk$status()#'loath_crustacean'
sim_seas_trough_bulk <- obj$task_bundle_get('loath_crustacean')
sim_seas_trough_bulk$tasks$e06d1ff3f15075dd10f5245f15fb7e35$log()
sim_seas_trough_result_list <- lapply(1:6, function(id){
  sim_seas_trough_bulk$tasks[[id]]$result()
})
View(sim_seas_trough_result_list[[2]]$history)
windows(10,7)
create_diag_figs(sim_seas_trough_result_list[[6]],'Seasonal','Trough')
sim_seas_peakplus3_bulk <- obj$enqueue_bulk(1:6, 
                                              function(x,obs_list,admin_list) {
                                                run_pmcmc(data = obs_list[[x]]$sim_obs_peakplus3,
                                                          n_particles = 200,
                                                          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                          max_EIR=1000,
                                                          max_steps = 1e7,
                                                          atol = 1e-5,
                                                          rtol = 1e-6,
                                                          n_steps = 1000,
                                                          n_threads = 4,
                                                          lag_rates = 10,
                                                          country = 'Burkina Faso',
                                                          admin_unit = admin_list[x],
                                                          preyears = 5,
                                                          seasonality_on = 1,
                                                          state_check = 0,
                                                          seasonality_check = 1)},obs_list=sim_seas_list,admin_list=admin_list)
sim_seas_peakplus3_bulk$status()#'wrinkly_triceratops'
sim_seas_peakplus3_bulk <- obj$task_bundle_get('wrinkly_triceratops')
sim_seas_peakplus3_bulk$tasks$`063aedf125e654e028b75380d511a62f`$log()
sim_seas_peakplus3_result_list <- lapply(2:6, function(id){
  sim_seas_peakplus3_bulk$tasks[[id]]$result()
})

sim_seas_troughplus3_bulk <- obj$enqueue_bulk(1:6, 
                                         function(x,obs_list,admin_list) {
                                           run_pmcmc(data = obs_list[[x]]$sim_obs_troughplus3,
                                                     n_particles = 200,
                                                     proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                     max_EIR=1000,
                                                     max_steps = 1e7,
                                                     atol = 1e-5,
                                                     rtol = 1e-6,
                                                     n_steps = 1000,
                                                     n_threads = 4,
                                                     lag_rates = 10,
                                                     country = 'Burkina Faso',
                                                     admin_unit = admin_list[x],
                                                     preyears = 5,
                                                     seasonality_on = 1,
                                                     state_check = 0,
                                                     seasonality_check = 1)},obs_list=sim_seas_list,admin_list=admin_list)
#'prosaic_herring'
sim_seas_troughplus3_bulk$status()
obj$login()
sim_seas_troughplus3_bulk <- obj$task_bundle_get('prosaic_herring')
sim_seas_troughplus3_result_list <- lapply(1:6, function(id){
  sim_seas_troughplus3_bulk$tasks[[id]]$result()
})


###Standard runs on same data as above
obj$login()
obj$cluster_load(TRUE)
obj$config
sim_std_peak_bulk <- obj$enqueue_bulk(1:6, function(x,obs_list,admin_list) 
{run_pmcmc(data = obs_list[[x]]$sim_obs_peak,
           n_particles = 200,
           proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
           max_EIR=1000,
           max_steps = 1e7,
           atol = 1e-5,
           rtol = 1e-6,
           n_steps = 1000,
           n_threads = 4,
           lag_rates = 10,
           country = 'Burkina Faso',
           admin_unit = admin_list[x],
           preyears = 5,
           seasonality_on = 0,
           state_check = 0,
           seasonality_check = 0)},obs_list=sim_seas_list,admin_list=admin_list)
sim_std_peak_bulk$status()#'fatherly_weevil'
sim_std_peak_bulk <- obj$task_bundle_get('fatherly_weevil')
sim_std_peak_result_list <- lapply(1:6, function(id){
  sim_std_peak_bulk$tasks[[id]]$result()
})

sim_std_trough_bulk <- obj$enqueue_bulk(1:6, 
                                         function(x,obs_list,admin_list) {
                                           run_pmcmc(data = obs_list[[x]]$sim_obs_trough,
                                                     n_particles = 200,
                                                     proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                     max_EIR=1000,
                                                     max_steps = 1e7,
                                                     atol = 1e-5,
                                                     rtol = 1e-6,
                                                     n_steps = 1000,
                                                     n_threads = 4,
                                                     lag_rates = 10,
                                                     country = 'Burkina Faso',
                                                     admin_unit = admin_list[x],
                                                     preyears = 5,
                                                     seasonality_on = 0,
                                                     state_check = 0,
                                                     seasonality_check = 0)},obs_list=sim_seas_list,admin_list=admin_list)
##'susceptible_noctilio'
sim_std_trough_bulk$status()
sim_std_trough_bulk <- obj$task_bundle_get('susceptible_noctilio')
sim_std_trough_result_list <- lapply(1:6, function(id){
  sim_std_trough_bulk$tasks[[id]]$result()
})

sim_std_peakplus3_bulk <- obj$enqueue_bulk(1:6, 
                                            function(x,obs_list,admin_list) {
                                              run_pmcmc(data = obs_list[[x]]$sim_obs_peakplus3,
                                                        n_particles = 200,
                                                        proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                        max_EIR=1000,
                                                        max_steps = 1e7,
                                                        atol = 1e-5,
                                                        rtol = 1e-6,
                                                        n_steps = 1000,
                                                        n_threads = 4,
                                                        lag_rates = 10,
                                                        country = 'Burkina Faso',
                                                        admin_unit = admin_list[x],
                                                        preyears = 5,
                                                        seasonality_on = 0,
                                                        state_check = 0,
                                                        seasonality_check = 0)},obs_list=sim_seas_list,admin_list=admin_list)
sim_std_peakplus3_bulk$status() #'contradictious_serval'
sim_std_peakplus3_bulk$tasks[[1]]$log()
sim_std_peakplus3_bulk <- obj$task_bundle_get('contradictious_serval')
sim_std_peakplus3_result_list <- lapply(2:6, function(id){
  sim_std_peakplus3_bulk$tasks[[id]]$result()
})

sim_std_troughplus3_bulk <- obj$enqueue_bulk(1:6, 
                                              function(x,obs_list,admin_list) {
                                                run_pmcmc(data = obs_list[[x]]$sim_obs_troughplus3,
                                                          n_particles = 200,
                                                          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                                          max_EIR=1000,
                                                          max_steps = 1e7,
                                                          atol = 1e-5,
                                                          rtol = 1e-6,
                                                          n_steps = 1000,
                                                          n_threads = 4,
                                                          lag_rates = 10,
                                                          country = 'Burkina Faso',
                                                          admin_unit = admin_list[x],
                                                          preyears = 5,
                                                          seasonality_on = 0,
                                                          state_check = 0,
                                                          seasonality_check = 0)},obs_list=sim_seas_list,admin_list=admin_list)
sim_std_troughplus3_bulk$status() #'underterrestrial_buckeyebutterfly'
sim_std_troughplus3_bulk <- obj$task_bundle_get('underterrestrial_buckeyebutterfly')
sim_std_troughplus3_result_list <- lapply(1:6, function(id){
  sim_std_troughplus3_bulk$tasks[[id]]$result()
})

##Summarize simulation runs
windows(10,7)
create_diag_figs(sim_seas_peak_result_list[[6]],'Seasonal','Peak')
create_diag_figs(sim_seas_trough_result_list[[6]],'Seasonal','Trough')
create_diag_figs(sim_seas_peakplus3_result_list[[5]],'Seasonal','Peak')
create_diag_figs(sim_seas_troughplus3_result_list[[6]],'Seasonal','Trough')

create_diag_figs(sim_std_peak_result_list[[6]],'Standard','Peak')
create_diag_figs(sim_std_trough_result_list[[6]],'Standard','Trough')
create_diag_figs(sim_std_peakplus3_result_list[[6]],'Standard','Peak')
create_diag_figs(sim_std_troughplus3_result_list[[6]],'Standard','Trough')


##Summarize diagnostics
create_diag_figs_bulk <- function(results){
  ar.df <- data.frame(sites = names(results),
                      acceptance_rate = sapply(1:length(results), function(x){
                        1 - coda::rejectionRate(as.mcmc(results[[x]]$mcmc))[[1]]})
  )
  
  ess.df <- bind_rows(lapply(1:length(results), 
                             function(x){
                               data.frame(sites = names(results[x]),
                                          t(coda::effectiveSize(as.mcmc(results[[x]]$mcmc))))
                             }))%>%
    melt(id.vars = 'sites')
  
  mcmc.df <- bind_rows(lapply(1:length(results), 
                              function(x){
                                df <- results[[x]]$mcmc
                                df$sites <- names(results[x])
                                df$step <- 1:nrow(df)
                                return(df)
                              }))
  
  ar.plot <- ggplot(ar.df,aes(x=sites,y=acceptance_rate))+
    geom_col()+
    labs(title='Acceptance Rate by Site',x='Site',y='Acceptance Rate')+
    geom_text(aes(label = round(acceptance_rate,digits=2)), vjust = -0.5)+
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  ess.plot <- ggplot(ess.df,aes(x=variable,y=value))+
    geom_col()+
    facet_wrap(vars(sites))+
    labs(title='ESS Values by Parameter and Site',x='Parameter',y='Effective Size')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  
  prior.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=log_prior))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='Log Prior Trace by Site',y='Log Prior')+
    theme(axis.title.x = element_blank())
  
  posterior.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=log_posterior))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='Log Posterior Trace by Site',y='Log Posterior')+
    theme(axis.title.x = element_blank())
  
  likelihood.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=log_likelihood))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='Log Likelihood Trace by Site',y='Log Likelihood')+
    theme(axis.title.x = element_blank())
  
  EIR_SD.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=EIR_SD))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='EIR_SD Trace by Site',y='EIR_SD')+
    theme(axis.title.x = element_blank())
  
  init_EIR.trace <- ggplot(mcmc.df)+
    geom_line(aes(x=step,y=log_init_EIR))+
    facet_wrap(vars(sites),scales='free_y')+
    labs(title='Log init_EIR Trace by Site',y='Log init_EIR')+
    theme(axis.title.x = element_blank())
  
  EIR_SD.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=EIR_SD),fill='darkgray')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='EIR_SD density by Site',x='Site',y='EIR_SD')
  
  log_init_EIR.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=log_init_EIR),fill='darkgray')+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='Log init_EIR density by Site',x='Site',y='Log init_EIR')
  
  init_EIR.density <- ggplot(mcmc.df)+
    geom_violin(aes(x=sites,y=exp(log_init_EIR)),fill='darkgray')+
    scale_y_log10()+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    labs(title='init_EIR density by Site',x='Site',y='init_EIR')
  
  diags <- list(ar.plot = ar.plot, ess.plot = ess.plot, 
                prior.trace = prior.trace, posterior.trace = posterior.trace, 
                likelihood.trace = likelihood.trace, EIR_SD.trace =EIR_SD.trace,
                init_EIR.trace = init_EIR.trace, EIR_SD.density = EIR_SD.density,
                log_init_EIR.density = log_init_EIR.density, init_EIR.density = init_EIR.density)
  
  
  return(diags)
}

sim_casc_med <- readRDS('sim/sim_datasets/sim_casc_med.rds')
sim_casc_hi <- readRDS('sim/sim_datasets/sim_casc_hi.rds')
sim_casc_low <- readRDS('sim/sim_datasets/sim_casc_low.rds')
sim_nord_med <- readRDS('sim/sim_datasets/sim_nord_med.rds')
sim_nord_hi <- readRDS('sim/sim_datasets/sim_nord_hi.rds')
sim_nord_low <- readRDS('sim/sim_datasets/sim_nord_low.rds')

sim_seas_list <- list(casc_low = sim_casc_low, casc_med = sim_casc_med, casc_hi = sim_casc_hi,
                      nord_low = sim_nord_low, nord_med = sim_nord_med, nord_hi = sim_nord_hi)

names(sim_seas_peak_result_list) <- names(sim_seas_list)
names(sim_seas_trough_result_list) <- names(sim_seas_list)
names(sim_std_peak_result_list) <- names(sim_seas_list)
names(sim_std_trough_result_list) <- names(sim_seas_list)
names(sim_seas_peakplus3_result_list) <- names(sim_seas_list)[2:6]
names(sim_seas_troughplus3_result_list) <- names(sim_seas_list)
names(sim_std_peakplus3_result_list) <- names(sim_seas_list)[2:6]
names(sim_std_troughplus3_result_list) <- names(sim_seas_list)

sim_seas_result_list <- list(peak = sim_seas_peak_result_list,
                             peakplus3 = sim_seas_peakplus3_result_list,
                             trough = sim_seas_trough_result_list,
                             troughplus3 = sim_seas_troughplus3_result_list)
sim_std_result_list <- list(peak = sim_std_peak_result_list,
                            peakplus3 = sim_std_peakplus3_result_list,
                            trough = sim_std_trough_result_list,
                            troughplus3 = sim_std_troughplus3_result_list)
sim_all_result_list <- list(standard = sim_std_result_list, seasonal = sim_seas_result_list)

##Show prevalence trajectories, comparing with seasonal deterministic model and simulated data

##Show posterior distributions of init_EIR and EIR_SD
  