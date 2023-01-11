##summarize results##
library(ggplot2)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
RColorBrewer::display.brewer.all()
source('addCIs.R')
source('plot_particle_filter.R')
source('plot_eir.R')
source('plot_inc.R')

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA)))

##Check a couple runs
result_64_10 <- run_64_10$result()

1 - coda::rejectionRate(as.mcmc(result_64_10$mcmc))
coda::effectiveSize(as.mcmc(result_64_10$mcmc))
proposal <-  cov(pmcmc_run2$pars)
summary(run2_mcmc)
windows(60,50)
plot(as.mcmc(result_64_10$mcmc))

result_64_50 <- run_64_50$result()

1 - coda::rejectionRate(as.mcmc(result_64_50$mcmc))
coda::effectiveSize(as.mcmc(result_64_50$mcmc))
proposal <-  cov(result_64_50$pars)
summary(as.mcmc(result_64_50$mcmc))
windows(60,50)
plot(as.mcmc(result_64_50$mcmc))

1 - coda::rejectionRate(as.mcmc(result_64_200$mcmc))
coda::effectiveSize(as.mcmc(result_64_200$mcmc))
cov(result_64_200$pars)
summary(as.mcmc(result_64_200$mcmc))
windows(60,50)
plot(as.mcmc(result_64_200$mcmc))

result_64_10 <- run_64_10$result()
result_64_50 <- run_64_50$result()
result_64_200 <- run_64_200$result()
result_32_10 <- run_32_10$result()
result_32_50 <- run_32_50$result()
result_32_200 <- run_32_200$result()
result_4_10 <- run_4_10$result()
result_4_50 <- run_4_50$result()
result_4_200 <- run_4_200$result()
result_1_10 <- run_1_10$result()
result_1_50 <- run_1_50$result()
result_1_200 <- run_1_200$result()

windows(60,50)
plot_particle_filter(result_32_50$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)
plot_eir(result_32_50$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)
plot_particle_filter(result_32_200$history,true_history=data_sim_comptest3,times=data_sim_comptest3$t)
plot_eir(result_32_200$history[,101:1000,],true_history=data_sim_comptest3,times=data_sim_comptest3$t)
plot_inc(result_32_200$history[,101:1000,],true_history=data_sim_comptest3,times=data_sim_comptest3$t)


result_32_50$history['EIR',,-1]

result_1_10$run_time
result_4_10$run_time
result_32_10$run_time
result_64_10$run_time
result_1_50$run_time
result_4_50$run_time
result_32_50$run_time
result_64_50$run_time
result_1_200$run_time
result_4_200$run_time
result_32_200$run_time
result_64_200$run_time

##Plot run times vs n_particles and n_threads
df_runtimes <- data.frame(run_time = as.numeric(c(result_1_10$run_time,
                                       result_4_10$run_time,
                                       result_32_10$run_time,
                                       result_64_10$run_time,
                                       result_1_50$run_time,
                                       result_4_50$run_time,
                                       result_32_50$run_time,
                                       result_64_50$run_time,
                                       result_1_200$run_time,
                                       result_4_200$run_time,
                                       result_32_200$run_time,
                                       result_64_200$run_time)),
                          particles = as.factor(c(result_1_10$particles,
                                        result_4_10$particles,
                                        result_32_10$particles,
                                        result_64_10$particles,
                                        result_1_50$particles,
                                        result_4_50$particles,
                                        result_32_50$particles,
                                        result_64_50$particles,
                                        result_1_200$particles,
                                        result_4_200$particles,
                                        result_32_200$particles,
                                        result_64_200$particles)),
                          threads = as.factor(c(result_1_10$threads,
                                      result_4_10$threads,
                                      result_32_10$threads,
                                      result_64_10$threads,
                                      result_1_50$threads,
                                      result_4_50$threads,
                                      result_32_50$threads,
                                      result_64_50$threads,
                                      result_1_200$threads,
                                      result_4_200$threads,
                                      result_32_200$threads,
                                      result_64_200$threads)))
ggplot(df_runtimes,aes(x=threads,y=run_time/60))+
  geom_point()+
  facet_grid(~particles)+
  scale_y_log10()+
  labs(x='Number of nodes',y='Run time in miuntes (log-scale)')
result <- result_1_10
rm(df)
var_name <- 'prev'
##Plot trajectories by threads and particles
make_df <- function(result,var_name=c('prev','EIR','inc')){
  index_list <- c(prev=1,EIR=2,inc=3)
  index <- index_list[[var_name]]

  df <- data.frame(t(result$history[index, , -1]))
  df %>% mutate(t=c(1:nrow(df)))%>%
    melt(id='t')%>%
    mutate(threads=result$threads,
           particles=result$particles,
           run_time=result$run_time/60)
}
result_list <- list(result_1_10,
                   result_4_10,
                   result_32_10,
                   result_1_50,
                   result_4_50,
                   result_32_50,
                   result_1_200,
                   result_4_200,
                   result_32_200)
df_prev_all <- lapply(1:9,function(i) make_df(result_list[[i]],var_name='prev'))
df_prev_all_final <- bind_rows(df_prev_all) %>%
  rename(prev=value) %>%
  mutate(threads = factor(threads),
         particles = factor(particles),
         t=t*30)
df_eir_all <- lapply(1:9,function(i) make_df(result_list[[i]],var_name='EIR'))
df_eir_all_final <- bind_rows(df_eir_all) %>%
  rename(EIR=value) %>%
  mutate(threads = factor(threads),
         particles = factor(particles),
         t=t*30-30)
df_inc_all <- lapply(1:9,function(i) make_df(result_list[[i]],var_name='inc'))
df_inc_all_final <- bind_rows(df_inc_all) %>%
  rename(inc=value) %>%
  mutate(threads = factor(threads),
         particles = factor(particles),
         t=t*30-30)

##Add prevalence confidence intervals
prev_cis <- addCIs(data_sim_comptest3,data_sim_comptest3$positive,data_sim_comptest3$tested)
##create a df to label run times on the facet grid
time_labels <- data.frame(run_time = as.factor(round(c(result_1_10$run_time/60,
                                                  result_4_10$run_time/60,
                                                  result_32_10$run_time/60,
                                                  result_1_50$run_time/60,
                                                  result_4_50$run_time/60,
                                                  result_32_50$run_time/60,
                                                  result_1_200$run_time/60,
                                                  result_4_200$run_time/60,
                                                  result_32_200$run_time/60),1)),
                          particles = as.factor(c(result_1_10$particles,
                                                  result_4_10$particles,
                                                  result_32_10$particles,
                                                  result_1_50$particles,
                                                  result_4_50$particles,
                                                  result_32_50$particles,
                                                  result_1_200$particles,
                                                  result_4_200$particles,
                                                  result_32_200$particles)),
                          threads = as.factor(c(result_1_10$threads,
                                                result_4_10$threads,
                                                result_32_10$threads,
                                                result_1_50$threads,
                                                result_4_50$threads,
                                                result_32_50$threads,
                                                result_1_200$threads,
                                                result_4_200$threads,
                                                result_32_200$threads)))
ggplot(df_prev_all_final)+
  geom_line(aes(x=t,y=prev,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=prev_cis,aes(x=t,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=prev_cis,aes(x=t,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  geom_label(data = time_labels, aes(label=paste0(run_time,' min')), 
             x = Inf, y = Inf, hjust=1, vjust=1,
             inherit.aes = FALSE)+
  scale_y_continuous(limits=c(0,max(c(prev_cis$upper,df_prev_all_final$prev))))+
  facet_grid(threads~particles, labeller=label_both)+
  labs(x='Time (days)',y='Prevalence')

ggplot(df_eir_all_final)+
  geom_line(aes(x=t,y=EIR,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=data_sim_comptest3,aes(x=t,y=EIR_true),pch = 19,
             col = "#666666")+
  geom_label(data = time_labels, aes(label=paste0(run_time,' min')), 
             x = Inf, y = Inf, hjust=1, vjust=1,
             inherit.aes = FALSE)+
  scale_y_continuous(limits=c(0,max(c(data_sim_comptest3$upper,df_eir_all_final$EIR))))+
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  facet_grid(threads~particles, labeller=label_both)+
  labs(x='Time (days)',y='EIR')

ggplot(df_inc_all_final)+
  geom_line(aes(x=t,y=inc,group=variable),col = "#A6CEE3",alpha=0.2)+
  geom_point(data=data_sim_comptest3,aes(x=t,y=inc_true),pch = 19,
             col = "#666666")+
  geom_label(data = time_labels, aes(label=paste0(run_time,' min')), 
             x = Inf, y = Inf, hjust=1, vjust=1,
             inherit.aes = FALSE)+
  scale_y_continuous(limits=c(0,max(c(df_inc_all_final$inc))))+
  facet_grid(threads~particles, labeller=label_both)+
  labs(x='Time (days)',y='Incidence (<5 years)')


##Create figures for 1 simulated data run
##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("run_pmcmc.R",
             "model_parameters.R","equilibrium-init-create-stripped.R")


ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))
config_32 <- didehpc::didehpc_config(cores = 32, parallel = TRUE)
obj_32 <- didehpc::queue_didehpc(ctx,config = config_32)
obj_32$task_list()

context::context_info("T:/jth/contexts")

task_32_200 <- obj_32$task_get('9c173380b6a53529d08b70ddcd892807')
task_32_200$log()
result_32_200 <- task_32_200$result()
##Add prevalence confidence intervals
prev_cis <- addCIs(data_sim_comptest3,data_sim_comptest3$positive,data_sim_comptest3$tested)
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
