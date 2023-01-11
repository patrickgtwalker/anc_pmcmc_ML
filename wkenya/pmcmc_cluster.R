#Run loop function
setwd('Q:/mode-test-pw')
# drat:::add("mrc-ide") # install.packages("drat") if this errors, then try again
# install.packages("didehpc")
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

data_raw <- readRDS('data_raw_120822.RDS')
root <- "contexts"
sources <- c("run_pmcmc.R",
             "MiP-given/model_parameters.R","MiP-given/equilibrium-init-create-stripped.R")
config <- didehpc::didehpc_config(cores = 4, parallel = TRUE)


ctx <- context::context_save("contexts", sources = sources,
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

obj <- didehpc::queue_didehpc(ctx,config = config)

obj$cluster_load(TRUE)
obj$config
obj$login()
pmcmc_sim2 <- obj$enqueue(run_pmcmc(data=data_raw,EIR_vol=0.3,n_particles = 100))
pmcmc_sim <- obj$task_get("8681bf6c1d3fde4486ffd630f0229273")
pmcmc_sim2 <- obj$task_get("77442db45ba2aa8ca1c44dc8028d5cfe")
pmcmc_sim$id #"8681bf6c1d3fde4486ffd630f0229273"
pmcmc_sim$status()
pmcmc_sim$times()
pmcmc_sim$log()
pmcmc_sim2$id #"77442db45ba2aa8ca1c44dc8028d5cfe"
pmcmc_sim2$status()
pmcmc_sim2$times()
pmcmc_sim2$log()
obj$unsubmit("8681bf6c1d3fde4486ffd630f0229273")

data_raw_cmis <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/data_raw_cmis.RDS')
data_raw_cmis_all <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/data_raw_cmis_all.RDS')
pmcmc_cmis2 <- obj$enqueue(run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 10,n_particles = 100))
pmcmc_desktop_3 <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 20,n_particles = 100)
pmcmc_cmis_first <- obj$task_get('7e4d21482a7cf2bb1803414ca0606710')
pmcmc_cmis1 <- obj$task_get("817207f4d579cf9ab67e72bda6beffc6")
pmcmc_cmis1$id #"817207f4d579cf9ab67e72bda6beffc6"
pmcmc_cmis1$status()
pmcmc_cmis1$times()
pmcmc_cmis1$log()
pmcmc_cmis2 <- obj$task_get("381937e3430b312e5f7c27fc2394b8f3")
pmcmc_cmis2$id #"381937e3430b312e5f7c27fc2394b8f3"
pmcmc_cmis2$status()
pmcmc_cmis2$times()
pmcmc_cmis2$log()

pars <- data.frame(max_steps = c(2e4,4e4,1e5), atol = c(5e-6,10e-5,10e-4), rtol = c(5e-6,10e-5,10e-4))
cmis_bulk <- obj$enqueue_bulk(pars, run_pmcmc, data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 10,n_particles = 100)
cmis_bulk$status()

na <- 21
age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
fd <- vector(length=20L)
for(i in 1:20){
  fd[i] <- 1-(1-state$fD0)/(1+((age[i]+age[i+1])/2/state$aD)^state$gammaD)
}
p_det <- state$d1 + (1-state$d1)/(1 + fd[]*(state$init_ID[,]/state$ID0)^state$kD)
prev0to59 <- state$init_T[1:state$age59,] + state$init_D[1:state$age59,]  + state$init_A[1:state$age59,]*p_det[1:state$age59,]
prev <- sum(prev0to59[,])/sum(state$den[1:state$age59])


source("run_pmcmc.R")
source("MiP-given/model_parameters.R")
source("MiP-given/equilibrium-init-create-stripped.R")
pmcmc_desktop <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,n_particles = 100)
saveRDS(pmcmc_desktop,'./pmcmc-run/pmcmc_desktop_1.RDS')

proposal_dist <- cov(pmcmc_desktop_1$pars)
pmcmc_desktop_2 <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,proposal_dist=proposal_dist,n_particles = 100)
saveRDS(pmcmc_desktop_2,'./pmcmc-run/pmcmc_desktop_2.RDS')

pmcmc_desktop_3 <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 10,n_particles = 100, max_steps = 1e7,
                             atol = 1e-3,
                             rtol = 1e-3)
saveRDS(pmcmc_desktop_3,'./pmcmc-run/pmcmc_desktop_3.RDS')

pmcmc_desktop_4 <- run_pmcmc(data=data_raw_cmis,EIR_vol=0.3,init_EIR_user = 5,n_particles = 100,
                             max_steps = 1e8,
                             atol = 1e-2,
                             rtol = 1e-2)
data_raw_cmis_all <- readRDS('./pmcmc-run/data_raw_cmis_all.RDS')
pmcmc_desktop_all <- run_pmcmc(data=data_raw_cmis_all,EIR_vol=0.3,init_EIR = 5,n_particles = 100,
                             max_steps = 1e8,
                             atol = 1e-2,
                             rtol = 1e-2)
saveRDS(pmcmc_desktop_all,'./pmcmc-run/pmcmc_desktop_all_EIR.RDS')
pmcmc_desktop_exp <- run_pmcmc(data=data_raw_cmis_all,EIR_vol=0.3,init_EIR = 5,n_particles = 100,
                               max_steps = 1e8,
                               atol = 1e-2,
                               rtol = 1e-2)
saveRDS(pmcmc_desktop_exp,'./pmcmc-run/pmcmc_desktop_all_exp.RDS')

pmcmc_desktop_exp <- run_pmcmc(data=data_raw_cmis_all,EIR_vol=0.1,init_EIR = 5,n_particles = 200,
                               max_steps = 1e6,
                               atol = 1e-6,
                               rtol = 1e-6)
saveRDS(pmcmc_desktop_exp,'./pmcmc-run/pmcmc_desktop_all_exp.RDS')
saveRDS(pmcmc_desktop_exp,'./pmcmc-run/cmis_good1.RDS')


data_raw_ng <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_NG_anc.RDS')
ng_test <- run_pmcmc(data=data_raw_ng,EIR_vol=0.3,init_EIR = 40,n_particles = 200,
                     max_steps = 1e6,
                     atol = 1e-4,
                     rtol = 1e-4)
saveRDS(ng_test,'ng_test.RDS')

data_raw_bf <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_BF_anc.RDS')
bf_test <- run_pmcmc(data=data_raw_bf,EIR_vol=0.1,init_EIR = 10,n_particles = 100,
                     max_steps = 1e6,
                     atol = 1e-4,
                     rtol = 1e-4)
saveRDS(bf_test,'bf_test.RDS')

data_sim <- readRDS('data_sim_180822.RDS')

data_raw <- data_gen(EIR_volatility = 0.6, init_EIR = 20)
saveRDS(data_raw,'data_sim2.RDS')

sim_test <- run_pmcmc(data=data_raw,EIR_vol=0.1,init_EIR = 20,n_particles = 100,
                     max_steps = 1e6,
                     atol = 1e-4,
                     rtol = 1e-4)
saveRDS(sim_test,'sim_test.RDS')
proposal <- matrix(c(0.0336,-0.000589,-0.000589,0.04942),ncol=2)
source("run_pmcmc.R")

pmcmc_eq <- run_pmcmc(data=data_raw_cmis_all,
                      # EIR_vol=0.3,
                      # init_EIR = 5,
                      n_particles = 200,
                      max_steps = 1e6,
                      atol = 1e-6,
                      rtol = 1e-6,
                      proposal_matrix = proposal,
                      n_steps = 1000)
saveRDS(pmcmc_eq,'cmis_all_1.RDS') #good run
eq_mcmc <- coda::as.mcmc(cbind(pmcmc_eq$probabilities, pmcmc_eq$pars))
1 - coda::rejectionRate(eq_mcmc)
coda::effectiveSize(eq_mcmc)
proposal <-  cov(pmcmc_eq$pars)
summary(eq_mcmc)
windows(60,50)
plot(eq_mcmc)
history2 <- pmcmc_eq$trajectories$state
plot_particle_filter(history2,true_history=data_raw_cmis_all,times=data_raw_cmis_all$t)
matplot(data_raw_cmis_all$t, t(history2[2, , -1]), type = "l",
        xlab = "Time", ylab = "EIR",
        col = "#A6CEE3", lty = 1, ylim = range(history2[2, , -1]))

matplot(data_raw_cmis_all$t, t(history[2, , -1]), type = "l",
        xlab = "Time", ylab = "EIR",
        col = "#A6CEE3", lty = 1, ylim = range(c(0,100)))

root <- "contexts"
sources <- c("run_pmcmc.R",
             "MiP-given/model_parameters.R","MiP-given/equilibrium-init-create-stripped.R")
config <- didehpc::didehpc_config(cores = 4, parallel = TRUE)


ctx <- context::context_save("contexts", sources = sources,
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

obj <- didehpc::queue_didehpc(ctx,config = config)

obj$cluster_load(TRUE)
obj$config
obj$login()
pmcmc_chain1 <- obj$enqueue(run_pmcmc(data=data_raw_cmis_all,
                                    # EIR_vol=0.3,
                                    # init_EIR = 5,
                                    n_particles = 200,
                                    max_steps = 1e6,
                                    atol = 1e-6,
                                    rtol = 1e-6,
                                    proposal_matrix = proposal,
                                    n_steps = 1000))
pmcmc_chain1$id
pmcmc_chain1$status()
pmcmc_chain1$times()
pmcmc_chain1$log()

pmcmc_chain2 <- obj$enqueue(run_pmcmc(data=data_raw_cmis_all,
                                      # EIR_vol=0.3,
                                      # init_EIR = 5,
                                      n_particles = 200,
                                      max_steps = 1e6,
                                      atol = 1e-6,
                                      rtol = 1e-6,
                                      proposal_matrix = proposal,
                                      n_steps = 1000))
pmcmc_chain2$id
pmcmc_chain2$status()
pmcmc_chain2$times()
pmcmc_chain2$log()

pmcmc_chain3 <- obj$enqueue(run_pmcmc(data=data_raw_cmis_all,
                                      # EIR_vol=0.3,
                                      # init_EIR = 5,
                                      n_particles = 200,
                                      max_steps = 1e6,
                                      atol = 1e-6,
                                      rtol = 1e-6,
                                      proposal_matrix = proposal,
                                      n_steps = 1000))
pmcmc_chain3$id
pmcmc_chain3$status()
pmcmc_chain3$times()
pmcmc_chain3$log()

##Get results from run 1 and 3
pmcmc_run1 <- pmcmc_chain1$result()
pmcmc_run2 <- pmcmc_chain2$result()
pmcmc_run3 <- pmcmc_chain3$result()

run1_mcmc <- coda::as.mcmc(cbind(pmcmc_run1$probabilities, pmcmc_run1$pars))
1 - coda::rejectionRate(run1_mcmc)
coda::effectiveSize(run1_mcmc)
proposal <-  cov(pmcmc_run1$pars)
summary(run1_mcmc)
windows(60,50)
plot(run1_mcmc)
history_run1 <- pmcmc_run1$trajectories$state
plot_particle_filter(history_run1,true_history=data_raw_cmis_all,times=data_raw_cmis_all$t)
matplot(data_raw_cmis_all$t, t(history_run1[2, , -1]), type = "l",
        xlab = "Time", ylab = "EIR",
        col = "#A6CEE3", lty = 1, ylim = range(history_run1[2, , -1]))
saveRDS(pmcmc_run1,'cmis_all_run1.RDS') #good run

run2_mcmc <- coda::as.mcmc(cbind(pmcmc_run2$probabilities, pmcmc_run2$pars))
1 - coda::rejectionRate(run2_mcmc)
coda::effectiveSize(run2_mcmc)
proposal <-  cov(pmcmc_run2$pars)
summary(run2_mcmc)
windows(60,50)
plot(run2_mcmc)
history_run2 <- pmcmc_run2$trajectories$state
plot_particle_filter(history_run2,true_history=data_raw_cmis_all,times=data_raw_cmis_all$t)
matplot(data_raw_cmis_all$t, t(history_run2[2, , -1]), type = "l",
        xlab = "Time", ylab = "EIR",
        col = "#A6CEE3", lty = 1, ylim = range(history_run2[2, , -1]))
saveRDS(pmcmc_run2,'cmis_all_run2.RDS') #good run


run3_mcmc <- coda::as.mcmc(cbind(pmcmc_run3$probabilities, pmcmc_run3$pars))
1 - coda::rejectionRate(run3_mcmc)
coda::effectiveSize(run3_mcmc)
proposal <-  cov(pmcmc_run3$pars)
summary(run3_mcmc)
windows(60,50)
plot(run3_mcmc)
history_run3 <- pmcmc_run3$trajectories$state
plot_particle_filter(history_run3,true_history=data_raw_cmis_all,times=data_raw_cmis_all$t)
matplot(data_raw_cmis_all$t, t(history_run3[2, , -1]), type = "l",
        xlab = "Time", ylab = "EIR",
        col = "#A6CEE3", lty = 1, ylim = range(history_run3[2, , -1]))
saveRDS(pmcmc_run3,'cmis_all_run3.RDS') #good run

all_runs <- as.mcmc(rbind(run1_mcmc[101:1000,],run2_mcmc[101:1000,],run3_mcmc[101:1000,]))
1 - coda::rejectionRate(all_runs)
coda::effectiveSize(all_runs)
plot(all_runs)
raft <- raftery.diag(run1_mcmc)
all_chains <- mcmc.list(run1_mcmc,run2_mcmc,run3_mcmc)
run1_mcmc_proc <- coda::as.mcmc(cbind(pmcmc_run1$probabilities[51:1000,], pmcmc_run1$pars[51:1000,]))
run2_mcmc_proc <- coda::as.mcmc(cbind(pmcmc_run2$probabilities[51:1000,], pmcmc_run2$pars[51:1000,]))
run3_mcmc_proc <- coda::as.mcmc(cbind(pmcmc_run3$probabilities[51:1000,], pmcmc_run3$pars[51:1000,]))
all_chains_proc <- mcmc.list(run1_mcmc_proc,run2_mcmc_proc,run3_mcmc_proc)

1 - coda::rejectionRate(all_chains_proc)
coda::effectiveSize(all_chains_proc)
plot(all_chains_proc)
summary(all_chains)
color_scheme_set("mix-pink-blue")
windows(60,50)
((bayesplot::mcmc_trace(all_chains_proc,pars = 'log_prior')+mcmc_dens_overlay(all_chains_proc,pars = 'log_prior'))/
  (bayesplot::mcmc_trace(all_chains_proc,pars = 'log_likelihood')+mcmc_dens_overlay(all_chains_proc,pars = 'log_likelihood'))/
  (bayesplot::mcmc_trace(all_chains_proc,pars = 'log_posterior')+mcmc_dens_overlay(all_chains_proc,pars = 'log_posterior'))/
  (bayesplot::mcmc_trace(all_chains_proc,pars = 'EIR_SD')+mcmc_dens_overlay(all_chains_proc,pars = 'EIR_SD'))/
  (bayesplot::mcmc_trace(all_chains_proc,pars = 'log_init_EIR')+mcmc_dens_overlay(all_chains_proc,pars = 'log_init_EIR'))) + plot_layout(guides = "collect") #& theme(legend.position = "bottom")
plot_particle_filter(history_run3,true_history=data_raw_cmis_all,times=data_raw_cmis_all$t)


ages_run1 <- history_run1[4:24,,-1]
#[agegroup,mcmc sample,time]
history_run1[4,101:1000,]
df_ages <- data.frame(agegroup=numeric(),
                      month=numeric(),
                      value=numeric())
for(age in 4:24){
  for(time in 2:62){
    temp <- data.frame(agegroup=age-3,
                       month=time,
                       value=history_run1[age,101:1000,time])
    df_ages <- rbind(df_ages,temp)
  }
}
windows(160,90)
ggplot(df_ages,aes(x=as.factor(agegroup),y=value))+
  geom_boxplot()+
  facet_wrap(~month)

ggplot(df_ages[df_ages$month<=13,],aes(x=as.factor(agegroup),y=value))+
  geom_boxplot()+
  facet_wrap(~month)


##PMCMC cluster for fitting prev by age
byage_1 <- obj$enqueue(run_pmcmc(data=data_raw_cmis_byage,proposal_matrix=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            n_particles=200,
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000))

byage_1$status()
byage_1$log()

byage_2 <- obj$enqueue(run_pmcmc(data=data_raw_cmis_byage,proposal_matrix=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                 n_particles=200,
                                 max_EIR=1000,
                                 max_steps = 1e7,
                                 atol = 1e-5,
                                 rtol = 1e-6,
                                 n_steps = 1000))
byage_2$status()
byage_2$log()

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("run_pmcmc.R",
             "MiP-given/model_parameters.R","MiP-given/equilibrium-init-create-stripped.R")


ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))
config_32 <- didehpc::didehpc_config(cores = 32, parallel = TRUE)
obj_32 <- didehpc::queue_didehpc(ctx,config = config_32)

obj_32$cluster_load(TRUE)
obj_32$config
obj_32$login()

byage_3 <- obj_32$enqueue(run_pmcmc(data=ANC_cMIS_for_pmcmc,proposal_matrix=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                 n_particles=200,
                                 max_EIR=1000,
                                 max_steps = 1e7,
                                 atol = 1e-5,
                                 rtol = 1e-6,
                                 n_steps = 1000,
                                 n_threads = 32))
byage_3$status()
byage_3$log()

byage_4 <- obj_32$enqueue(run_pmcmc(data=ANC_cMIS_for_pmcmc,proposal_matrix=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                    n_particles=200,
                                    max_EIR=1000,
                                    max_steps = 1e7,
                                    atol = 1e-5,
                                    rtol = 1e-6,
                                    n_steps = 1000,
                                    n_threads = 32))
byage_4$status()
byage_4$log()
byage_4$times()

pmcmc_byage3 <- byage_3$result()
saveRDS(pmcmc_byage3,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/pmcmc/pmcmc_byage3_ANCcMIS_270922.rds')

byage3_mcmc <- coda::as.mcmc(cbind(pmcmc_byage3$probabilities, pmcmc_byage3$pars))


windows(60,50)
1 - coda::rejectionRate(byage3_mcmc)
coda::effectiveSize(byage3_mcmc)
plot(byage3_mcmc)
summary(byage3_mcmc)
history_byage3 <- pmcmc_byage3$trajectories$state

pmcmc_byage4 <- byage_4$result()
saveRDS(pmcmc_byage4,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/pmcmc/pmcmc_byage4_ANCcMIS_270922.rds')

byage4_mcmc <- coda::as.mcmc(cbind(pmcmc_byage4$probabilities, pmcmc_byage4$pars))

windows(60,50)
1 - coda::rejectionRate(byage4_mcmc)
coda::effectiveSize(byage4_mcmc)
plot(byage4_mcmc)
summary(byage4_mcmc)
history_byage4 <- pmcmc_byage4$trajectories$state

color_scheme_set("mix-pink-blue")

##Extract run trajectories and remove burnin (50 steps)
eir_history_byage3 <- data.frame(t(history_byage3['EIR', , -1]))
long_eir_history_byage3 <- eir_history_byage3%>% mutate(t=c(1:nrow(eir_history_byage3)))%>%
  melt(id='t')

ggplot(long_eir_history_byage3)+
  geom_line(aes(x=t,y=value,group=variable),col = "#A6CEE3",alpha=1)+
  # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')

eir_history_byage3 <- data.frame(t(history_byage3['EIR', 51:1000, -1]))
long_eir_history_byage3 <- eir_history_byage3%>% mutate(t=c(1:nrow(eir_history_byage3)))%>%
  melt(id='t')

ggplot(long_eir_history_byage3)+
  geom_line(aes(x=t,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')

eir_history_byage4 <- data.frame(t(history_byage4['EIR', 51:1000, -1]))
long_eir_history_byage4 <- eir_history_byage4%>% mutate(t=c(1:nrow(eir_history_byage4)))%>%
  melt(id='t')

ggplot(long_eir_history_byage4)+
  geom_line(aes(x=t,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')

history_byage4_trunc <- history_byage3[,,]
