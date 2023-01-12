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


source("run_pmcmc.R")
source("shared/model_parameters.R")
source("shared/equilibrium-init-create-stripped.R")

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("run_pmcmc.R",
             "shared/model_parameters.R","shared/equilibrium-init-create-stripped.R")

ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))
config_32 <- didehpc::didehpc_config(cores = 32, parallel = TRUE)
obj_32 <- didehpc::queue_didehpc(ctx,config = config_32)

obj_32$cluster_load(TRUE)
obj_32$config
obj_32$login()

byage <- obj_32$enqueue(run_pmcmc(data=ANC_cMIS_for_pmcmc,proposal_matrix=matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                 n_particles=200,
                                 max_EIR=1000,
                                 max_steps = 1e7,
                                 atol = 1e-5,
                                 rtol = 1e-6,
                                 n_steps = 1000,
                                 n_threads = 32))
byage$status()
byage$log()

pmcmc_byage <- byage$result()
saveRDS(pmcmc_byage,'C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/pmcmc/pmcmc_byage_ANCcMIS_270922.rds')

byage_mcmc <- coda::as.mcmc(cbind(pmcmc_byage$probabilities, pmcmc_byage$pars))


windows(60,50)
color_scheme_set("mix-pink-blue")
1 - coda::rejectionRate(byage_mcmc)
coda::effectiveSize(byage_mcmc)
plot(byage_mcmc)
summary(byage_mcmc)
history_byage <- pmcmc_byage$trajectories$state


##Extract run trajectories and remove burnin (50 steps)
eir_history_byage <- data.frame(t(history_byage['EIR', , -1]))
long_eir_history_byage <- eir_history_byage%>% mutate(t=c(1:nrow(eir_history_byage)))%>%
  melt(id='t')

ggplot(long_eir_history_byage)+
  geom_line(aes(x=t,y=value,group=variable),col = "#A6CEE3",alpha=1)+
  # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')

eir_history_byage <- data.frame(t(history_byage['EIR', 51:1000, -1]))
long_eir_history_byage <- eir_history_byage%>% mutate(t=c(1:nrow(eir_history_byage)))%>%
  melt(id='t')

ggplot(long_eir_history_byage)+
  geom_line(aes(x=t,y=value,group=variable),col = "#A6CEE3",alpha=0.2)+
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  labs(x='Time (days)',y='EIR')

