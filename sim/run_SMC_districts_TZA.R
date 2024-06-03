library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library('mcstate')
#library(didehpc)
library(tidyverse)
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
source('shared/SMC_functions.R')
### RUN MODEL WITH SMC ####

itigi_draws<-as.data.frame(readRDS("./smc_stuff/itigi_betaa_4pw_corrected.rds"))
itigi_sims_step<-get_smc_sims(0.2840339,itigi_draws,1840-365,4,30,4,7,40)

itigi_summary<-itigi_sims_step%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))
ggplot(itigi_summary,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  geom_ribbon(fill="#FF410DFF",lwd=2,aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30))+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  geom_ribbon(fill="#6EE2FFFF",lwd=2,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30))+
  # geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per month")

Malinyi_draws<-as.data.frame(readRDS("./smc_stuff/malinyi_betaa_4pw.rds"))
Malinyi_sims_step<-get_smc_sims(2.207521,Malinyi_draws,1840-365,4,30,4,7,40)

Malinyi_summary<-Malinyi_sims_step%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))
date_zero<-2016
ggplot(Malinyi_summary,aes(x=date_zero+t-30,y=incu5_med*1000*30))+
  geom_line(col="#FF410DFF",lwd=2)+
  geom_ribbon(fill="#FF410DFF",lwd=2,aes(ymin=incu5_lower*1000*30,ymax=incu5_upper*1000*30))+
  # geom_line(data=median_inc_05,aes(y=med_inc05*1000*30),lwd=3,col="#FF410DFF")+
  geom_vline(xintercept=date_zero,lty=2)+
  geom_line(aes(y=incu5_smc_med*1000*30),col="#6EE2FFFF",lwd=2)+
  geom_ribbon(fill="#6EE2FFFF",lwd=2,aes(ymin=incu5_smc_lower*1000*30,ymax=incu5_smc_upper*1000*30))+
 # geom_line(data=median_inc_05,aes(y=med_inc05_smc*1000*30),lwd=3,col="#6EE2FFFF")+
  ylab("Cases per 1000 children under 5 per month")


write_csv(Malinyi_summary,"./smc_stuff/malinyi_smc.csv")
write_csv(itigi_summary,"./smc_stuff/itigi_smc.csv")



ghana_draws<-as.data.frame(readRDS("./smc_stuff/ghana_betaa_flat3.rds"))

#%>%
 # pivot_longer(!c(date,t),names_to = "draw",values_to="Emergence_rate")



ghana_draws$mean=rowMeans(ghana_draws[,-c(1:2)])
ghana_draws$upper=row(ghana_draws[,-c(1:2)])
ghana_quants=as.data.frame(t(apply(ghana_draws[,-c(1:2)], 1, quantile, probs = c(0.2, 0.5, 0.8))))

ghana_quants$date<-ghana_draws$date
names(ghana_quants)[1:3]<-c("lower","med","upper")
tiff("./smc_stuff/Navrongo_mosq.tif",height=7,width=15,unit="cm",res=300)
ggplot(ghana_draws,aes(x=date,y=mean))+
  geom_line(size=2,color="purple")+ 
  #geom_ribbon(aes(ymin=lower,ymax=upper),fill="purple",alpha=0.2)+
  ylab("Mosquito emergence rate")+
  scale_x_date(labels = scales::date_format("%b-%Y"),breaks = seq(min(ghana_draws$date), max(ghana_draws$date), by = "2 months"))
dev.off()

row
ghana_summary_step_2yr<-ghana_draws%>%
  group_by(t)%>%
  summarise(incu5_med=quantile(inc_target,0.5),incu5_lower=quantile(inc_target,0.025),incu5_upper=quantile(inc_target,0.975),
            incu5_smc_med=quantile(inc_target_smc,0.5),incu5_smc_lower=quantile(inc_target_smc,0.025),incu5_smc_upper=quantile(inc_target_smc,0.975))