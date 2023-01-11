########################################
####### Function to plot particle trajectories
addCIs<-function(df,Ys,Ns){
  df$mean<-NA
  df$upper<-NA
  df$lower<-NA
  CIs<-binom.confint(Ys,Ns,method="exact")
  df$mean[Ns>0]<-CIs$mean[Ns>0]
  df$upper[Ns>0]<-CIs$upper[Ns>0]
  df$lower[Ns>0]<-CIs$lower[Ns>0]
  return(df)
}

##Figures##
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(zoo)
library(patchwork)
library(coda)
library(dplyr)
library(binom)
RColorBrewer::display.brewer.all(colorblindFriendly = TRUE)
RColorBrewer::display.brewer.pal(n=3,name='RdYlBu')
brewer.pal(n = 2, name = "RdYlBu")

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))


pmcmc_byage3 <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/pmcmc/pmcmc_byage3_ANCcMIS_270922.rds')
pmcmc_byage4 <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/WesternKenya/pmcmc/pmcmc_byage4_ANCcMIS_270922.rds')

byage3_mcmc <- coda::as.mcmc(cbind(pmcmc_byage3$probabilities, pmcmc_byage3$pars))

windows(60,50)
##Run 3 diagnostics
1 - coda::rejectionRate(byage3_mcmc)
coda::effectiveSize(byage3_mcmc)
plot(byage3_mcmc)
summary(byage3_mcmc)
history_byage3 <- pmcmc_byage3$trajectories$state

#Run 4 diagnostics
1 - coda::rejectionRate(byage4_mcmc)
coda::effectiveSize(byage4_mcmc)
plot(byage4_mcmc)
summary(byage4_mcmc)
history_byage4 <- pmcmc_byage4$trajectories$state

##Combine runs with a 50 step burnin
pmcmc_byage3_proc <- coda::as.mcmc(cbind(pmcmc_byage3$probabilities[51:1000,], pmcmc_byage3$pars[51:1000,]))
pmcmc_byage4_proc <- coda::as.mcmc(cbind(pmcmc_byage4$probabilities[51:1000,], pmcmc_byage4$pars[51:1000,]))
all_chains_proc <- mcmc.list(pmcmc_byage3_proc,pmcmc_byage4_proc)

##Diagnostics of combined runs
1 - coda::rejectionRate(all_chains_proc)
coda::effectiveSize(all_chains_proc)
plot(all_chains_proc)
summary(all_chains)
color_scheme_set("mix-pink-blue")
((bayesplot::mcmc_trace(all_chains_proc,pars = 'log_prior')+mcmc_dens_overlay(all_chains_proc,pars = 'log_prior'))/
    (bayesplot::mcmc_trace(all_chains_proc,pars = 'log_likelihood')+mcmc_dens_overlay(all_chains_proc,pars = 'log_likelihood'))/
    (bayesplot::mcmc_trace(all_chains_proc,pars = 'log_posterior')+mcmc_dens_overlay(all_chains_proc,pars = 'log_posterior'))/
    (bayesplot::mcmc_trace(all_chains_proc,pars = 'EIR_SD')+mcmc_dens_overlay(all_chains_proc,pars = 'EIR_SD'))/
    (bayesplot::mcmc_trace(all_chains_proc,pars = 'log_init_EIR')+mcmc_dens_overlay(all_chains_proc,pars = 'log_init_EIR'))) + plot_layout(guides = "collect") #& theme(legend.position = "bottom")


##Plot EIR
eir_history_byage3 <- data.frame(t(history_byage3['EIR', 51:1000, -1]))
long_eir_history_byage3 <- eir_history_byage3%>% mutate(t=c(1:nrow(eir_history_byage3)))%>%
  melt(id='t')%>%
  mutate(chain=1,
         time=t)
eir_history_byage4 <- data.frame(t(history_byage4['EIR', 51:1000, -1]))
long_eir_history_byage4 <- eir_history_byage4%>% mutate(t=c(1:nrow(eir_history_byage4)))%>%
  melt(id='t')%>%
  mutate(chain=2,
         time=t)

##Create upper and lower bounds by time to make a ribbon
dates <- seq(as.Date("2015-4-1"), as.Date("2022-2-1"), by = "months")
eir_history_both_notgrouped <- rbind(long_eir_history_byage3,long_eir_history_byage4)%>%
  mutate(epoch=ifelse(time>61,'Without cMIS','With cMIS'))

eir_history_both <- rbind(long_eir_history_byage3,long_eir_history_byage4)%>%
  group_by(time)%>%
  summarise(median=median(value),
            mean=mean(value),
            upper=quantile(value,probs=0.975),
            lower=quantile(value,probs=0.025))%>%
  mutate(epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates)
# eir_history_chains <- rbind(long_eir_history_byage3,long_eir_history_byage4)%>%
#   group_by(time,chain)%>%
#   summarise(median=median(value),
#             mean=mean(value),
#             upper=quantile(value,probs=0.975),
#             lower=quantile(value,probs=0.025))%>%
#   mutate(epoch=ifelse(time>61,'Without cMIS','With cMIS'))

# (ggplot(eir_history_both_notgrouped)+
#   geom_line(aes(x=t,y=value,group=c(variable,chain)),col = "#A6CEE3",alpha=0.2)+
#   scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
#   labs(x='Time (days)',y='EIR'))
colors <- c('#0571b0','#ca0020')
datebreaks <- as.Date(c('2015-4-1','2016-1-1','2017-1-1','2018-1-1','2019-1-1','2020-1-1','2021-1-1','2022-1-1'))
windows(7,5)
eir_plot <- ggplot(eir_history_both)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='EIR')+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1))
eir_plot
eir_plot_4comp <- ggplot(eir_history_both)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='EIR')+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
eir_plot_4comp
eir_plot_4comp_nollin <- ggplot(eir_history_both)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='EIR')+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
eir_plot_4comp_nollin

eir_plot_4comp_onlyllin <- ggplot(eir_history_both)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0)+
  geom_line(aes(x=month,y=median,color=epoch),size=1,alpha=0)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='EIR')+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
eir_plot_4comp_onlyllin

# ggplot(eir_history_chains)+
#   geom_ribbon(aes(x=time,ymin=lower,ymax=upper,fill=epoch),alpha=0.5)+
#   geom_line(aes(x=time,y=median,color=epoch),size=1)+
#   scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
#   facet_wrap(~as.factor(chain))

###Plot Incidence
inc_history_byage3 <- data.frame(t(history_byage3['inc', 51:1000, -1]))
long_inc_history_byage3 <- inc_history_byage3%>% mutate(t=c(1:nrow(inc_history_byage3)))%>%
  melt(id='t')%>%
  mutate(chain=1,
         time=t)
inc_history_byage4 <- data.frame(t(history_byage4['inc', 51:1000, -1]))
long_inc_history_byage4 <- inc_history_byage4%>% mutate(t=c(1:nrow(inc_history_byage4)))%>%
  melt(id='t')%>%
  mutate(chain=2,
         time=t)
inc_history_both <- rbind(long_inc_history_byage3,long_inc_history_byage4)%>%
  group_by(time)%>%
  summarise(median=median(value),
            mean=mean(value),
            upper=quantile(value,probs=0.975),
            lower=quantile(value,probs=0.025))%>%
  mutate(epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates)
inc_plot <- ggplot(inc_history_both)+
  geom_ribbon(aes(x=month,ymin=lower*10000,ymax=upper*10000,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median*10000,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='Incidence in <5 year olds\n(per 10,000 persons)')+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
inc_plot
inc_plot_nollin <- ggplot(inc_history_both)+
  geom_ribbon(aes(x=month,ymin=lower*10000,ymax=upper*10000,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median*10000,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='Incidence in <5 year olds\n(per 10,000 persons)')+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
inc_plot_nollin

inc_plot_onlyllin <- ggplot(inc_history_both)+
  geom_ribbon(aes(x=month,ymin=lower*10000,ymax=upper*10000,fill=epoch),alpha=0)+
  geom_line(aes(x=month,y=median*10000,color=epoch),size=1,alpha=0)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='Incidence in <5 year olds\n(per 10,000 persons)')+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
inc_plot_onlyllin

eir_plot_4comp/inc_plot + plot_layout(guides = "collect") 
eir_plot_4comp_onlyllin/inc_plot_onlyllin + plot_layout(guides = "collect") 
eir_plot_4comp_nollin/inc_plot_nollin + plot_layout(guides = "collect") 

##Plot fitted cMIS prevalence all ages
data_raw_cmis_all <- readRDS('./pmcmc-run/data_raw_cmis_all.RDS')
cmis_data_cis <- addCIs(data_raw_cmis_all,data_raw_cmis_all$positive,data_raw_cmis_all$tested)
dates <- seq(as.Date("2015-4-1"), as.Date("2022-2-1"), by = "months")
cmis_data_cis$month <- dates[1:61]

prev_history_byage3 <- data.frame(t(history_byage3['prev', 51:1000, -1]))
long_prev_history_byage3 <- prev_history_byage3%>% mutate(t=c(1:nrow(prev_history_byage3)))%>%
  melt(id='t')%>%
  mutate(chain=1,
         time=t)
prev_history_byage4 <- data.frame(t(history_byage4['prev', 51:1000, -1]))
long_prev_history_byage4 <- prev_history_byage4%>% mutate(t=c(1:nrow(prev_history_byage4)))%>%
  melt(id='t')%>%
  mutate(chain=2,
         time=t)
prev_history_both <- rbind(long_prev_history_byage3,long_prev_history_byage4)%>%
  group_by(time)%>%
  summarise(median=median(value),
            mean=mean(value),
            upper=quantile(value,probs=0.975),
            lower=quantile(value,probs=0.025))%>%
  mutate(epoch=ifelse(time>61,'Without cMIS','With cMIS'),
         month=dates)

prev_plot <- ggplot(prev_history_both)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_data_cis,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_data_cis,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='Prevalence in cMIS\n(All Ages)')+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
prev_plot
prev_plot_4comp <- ggplot(prev_history_both)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_data_cis,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_data_cis,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='Prevalence in cMIS\n(All Ages)')+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
prev_plot_4comp
prev_plot_4comp/eir_plot_4comp/inc_plot + plot_layout(guides = "collect") 

# plot_particle_filter(history_byage3[, 51:1000,],true_history = data_raw_cmis_all, times=data_raw_cmis_all$t, obs_end = 83*30)

##Plot prevalence by age groups
cmisage_history_byage3_age1 <- data.frame(t(history_byage3[4, 51:1000, -1]))
long_cmisage_history_byage3_age1 <- cmisage_history_byage3_age1%>% mutate(t=c(1:nrow(cmisage_history_byage3_age1)))%>%
  melt(id='t')%>%
  mutate(chain=1,
         time=t)%>%
  group_by(time)%>%
  summarise(median=median(value),
            mean=mean(value),
            upper=quantile(value,probs=0.975),
            lower=quantile(value,probs=0.025))%>%
  mutate(epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates)

# cmisage_history_byage4_age1 <- data.frame(t(history_byage4[4, 51:1000, -1]))
# long_cmisage_history_byage4_age1 <- cmisage_history_byage4_age1%>% mutate(t=c(1:nrow(cmisage_history_byage4_age1)))%>%
#   melt(id='t')%>%
#   mutate(chain=1,
#          time=t)
# cmis_age1_history_both <- rbind(long_cmisage_history_byage3_age1,long_cmisage_history_byage4_age1)%>%
#   group_by(time)%>%
#   summarise(median=median(value),
#             mean=mean(value),
#             upper=quantile(value,probs=0.975),
#             lower=quantile(value,probs=0.025))%>%
#   mutate(epoch=ifelse(time>61,'Without cMIS','With cMIS'),
#          month=dates)
data_raw_cmis_byage <- readRDS('data_raw_cmis_byage.RDS')
cmis_data_cis_age1 <- addCIs(data_raw_cmis_byage[,c('t','tested_1','positive_1')],data_raw_cmis_byage$positive_1,data_raw_cmis_byage$tested_1)
dates <- seq(as.Date("2015-4-1"), as.Date("2022-2-1"), by = "months")
cmis_data_cis_age1$month <- dates[1:61]

cmis_age1prev_plot <- ggplot(long_cmisage_history_byage3_age1)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_data_cis_age1,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_data_cis_age1,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='<5 years\nPrevalence')+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
cmis_age1prev_plot
source('re_cat_cmis.R')
cmis_age2_recat <- re_cat_cmis(history_byage3,anc_index_min = 2, anc_index_max = 4)
cmis_age2_history <- cmis_age2_recat[[1]]%>%
  group_by(time)%>%
  dplyr::summarise(median=median(value),
                   mean=mean(value),
                   upper=quantile(value,probs=0.975,na.rm=TRUE),
                   lower=quantile(value,probs=0.025,na.rm=TRUE))%>%
  mutate(epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates)

cmis_age2_observed <- cmis_age2_recat[[2]]%>%
  mutate(time=1:61,
         epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates[1:61])


cmis_age2prev_plot <- ggplot(cmis_age2_history)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_age2_observed,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_age2_observed,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='5-12 years\nPrevalence')+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
cmis_age2prev_plot

# cmisage_history_byage3_age2 <- data.frame(t(history_byage3[5, 51:1000, -1]))
# long_cmisage_history_byage3_age2 <- cmisage_history_byage3_age2%>% mutate(t=c(1:nrow(cmisage_history_byage3_age2)))%>%
#   melt(id='t')%>%
#   mutate(chain=1,
#          time=t)
# cmisage_history_byage4_age2 <- data.frame(t(history_byage4[5, 51:1000, -1]))
# long_cmisage_history_byage4_age2 <- cmisage_history_byage4_age2%>% mutate(t=c(1:nrow(cmisage_history_byage4_age2)))%>%
#   melt(id='t')%>%
#   mutate(chain=1,
#          time=t)
# cmis_age2_history_both <- rbind(long_cmisage_history_byage3_age2,long_cmisage_history_byage4_age2)%>%
#   group_by(time)%>%
#   summarise(median=median(value),
#             mean=mean(value),
#             upper=quantile(value,probs=0.975),
#             lower=quantile(value,probs=0.025))%>%
#   mutate(epoch=ifelse(time>61,'Without cMIS','With cMIS'),
#          month=dates)
# 
# cmis_data_cis_age2 <- addCIs(data_raw_cmis_byage[,c('t','tested_2','positive_2')],data_raw_cmis_byage$positive_2,data_raw_cmis_byage$tested_2)
# dates <- seq(as.Date("2015-4-1"), as.Date("2022-2-1"), by = "months")
# cmis_data_cis_age2$month <- dates[1:61]
# 
# cmis_age2prev_plot <- ggplot(cmis_age2_history_both)+
#   geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
#   geom_line(aes(x=month,y=median,color=epoch),size=1)+
#   geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
#   geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
#   geom_point(data=cmis_data_cis_age2,aes(x=month,y=mean),pch = 19,
#              col = "#666666")+
#   geom_errorbar(data=cmis_data_cis_age2,aes(x=month,ymin=lower,ymax=upper),width = 0,
#                 col = "#666666")+
#   scale_color_manual(values=colors)+
#   scale_fill_manual(values=colors)+
#   scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
#   scale_y_continuous(limits=c(0,1))+
#   labs(x='Date',y='Prevalence in cMIS\n(5-10 years old)')+
#   theme(legend.title = element_blank(),
#         axis.text.x=element_blank(),
#         axis.title.x = element_blank(),
#         axis.ticks.x = element_line(size = 0.5), 
#         axis.ticks.length = unit(3, "pt"))
# cmis_age2prev_plot

##Recategorize middle age groups in cmis data and plot by age group
cmis_cba_recat <- re_cat_cmis(history_byage3)

cmis_cba_history <- cmis_cba_recat[[1]]%>%
  group_by(time)%>%
  dplyr::summarise(median=median(value),
                   mean=mean(value),
                   upper=quantile(value,probs=0.975,na.rm=TRUE),
                   lower=quantile(value,probs=0.025,na.rm=TRUE))%>%
  mutate(epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates)

cmis_cba_observed <- cmis_cba_recat[[2]]%>%
  mutate(time=1:61,
         epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates[1:61])


cmis_agecbaprev_plot <- ggplot(cmis_cba_history)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_cba_observed,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_cba_observed,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='13-44 years\nPrevalence')+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
cmis_agecbaprev_plot

cmisage_history_byage3_agelast <- data.frame(t(history_byage3[28, 51:1000, -1]))
long_cmisage_history_byage3_agelast <- cmisage_history_byage3_agelast%>% mutate(t=c(1:nrow(cmisage_history_byage3_agelast)))%>%
  melt(id='t')%>%
  mutate(chain=1,
         time=t)%>%
  group_by(time)%>%
  summarise(median=median(value),
            mean=mean(value),
            upper=quantile(value,probs=0.975),
            lower=quantile(value,probs=0.025))%>%
  mutate(epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates)
# cmisage_history_byage4_agelast <- data.frame(t(history_byage4[28, 51:1000, -1]))
# long_cmisage_history_byage4_agelast <- cmisage_history_byage4_agelast%>% mutate(t=c(1:nrow(cmisage_history_byage4_agelast)))%>%
#   melt(id='t')%>%
#   mutate(chain=1,
#          time=t)
# cmis_agelast_history_both <- rbind(long_cmisage_history_byage3_agelast,long_cmisage_history_byage4_agelast)%>%
#   group_by(time)%>%
#   summarise(median=median(value),
#             mean=mean(value),
#             upper=quantile(value,probs=0.975),
#             lower=quantile(value,probs=0.025))%>%
#   mutate(epoch=ifelse(time>61,'Without cMIS','With cMIS'),
#          month=dates)

cmis_data_cis_agelast <- addCIs(data_raw_cmis_byage[,c('t','tested_25','positive_25')],data_raw_cmis_byage$positive_25,data_raw_cmis_byage$tested_25)
dates <- seq(as.Date("2015-4-1"), as.Date("2022-2-1"), by = "months")
cmis_data_cis_agelast$month <- dates[1:61]

cmis_agelastprev_plot <- ggplot(long_cmisage_history_byage3_agelast)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_data_cis_agelast,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_data_cis_agelast,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='>44 years\nPrevalence')+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
cmis_agelastprev_plot
cmis_age1prev_plot/cmis_age2prev_plot/cmis_agecbaprev_plot/cmis_agelastprev_plot+ plot_layout(guides = "collect") 

##Replot but without LLIn dates and without ANC only
cmis_age1prev_plot_nollin <- ggplot(long_cmisage_history_byage3_age1)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_data_cis_age1,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_data_cis_age1,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='<5 years\nPrevalence')+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
cmis_age1prev_plot_nollin
cmis_age2prev_plot_nollin <- ggplot(cmis_age2_history)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_age2_observed,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_age2_observed,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='5-12 years\nPrevalence')+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
cmis_age2prev_plot_nollin
cmis_agecbaprev_plot_nollin <- ggplot(cmis_cba_history)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_cba_observed,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_cba_observed,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='13-44 years\nPrevalence')+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
cmis_agecbaprev_plot_nollin
cmis_agelastprev_plot_nollin <- ggplot(long_cmisage_history_byage3_agelast)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=cmis_data_cis_agelast,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=cmis_data_cis_agelast,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='>44 years\nPrevalence')+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.title.y = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
cmis_agelastprev_plot_nollin
nollin <- cmis_age1prev_plot_nollin/cmis_age2prev_plot_nollin/cmis_agecbaprev_plot_nollin/cmis_agelastprev_plot_nollin

##ANC prevalence by age group
change_prev_by_OR<-function(prev,OR){
  get_odds<-prev/(1-prev)*OR
  return(get_odds/(1+get_odds))
}
##Gravida 1
ancage_history_byage3_age22 <- data.frame(t(history_byage3[16, 51:1000, -1]))
history_byage3[,1,1]
long_ancg1_history_byage3_age22 <- ancage_history_byage3_age22%>% mutate(t=c(1:nrow(ancage_history_byage3_age22)))%>%
  melt(id='t')%>%
  mutate(chain=1,
         time=t,
         prev=change_prev_by_OR(value,1.1422))
ancage_history_byage4_age22 <- data.frame(t(history_byage4[16, 51:1000, -1]))
long_ancg1_history_byage4_age22 <- ancage_history_byage4_age22%>% mutate(t=c(1:nrow(ancage_history_byage4_age22)))%>%
  melt(id='t')%>%
  mutate(chain=1,
         time=t,
         prev=change_prev_by_OR(value,1.1422))
ancg1_age22_history_both <- rbind(long_ancg1_history_byage3_age22,long_ancg1_history_byage4_age22)%>%
  group_by(time)%>%
  summarise(median=median(value),
            mean=mean(value),
            upper=quantile(value,probs=0.975),
            lower=quantile(value,probs=0.025))%>%
  mutate(epoch=ifelse(time>61,'Without cMIS','With cMIS'),
         month=dates)

readRDS(ANC_cMIS_for_pmcmc,'ANC_cMIS_for_pmcmc (1).RDS')
ancg1_data_cis_age22 <- addCIs(ANC_cMIS_for_pmcmc[,c('t','tested_g1_09','positive_g1_09')],ANC_cMIS_for_pmcmc$positive_g1_09,ANC_cMIS_for_pmcmc$tested_g1_09)
dates <- seq(as.Date("2015-4-1"), as.Date("2022-2-1"), by = "months")
cmis_data_cis_agelast$month <- dates[1:61]



##Recategorize age groups and plot by gravidity
grav1_history <- re_cat_anc(history_byage3,grav_code = 1)
grav23_history <- re_cat_anc(history_byage3,grav_code = 2)
grav4p_history <- re_cat_anc(history_byage3,grav_code = 3)

all_history <- rbind(grav1_history[[1]],grav23_history[[1]],grav4p_history[[1]])%>%
  group_by(grav,time)%>%
  dplyr::summarise(median=median(value),
            mean=mean(value),
            upper=quantile(value,probs=0.975,na.rm=TRUE),
            lower=quantile(value,probs=0.025,na.rm=TRUE))%>%
  mutate(epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=dates)

all_observed <- rbind(grav1_history[[2]],grav23_history[[2]],grav4p_history[[2]])%>%
  mutate(time=rep(1:83,3),
         epoch=factor(ifelse(time>61,'ANC Only','cMIS + ANC'),levels=c('cMIS + ANC','ANC Only')),
         month=rep(dates,3))


prev_plot <- ggplot(all_history)+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=all_observed,aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=all_observed,aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  labs(x='Date',y='ANC Prevalence')+
  facet_grid(rows = vars(as.factor(grav)),labeller = labeller(`as.factor(grav)` = c(`1`='Gravida 1',`2`='Gravida 2-3', `3`='Gravida 4+')))+
  theme(legend.title = element_blank(),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
prev_plot

prev_grav1_plot <- ggplot(all_history[all_history$grav==1,])+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=all_observed[all_observed$grav==1,],aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=all_observed[all_observed$grav==1,],aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='Gravida 1\nPrevalence')+
  # facet_grid(rows = vars(as.factor(grav)),labeller = labeller(`as.factor(grav)` = c(`1`='Gravida 1',`2`='Gravida 2-3', `3`='Gravida 4+')))+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
prev_grav1_plot

prev_grav23_plot <- ggplot(all_history[all_history$grav==2,])+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=all_observed[all_observed$grav==2,],aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=all_observed[all_observed$grav==2,],aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='Gravida 2-3\nPrevalence')+
  # facet_grid(rows = vars(as.factor(grav)),labeller = labeller(`as.factor(grav)` = c(`1`='Gravida 1',`2`='Gravida 2-3', `3`='Gravida 4+')))+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
prev_grav23_plot

prev_grav4p_plot <- ggplot(all_history[all_history$grav==3,])+
  geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=epoch),alpha=0.2)+
  geom_line(aes(x=month,y=median,color=epoch),size=1)+
  # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # geom_vline(xintercept = as.Date('2020-7-1'),size=1,linetype='dashed',color='gray') +
  geom_point(data=all_observed[all_observed$grav==3,],aes(x=month,y=mean),pch = 19,
             col = "#666666")+
  geom_errorbar(data=all_observed[all_observed$grav==3,],aes(x=month,ymin=lower,ymax=upper),width = 0,
                col = "#666666")+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(breaks= datebreaks, date_labels = "%b %Y")+
  scale_y_continuous(limits=c(0,1))+
  labs(x='Date',y='Gravida 4+\nPrevalence')+
  # facet_grid(rows = vars(as.factor(grav)),labeller = labeller(`as.factor(grav)` = c(`1`='Gravida 1',`2`='Gravida 2-3', `3`='Gravida 4+')))+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"))
prev_grav4p_plot
windows(20,10)
cmis_plots <- cmis_age1prev_plot/cmis_age2prev_plot/cmis_agecbaprev_plot/cmis_agelastprev_plot
grav_plots <- prev_grav1_plot/prev_grav23_plot/prev_grav4p_plot
(cmis_age1prev_plot/cmis_age2prev_plot/cmis_agecbaprev_plot/cmis_agelastprev_plot) + 
  (prev_grav1_plot/prev_grav23_plot/prev_grav4p_plot)+plot_layout(guides = "collect") 
(cmis_plots | grav_plots) +plot_layout(guides = "collect")

(nollin | grav_plots) +plot_layout(guides = "collect")

cmis_age1prev_plot + plot_spacer() + 
  cmis_age2prev_plot + prev_grav1_plot +
  cmis_agecbaprev_plot + prev_grav23_plot +
  cmis_agelastprev_plot + prev_grav4p_plot +plot_layout(ncol = 2, guides = "collect")

cmis_age1prev_plot + plot_spacer() + 
  cmis_age2prev_plot + prev_grav1_plot +
  cmis_agecbaprev_plot + prev_grav23_plot +
  cmis_agelastprev_plot + prev_grav4p_plot +plot_layout(ncol = 2, guides = "collect")
