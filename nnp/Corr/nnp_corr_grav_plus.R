library(R2OpenBUGS)
library(coda)
library(car)
library(MASS)
library(plotrix)
library(binom)
library(scales)
library(zoo)
library(plyr)
library(dplyr)
library(ggplot2)
library(viridis)
library(ggpubr)
library(readxl)

## fn to return prevalence from log_odds
get_prev_from_log_odds<-function(log_odds){
  return(exp(log_odds)/(1+exp(log_odds)))
}

## fn to return odds from prevalence
get_odds_from_prev<-function(prev){
  return(prev/(1-prev))
}

addCIs<-function(df,Ys.cs,Ns.cs,Ys.anc,Ns.anc){
  df$mean.cs<-NA
  df$upper.cs<-NA
  df$lower.cs<-NA
  CIs.cs<-binom.confint(Ys.cs,Ns.cs,method="exact")
  df$mean.cs[Ns.cs>0]<-CIs.cs$mean[Ns.cs>0]
  df$upper.cs[Ns.cs>0]<-CIs.cs$upper[Ns.cs>0]
  df$lower.cs[Ns.cs>0]<-CIs.cs$lower[Ns.cs>0]
  df$mean.anc<-NA
  df$upper.anc<-NA
  df$lower.anc<-NA
  CIs.anc<-binom.confint(Ys.anc,Ns.anc,method="exact")
  df$mean.anc[Ns.anc>0]<-CIs.anc$mean[Ns.anc>0]
  df$upper.anc[Ns.anc>0]<-CIs.anc$upper[Ns.anc>0]
  df$lower.anc[Ns.anc>0]<-CIs.anc$lower[Ns.anc>0]
  return(df)
}
##Load in data data grouped by site and month
NG_cs <- readRDS('nnp/data/NG_CS_all_grouped_site.rds')
NG_anc <- readRDS('nnp/data/NG_ANC_mother_grouped_sitegrav.rds')
MZ_cs <- readRDS('nnp/data/MZ_CS_all_grouped_site_180123.rds')
MZ_anc <- readRDS('nnp/data/MZ_ANC_mother_grouped_sitegrav_0822.rds')
BF_cs <- readRDS('nnp/data/BF_CS_all_grouped_site_0822.rds')
BF_anc <- readRDS('nnp/data/BF_ANC_mother_grouped_sitegrav.rds')

NG_cs$country <- 'Nigeria'
NG_anc$country <- 'Nigeria'
MZ_cs$country <- 'Mozambique'
MZ_anc$country <- 'Mozambique'
BF_cs$country <- 'Burkina Faso'
BF_anc$country <- 'Burkina Faso'

all_cs <- rbind(NG_cs,MZ_cs,BF_cs)
all_anc <- rbind(NG_anc,MZ_anc,BF_anc)

all_both <- merge(all_cs,all_anc, by = c('site','month','country'), suffixes = c('.cs','.anc'))
all_both_total <- all_both %>%
  group_by(country,site,month,.drop=FALSE)%>%
  dplyr::summarise(positive.anc=sum(positive.anc),total.anc=sum(total.anc),positive.cs=mean(positive.cs),total.cs=mean(total.cs),
                   mean.cs=mean(mean.cs),upper.cs=mean(upper.cs),lower.cs=mean(lower.cs),
                   mean.anc=mean(mean.anc),upper.anc=mean(upper.anc),lower.anc=mean(lower.anc))%>%
  mutate(grav_cat = 'All pregnancies')
all_both <- rbind(all_both,all_both_total)


vaneijk <- read_excel('nnp/data/paper_data_strat_G.xlsx')
primi_ve <- vaneijk[vaneijk$subG == 'P',3:6]%>%
  dplyr::rename(positive.cs=child_Y,
         total.cs=child_N,
         positive.anc=preg_Y,
         total.anc=preg_N)
multi_ve <- vaneijk[vaneijk$subG == 'M',3:6]%>%
  dplyr::rename(positive.cs=child_Y,
         total.cs=child_N,
         positive.anc=preg_Y,
         total.anc=preg_N)
all_ve <- vaneijk %>%
  subset(subG %in% c('P','M')) %>%
  group_by(Site)%>%
  dplyr::summarise(inf_prev_n=mean(child_Y),
            inf_N=mean(child_N),
            ANC_prev_n=sum(preg_Y),
            ANC_N=sum(preg_N))%>%
  dplyr::select(-Site) %>%
  dplyr::rename(positive.cs=inf_prev_n,
         total.cs=inf_N,
         positive.anc=ANC_prev_n,
         total.anc=ANC_N)
primi_ve <- addCIs(primi_ve,primi_ve$positive.cs,primi_ve$total.cs,primi_ve$positive.anc,primi_ve$total.anc)%>%
  mutate(country='Van Eijk')
multi_ve <- addCIs(multi_ve,multi_ve$positive.cs,multi_ve$total.cs,multi_ve$positive.anc,multi_ve$total.anc)%>%
  mutate(country='Van Eijk')
all_ve <- addCIs(all_ve,all_ve$positive.cs,all_ve$total.cs,all_ve$positive.anc,all_ve$total.anc)%>%
  mutate(country='Van Eijk')

all_nnp_total <- all_both_total %>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
         `Burkina Faso`='NNP - Burkina Faso',
         `Mozambique`='NNP - Mozambique'))
all_nnp_pg <- all_both[all_both$grav_cat=='Gravidities 1',]%>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                        `Burkina Faso`='NNP - Burkina Faso',
                        `Mozambique`='NNP - Mozambique'))
all_nnp_mg <- all_both[all_both$grav_cat %in% c('Gravidities 2-3','Gravidities 4+'),]%>%
  group_by(country,site,month)%>%
  dplyr::summarise(positive.cs=mean(positive.cs),
            total.cs=mean(total.cs),
            positive.anc=sum(positive.anc),
            total.anc=sum(total.anc))%>%
  dplyr::mutate(country=dplyr::recode(country,`Nigeria`='NNP - Nigeria',
                        `Burkina Faso`='NNP - Burkina Faso',
                        `Mozambique`='NNP - Mozambique'))
all_nnp_mg <- addCIs(all_nnp_mg,all_nnp_mg$positive.cs,all_nnp_mg$total.cs,all_nnp_mg$positive.anc,all_nnp_mg$total.anc)

names(all_nnp)
names(all_mipmon)
mipmon <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/Analysis/mipmon/mipmon4corr.RDS')
all_mipmon_total <- mipmon[mipmon$grav=='All',]%>%
  group_by(site,month)%>%
  dplyr::summarise(positive.cs=mean(positive.cs),
            total.cs=mean(total.cs),
            positive.anc=sum(positive.anc),
            total.anc=sum(total.anc))
all_mipmon_pg <- mipmon[mipmon$grav=='primi',]
all_mipmon_mg <- mipmon[mipmon$grav=='multi',]
all_mipmon_total <- addCIs(all_mipmon_total,all_mipmon_total$positive.cs,all_mipmon_total$total.cs,all_mipmon_total$positive.anc,all_mipmon_total$total.anc)%>%
  mutate(country='MiPMon')
all_mipmon_pg <- addCIs(all_mipmon_pg,all_mipmon_pg$positive.cs,all_mipmon_pg$total.cs,all_mipmon_pg$positive.anc,all_mipmon_pg$total.anc)%>%
  mutate(country='MiPMon')
all_mipmon_mg <- addCIs(all_mipmon_mg,all_mipmon_mg$positive.cs,all_mipmon_mg$total.cs,all_mipmon_mg$positive.anc,all_mipmon_mg$total.anc)%>%
  mutate(country='MiPMon')

all_data_total <- plyr::rbind.fill(all_nnp_total,all_ve,all_mipmon_total)
all_data_pg <- plyr::rbind.fill(all_nnp_pg,primi_ve,all_mipmon_pg)
all_data_mg <- plyr::rbind.fill(all_nnp_mg,multi_ve,all_mipmon_mg)

all_data_total4model <- all_data_total %>%
  filter(country != 'NNP - Burkina Faso') %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')))
all_data_pg4model <- all_data_pg %>%
  filter(country != 'NNP - Burkina Faso') %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')))
all_data_mg4model <- all_data_mg %>%
  filter(country != 'NNP - Burkina Faso') %>%
  filter(!(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')))

country_levels <- c("NNP - Mozambique","NNP - Nigeria","MiPMon","Van Eijk",'NNP - SMC implemented')
country_levels_nosmc <- c("NNP - Mozambique","NNP - Nigeria","MiPMon","Van Eijk")
all_data_total4fig <- all_data_total %>%
  mutate(country = factor(ifelse((country=='NNP - Burkina Faso')|(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')),'NNP - SMC implemented',country),
         levels = country_levels))

all_data_pg4fig <- all_data_pg %>%
  mutate(country = factor(ifelse((country=='NNP - Burkina Faso')|(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')),'NNP - SMC implemented',country),
                          levels = country_levels))
all_data_mg4fig <- all_data_mg %>%
  mutate(country = factor(ifelse((country=='NNP - Burkina Faso')|(country == 'NNP - Nigeria' & (site == 'Asa' | site == 'Moro') & month == as.yearmon('Nov 2021')),'NNP - SMC implemented',country),
                          levels = country_levels))

N_sites=nrow(all_data_total4model)

####All pregnancies####
windows(height=5,width=10)
windows(height=5,width=3.3)
par(mar=c(5,5,1,1))

### now define the model in bugs terms to fit to the data
preg_child_model<-function(){
  for (i in 1:N) {
    pos_child[i] ~ dbin(p_child[i],total_child[i])
    pos_preg[i] ~ dbin(p_preg[i],total_preg[i])
    logit(p_child[i]) <- log_odds_child[i]
    logit(p_preg[i]) <- log_odds_child[i]+log_OR_p_v_c[i]
    log_odds_child[i]~dnorm(av_lo_child,sigma_c_inv)
    log_OR_p_v_c[i]<-RE_intercept[i]+gradient*(log_odds_child[i]-av_lo_child)
    RE_intercept[i]~dnorm(intercept,sigma_int_inv)
  }
  intercept~dnorm(0,0.001)
  gradient ~ dnorm(0,0.001)
  av_lo_child~dnorm(0,0.001)
  sigma_c_inv~dgamma(0.001,0.001)
  sigma_int_inv~dgamma(0.001,0.001)
  sigma_child<-1/sqrt(sigma_c_inv)
  sigma_intercept<-1/sqrt(sigma_int_inv)
}

## write to directory
model.file <- file.path(tempdir(),"model.txt")
write.model(preg_child_model, model.file)

## put in bugs data format
data_list<-list(pos_child=all_data_total4model$positive.cs,
                total_child=all_data_total4model$total.cs,
                pos_preg=all_data_total4model$positive.anc,
                total_preg=all_data_total4model$total.anc,
                N=N_sites)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept=rnorm(1),gradient=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept","gradient","sigma_child","sigma_intercept")

## run the model
run_model <- bugs(data_list, inits, params,
                model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient,c(0.025,0.5,0.975))
quantile(intercept,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient)*(logodds_child-median(av_lo_child))+median(intercept))
pred <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

##Try to make similar graph in ggplot
colors <- c(viridis(4,begin=0.2,end=0.9),"#D95F02")
colors_nosmc <- c(viridis(4,begin=0.2,end=0.9))

grav_all <- ggplot(all_data_total4fig)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='All pregnancies',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
grav_all
##plot without smc
grav_all_nosmc <- ggplot(all_data_total4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='All pregnancies',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
grav_all_nosmc

###########  
##Primigrav
###########
##Model has been defined above
## put in bugs data format
data_list_pg<-list(pos_child=all_data_pg4model$positive.cs,
                total_child=all_data_pg4model$total.cs,
                pos_preg=all_data_pg4model$positive.anc,
                total_preg=all_data_pg4model$total.anc,
                N=N_sites)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept=rnorm(1),gradient=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept","gradient","sigma_child","sigma_intercept")

## run the model
run_model_pg <- bugs(data_list_pg, inits, params,
                  model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model_pg)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient,c(0.025,0.5,0.975))
quantile(intercept,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient)*(logodds_child-median(av_lo_child))+median(intercept))
pred_pg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

##save a sample of 1000 simulations for use in pmcmc model fitting:
pg_corr_sample <- sample_n(as.data.frame(run_model_pg$sims.matrix),1000)
saveRDS(pg_corr_sample,'pg_corr_sample.RDS')
##Try to make similar graph in ggplot
grav_pg <- ggplot(all_data_pg4fig)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_pg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_pg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Primigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
grav_pg

##Plot points without smc
grav_pg <- ggplot(all_data_pg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_pg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_pg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Primigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
grav_pg
###########  
##Multi Grav
###########
##Model has been defined above
## put in bugs data format
data_list_mg<-list(pos_child=all_data_mg4model$positive.cs,
                   total_child=all_data_mg4model$total.cs,
                   pos_preg=all_data_mg4model$positive.anc,
                   total_preg=all_data_mg4model$total.anc,
                   N=N_sites-1)

## fn for random initial conditions
inits<-function(){
  list(av_lo_child=rnorm(1),intercept=rnorm(1),gradient=rnorm(1),sigma_c_inv=runif(1),sigma_int_inv=runif(1))
}
## params of interest
params<-c("av_lo_child","intercept","gradient","sigma_child","sigma_intercept")

## run the model
run_model_mg <- bugs(data_list_mg, inits, params,
                     model.file, n.iter=50000,n.burnin=10000)

##attach the output
attach.bugs(run_model_mg)
### should recapture params - something slightly funny with gradient- need to check SD formats and priors
quantile(gradient,c(0.025,0.5,0.975))
quantile(intercept,c(0.025,0.5,0.975))
quantile(av_lo_child,c(0.025,0.5,0.975))
get_prev_from_log_odds(quantile(av_lo_child,c(0.025,0.5,0.975)))
## lets plot the relationship (apols for horrible base R)
prev_child=seq(0.01,0.99,by=0.01)
logodds_child=log(get_odds_from_prev(prev_child))
prev_preg_lower=array(dim=length(logodds_child))
prev_preg_upper=array(dim=length(logodds_child))
for(i in 1:length(logodds_child)){
  prev_preg_lower[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.025))
  prev_preg_upper[i]=get_prev_from_log_odds(quantile(logodds_child[i]+intercept+gradient*(logodds_child[i]-av_lo_child),0.975))
}
prev_preg_median=get_prev_from_log_odds(logodds_child+median(gradient)*(logodds_child-median(av_lo_child))+median(intercept))
pred_mg <- data.frame(prev_child,prev_preg_lower,prev_preg_upper,prev_preg_median)

##Try to make similar graph in ggplot
grav_mg <- ggplot(all_data_mg4fig)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_mg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_mg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Multigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
grav_mg

##No smc
grav_mg <- ggplot(all_data_mg4model)+
  geom_point(aes(x=mean.cs*100,y=mean.anc*100,col=country),size=3)+
  geom_errorbar(aes(x=mean.cs*100,ymin=lower.anc*100,ymax=upper.anc*100,col=country),width=0)+
  geom_errorbarh(aes(y=mean.anc*100,xmin=lower.cs*100,xmax=upper.cs*100,col=country),height=0)+
  scale_color_manual(name="Study",values=colors_nosmc)+
  scale_y_continuous(limits = c(0,100))+
  scale_x_continuous(limits = c(0,100))+
  geom_abline(size=0.8,linetype='dashed')+
  geom_ribbon(data=pred_mg,aes(x=prev_child*100,ymin=prev_preg_lower*100,ymax=prev_preg_upper*100),alpha=0.2)+
  geom_line(data=pred_mg,aes(x=prev_child*100,y=prev_preg_median*100),size=1)+
  theme_bw()+
  theme(legend.position = 'bottom')+
  labs(title='Multigravida',x='Cross-section Prevalence (<5 yo)',y='ANC Prevalence')
grav_mg

windows(height=5,width=13)

ggarrange(grav_all,grav_pg,grav_mg,
          ncol = 3, nrow = 1,
          common.legend = TRUE, legend = "bottom")
