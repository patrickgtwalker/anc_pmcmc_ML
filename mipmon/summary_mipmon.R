library(tidyverse)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
library(bayesplot)
library(binom)
source('shared/addCIs.R')
library(lubridate)
library(matrixStats)
library(epitools)
library(coda)

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))
mipmon_data_pg
View(mipmon_data_mg)
##Arrange data for plotting
df_data_pg <- data.frame(t=numeric(),
                         month=Date(),
                         tested=integer(),
                         positive=numeric(),
                         mean=numeric(),
                         upper=numeric(),
                         lower=numeric(),
                         district=character(),
                         grav=character())

for(i in 1:length(mipmon_data_pg)){
  print(i)
  data_cis <- addCIs(mipmon_data_pg[[i]],mipmon_data_pg[[i]]$positive,mipmon_data_pg[[i]]$tested)%>%
    mutate(district=names(mipmon_data_pg[i]),
           grav = 'pg')
  df_data_pg <- rbind(df_data_pg,data_cis)
}
df_data_pg%>%group_by(district)%>%summarise(min=min(month),max=max(month))

df_data_mg <- data.frame(t=numeric(),
                         month=Date(),
                         tested=integer(),
                         positive=numeric(),
                         mean=numeric(),
                         upper=numeric(),
                         lower=numeric(),
                         district=character(),
                         grav=character())

for(i in 1:length(mipmon_data_mg)){
  print(i)
  data_cis <- addCIs(mipmon_data_mg[[i]],mipmon_data_mg[[i]]$positive,mipmon_data_mg[[i]]$tested)%>%
    mutate(district=names(mipmon_data_mg[i]),
           grav = 'mg')
  df_data_mg <- rbind(df_data_mg,data_cis)
}

df_data_all <- rbind(df_data_pg,df_data_mg)
colors <- c(`Ilha Josina` = "#1B9E77", `Magude Sede` = "#D95F02", `Manhica` = "#377EB8", `All Sites` = "#999999")

##Plot prevalence by gravidity and district
prev_plot <- ggplot(df_data_all)+
  geom_point(aes(x=as.Date(month),y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
  geom_errorbar(aes(x=as.Date(month),ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors)+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC Prevalence')+
  facet_grid(grav~district)+
  theme(legend.title = element_blank(),
        # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"),
  )
windows(10,10)
prev_plot

##Create diagnostic plots
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
names(mipmon_mg_bulk_seas_results) <- c('Ilha Josina', 'Magude Sede', 'Manhica', 'All Sites')

mipmon_diag_plots <- create_diag_figs_bulk(mipmon_mg_bulk_seas_results)

windows(10,10)
mipmon_diag_plots$ar.plot
mipmon_diag_plots$ess.plot
mipmon_diag_plots$posterior.trace
mipmon_diag_plots$EIR_SD.density
mipmon_diag_plots$init_EIR.density

create_prev_plots <- function(results,data_list=nnp_list){
  districts <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))

  colors <- c(`Ilha Josina` = "#1B9E77", `Magude Sede` = "#D95F02", `Manhica` = "#377EB8", `All Sites` = "#999999")
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in 1:4){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', 51:1000, -1]))
    long_prev_sum <- prev_history%>% mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df <- rbind(df,long_prev_sum)
    
  }
  df_data <- data.frame(t=numeric(),
                        tested=integer(),
                        positive=numeric(),
                        mean=numeric(),
                        upper=numeric(),
                        lower=numeric(),
                        district=character())
  
  for(i in start:(start+number-1)){
    data_cis <- addCIs(data_list[[i]],data_list[[i]]$positive,data_list[[i]]$tested)%>%
      mutate(district=names(data_list[i]),
             month = dates_list[[country]])
    df_data <- rbind(df_data,data_cis)
  }
  
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-1-1'), xmax = as.Date('2022-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2022-7-1'), xmax = as.Date('2022-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
  )
  prev_plot <- annotations[[country]]+
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    labs(x='Date',y='ANC Prevalence')+
    theme(legend.title = element_blank(),
          # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"),
    )
  prev_plot
}

create_dashboard_pgmg_plots <- function(results,prev_pg,prev_mg,prev_all=NULL,coefs_pg_df,coefs_mg_df,incidence){
  ## fn to return prevalence from log_odds
  get_prev_from_log_odds<-function(log_odds){
    return(exp(log_odds)/(1+exp(log_odds)))
  }
  
  ## fn to return odds from prevalence
  get_odds_from_prev<-function(prev){
    return(prev/(1-prev))
  }
  
  ##District names
  districts <- c('Ilha Josina', 'Magude Sede', 'Manhiça', 'All Sites')
  
  ##Color legend for districts
  colors <- c(`Ilha Josina` = "#1B9E77", `Magude Sede` = "#D95F02", `Manhiça` = "#377EB8", `All Sites` = "#999999")
  
  ##Cross-sectional data (Don't think we have for MipMon sites?)
  # cs_data <- cs_data_list%>%
  #   dplyr::rename(district=site)%>%
  #   mutate(month=as.Date(month))
  
  df_prev <- data.frame(time = integer(),
                        median = numeric(),
                        mean = numeric(),
                        upper = numeric(),
                        lower = numeric(),
                        district = character(),
                        month = character())
  df_prev_sample <- data.frame(time = integer(),
                               value = numeric(),
                               variable = character(),
                               district = character(),
                               month = character())
  for(i in 1:4){
    prev_history <- data.frame(t(results[[i]]$history['prev_05', 501:1000, -1]))

    long_prev_sum <- prev_history%>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = rowMedians(as.matrix(sapply(1:nrow(coefs_mg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      group_by(time)%>%
      dplyr::summarise(median_pg=median(prev_pg),
                       mean_pg=mean(prev_pg),
                       upper_pg=quantile(prev_pg,probs=0.975),
                       lower_pg=quantile(prev_pg,probs=0.025),
                       median_mg=median(prev_mg),
                       mean_mg=mean(prev_mg),
                       upper_mg=quantile(prev_mg,probs=0.975),
                       lower_mg=quantile(prev_mg,probs=0.025))%>%
      mutate(district = districts[i])
    long_prev_sum$month <- as.Date(prev_pg[[i]]$month)
    df_prev <- rbind(df_prev,long_prev_sum)
    
    prev_sample <- prev_history[, sample(ncol(prev_history), 100)] %>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      mutate(logodds_child = log(get_odds_from_prev(value)),
             prev_pg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_pg_df$gradient[x]*(logodds_child-coefs_pg_df$av_lo_child[x])+coefs_pg_df$intercept[x])
             }))),
             prev_mg = rowMedians(as.matrix(sapply(1:nrow(coefs_pg_df), function(x){
               prev_preg <- get_prev_from_log_odds(logodds_child+coefs_mg_df$gradient[x]*(logodds_child-coefs_mg_df$av_lo_child[x])+coefs_mg_df$intercept[x])
             }))))%>%
      mutate(district = districts[i])
    prev_sample$month <- as.Date(as.yearmon(rep(as.character(prev_pg[[i]]$month),100)))
    df_prev_sample <- rbind(df_prev_sample,prev_sample)
  }
  
  df_data_pg <- data.frame(t=numeric(),
                           month=Date(),
                           tested=integer(),
                           positive=numeric(),
                           mean=numeric(),
                           upper=numeric(),
                           lower=numeric(),
                           district=character())
  
  for(i in 1:4){
    data_cis <- addCIs(prev_pg[[i]],prev_pg[[i]]$positive,prev_pg[[i]]$tested)%>%
      mutate(district=districts[i])
    df_data_pg <- rbind(df_data_pg,data_cis)
  }
  df_data_pg$month <- as.Date(df_data_pg$month)
  
  df_data_mg <- data.frame(t=numeric(),
                           month=Date(),
                           tested=integer(),
                           positive=numeric(),
                           mean=numeric(),
                           upper=numeric(),
                           lower=numeric(),
                           district=character())
  
  for(i in 1:4){
    data_cis <- addCIs(prev_mg[[i]],prev_mg[[i]]$positive,prev_mg[[i]]$tested)%>%
      mutate(district=districts[i])
    df_data_mg <- rbind(df_data_mg,data_cis)
  }
  df_data_mg$month <- as.Date(df_data_mg$month)
  
  # df_data_all <- data.frame(t=numeric(),
  #                          tested=integer(),
  #                          positive=numeric(),
  #                          mean=numeric(),
  #                          upper=numeric(),
  #                          lower=numeric(),
  #                          district=character())
  
  # for(i in start:(start+number-1)){
  #   data_cis <- addCIs(prev_all[[i]],prev_all[[i]]$positive,prev_all[[i]]$tested)%>%
  #     mutate(district=names(prev_all[i]),
  #            month = dates_list[[country]])
  #   df_data_all <- rbind(df_data_all,data_cis)
  # }
  
  df_inc <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       district = character(),
                       month = character())
  df_inc_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              district = character(),
                              month = character())
  for(i in 1:4){
    inc_history <- data.frame(t(results[[i]]$history['clininc_05', 501:1000, -1]))
    
    long_inc_sum <- inc_history%>%
      dplyr::mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value)*10000*30,
                       mean=mean(value)*10000*30,
                       upper=quantile(value,probs=0.975)*10000*30,
                       lower=quantile(value,probs=0.025)*10000*30)%>%
      dplyr::mutate(district = districts[i])
    long_inc_sum$month <- as.Date(prev_pg[[i]]$month)
    df_inc <- rbind(df_inc,long_inc_sum)
    
    inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
      mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value*10000*30,
                    district = districts[i])
    inc_sample$month <-  as.Date(as.yearmon(rep(as.character(prev_pg[[i]]$month),100)))
    df_inc_sample <- rbind(df_inc_sample,inc_sample)
  }
  
  df_eir <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       district = character(),
                       month = character())
  df_eir_sample <- data.frame(time = integer(),
                              value = numeric(),
                              variable = character(),
                              district = character(),
                              month = character())
  for(i in 1:4){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 501:1000, -1]))
    
    long_eir_sum <- eir_history%>%
      dplyr::mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      group_by(time)%>%
      dplyr::summarise(median=median(value),
                       mean=mean(value),
                       upper=quantile(value,probs=0.975),
                       lower=quantile(value,probs=0.025))%>%
      dplyr::mutate(district = districts[i])
    long_eir_sum$month <- as.Date(prev_pg[[i]]$month)
    df_eir <- rbind(df_eir,long_eir_sum)
    
    eir_sample <- eir_history[, sample(ncol(eir_history), 100)] %>%
      mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      dplyr::rename(time=t)%>%
      dplyr::mutate(value = value,
                    district = districts[i])
    eir_sample$month <- as.Date(as.yearmon(rep(as.character(prev_pg[[i]]$month),100)))
    df_eir_sample <- rbind(df_eir_sample,eir_sample)
  }
  
  est_inc_plot <- ggplot()+
    geom_line(data=df_inc_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "'%y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(district~.)+
    labs(title = 'Estimated Incidence\nper 10,000 person-months')+
    coord_cartesian(ylim=c(0,1200))+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  blank_inc_plot <- ggplot()+
    geom_line(data=df_inc_sample,aes(x=month,y=value,color=district,group=variable),alpha=0)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_inc,aes(x=month,y=median,color=district,group=district),size=1,alpha=0)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "'%y")+
    # coord_cartesian(ylim=c(0, 1000))+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    facet_grid(district~.)+
    labs(title = 'Estimated Incidence\nper 10,000 person-months')+
    coord_cartesian(ylim=c(0,1200))+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt")
    )
  
  # print('est_inc_plot')
  obs_prev_plot_mg <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    # geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    # geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    geom_point(data=df_data_mg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_mg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "'%y")+
    facet_grid(district~.)+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(df_data_mg$month))+
    labs(title = 'ANC Prevalence\nMultigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  est_prev_plot_mg <- ggplot()+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    # geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),pch= 2, size=1.5,color="black")+
    # geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color='black')+
    geom_line(data=df_prev_sample,aes(x=month,y=prev_mg,color=district,group=variable),alpha=0.2)+
    geom_line(data=df_prev,aes(x=month,y=median_mg,color=district,group=district),linewidth=1)+
    geom_point(data=df_data_mg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_mg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "'%y")+
    facet_grid(district~.)+
    scale_y_continuous(limits = c(0,1))+
    coord_cartesian(ylim=c(0, 1),xlim = range(df_data_mg$month))+
    labs(title = 'ANC Prevalence\nMultigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # print('obs_prev_plot_mg')
  # # obs_prev_plot_all <- annotations[[country]]+
  # #   # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
  # #   # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
  # #   # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
  # #   geom_point(data=df_data_all,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
  # #   geom_errorbar(data=df_data_all,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
  # #   scale_color_manual(values=colors)+
  # #   scale_fill_manual(values=colors)+
  # #   scale_x_date(date_labels = "%b %Y")+
  # #   facet_grid(district~.)+
  # #   scale_y_continuous(limits = c(0,1))+
  # #   coord_cartesian(ylim=c(0, 1))+
  # #   labs(title = 'ANC Prevalence\nAll gravidities')+
  # #   theme(legend.title = element_blank(),
  # #         axis.title.x = element_blank(),
  # #         axis.title.y = element_blank(),
  # #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  # #         axis.ticks.x = element_line(linewidth = 0.5),
  # #         axis.ticks.length = unit(3, "pt"),
  # #         legend.position = 'none'
  # #   )
  est_prev_plot_pg <- ggplot()+
    # geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    # geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
    geom_line(data=df_prev_sample,aes(x=month,y=prev_pg,color=district,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(data=df_prev,aes(x=month,y=median_pg,color=district,group=district),linewidth=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data_pg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_pg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "'%y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(df_data_pg$month))+
    facet_grid(district~.)+
    labs(title = 'ANC Prevalence\nPrimigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  obs_prev_plot_pg <- ggplot()+
    # geom_point(data=cs_data,aes(x=month,y=mean,color=district,group=district),size=1.5,pch=2,color="black")+
    # geom_errorbar(data=cs_data,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,color=district,group=district),color="black")+
    geom_point(data=df_data_pg,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data_pg,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "'%y")+
    scale_y_continuous(limits=c(0,1))+
    coord_cartesian(xlim = range(df_data_pg$month))+
    facet_grid(district~.)+
    labs(title = 'ANC Prevalence\nPrimigrav')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  # print('est_prev_plot')
  incidence <- plyr::rbind.fill(incidence,data.frame(district='All Sites'))
  obs_inc_plot <- ggplot()+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=inc,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "'%y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    facet_grid(district~.)+
    labs(title='Observed Incidence\nper 10,000 person-months')+
    # coord_cartesian(ylim=c(0, ifelse(country=='NG',1000,max(incidence$inc)*2)))+
    # coord_cartesian(ylim=c(0, max(incidence$inc)*2))+
    coord_cartesian(ylim=c(0,1200))+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt"))
  # print('obs_inc_plot')
  obs_pos_plot <- ggplot()+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(data=incidence,aes(x=date_ex,y=prop_pos,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "'%y")+
    scale_y_continuous(limits=c(0,1))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    # labs(x='Date',y='EIR')+
    labs(title = 'Test positivity\nproportion')+
    facet_grid(district~.)+
    theme(legend.title = element_blank(),
          legend.position = 'none',
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5),
          axis.ticks.length = unit(3, "pt"))
  # print('obs_pos_plot')
  # # obs_test_plot <- annotations[[country]]+
  # #   # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
  # #   geom_line(aes(x=date_ex,y=inc,color=district,group=district),size=1)+
  # #   # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
  # #   scale_color_manual(values=colors)+
  # #   scale_fill_manual(values=colors)+
  # #   scale_x_date(date_labels = "%b %Y")+
  # #   # scale_y_continuous(limits=c(0,500))+
  # #   # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
  # #   # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
  # #   labs(x='Date',y='Clinical Incidence\nper 10,000 person-months')+
  # #   # labs(x='Date',y='EIR')+
  # #   facet_grid(~district)+
  # #   theme(legend.title = element_blank(),
  # #         axis.title.x = element_blank(),
  # #         axis.text.x=element_text(angle=45, hjust=1, vjust=1),
  # #         axis.ticks.x = element_line(size = 0.5), 
  # #         axis.ticks.length = unit(3, "pt"))
  full_dash <- est_prev_plot_pg+est_prev_plot_mg+est_inc_plot+obs_inc_plot+obs_pos_plot+ plot_layout(guides = "collect",ncol=5)
  obs_data_dash <- obs_prev_plot_pg+obs_prev_plot_mg+est_inc_plot+obs_inc_plot+ plot_layout(guides = "collect",ncol=4)
  pres_dash_1 <- obs_prev_plot_pg+obs_prev_plot_mg+blank_inc_plot+obs_inc_plot+ plot_layout(guides = "collect",ncol=4)
  pres_dash_2 <- est_prev_plot_pg+est_prev_plot_mg+blank_inc_plot+obs_inc_plot+ plot_layout(guides = "collect",ncol=4)
  pres_dash_fin <- est_prev_plot_pg+est_prev_plot_mg+est_inc_plot+obs_inc_plot+ plot_layout(guides = "collect",ncol=4)
  
  return(list(full_dash = full_dash,
              obs_data_dash = obs_data_dash,
              pres_dash_1 = pres_dash_1,
              pres_dash_2 = pres_dash_2,
              pres_dash_fin = pres_dash_fin))
  # obs_prev_plot_mg
}

##Read in correlation coefficients
coefs_pg_df <- as.data.frame(readRDS('./nnp/Corr/pg_corr_sample.RDS'))
coefs_mg_df <- as.data.frame(readRDS('./nnp/Corr/mg_corr_sample.RDS'))

##read in and process clinical incidence data
manhica_inc <- read_csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/isglobal_cism_data/weekly_OPD_cases_2014_2019_6_posts_age_5.csv')
magude_inc <- read_csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/isglobal_cism_data/RRS_data_age.csv')


manhica_inc_filtered <- manhica_inc[manhica_inc$place %in% c('Manhiça','Ilha Josina')&!is.na(manhica_inc$malaria),]%>%
  group_by(place, week, year)%>%
  dplyr::summarise(Tested=sum(testdone),
                   Positive=sum(malaria),
                   Total=sum(visits),
                   Population=ifelse(place == 'Manhiça',round(78479*0.07),round(8288*0.07)))#population from Nhacolo, et al. 2021. Int J Epidemiology. Took population of sub-district in 2019 times the percentage of the population <5yo
manhica_inc_forplot <- addCIs_inc(manhica_inc_filtered,manhica_inc_filtered$Positive,manhica_inc_filtered$Population)%>%
  mutate(date_ex=as.Date(paste(year, week, '1'), "%Y %U %u"),
         prop_pos = Positive/Tested,
         week = formatC(week, width = 2, format = "d", flag = "0"))%>%
  filter(date_ex>=as.Date('2016-10-31')&date_ex<=as.Date('2019-10-31'))%>%
  dplyr::rename(district = place,
                inc = mean)

library(ISOweek)
source('shared/addCIs_inc.R')
magude_inc_filtered <- magude_inc[!is.na(magude_inc$mal_a),]%>%
  mutate(week = formatC(week, width = 2, format = "d", flag = "0"))%>%
  group_by(week, yr)%>%
  dplyr::summarise(Tested=sum(tot_test_a),
                   Positive=sum(mal_a),
                   Total=sum(tot_a),
                   Population=round(44203*0.08))#population from Galatas, et al. 2020. BMJ. Took population of Magude Sede in 2016 times the percentage of the population <5yo
magude_inc_forplot <- addCIs_inc(magude_inc_filtered,magude_inc_filtered$Positive,magude_inc_filtered$Population)%>%
  mutate(date_ex=ISOweek2date(paste0(yr,'-W', week,'-1')),
         prop_pos = Positive/Tested)%>%
  filter(date_ex>=as.Date('2016-10-31')&date_ex<=as.Date('2019-10-31'))%>%
  mutate(district = 'Magude Sede')%>%
  dplyr::rename(inc = mean)
all_inc_forplot <- rbind(manhica_inc_forplot,magude_inc_forplot)
total_inc <- all_inc_forplot %>%
  mutate(year = ifelse(district=='Magude Sede',yr,year))%>%
  group_by(week,year)%>%
  summarise(Tested = sum(Tested),
            Positive = sum(Positive),
            Population = sum(Population),
            Total = sum(Total))%>%
  mutate(date_ex = as.Date(paste(year, week, '1'), "%Y %U %u"))
total_inc <- addCIs_inc(total_inc,total_inc$Positive,total_inc$Population)%>%
  rename(inc = mean)%>%
  mutate(district = 'All Sites',
         prop_pos = Positive/Tested)
all_inc_forplot <- rbind(all_inc_forplot,total_inc)
mipmon_dash <- create_dashboard_pgmg_plots(results=mipmon_mg_bulk_seas_results,
                            prev_pg=mipmon_data_pg,
                            prev_mg=mipmon_data_mg,
                            coefs_pg_df=coefs_pg_df,
                            coefs_mg_df=coefs_mg_df,
                            incidence=all_inc_forplot)
rep(mipmon_data_pg[[1]]$month,2)
windows(30,15)
mipmon_dash$full_dash
mipmon_dash$obs_data_dash
mipmon_dash$pres_dash
mipmon_dash$pres_dash_1
mipmon_dash$pres_dash_2
mipmon_dash$pres_dash_fin
ggsave(paste0('Q:/anc_pmcmc/mipmon/dash_mipmon_pgmg_seas_300523.pdf'), plot = mipmon_dash$full_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/mipmon/dash_mipmon_pgmg_seas_300523_obs.pdf'), plot = mipmon_dash$obs_data_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/mipmon/dash_mipmon_pgmg_seas_300523_pres.pdf'), plot = mipmon_dash$pres_dash, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/mipmon/dash_mipmon_pgmg_seas_300523_pres_1.pdf'), plot = mipmon_dash$pres_dash_1, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/mipmon/dash_mipmon_pgmg_seas_300523_pres_2.pdf'), plot = mipmon_dash$pres_dash_2, width = 12, height = 6)
ggsave(paste0('Q:/anc_pmcmc/mipmon/dash_mipmon_pgmg_seas_300523_pres_fin.pdf'), plot = mipmon_dash$pres_dash_fin, width = 12, height = 6)
