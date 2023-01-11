##summarize results##
library(ggplot2)
library(reshape2)
library(ggpubr)
library(zoo)
library(patchwork)
library(RColorBrewer)
source('addCIs.R')

theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

result <- mz_chemba
result <- ng_asa
result <- bf_banfora
create_diag_figs <- function(result){
  print('acceptance rate')
  print(1 - coda::rejectionRate(as.mcmc(result$mcmc)))
  print('effective size')
  print(coda::effectiveSize(as.mcmc(result$mcmc)))

  diag <- ((bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'log_prior')+mcmc_dens(result$mcmc[51:1000,],pars = 'log_prior'))/
      (bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'log_likelihood')+mcmc_dens(result$mcmc[51:1000,],pars = 'log_likelihood'))/
      (bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'log_posterior')+mcmc_dens(result$mcmc[51:1000,],pars = 'log_posterior'))/
      (bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'EIR_SD')+mcmc_dens(result$mcmc[51:1000,],pars = 'EIR_SD'))/
      (bayesplot::mcmc_trace(result$mcmc[51:1000,],pars = 'log_init_EIR')+mcmc_dens(result$mcmc[51:1000,],pars = 'log_init_EIR'))) + plot_layout(guides = "collect") #& theme(legend.position = "bottom")
  
  return(diag)
}

diag_plots <- lapply(1:10,function(i) create_diag_figs(nnp_result_list[[i]]))
diag_pg_plots <- lapply(1:10,function(i) create_diag_figs(nnp_pg_result_list[[i]]))
diag_pgcorr_plots <- lapply(1:10,function(i) create_diag_figs(nnp_pgcorr_result_list[[i]]))
diag_pgvol_plots <- lapply(1:10,function(i) create_diag_figs(nnp_pgvol_result_list[[i]]))
diag_pgseas_plots <- lapply(1:10,function(i) create_diag_figs(nnp_pgseas_result_list[[i]]))
diag_pgorig_plots <- lapply(1:10,function(i) create_diag_figs(nnp_pgorig_result_list[[i]]))

windows(10,7)
diag_pg_plots[10]
diag_pgcorr_plots[9]
diag_pgvol_plots[9]
diag_pgseas_plots[10]
diag_pgorig_plots[10]

country <- 'BF'
country <- 'MZ'
country <- 'NG'

results <- nnp_result_list
nnp_list <- list(Banfora = data_raw_bf_banfora,Gaoua = data_raw_bf_gaoua,Orodara = data_raw_bf_orodara,
                 Changara = data_raw_mz_changara, Chemba = data_raw_mz_chemba,Guro = data_raw_mz_guro,
                 Asa = data_raw_ng_asa,Ejigbo = data_raw_ng_ejigbo,`Ife North` = data_raw_ng_ifenorth, Moro = data_raw_ng_moro)
nnp_pg_list <- list(Banfora = data_raw_bf_pg_banfora,Gaoua = data_raw_bf_pg_gaoua,Orodara = data_raw_bf_pg_orodara,
                 Changara = data_raw_mz_pg_changara, Chemba = data_raw_mz_pg_chemba,Guro = data_raw_mz_pg_guro,
                 Asa = data_raw_ng_pg_asa,Ejigbo = data_raw_ng_pg_ejigbo,`Ife North` = data_raw_ng_pg_ifenorth, Moro = data_raw_ng_pg_moro)

create_prev_plots <- function(results,data_list=nnp_list,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                 MZ = c('Changara','Chemba','Guro'),
                 NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -1]))
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
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
      
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

##Prev plots
BF_prev <- create_prev_plots(results = nnp_result_list,country='BF')
MZ_prev <- create_prev_plots(results = nnp_result_list,country='MZ')
NG_prev <- create_prev_plots(results = nnp_result_list,country='NG')
BF_prev_aa <- create_prev_plots(results = nnp_result_list_aa,country='BF')
MZ_prev_aa <- create_prev_plots(results = nnp_result_list_aa,country='MZ')
NG_prev_aa <- create_prev_plots(results = nnp_result_list_aa,country='NG')

BF_prev_pg <- create_prev_plots(results = nnp_pg_result_list,data_list=nnp_pg_list,country='BF')
MZ_prev_pg <- create_prev_plots(results = nnp_pg_result_list,data_list=nnp_pg_list,country='MZ')
NG_prev_pg <- create_prev_plots(results = nnp_pg_result_list,data_list=nnp_pg_list,country='NG')

BF_prev_pgcorr <- create_prev_plots(results = nnp_pgcorr_result_list,data_list=nnp_pg_list,country='BF')
MZ_prev_pgcorr <- create_prev_plots(results = nnp_pgcorr_result_list,data_list=nnp_pg_list,country='MZ')
NG_prev_pgcorr <- create_prev_plots(results = nnp_pgcorr_result_list,data_list=nnp_pg_list,country='NG')

BF_prev_pgvol <- create_prev_plots(results = nnp_pgvol_result_list,data_list=nnp_pg_list,country='BF')
MZ_prev_pgvol <- create_prev_plots(results = nnp_pgvol_result_list,data_list=nnp_pg_list,country='MZ')
NG_prev_pgvol <- create_prev_plots(results = nnp_pgvol_result_list,data_list=nnp_pg_list,country='NG')

BF_prev_pgvol
NG_prev_pgvol
BF_prev_pgseas <- create_prev_plots(results = nnp_pgseas_result_list,data_list=nnp_pg_list,country='BF')
MZ_prev_pgseas <- create_prev_plots(results = nnp_pgseas_result_list,data_list=nnp_pg_list,country='MZ')
NG_prev_pgseas <- create_prev_plots(results = nnp_pgseas_result_list,data_list=nnp_pg_list,country='NG')

windows(7,5)
windows(10,7)
BF_prev
MZ_prev
NG_prev
BF_prev_pg
MZ_prev_pg
NG_prev_pg
BF_prev_pgcorr
MZ_prev_pgcorr
NG_prev_pgcorr

windows(10,5)
BF_prev_aa + BF_prev
MZ_prev_aa + MZ_prev
NG_prev_aa + NG_prev
BF_prev_pg+BF_prev_pgcorr
MZ_prev_pg+MZ_prev_pgcorr
NG_prev_pg+NG_prev_pgcorr

BF_prev+BF_prev_pgvol
MZ_prev+MZ_prev_pgvol
NG_prev+NG_prev_pgvol

##Inc and EIR plots
create_inc_plots <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    inc_history <- data.frame(t(results[[i]]$history['inc', 51:1000, -1]))
    long_inc_sum <- inc_history%>% mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value)*10000*30,
                mean=mean(value)*10000*30,
                upper=quantile(value,probs=0.975)*10000*30,
                lower=quantile(value,probs=0.025)*10000*30)%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df <- rbind(df,long_inc_sum)
    
  }

  df_eir <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))
    long_eir_sum <- eir_history%>% mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df_eir <- rbind(df_eir,long_eir_sum)
  }
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits = c(0,3750)),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-1'), ymin = 0, ymax = 400,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,400)),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,3750))

  )
  inc_plot <- annotations[[country]]+
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='Estimated Incidence (<5 year old)\nper 10,000 person-months')+
    # labs(x='Date',y='EIR')+
    theme(legend.title = element_blank(),
          # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  # +
  #   coord_cartesian(ylim=c(0, 40))
  inc_plot
}

BF_inc_aa <- create_inc_plots(results = nnp_result_list_aa,country='BF')
MZ_inc_aa <- create_inc_plots(results = nnp_result_list_aa,country='MZ')
NG_inc_aa <- create_inc_plots(results = nnp_result_list_aa,country='NG')
BF_inc <- create_inc_plots(results = nnp_result_list,country='BF')
MZ_inc <- create_inc_plots(results = nnp_result_list,country='MZ')
NG_inc <- create_inc_plots(results = nnp_result_list,country='NG')

BF_inc_pg <- create_inc_plots(results = nnp_pg_result_list,country='BF')
MZ_inc_pg <- create_inc_plots(results = nnp_pg_result_list,country='MZ')
NG_inc_pg <- create_inc_plots(results = nnp_pg_result_list,country='NG')
BF_inc_pgcorr <- create_inc_plots(results = nnp_pgcorr_result_list,country='BF')
MZ_inc_pgcorr <- create_inc_plots(results = nnp_pgcorr_result_list,country='MZ')
NG_inc_pgcorr <- create_inc_plots(results = nnp_pgcorr_result_list,country='NG')
BF_inc_pgvol <- create_inc_plots(results = nnp_pgvol_result_list,country='BF')
MZ_inc_pgvol <- create_inc_plots(results = nnp_pgvol_result_list,country='MZ')
NG_inc_pgvol <- create_inc_plots(results = nnp_pgvol_result_list,country='NG')

windows(7,5)
BF_inc
MZ_inc
NG_inc
BF_inc_aa + BF_inc
MZ_inc_aa + MZ_inc
NG_inc_aa + NG_inc
BF_inc_pg + BF_inc_pgcorr
MZ_inc_pg + MZ_inc_pgcorr
NG_inc_pg + NG_inc_pgcorr

BF_inc + BF_inc_pgvol
MZ_inc + MZ_inc_pgvol
NG_inc + NG_inc_pgvol

country <- 'MZ'

create_relinc_plots <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    inc_history <- data.frame(t(results[[i]]$history['inc', 51:1000, -1]))
    long_inc_sum <- inc_history%>% mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value)*10000,
                mean=mean(value)*10000,
                upper=quantile(value,probs=0.975)*10000,
                lower=quantile(value,probs=0.025)*10000)%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    long_relinc_sum <- long_inc_sum %>%
      mutate(median=median/max(long_inc_sum$median),
             lower=lower/max(long_inc_sum$upper),
             upper=upper/max(long_inc_sum$upper),
      )
    df <- rbind(df,long_relinc_sum)
    
  }
  
  
  # df_eir <- data.frame(time = integer(),
  #                      median = numeric(),
  #                      mean = numeric(),
  #                      upper = numeric(),
  #                      lower = numeric(),
  #                      distrit = character())
  # 
  # for(i in start:(start+number-1)){
  #   eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))
  #   long_eir_sum <- eir_history%>% mutate(t=c(1:nrow(eir_history)))%>%
  #     melt(id='t')%>%
  #     rename(time=t)%>%
  #     group_by(time)%>%
  #     summarise(median=median(value),
  #               mean=mean(value),
  #               upper=quantile(value,probs=0.975),
  #               lower=quantile(value,probs=0.025))%>%
  #     mutate(district = districts[[i-start+1]],
  #            month = dates_list[[country]])
  #   df_eir <- rbind(df_eir,long_eir_sum)
  # }
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-311'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-31'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
    
  )
  inc_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='Normalised Median Incidence\n(<5 year old)')+
    # labs(x='Date',y='EIR')+
    theme(legend.title = element_blank(),
          # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  # +
  #   coord_cartesian(ylim=c(0, 40))
  inc_plot
}
create_relinc_plots(nnp_result_list,country='BF')
create_relinc_plots(nnp_result_list,country='MZ')
create_relinc_plots(nnp_result_list,country='NG')


##EIR plots
create_eir_plots <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df_eir <- data.frame(time = integer(),
                       median = numeric(),
                       mean = numeric(),
                       upper = numeric(),
                       lower = numeric(),
                       distrit = character())
  
  for(i in start:(start+number-1)){
    eir_history <- data.frame(t(results[[i]]$history['EIR', 51:1000, -1]))
    long_eir_sum <- eir_history%>% mutate(t=c(1:nrow(eir_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value),
                mean=mean(value),
                upper=quantile(value,probs=0.975),
                lower=quantile(value,probs=0.025))%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df_eir <- rbind(df_eir,long_eir_sum)
  }
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(df_eir)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = 1000,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = 1000,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits = c(0,105)),
    MZ = ggplot(df_eir)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-1'), ymin = 0, ymax = 1000,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,400)),
    NG = ggplot(df_eir)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = 1000,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,45))
    
  )
  eir_plot <- annotations[[country]]+
    geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='EIR')+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          # axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  # +
  #   coord_cartesian(ylim=c(0, 40))
  eir_plot
}

BF_eir_pg <- create_eir_plots(results = nnp_pg_result_list,country='BF')
MZ_eir_pg <- create_eir_plots(results = nnp_pg_result_list,country='MZ')
NG_eir_pg <- create_eir_plots(results = nnp_pg_result_list,country='NG')

BF_eir_pgcorr <- create_eir_plots(results = nnp_pgcorr_result_list,country='BF')
MZ_eir_pgcorr <- create_eir_plots(results = nnp_pgcorr_result_list,country='MZ')
NG_eir_pgcorr <- create_eir_plots(results = nnp_pgcorr_result_list,country='NG')

BF_eir_pgvol <- create_eir_plots(results = nnp_pgvol_result_list,country='BF')
NG_eir_pgvol <- create_eir_plots(results = nnp_pgvol_result_list,country='NG')

windows(15,5)
BF_eir_pg + BF_eir_pgcorr
MZ_eir_pg + MZ_eir_pgcorr
NG_eir_pg + NG_eir_pgcorr

BF_eir_pg + BF_eir_pgvol

##Compare posterior distributions
# Each country separately
# Grab posterior distributions for each run type
# Melt dataframe
# Add variable, country name, and run type as new variable
# rbind 4 different run types
# Plot as violin plot (maybe ridge plot), variables as facets, run types as colors

#Burkina Faso
nnp_result_list_aa[1]
names(nnp_result_list_aa)

create_comp_plots <- function(aa,u5,pg,pgcorr,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors <- c(`All ages` = "#1B9E77", `Under 5 years old` = "#999999", `Primigrav` = "#D95F02", `Primigrav Correlation` = "#377EB8")
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]

  df <- data.frame(value = numeric(),
                   variable = character(),
                   run = character(),
                   distrit = character())
  aa[[1]]$mcmc
##All ages results  
  for(i in start:(start+number-1)){
    mcmc <- data.frame(aa[[i]]$mcmc[51:1000,])
    long_mcmc <- mcmc%>% mutate(step=c(1:nrow(mcmc)))%>%
      melt(id='step')%>%
      mutate(district = districts[[i-start+1]],
             run = 'All ages')
    df <- rbind(df,long_mcmc)
    
  }
  ##Under 5 results  
  for(i in start:(start+number-1)){
    mcmc <- data.frame(u5[[i]]$mcmc[51:1000,])
    long_mcmc <- mcmc%>% mutate(step=c(1:nrow(mcmc)))%>%
      melt(id='step')%>%
      mutate(district = districts[[i-start+1]],
             run = 'Under 5 years old')
    df <- rbind(df,long_mcmc)
  }
  ##Primigravs 
  for(i in start:(start+number-1)){
    mcmc <- data.frame(pg[[i]]$mcmc[51:1000,])
    long_mcmc <- mcmc%>% mutate(step=c(1:nrow(mcmc)))%>%
      melt(id='step')%>%
      mutate(district = districts[[i-start+1]],
             run = 'Primigrav')
    df <- rbind(df,long_mcmc)
  }
  ##Primigravs 
  for(i in start:(start+number-1)){
    mcmc <- data.frame(pgcorr[[i]]$mcmc[51:1000,])
    long_mcmc <- mcmc%>% mutate(step=c(1:nrow(mcmc)))%>%
      melt(id='step')%>%
      mutate(district = districts[[i-start+1]],
             run = 'Primigrav Correlation')
    df <- rbind(df,long_mcmc)
  }
  
  df_out <- df%>%
    mutate(run=factor(run,levels=c('All ages','Under 5 years old','Primigrav','Primigrav Correlation')),
           value = ifelse(variable=='log_init_EIR',exp(value),value))
  levels(df_out$variable)[levels(df_out$variable)=='log_init_EIR'] <- 'init_EIR'

  comp_prior_plot <- ggplot(data=df_out[df_out$variable=='log_prior',])+
    geom_violin(aes(x=district,y=value,fill=run))+
    scale_fill_manual(values=colors)+
    labs(y='Log Prior')+
    facet_wrap(~district,scales='free_x')+
    theme(axis.title.x = element_blank())
  comp_likelihood_plot <- ggplot(data=df_out[df_out$variable=='log_likelihood',])+
    geom_violin(aes(x=district,y=value,fill=run))+
    scale_fill_manual(values=colors)+
    labs(y='Log Likelihood')+
    facet_wrap(~district,scales='free_x')+
    theme(axis.title.x = element_blank())
  comp_posterior_plot <- ggplot(data=df_out[df_out$variable=='log_posterior',])+
    geom_violin(aes(x=district,y=value,fill=run))+
    scale_fill_manual(values=colors)+
    labs(y='Log Posterior')+
    facet_wrap(~district,scales='free_x')+
    theme(axis.title.x = element_blank())
  
  probs_plot <- comp_prior_plot/comp_likelihood_plot/comp_posterior_plot +  plot_layout(guides = "collect")
  
  comp_eirsd_plot <- ggplot(data=df_out[df_out$variable=='EIR_SD',])+
    geom_violin(aes(x=district,y=value,fill=run))+
    scale_fill_manual(values=colors)+
    labs(y='EIR Volatility')+
    facet_wrap(~district,scales='free_x')+
    theme(axis.title.x = element_blank())
  comp_eir_plot <- ggplot(data=df_out[df_out$variable=='init_EIR',])+
    geom_violin(aes(x=district,y=value,fill=run))+
    scale_fill_manual(values=colors)+
    labs(y='Initial EIR')+
    facet_wrap(~district,scales='free_x')+
    scale_y_log10()+
    theme(axis.title.x = element_blank())
  params_plot <- comp_eirsd_plot/comp_eir_plot +  plot_layout(guides = "collect")
}
windows(10,15)
bf_comp_plot <- create_comp_plots(aa=nnp_result_list_aa,u5=nnp_result_list,
                                  pg=nnp_pg_result_list,pgcorr=nnp_pgcorr_result_list,
                                  country = 'BF')
##For presentation
create_obsprev_plots <- function(results,data_list=nnp_list,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   distrit = character())
  
  for(i in start:(start+number-1)){
    prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -1]))
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
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
    
  )
  prev_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    # geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    labs(y='ANC Prevalence')+
    facet_grid(~district)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  prev_plot
}
BF_obsprev_pgvol <- create_obsprev_plots(results = nnp_pgvol_result_list,data_list=nnp_pg_list,country='BF')
NG_obsprev_pgvol <- create_obsprev_plots(results = nnp_pgvol_result_list,data_list=nnp_pg_list,country='NG')
BF_obsprev_pgvol / plot_spacer()
NG_obsprev_pgvol / plot_spacer()

create_estprev_plots <- function(results,data_list=nnp_list,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   district = character(),
                   month = character())
  df_sample <- data.frame(time = integer(),
                          value = numeric(),
                          variable = character(),
                          district = character(),
                          month = character())
  for(i in start:(start+number-1)){
    prev_history <- data.frame(t(results[[i]]$history['prev', 51:1000, -1]))
      
    
    long_prev_sum <- prev_history%>%
      mutate(t=c(1:nrow(prev_history)))%>%
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
    
    prev_sample <- prev_history[, sample(ncol(prev_history), 100)] %>%
      mutate(t=c(1:nrow(prev_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      mutate(district = districts[[i-start+1]],
             month = rep(dates_list[[country]],100))
    df_sample <- rbind(df_sample,prev_sample)
    
    
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
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999"),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-30'), ymin = 0, ymax = 1,alpha = .1,fill = "#999999")
    
  )
  prev_plot <- annotations[[country]]+
    geom_line(data=df_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.2)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_vline(xintercept = c(as.Date('2017-7-1'),as.Date('2021-5-1')),size=1,linetype='dashed') +
    geom_point(data=df_data,aes(x=month,y=mean,color=district,group=district),pch = 19,position=position_dodge(width=10))+
    geom_errorbar(data=df_data,aes(x=month,ymin=lower,ymax=upper,color=district,group=district),width = 0,position=position_dodge(width=10))+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    labs(y='ANC Prevalence')+
    facet_grid(~district)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"),
          legend.position = 'none'
    )
  prev_plot
}
windows(10,10)

BF_estprev_pgvol <- create_estprev_plots(results = nnp_pgvol_result_list,data_list=nnp_pg_list,country='BF')
NG_estprev_pgvol <- create_estprev_plots(results = nnp_pgvol_result_list,data_list=nnp_pg_list,country='NG')
BF_estprev_pgvol / plot_spacer()
NG_estprev_pgvol / plot_spacer()
windows(6,7)
BF_estprev_pgvol/BF_inc_pgvol
NG_estprev_pgvol/NG_inc_pgvol

bf_hmis <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Burkina Faso/Routine HMIS/BF_Routine_Mar2022.xlsx')
ng_hmis <- readxl::read_excel('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Nigeria/NNP Nigeria HMIS Data 2019-2022 - LGA.xlsx')
bf_hmis$inc <- bf_hmis$Confirmed/bf_hmis$Population*10000

bf_hmis_sum <- bf_hmis %>%
  group_by(Month,Year,`Net type`)%>%
  summarise(confirmed = sum(Confirmed),
            population = sum(Population))%>%
  mutate(inc = confirmed*10000/population)

bf_hmis_sum$date = as.yearmon(paste(bf_hmis_sum$Year, bf_hmis_sum$Month), "%Y %b")

colors_nets <- c(IG2 = "#1B9E77", Stdr = "#999999", PBO = "#D95F02")
windows(10,4)
ggplot(bf_hmis_sum)+
  geom_line(aes(x=date,y=inc,col=`Net type`),size=1)+
  scale_y_continuous(limits=c(0,900))+
  scale_color_manual(values=colors_nets)

create_inc_plots_threads <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  df <- data.frame(time = integer(),
                   median = numeric(),
                   mean = numeric(),
                   upper = numeric(),
                   lower = numeric(),
                   district = character(),
                   month = character())
  df_sample <- data.frame(time = integer(),
                          value = numeric(),
                          variable = character(),
                          district = character(),
                          month = character())
  for(i in start:(start+number-1)){
    inc_history <- data.frame(t(results[[i]]$history['inc', 51:1000, -1]))
    
    long_inc_sum <- inc_history%>%
      mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      group_by(time)%>%
      summarise(median=median(value)*10000*30,
                mean=mean(value)*10000*30,
                upper=quantile(value,probs=0.975)*10000*30,
                lower=quantile(value,probs=0.025)*10000*30)%>%
      mutate(district = districts[[i-start+1]],
             month = dates_list[[country]])
    df <- rbind(df,long_inc_sum)
    
    inc_sample <- inc_history[, sample(ncol(inc_history), 100)] %>%
      mutate(t=c(1:nrow(inc_history)))%>%
      melt(id='t')%>%
      rename(time=t)%>%
      mutate(value = value*10000*30,
             district = districts[[i-start+1]],
             month = rep(dates_list[[country]],100))
    df_sample <- rbind(df_sample,inc_sample)
    
    
  }
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(df)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits = c(0,3750)),
    MZ = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-1'), ymin = 0, ymax = 400,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,400)),
    NG = ggplot(df)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,3750))+
      coord_cartesian(ylim=c(0, 1000))
    
  )
  inc_plot <- annotations[[country]]+
    geom_line(data=df_sample,aes(x=month,y=value,color=district,group=variable),alpha=0.1)+
    # geom_ribbon(aes(x=month,ymin=lower,ymax=upper,fill=district,group=district),alpha=0.2)+
    geom_line(aes(x=month,y=median,color=district,group=district),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='Estimated Incidence\nper 10,000 person-months')+
    # labs(x='Date',y='EIR')+
    facet_grid(~district)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  # +
  #   coord_cartesian(ylim=c(0, 40))
  inc_plot
}
windows(7,5)

BF_inc_pgvol_threads <- create_inc_plots_threads(results = nnp_pgvol_result_list,country='BF')
MZ_inc_pgvol_threads <- create_inc_plots_threads(results = nnp_pgvol_result_list,country='MZ')
NG_inc_pgvol_threads <- create_inc_plots_threads(results = nnp_pgvol_result_list,country='NG')



create_obsinc_plots <- function(results,country=c('BF','MZ','NG')){
  district_list <- list(BF = c('Banfora','Gaoua','Orodara'),
                        MZ = c('Changara','Chemba','Guro'),
                        NG = c('Asa','Ejigbo','Ife North','Moro'))
  dates_list <- list(BF = seq(as.Date('2020-9-1'),as.Date('2022-5-1'),by='months'),
                     MZ = seq(as.Date('2020-12-1'),as.Date('2021-9-1'),by='months'),
                     NG = seq(as.Date('2020-11-1'),as.Date('2021-12-1'),by='months'))
  colors_list <- list(BF = c(Banfora = "#1B9E77", Gaoua = "#999999", Orodara = "#D95F02"),
                      MZ = c(Changara = "#D95F02", Chemba = "#999999", Guro = "#1B9E77"),
                      NG = c(Asa = "#1B9E77", Ejigbo = "#999999", `Ife North` = "#D95F02", Moro = "#377EB8"))
  
  start_list <- c(BF = 1, MZ = 4, NG = 7)
  number_list <- c(BF = 3, MZ = 3, NG = 4)
  
  districts <- district_list[[country]]
  start <- start_list[[country]]
  number <- number_list[[country]]
  colors <- colors_list[[country]]
  
  # ratio <- 1.5 * max(df$upper)/max(df_eir$median)
  annotations <- list(
    BF = ggplot(results)+
      annotate("rect", xmin = as.Date('2020-9-1'), xmax = as.Date('2020-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      annotate("rect", xmin = as.Date('2021-6-1'), xmax = as.Date('2021-10-1'), ymin = 0, ymax = 3750,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits = c(0,3750))+
      coord_cartesian(ylim=c(0, 2000)),
    MZ = ggplot(results)+
      annotate("rect", xmin = as.Date('2021-1-1'), xmax = as.Date('2021-6-1'), ymin = 0, ymax = 400,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,900)),
    NG = ggplot(results)+
      annotate("rect", xmin = as.Date('2021-7-1'), xmax = as.Date('2021-11-1'), ymin = 0, ymax = 1500,alpha = .1,fill = "#999999")+
      scale_y_continuous(limits=c(0,1500))+
      coord_cartesian(ylim=c(0, 1000))
    
  )
  inc_plot <- annotations[[country]]+
    # geom_ribbon(aes(x=date_ex,ymin=lower,ymax=upper,fill=Distrist,group=Distrist),alpha=0.2)+
    geom_line(aes(x=date_ex,y=inc,color=Distrist,group=Distrist),size=1)+
    # geom_line(data=df_eir,aes(x=month,y=median,color=district,group=district),size=1,linetype='dashed')+
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors)+
    scale_x_date(date_labels = "%b %Y")+
    # scale_y_continuous(limits=c(0,500))+
    # scale_y_log10(breaks=c(.001,0.01,.1,1,10,100,1000),labels=c(.001,0.01,.1,1,10,100,1000))+
    # scale_y_continuous(sec.axis = sec_axis(~ . /ratio, name = "EIR"))+
    labs(x='Date',y='Clinical Incidence\nper 10,000 person-months')+
    # labs(x='Date',y='EIR')+
    facet_grid(~Distrist)+
    theme(legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x=element_text(angle=45, hjust=1, vjust=1),
          axis.ticks.x = element_line(size = 0.5), 
          axis.ticks.length = unit(3, "pt"))
  inc_plot
}
library(epitools)
source('addCIs_inc.R')
bf_hmis$date = as.yearmon(paste(bf_hmis$Year, bf_hmis$Month), "%Y %b")
bf_hmis_forplot <- bf_hmis[bf_hmis$Distrist %in% c('Banfora','Gaoua','Orodara')&bf_hmis$date>=as.yearmon('Sep 2020')&
                     bf_hmis$date<=as.yearmon('May 2022')&!is.na(bf_hmis$Confirmed),]
bf_hmis_forplot <- addCIs_inc(bf_hmis_forplot,bf_hmis_forplot$Confirmed,bf_hmis_forplot$Population)
bf_hmis_forplot <- mutate(bf_hmis_forplot, date_ex = as.Date(date, frac = 0))
BF_obsinc_plot <- create_obsinc_plots(bf_hmis_forplot,'BF')

ng_hmis_forplot <- ng_hmis%>%
  mutate(date=as.Date(period))%>%
  filter(date>=as.Date('2020-11-1')&date<=as.Date('2021-12-1'))%>%
  rename(date_ex = date,
         Distrist = lga,
         inc = incidence)

NG_obsinc_plot <- create_obsinc_plots(ng_hmis_forplot,'NG')

pois.daly(bf_hmis_forplot$Confirmed[1],bf_hmis_forplot$Population[1])
windows(10,10)
BF_inc_pgvol_threads/BF_obsinc_plot + plot_layout(guides = "collect")
NG_inc_pgvol_threads / NG_obsinc_plot + plot_layout(guides = "collect")
