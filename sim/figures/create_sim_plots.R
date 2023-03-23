##Show prevalence trajectories, comparing with seasonal deterministic model and simulated data
results <- sim_all_result_list
data_list <- sim_seas_list
create_sim_plots <- function(results,data_list,rainfall){
    mcmc.df <- bind_rows(lapply(1:length(results), 
                                function(x){
                                  print(names(results[x]))
                                  df <- bind_rows(lapply(1:length(results[[x]]), 
                                                         function(y,subresults=results[[x]]){
                                                           print(names(subresults[y]))
                                                           df <- bind_rows(lapply(1:length(subresults[[y]]), 
                                                                                  function(z,subsubresults=subresults[[y]]){
                                                                                    print('We made it!')
                                                                                    print(names(subsubresults[z]))
                                                                                    df <- subsubresults[[z]]$mcmc[251:1000,]
                                                                                    nrow(df)
                                                                                    df$site <- names(subsubresults[z])
                                                                                    df$step <- 1:nrow(df)
                                                                                    
                                                                                    return(df)
                                                                                  }))
                                                           
                                                           df$start <- names(subresults[y])
                                                           return(df)
                                                         }))
                                  df$model <- names(results[x])
                                  return(df)
                                }))
    mcmc.df <- mcmc.df %>% 
      separate_wider_delim(site, "_", names = c("site", "transmission"))
    
    get_history <- function(measure){
      bind_rows(lapply(1:length(results), 
                       function(x){
                         df <- bind_rows(lapply(1:length(results[[x]]), 
                                                function(y,subresults=results[[x]]){
                                                  df <- bind_rows(lapply(1:length(subresults[[y]]), 
                                                                         function(z,subsubresults=subresults[[y]],start=names(subresults[y])){
                                                                           history.df <- as.data.frame(t(subsubresults[[z]]$history[measure, 251:1000, -length(data_list[[z]]$sim_obs_peak$month)]))
                                                                           prev_history <- history.df%>%
                                                                             mutate(t=c(1:nrow(history.df)))%>%
                                                                             melt(id='t')%>%
                                                                             dplyr::rename(time=t)%>%
                                                                             group_by(time)%>%
                                                                             dplyr::summarise(median=median(value),
                                                                                              mean=mean(value),
                                                                                              upper=quantile(value,probs=0.975),
                                                                                              lower=quantile(value,probs=0.025))%>%
                                                                             mutate(site = names(subsubresults[z]))
                                                                           if(start=='peak'){
                                                                             prev_history$month <- as.Date(data_list[[z]]$sim_obs_peak$month)
                                                                             
                                                                           }
                                                                           if(start=='peakplus3'){
                                                                             prev_history$month <- as.Date(data_list[[z]]$sim_obs_peakplus3$month)
                                                                             
                                                                           }
                                                                           if(start=='trough'){
                                                                             prev_history$month <- as.Date(data_list[[z]]$sim_obs_trough$month)
                                                                             
                                                                           }
                                                                           else{
                                                                             prev_history$month <- as.Date(data_list[[z]]$sim_obs_troughplus3$month)
                                                                           }
                                                                           return(prev_history)
                                                                         }))
                                                  
                                                  df$start <- names(subresults[y])
                                                  return(df)
                                                }))
                         df$model <- names(results[x])
                         return(df)
                       }))
    }
    prev.df <- get_history('prev')
    inc.df <- get_history('inc')
    true_val <- bind_rows(lapply(1:length(data_list), 
                                       function(x){
                                         df <- data_list[[x]]$true_val
                                         df$site <- names(data_list[x])
                                         return(df)
                                       }))
    true_val$date <- as.Date(true_val$date)
    true_val$site <- factor(true_val$site,levels = c('casc_low','casc_med','casc_hi','nord_low','nord_med','nord_hi'))
    sim_obs <- bind_rows(lapply(1:length(data_list), 
                                function(x){
                                  df.peak <- data_list[[x]]$sim_obs_peak
                                  df.peak <- addCIs(df.peak,df.peak$positive,df.peak$tested)
                                  df.peak$start <- 'peak'
                                  df.trough <- data_list[[x]]$sim_obs_trough
                                  df.trough <- addCIs(df.trough,df.trough$positive,df.trough$tested)
                                  df.trough$start <- 'trough'
                                  df <- bind_rows(df.peak,df.trough)
                                  df$site <- names(data_list[x])
                                  return(df)
                                }))
                                       
    sim_obs$site <- factor(sim_obs$site,levels = c('casc_low','casc_med','casc_hi','nord_low','nord_med','nord_hi'))
    mcmc.df$site <- factor(mcmc.df$site,levels = c('casc_low','casc_med','casc_hi','nord_low','nord_med','nord_hi'))
    prev.df$site <- factor(prev.df$site,levels = c('casc_low','casc_med','casc_hi','nord_low','nord_med','nord_hi'))
    inc.df$site <- factor(inc.df$site,levels = c('casc_low','casc_med','casc_hi','nord_low','nord_med','nord_hi'))
    mcmc.df$model <- factor(mcmc.df$model,levels = c('standard','seasonal'))
    prev.df$model <- factor(prev.df$model,levels = c('standard','seasonal'))
    inc.df$model <- factor(inc.df$model,levels = c('standard','seasonal'))
    
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
      facet_grid(.~site) +
      theme(
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"),
        axis.title.x = element_blank()
      )
    true_inc_plot <- ggplot(true_val)+
      geom_line(aes(x=date,y=inc*1000),linewidth=1)+
      scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
      coord_cartesian(ylim = c(0,10),xlim=as.Date(c('2014-01-01','2019-12-01')))+
      labs(x='Date',y='Incidence\nper 1,000 people')+
      facet_grid(.~site) +
      theme(
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"),
        axis.title.x = element_blank()
      )
    obs_prev_plot <- ggplot(true_val)+
      geom_line(aes(x=date,y=prev),linewidth=1)+
      geom_point(data=sim_obs,aes(x=date,y=mean),pch = 19,color="#999999")+
      geom_errorbar(data=sim_obs,aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
      scale_x_date(date_labels = "%b %Y",date_breaks = "1 year")+
      coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
      labs(x='Date',y='Prevalence')+
      facet_grid(start~site) +
      theme(
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"),
        axis.title.x = element_blank()
      )

    est_prev_plot_casc_hi <- ggplot(true_val[true_val$site=='casc_hi',])+
      geom_line(aes(x=date,y=prev),linewidth=1)+
      geom_point(data=sim_obs[sim_obs$site=='casc_hi',],aes(x=date,y=mean),pch = 19,color="#999999")+
      geom_errorbar(data=sim_obs[sim_obs$site=='casc_hi',],aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
      # geom_line(data=df_sample[df_sample$model=='Standard'&df_sample$start=='Peak',],aes(x=date,y=value,group=variable),alpha=0.1,color="#377EB8")+
      geom_ribbon(data=prev.df[prev.df$site=='casc_hi',],aes(x=month,ymin=lower,ymax=upper),alpha=0.2,fill="#377EB8")+
      geom_line(data=prev.df[prev.df$site=='casc_hi',],aes(x=month,y=median),linewidth=1,color="#377EB8")+
      scale_x_date(date_labels = "%Y",date_breaks = "1 year")+
      coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
      labs(x='Date',y='Prevalence')+
      facet_grid(start~model) +
      theme(
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"),
        axis.title.x = element_blank()
      )
    est_inc_plot_casc_hi <- ggplot(true_val[true_val$site=='casc_hi',])+
      geom_line(aes(x=date,y=inc*1000),linewidth=1)+
      # geom_line(data=inc.df[inc.df$site=='casc_hi',],aes(x=date,y=value*1000),alpha=0.1,color="#1B9E77")+
      geom_ribbon(data=inc.df[inc.df$site=='casc_hi',],aes(x=month,ymin=lower,ymax=upper),alpha=0.2)+
      geom_line(data=inc.df[inc.df$site=='casc_hi',],aes(x=month,y=median*1000),linewidth=1,color="#1B9E77")+
      scale_x_date(date_labels = "%Y",date_breaks = "1 year")+
      coord_cartesian(ylim = c(0,10),xlim=as.Date(c('2014-01-01','2019-12-01')))+
      labs(x='Date',y='Incidence\nper 1,000 people')+
      facet_grid(start~model) +
      theme(
        axis.ticks.x = element_line(size = 0.5), 
        axis.ticks.length = unit(3, "pt"),
        axis.title.x = element_blank())
    
    
}
prev.df.casc_hi_tr <- prev.df[prev.df$site=='casc_hi'&prev.df$start=='trough',]
inc.df.casc_hi_tr <- inc.df[inc.df$site=='casc_hi'&inc.df$start=='trough',]

sim_std_trough_result_list$casc_hi
prev_history_std <- data.frame(t(sim_std_trough_result_list$casc_hi$history['prev', 501:1000, -length(sim_seas_list$casc_hi$sim_obs_trough$month)]))
prev_history_seas <- data.frame(t(sim_seas_trough_result_list$casc_hi$history['prev', 501:1000, -length(sim_seas_list$casc_hi$sim_obs_trough$month)]))

sample_seas <- sample(ncol(prev_history_seas),100)
matplot(as.matrix(prev_history_std[,sample_seas]))
prev_sample_std <- prev_history_std[,sample_seas] %>%
  mutate(t=c(1:nrow(prev_history_std)))%>%
  melt(id='t')%>%
  rename(time=t)%>%
  mutate(model = 'standard',
         start = 'trough',
         date = rep(as.Date(sim_seas_list$casc_hi$sim_obs_trough$month),100))

prev_sample_seas <- prev_history_seas[, sample_seas] %>%
  mutate(t=c(1:nrow(prev_history_seas)))%>%
  melt(id='t')%>%
  rename(time=t)%>%
  mutate(model = 'seasonal',
         start = 'trough',
         date = rep(as.Date(sim_seas_list$casc_hi$sim_obs_trough$month),100))
prev_sample <- rbind(prev_sample_std,prev_sample_seas)
prev_sample$group <- paste0(prev_sample$variable,prev_sample$model)

inc_history_std <- data.frame(t(sim_std_trough_result_list$casc_hi$history['inc', 501:1000, -length(sim_seas_list$casc_hi$sim_obs_trough$month)]))
inc_history_seas <- data.frame(t(sim_seas_trough_result_list$casc_hi$history['inc', 501:1000, -length(sim_seas_list$casc_hi$sim_obs_trough$month)]))

inc_sample_std <- inc_history_std[, sample_seas] %>%
  mutate(t=c(1:nrow(inc_history_std)))%>%
  melt(id='t')%>%
  rename(time=t)%>%
  mutate(model = 'standard',
         start = 'trough',
         date = rep(as.Date(sim_seas_list$casc_hi$sim_obs_trough$month),100))

inc_sample_seas <- inc_history_seas[, sample_seas] %>%
  mutate(t=c(1:nrow(inc_history_seas)))%>%
  melt(id='t')%>%
  rename(time=t)%>%
  mutate(model = 'seasonal',
         start = 'trough',
         date = rep(as.Date(sim_seas_list$casc_hi$sim_obs_trough$month),100))
inc_sample <- rbind(inc_sample_std,inc_sample_seas)
inc_sample$group <- paste0(inc_sample$variable,inc_sample$model)

prehistory_seas <- data.frame(bind_rows(sim_seas_trough_result_list$casc_hi$seas_history[sample_seas+500],.id='id'))
prehistory_seas$date <- as.Date(min(as.Date(inc_sample$date))-(max(prehistory_seas$t)-prehistory_seas$t))
prehistory_seas$model <- 'seasonal'
prehistory_std <- data.frame(date=prehistory_seas$date,
                             id=prehistory_seas$id)

prehistory_std$inc <- bind_rows(lapply(1:100,function(x){
  inc <- inc_sample_std[inc_sample_std$date==min(inc_sample_std$date),]$value[x]
  return(data.frame(inc=rep(inc,100)))
}))
prehistory_std$inc <- prehistory_std$inc$inc
prehistory_std$prev <- bind_rows(lapply(1:100,function(x){
  prev <- prev_sample_std[prev_sample_std$date==min(prev_sample_std$date),]$value[x]
  return(data.frame(prev=rep(prev,100)))
}))
prehistory_std$prev <- prehistory_std$prev$prev
prehistory_std$model <- 'standard'
prehistory <- bind_rows(prehistory_std,prehistory_seas)
prehistory$group <- paste0(prehistory$id,prehistory$model)
brewer.pal(8,'Paired')
display.brewer.pal('Paired')
windows(10,7)
"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00"
colors <- c(standard = "#E31A1C", seasonal = "#33A02C")
est_prev_plot_casc_hi <- ggplot(true_val[true_val$site=='casc_hi',])+
  geom_line(aes(x=date,y=prev),linewidth=1)+
  geom_point(data=sim_obs[sim_obs$site=='casc_hi',],aes(x=date,y=mean),pch = 19,color="#999999")+
  geom_errorbar(data=sim_obs[sim_obs$site=='casc_hi',],aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
  geom_line(data=prev_sample[prev_sample$model=='standard',],aes(x=date,y=value,group=group,color=model),alpha=0.1)+
  geom_line(data=prehistory[prehistory$model=='standard',],aes(x=date,y=prev,group=group,color=model),alpha=0.1)+
  # geom_ribbon(data=prev.df[prev.df$site=='casc_hi',],aes(x=month,ymin=lower,ymax=upper),alpha=0.2,fill="#377EB8")+
  geom_line(data=prev.df.casc_hi_tr[prev.df.casc_hi_tr$model=='standard',],aes(x=month,y=median,color=as.factor(model)),linewidth=1)+
  scale_x_date(date_labels = "%Y",date_breaks = "1 year")+
  scale_color_manual(values = colors)+
  coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Prevalence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank()
  )
est_prev_plot_casc_hi <- ggplot(true_val[true_val$site=='casc_hi',])+
  geom_line(aes(x=date,y=prev),linewidth=1)+
  geom_point(data=sim_obs[sim_obs$site=='casc_hi',],aes(x=date,y=mean),pch = 19,color="#999999")+
  geom_errorbar(data=sim_obs[sim_obs$site=='casc_hi',],aes(x=date,ymin=lower,ymax=upper),width = 0,color="#999999")+
  geom_line(data=prev_sample,aes(x=date,y=value,group=group,color=model),alpha=0.1)+
  geom_line(data=prehistory,aes(x=date,y=prev,group=group,color=model),alpha=0.1)+
  # geom_ribbon(data=prev.df[prev.df$site=='casc_hi',],aes(x=month,ymin=lower,ymax=upper),alpha=0.2,fill="#377EB8")+
  geom_line(data=prev.df.casc_hi_tr,aes(x=month,y=median,color=as.factor(model)),linewidth=1)+
  scale_x_date(date_labels = "%Y",date_breaks = "1 year")+
  scale_color_manual(values = colors)+
  coord_cartesian(ylim = c(0,0.8),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  labs(x='Date',y='Prevalence')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank()
  )

est_inc_plot_casc_hi <- ggplot(true_val[true_val$site=='casc_hi',])+
  geom_line(aes(x=date,y=inc*1000),linewidth=1)+
  geom_line(data=inc_sample[inc_sample$model=='standard',],aes(x=date,y=value*1000,group=group,color=model),alpha=0.1)+
  # geom_ribbon(data=inc.df[inc.df$site=='casc_hi',],aes(x=month,ymin=lower,ymax=upper),alpha=0.2)+
  geom_line(data=prehistory[prehistory$model=='standard',],aes(x=date,y=inc*1000,group=group,color=model),alpha=0.1)+
  geom_line(data=inc.df.casc_hi_tr[inc.df.casc_hi_tr$model=='standard',],aes(x=month,y=median*1000,color=model),linewidth=1)+
  scale_x_date(date_labels = "%Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,25),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  scale_color_manual(values = colors)+
  labs(x='Date',y='Incidence\nper 1,000 people')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank())
est_inc_plot_casc_hi <- ggplot(true_val[true_val$site=='casc_hi',])+
  geom_line(aes(x=date,y=inc*1000),linewidth=1)+
  geom_line(data=inc_sample,aes(x=date,y=value*1000,group=group,color=model),alpha=0.1)+
  # geom_ribbon(data=inc.df[inc.df$site=='casc_hi',],aes(x=month,ymin=lower,ymax=upper),alpha=0.2)+
  geom_line(data=prehistory,aes(x=date,y=inc*1000,group=group,color=model),alpha=0.1)+
  geom_line(data=inc.df.casc_hi_tr,aes(x=month,y=median*1000,color=model),linewidth=1)+
  scale_x_date(date_labels = "%Y",date_breaks = "1 year")+
  coord_cartesian(ylim = c(0,25),xlim=as.Date(c('2014-01-01','2019-12-01')))+
  scale_color_manual(values = colors)+
  labs(x='Date',y='Incidence\nper 1,000 people')+
  theme(
    axis.ticks.x = element_line(size = 0.5), 
    axis.ticks.length = unit(3, "pt"),
    axis.title.x = element_blank())

windows(10,5)
both_final <- est_prev_plot_casc_hi + est_inc_plot_casc_hi + plot_layout(guides='collect')

table(is.na(mcmc.df$site))
table(is.na(mcmc.df$transmission))
df.na <- mcmc.df[is.na(mcmc.df$site)|is.na(mcmc.df$transmission),]
init_eir_comparison <- ggplot(mcmc.df)+
  geom_violin(aes(x=start,y=exp(log_init_EIR),fill=model))+
  facet_grid(factor(transmission,levels=c('low','med','hi'))~site)+
  scale_y_log10()
vol_eir_comparison <- ggplot(mcmc.df)+
  geom_violin(aes(x=start,y=EIR_SD,fill=model))+
  facet_grid(factor(transmission,levels=c('low','med','hi'))~site)
post_comparison <- ggplot(mcmc.df)+
  geom_violin(aes(x=start,y=log_posterior,fill=model))+
  facet_grid(factor(transmission,levels=c('low','med','hi'))~site)
likelihood_comparison <- ggplot(mcmc.df)+
  geom_violin(aes(x=start,y=log_likelihood,fill=model))+
  facet_grid(factor(transmission,levels=c('low','med','hi'))~site)
