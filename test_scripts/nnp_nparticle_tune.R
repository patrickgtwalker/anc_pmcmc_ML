source('./test_scripts/run_pmcmc_mg_particle.R')


1,4,8
country <- c('Burkina Faso','Burkina Faso','Burkina Faso',
             'Mozambique','Mozambique','Mozambique',
             'Nigeria','Nigeria','Nigeria','Nigeria')
admin <- c('Cascades','Sud-Ouest','Haut-Bassins',
           'Tete','Sofala','Manica',
           'Kwara','Osun','Osun','Kwara')
i <- 1
x<- 2
View(nnp_mgcorr_bulk_std_results)
fixed_pars <- as.data.frame(t(sapply(c(1,4,8), function(x){
  return(nnp_mgcorr_bulk_std_results[[x]]$mcmc[1000,c('EIR_SD','log_init_EIR')])
})))
particle_test <- lapply(1:3, function(i) {
  index <- c(1,4,8)
  list <- lapply(c(2,4,16,64,128),function(x,data_pg,data_mg){
    test <- run_pmcmc_mg_particle(data_raw_pg = data_pg[[index[i]]],
             data_raw_mg = data_mg[[index[i]]],
             n_particles = x,
             proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
             max_EIR=1000,
             max_steps = 1e7,
             atol = 1e-5,
             rtol = 1e-6,
             n_steps = 100,
             n_threads = 1,
             lag_rates = 10,
             seasonality_on = 0,
             state_check = 0,
             fixed_EIR_SD = fixed_pars[[i,1]],
             fixed_log_init_EIR = fixed_pars[[i,2]])
    return(var(test$probs$log_likelihood))
    },data_pg=nnp_pg_list,data_mg=nnp_mg_list)
  names(list) <- c(2,4,16,64,128)
  return(list)
  })
saveRDS(particle_test,'./test_scripts/particle_tune_010323.RDS')
saveRDS(particle_test,'./test_scripts/particle_tune_020323.RDS')

banfora.particle.var <- tibble::rownames_to_column(as.data.frame(t(bind_rows(particle_test[[1]]))),'n')%>%
  mutate(site='Banfora')
ggplot(banfora.particle.var,aes(x=as.numeric(n),y=V1))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  geom_hline(yintercept = 1)+
  scale_x_log10()
changara.particle.var <- tibble::rownames_to_column(as.data.frame(t(bind_rows(particle_test[[2]]))),'n')%>%
  mutate(site='Changara')
ejigbo.particle.var <- tibble::rownames_to_column(as.data.frame(t(bind_rows(particle_test[[3]]))),'n')%>%
  mutate(site='Ejigbo')
all.particle.var <- bind_rows(banfora.particle.var,changara.particle.var,ejigbo.particle.var)
ggplot(changara.particle.var,aes(x=as.numeric(n),y=V1))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  geom_hline(yintercept = 1)+
  scale_x_log10()
ggplot(ejigbo.particle.var,aes(x=as.numeric(n),y=V1))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE)+
  geom_hline(yintercept = 1)+
  scale_x_log10()+
  scale_y_log10()

windows(10,5)
ggplot(all.particle.var,aes(x=as.numeric(n),y=V1))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE, fullrange=TRUE)+
  geom_hline(yintercept = 1)+
  scale_x_log10(limits=c(1,500),breaks = c(1,2,4,8,16,32,64,128,256,512))+
  scale_y_log10(breaks = c(0.1,1,10,100,1000),labels = c('0.1','1','10','100','1000'))+
  facet_grid(.~site)+
  labs(x='Number of Particles',y='Var of Log Likelihood')
