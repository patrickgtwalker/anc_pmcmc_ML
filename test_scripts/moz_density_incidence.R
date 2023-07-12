source('shared/equilibrium-init-create-stripped.R')
source('shared/model_parameters.R')

#Set default ggplot theme
theme_set(theme_minimal()+
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  strip.background = element_blank(),
                  panel.border = element_rect(colour = "black",fill=NA),
                  legend.position = 'bottom'))

##Set initial conditions and parameters to match usual runs
ages <-  c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
prop_treated <- 0.4
het_brackets <- 5

mpl_pf <- model_param_list_create(init_age = ages,
                                  prop_treated = prop_treated,
                                  het_brackets = het_brackets,
                                  max_EIR = 1000,
                                  state_check = 0,
                                  lag_rates = 10,
                                  country = 'Burkina Faso',
                                  admin_unit = 'Cascades',
                                  start_stoch = as.Date('2015-01-01'),
                                  time_origin = as.Date('2010-01-01'),
                                  seasonality_on = 0,
                                  EIR_SD = 1)

##Run equilibrium function at different EIR values and compare:
## EIR, mosquito density, incidence, and prevalence

EIRs <- c(c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000),c(200,300,400,600,700,800,900))

#Function to get relevant values from equilibrium solution
get_eq <- function(EIR,mpl = mpl_pf){
  state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                            init_EIR = EIR,
                                            ft = mpl$prop_treated,
                                            model_param_list = mpl,
                                            het_brackets = mpl$het_brackets,
                                            state_check = mpl$state_check)
  data.frame(init_EIR = EIR,
             mv0 = state$mv0,
             prev05 = state$prev05,
             inc05 = state$inc05)
}

##Loop equilibrium function over range of EIR values
eq_list <- lapply(EIRs, get_eq)
##Create dataframe from list of equilibrium output values
eq_df <- do.call('rbind', eq_list)

windows(10,10)
ylim.prim <- c(0, 1)   # prev limits
ylim.sec <- c(0, 0.01)    # incidence limits

##Transform incidence to same scale as prevalence
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

ggplot(eq_df)+
  # geom_line(aes(x=init_EIR,y=mv0),color="#1582AD",size=1)+
  geom_line(aes(x=init_EIR,y=prev05),color="#EFBB12",size=1)+
  geom_line(aes(x=init_EIR,y=inc05*b+a),color="#CE5126",size=1)+
  scale_y_continuous("Prevalence <5yo", sec.axis = sec_axis(~ (. - a)/b, name = "Incidence <5yo"),limits = c(0,NA)) 
  
###Run seasonal model at different levels of EIR
season_model <- odin::odin("shared/odin_model_stripped_seasonal.R")
init_EIR <- 10

## Run equilibrium function
mpl <- mpl_pf
get_seas_eq <- function(EIR,mpl = mpl_pf){
  state <- equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                            init_EIR = EIR,
                                            ft = mpl$prop_treated,
                                            model_param_list = mpl,
                                            het_brackets = mpl$het_brackets,
                                            state_check = mpl$state_check)
  state_use <- state[names(state) %in% coef(season_model)$name]
  # create model with initial values
  mod <- season_model$new(user = state_use, use_dde = TRUE)
  
  # tt <- c(0, preyears*365+as.integer(difftime(mpl$start_stoch,mpl$time_origin,units="days")))
  tt <- seq(0, 10*365,length.out=5*365)
  
  # run seasonality model
  mod_run <- mod$run(tt, verbose=FALSE,step_size_max=9)
  out <- mod$transform_variables(mod_run)

  data.frame(t = out$t,
             init_EIR = rep(EIR,length(out$t)),
             EIR = out$EIR_init[,1,1],
             EIR_rel = out$EIR_init[,1,1]/max(out$EIR_init[,1,1]),
             mv0 = out$mv_init,
             mv0_rel = out$mv_init/max(out$mv_init),
             prev05 = out$prev,
             prev_all = out$prev_all,
             inc05 = out$inc05,
             inc_all = out$inc)
}

EIRs <- c(0.01,0.1,1,5,10,50,100,500,1000)

seas_eq_list <- lapply(EIRs, get_seas_eq)
seas_eq_df <- do.call('rbind', seas_eq_list)

table(seas_eq_df$init_EIR)
ggplot(seas_eq_df)+
  geom_line(aes(x=t,y=prev05),color="#EFBB12",size=1)+
  facet_wrap(.~init_EIR)
ggplot(seas_eq_df)+
  geom_line(aes(x=t,y=inc05),color="#CE5126",size=1)+
  facet_wrap(.~init_EIR)

ylim.prim <- c(0, 1)   # prev limits
ylim.sec <- c(0, 0.015)    # incidence limits

##Transform incidence to same scale as prevalence
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

ggplot(seas_eq_df)+
  geom_line(aes(x=t,y=inc05*b+a),color="#CE5126",size=1)+
  geom_line(aes(x=t,y=prev05),color="#999999",size=1)+
  facet_wrap(.~init_EIR)+
  scale_x_continuous(limits = c(10*365-2*365,NA))+
  scale_y_continuous("Prevalence <5yo", sec.axis = sec_axis(~ (. - a)/b, name = "Incidence <5yo"),limits = c(0,NA)) 

ggplot(seas_eq_df)+
  geom_line(aes(x=t,y=mv0),color="#999999",size=1)+
  facet_wrap(.~init_EIR)+
  scale_x_continuous(limits = c(10*365-2*365,NA))
ggplot(seas_eq_df)+
  geom_line(aes(x=t,y=EIR),color="#999999",size=1)+
  facet_wrap(.~init_EIR)+
  scale_x_continuous(limits = c(10*365-2*365,NA))

##just look at shapes.
ylim.prim <- c(0, 1)   # prev limits
ylim.sec <- c(0, 1)    # incidence limits

##Transform incidence to same scale as prevalence
b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

windows(10,10)
ggplot(seas_eq_df[seas_eq_df$init_EIR==100,])+
  geom_line(aes(x=t,y=EIR),color="#999999",size=1)+
  geom_line(aes(x=t,y=mv0*b+a),color="#1F78B4",size=1)+
  # facet_wrap(.~init_EIR,scales='free_y')+
  scale_x_continuous(limits = c(10*365-2*365,NA))+
  scale_y_continuous("EIR", sec.axis = sec_axis(~ (. - a)/b, name = "Mosquito Density"),limits = c(0,NA)) 

ggplot(seas_eq_df[seas_eq_df$init_EIR>=5,])+
  geom_line(aes(x=t,y=EIR_rel),color="#999999",size=1)+
  geom_line(aes(x=t,y=mv0_rel*b+a),color="#1F78B4",size=1)+
  facet_wrap(.~init_EIR)+
  scale_x_continuous(limits = c(10*365-2*365,NA))+
  scale_y_continuous("Relative EIR", sec.axis = sec_axis(~ (. - a)/b, name = "Relative Mosquito Density"),limits = c(0,NA)) 


ggplot(seas_eq_df)+
  geom_line(aes(x=t,y=mv0_rel,color=as.factor(init_EIR)),size=1)+
  scale_x_continuous(limits = c(10*365-2*365,NA))+
  scale_y_continuous("Relative Mosquito Density",limits = c(0,NA)) 

last_two <- seas_eq_df%>%
  filter(t>=10*365-2*365)%>%
  group_by(init_EIR)%>%
  mutate(mv0_rel = mv0/max(mv0),
         EIR_rel = EIR/max(EIR))
ggplot(last_two)+
  geom_line(aes(x=t,y=EIR_rel,color=as.factor(init_EIR)),size=1)+
  scale_x_continuous(limits = c(10*365-2*365,NA))+
  scale_y_continuous("Relative Mosquito Density",limits = c(0,NA)) 

"#CE5126"