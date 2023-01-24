data_gen <- function(EIR_volatility,
                     init_EIR=100,
                     max_EIR=1000,
                     model_file="shared/odin_model_stripped_matched.R"){
cat('EIR_vol = ',EIR_volatility,' init_EIR = ',init_EIR,'\n')
init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)

prop_treated <- 0.4
rA_preg <- 0.00512821
rU_preg <- 0.00906627
het_brackets <- 5

################## generate the data ######################
# generate random walk of EIR (recursive fn)
genRandWalk <- function(x,vol,randWalk) {
  if (x == 0)    return (randWalk)
  else return(genRandWalk(x-1,vol,c(randWalk,min(exp(log(randWalk[length(randWalk)])+rnorm(1)*vol),max_EIR))))
}
EIR_times<-seq(0,1800,by=30)

### just a random walk on logscale
EIR_vals=genRandWalk(length(EIR_times)-1,EIR_volatility,init_EIR)

##set up the simulation for the simualted data 
time<- 5*365
out_step=30

mpl <- model_param_list_create(init_EIR = init_EIR,
                               init_ft = prop_treated,
                               EIR_times=EIR_times,
                               EIR_vals=EIR_vals
                               )

pars <- equilibrium_init_create_stripped(age_vector = init_age,
                                init_EIR = init_EIR,
                                ft = prop_treated,
                                model_param_list = mpl,
                                het_brackets = het_brackets)

##The malaria model but only on human side (i.e. no mosquitoes to worry about)
generator <- odin(model_file)
state_use <- pars[names(pars) %in% coef(generator)$name]

# create model with initial values
mod <- generator(user = state_use, use_dde = TRUE)
tt <- seq(0, time, out_step)

# run the simulation to base the data
start.time <- Sys.time()
mod_run <- mod$run(tt, step_max_n = 1e7,
                   atol = 1e-5,
                   rtol = 1e-5)
print(Sys.time()-start.time)

# shape output
out <- mod$transform_variables(mod_run)

# plot data and generate data
plot(out$t,out$prev,col="white")
lines(out$t,out$prev,col="blue",lwd=4)
tested<-round(rnorm(length(out$prev_all),220,40))
positive<-rbinom(length(out$prev_all),tested,out$prev_all)
month <- seq.Date(from = as.Date('2015-01-01'),by = 'month',length.out = length(tt))
data_raw<-data.frame(t=out$t+30,
                     month=as.yearmon(month),
                     tested=tested,
                     positive=positive,
                     prev_true=out$prev_all,
                     EIR_true=EIR_vals,
                     vol_true=EIR_volatility,
                     inc_true=out$incunder5)
return(data_raw)
}
