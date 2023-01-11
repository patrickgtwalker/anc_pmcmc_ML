
run_pmcmc <- function(data_raw,
                      n_particles=200,
                      proposal_matrix,
                      max_EIR=1000,
                      # EIR_vol,
                      # proposal_dist,
                      # init_EIR = 100,
                      max_steps = 1e7,
                      atol = 1e-3,
                      rtol = 1e-6,
                      n_steps = 500,
                      n_threads = 4){
  ######## run pMCMC with same model with same log(EIR) random walk but within odin.dust

  data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)
  
  # compare <- function(state, observed, pars = NULL) {
  #   positive_by_age <- unlist(observed[grep('positive_',names(observed))])
  #   tested_by_age <- unlist(observed[grep('tested_',names(observed))])
  # 
  #   sapply(1:ncol(state), function(i) sum(dbinom(x=positive_by_age,size=tested_by_age,prob=state[,i],log=TRUE),na.rm = TRUE))
  # }
  # 
  compare <- function(state, observed, pars = NULL) {
    positive_by_age <- unlist(observed[grep('positive_cmis',names(observed))])
    tested_by_age <- unlist(observed[grep('tested_cmis',names(observed))])
    
    positive_by_age_g1 <- unlist(observed[grep('positive_g1_',names(observed))])
    tested_by_age_g1 <- unlist(observed[grep('tested_g1_',names(observed))])
    positive_by_age_g2_3 <- unlist(observed[grep('positive_g2_',names(observed))])
    tested_by_age_g2_3 <- unlist(observed[grep('tested_g2_',names(observed))])
    positive_by_age_g3p <- unlist(observed[grep('positive_g3_',names(observed))])
    tested_by_age_g3p <- unlist(observed[grep('tested_g3_',names(observed))])
    
    # CBA_mat<-state[first_cat:second_cat,]
    g1_OR <- 1.1422
    g2_3_OR <- 0.8599
    g3p_OR <- 0.9137
    g1_mat<-change_prev_by_OR(state[5:24,],g1_OR)
    g2_3_mat<-change_prev_by_OR(state[5:24,],g2_3_OR)
    g3p_mat<-change_prev_by_OR(state[5:24,],g3p_OR)
    
    sapply(1:ncol(state), function(i){
      sum(dbinom(x=positive_by_age,size=tested_by_age,prob=state[,i],log=TRUE),na.rm = TRUE)+
        sum(dbinom(x=positive_by_age_g1,size=tested_by_age_g1,prob=g1_mat[,i],log=TRUE),na.rm = TRUE)+
        sum(dbinom(x=positive_by_age_g2_3,size=tested_by_age_g2_3,prob=g2_3_mat[,i],log=TRUE),na.rm = TRUE)+
        sum(dbinom(x=positive_by_age_g3p,size=tested_by_age_g3p,prob=g3p_mat[,i],log=TRUE),na.rm = TRUE)
    })
  }
  
  change_prev_by_OR<-function(prev,OR){
    get_odds<-prev/(1-prev)*OR
    return(get_odds/(1+get_odds))
  }
  ##Simple compare function for troubleshooting
  # compare <- function(state, observed, pars = NULL) {
  #   dbinom(x = observed$positive_1,
  #          size = observed$tested_1,
  #          prob = state[1,],
  #          log = TRUE)
  # }
  
  
  index <- function(info) {
    list(run = c(prev = info$index$prevout),
         state = c(prev = info$index$prev_all,
                   EIR = info$index$EIR_out,
                   inc = info$index$inc,
                   prevage = info$index$prevout))
  }
  
  stochastic_schedule <- seq(from = 30, by = 30, to = max(data$t_end))
  
  
  transform <- function(theta) {
    init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 4, 5, 7.5, 10, 11, 12, 13, 14, 15,
                  16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 50, 60, 70, 80)
    age_groups <- c(0, 4, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 80)
    ag_dim <- as.double(length(age_groups))

    age_convert <- sapply(1:length(age_groups),function(i) as.numeric(which(init_age == age_groups[i])[1]))
    # age_convert <- c(0, 4, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 80)
    prop_treated <- 0.4
    het_brackets <- 5
    
    init_EIR <- exp(theta[["log_init_EIR"]])
    EIR_vol <- theta[["EIR_SD"]]
    mpl_pf <- model_param_list_create(EIR_SD=EIR_vol,max_EIR=max_EIR,ag_dim=ag_dim,age_convert=age_convert)
    equilibrium_init_create_stripped(age_vector = init_age,
                                     EIR = init_EIR,
                                     ft = prop_treated,
                                     model_param_list = mpl_pf,
                                     het_brackets = het_brackets)
  }
  
  #### NB the volatility and initial EIR is hard-coded in the odinmodelmatchedstoch bw lines 230 and 234###
  model <- odin.dust::odin_dust("original_malaria/odinmodelmatchedstoch.R")
  
  #set.seed(1)
  
  ### Set particle filter
  pf <- mcstate::particle_filter$new(data, model, n_particles, compare,
                                     index = index, seed = 1L,
                                     stochastic_schedule = stochastic_schedule,
                                     ode_control = mode::mode_control(max_steps = max_steps, atol = atol, rtol = rtol),
                                     n_threads = n_threads)
  
  ### Set pmcmc control
  control <- mcstate::pmcmc_control(
    n_steps,
    save_state = TRUE,
    save_trajectories = TRUE,
    progress = TRUE,
    n_chains = 1,
    n_workers = 1,
    n_threads_total = n_threads,
    rerun_every = 50,
    rerun_random = TRUE)
  
  ### Set pmcmc parameters
  EIR_SD <- mcstate::pmcmc_parameter("EIR_SD", 0.3, min = 0,max=2.5,
                                     prior = function(p) dexp(p, rate = 5, log = TRUE))
  log_init_EIR <- mcstate::pmcmc_parameter("log_init_EIR", 1.5, min = -8.5, max = 8.5,
                                       prior = function(p) dnorm(p, mean = 0, sd = 10, log = TRUE) + p) #Add p to adjust for sampling on log scale
  
  pars = list(EIR_SD = EIR_SD, log_init_EIR = log_init_EIR)

  mcmc_pars <- mcstate::pmcmc_parameters$new(pars,
                                             proposal_matrix,
                                             transform = transform)

  ### Run pMCMC
  start.time <- Sys.time()
  pmcmc_run <- mcstate::pmcmc(mcmc_pars, pf, control = control)
  print(Sys.time()-start.time)
  return(pmcmc_run)
}