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
                      n_threads = 4,
                      lag_rates = 10,
                      state_check = 0){
  ######## run pMCMC with same model with same log(EIR) random walk but within odin.dust
  
  data <- mcstate::particle_filter_data(data_raw, time = "t", rate = NULL, initial_time = 0)
  
  compare <- function(state, observed, pars = NULL) {
    dbinom(x = observed$positive,
           size = observed$tested,
           prob = state[1,],
           log = TRUE)
  }
  ##in state have a prevalence by age 
  index <- function(info) {
    list(run = c(prev = info$index$prev),
         state = c(prev = info$index$prev,
                   EIR = info$index$EIR_out,
                   inc = info$index$inc,
                   Sout = info$index$Sout,
                   Tout = info$index$Tout,
                   Dout = info$index$Dout,
                   Aout = info$index$Aout,
                   Uout = info$index$Uout,
                   Pout = info$index$Pout))
  }
  
  stochastic_schedule <- seq(from = 30, by = 30, to = 1830)
  
  init_age <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80)
  prop_treated <- 0.4
  het_brackets <- 5
  
  mpl_pf <- model_param_list_create(init_age = init_age,
                                    pro_treated = prop_treated,
                                    het_brackets = het_brackets,
                                    max_EIR = max_EIR,
                                    state_check = state_check,
                                    lag_rates = lag_rates)
  # print(mpl_pf$state_check)
  
  transform <- function(mpl){
    function(theta) {
      
      init_EIR <- exp(theta[["log_init_EIR"]])
      EIR_vol <- theta[["EIR_SD"]]
      mpl <- append(mpl,list(EIR_SD = theta[["EIR_SD"]]))
      
      equilibrium_init_create_stripped(age_vector = mpl$init_age,
                                       EIR = init_EIR,
                                       ft = prop_treated,
                                       model_param_list = mpl,
                                       het_brackets = het_brackets,
                                       state_check = mpl$state_check)
    }
  }
  
  #### NB the volatility and initial EIR is hard-coded in the odinmodelmatchedstoch bw lines 230 and 234###
  model <- odin.dust::odin_dust("shared/odinmodelmatchedstoch.R")
  
  set.seed(1)
  
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
  # EIR_SD <- mcstate::pmcmc_parameter("EIR_SD", 0.3, min = 0,max=2.5,
  #                                    prior = function(p) dexp(p, rate = 5, log = TRUE))
  EIR_SD <- mcstate::pmcmc_parameter("EIR_SD", 0.3, min = 0,max=2.5,
                                     prior = function(p) dlnorm(p, meanlog = -.2, sdlog = 0.5, log = TRUE))
  log_init_EIR <- mcstate::pmcmc_parameter("log_init_EIR", 1.5, min = -8.5, max = 8.5,
                                           prior = function(p) dnorm(p, mean = 0, sd = 10, log = TRUE) + p) #Add p to adjust for sampling on log scale
  
  pars = list(EIR_SD = EIR_SD, log_init_EIR = log_init_EIR)
  
  mcmc_pars <- mcstate::pmcmc_parameters$new(pars,
                                             proposal_matrix,
                                             transform = transform(mpl_pf))
  
  ### Run pMCMC
  start.time <- Sys.time()
  pmcmc_run <- mcstate::pmcmc(mcmc_pars, pf, control = control)
  run_time <- difftime(Sys.time(),start.time,units = 'secs')
  print(run_time)
  pars <- pmcmc_run$pars
  probs <- pmcmc_run$probabilities
  mcmc <- coda::as.mcmc(cbind(probs, pars))
  
  to_return <- list(threads = n_threads,
                    particles = n_particles,
                    run_time = run_time,
                    mcmc = as.data.frame(mcmc),
                    pars = as.data.frame(pars),
                    probs = as.data.frame(probs),
                    history = pmcmc_run$trajectories$state)
  
  return(to_return)
}