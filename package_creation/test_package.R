library('devtools')
##Install these from github, if needed
remotes::install_github("mrc-ide/mcstate", upgrade = TRUE)
remotes::install_github("mrc-ide/dust", upgrade = TRUE)
remotes::install_github("mrc-ide/mode", upgrade = TRUE)
remotes::install_github("mrc-ide/odin.dust", upgrade = TRUE)
remotes::install_github("mrc-ide/didehpc", upgrade = TRUE)
#devtools::install_github('klutometis/roxygen')
#library(roxygen2)

## Load package from github
devtools::install_github('jt-hicks/sifter@issue-5',force = TRUE)
library(sifter)

test_run_sifter <- sifter::run_pmcmc(data_raw = sifter::data_sim, #I've added data_sim to the package for an easy test
                                     n_particles = 200,
                                     proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                     max_EIR=1000,
                                     max_steps = 1e7,
                                     atol = 1e-6,
                                     rtol = 1e-6,
                                     n_steps = 100,
                                     n_threads = 2,
                                     lag_rates = 10,
                                     country = 'Burkina Faso',
                                     admin_unit = 'Cascades',
                                     seasonality_on = 1,
                                     state_check = 0,
                                     seasonality_check = 0,
                                     stoch_param = 'betaa')
