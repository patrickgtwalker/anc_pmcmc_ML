## MODEL VARIABLES
##------------------------------------------------------------------------------
smc_cov<-user()
dim(cov_vect) <- c(2)
cov_vect[1]<-1-smc_cov
cov_vect[2]<-smc_cov
dim(on_vect) <- c(na)
on_vect[2:age05]<-1
on_vect[1]<-0
on_vect[(1+age05):na]<-0
na <- user() # number of age categories
nh <- user() # number of biting heterogeneity categories
ft <- user() # proportion of cases treated

##------------------------------------------------------------------------------
##################
## HUMAN STATES ##
##################
##------------------------------------------------------------------------------

# Human states as specified in full transmission model
# http://journals.plos.org/plosmedicine/article?id=10.1371%2Fjournal.pmed.1000324
# http://www.nature.com/articles/ncomms4136

# fitted parameters for human compartments:
eta <- user() # death rate for exponential population distribution
age_rate[] <- user() # rate at which humans move through age categories
dim(age_rate) <- na
het_wt[] <- user() # weights of the different heterogeneous biting categories
dim(het_wt) <- nh
rA <- user() # rate of movement from A -> U
rT <- user() # rate of treatment working: T -> P
rD <- user() #  rate from D -> A
rU <- user() # rate of clearance of subpatent infection U -> S
rP <- user() # rate at which prophylaxis wears off P -> S

# S - SUSCEPTIBLE
init_S[,] <- user()
dim(init_S) <- c(na,nh)
initial(S[,,]) <- init_S[i,j]*cov_vect[k]
dim(S) <- c(na,nh,2)

deriv(S[1, 1:nh,1:2]) <- -effective_FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] +
  eta*H*het_wt[j]*cov_vect[k] - (eta+age_rate[i])*S[i,j,k]
deriv(S[2:na, 1:nh,1:2]) <- -effective_FOI[i,j,k]*S[i,j,k] + rP*P[i,j,k] + rU*U[i,j,k] -
  (eta+age_rate[i])*S[i,j,k] + age_rate[i-1]*S[i-1,j,k]

# T- SUCCESSFULLY TREATED
init_T[,] <- user()
dim(init_T) <- c(na,nh)
initial(T[,,]) <- init_T[i,j]*cov_vect[k]
dim(T) <- c(na,nh,2)
deriv(T[1, 1:nh,1:2]) <- ft*clin_inc[i,j,k] - rT*T[i,j,k] -
  (eta+age_rate[i])*T[i,j,k]
deriv(T[2:na, 1:nh,1:2]) <- ft*clin_inc[i,j,k] - rT*T[i,j,k] -
  (eta+age_rate[i])*T[i,j,k] + age_rate[i-1]*T[i-1,j,k]

# D - CLEAR DISEASE
init_D[,] <- user()
dim(init_D) <- c(na,nh)
initial(D[,,]) <- init_D[i,j]*cov_vect[k]
dim(D) <- c(na,nh,2)

deriv(D[1, 1:nh,1:2]) <- (1-ft)*clin_inc[i,j,k] - rD*D[i,j,k] -
  (eta+age_rate[i])*D[i,j,k]
deriv(D[2:na, 1:nh,1:2]) <- (1-ft)*clin_inc[i,j,k] - rD*D[i,j,k] -
  (eta+age_rate[i])*D[i,j,k] + age_rate[i-1]*D[i-1,j,k]

# A - ASYMPTOMATIC DISEASE
init_A[,] <- user()
dim(init_A) <- c(na,nh)
initial(A[,,]) <- init_A[i,j]*cov_vect[k]
dim(A) <- c(na,nh,2)

deriv(A[1, 1:nh,1:2]) <- (1-phi[i,j,k])*effective_FOI[i,j,k]*Y[i,j,k] - effective_FOI[i,j,k]*A[i,j,k] +
  rD*D[i,j,k] - rA*A[i,j,k] - (eta+age_rate[i])*A[i,j,k]
deriv(A[2:na, 1:nh,1:2]) <- (1-phi[i,j,k])*effective_FOI[i,j,k]*Y[i,j,k] - effective_FOI[i,j,k]*A[i,j,k] +
  rD*D[i,j,k] - rA*A[i,j,k] - (eta+age_rate[i])*A[i,j,k] + age_rate[i-1]*A[i-1,j,k]

# U - SUBPATENT DISEASE
init_U[,] <- user()
dim(init_U) <- c(na,nh)
initial(U[,,]) <- init_U[i,j]*cov_vect[k]
dim(U) <- c(na,nh,2)

deriv(U[1, 1:nh,1:2]) <- rA*A[i,j,k] - effective_FOI[i,j,k]*U[i,j,k] - rU*U[i,j,k] -
  (eta+age_rate[i])*U[i,j,k]
deriv(U[2:na, 1:nh,1:2]) <- rA*A[i,j,k] - effective_FOI[i,j,k]*U[i,j,k] - rU*U[i,j,k] -
  (eta+age_rate[i])*U[i,j,k] + age_rate[i-1]*U[i-1,j,k]

# P - PROPHYLAXIS
init_P[,] <- user()
dim(init_P) <- c(na,nh)
initial(P[,,]) <- init_P[i,j]*cov_vect[k]
dim(P) <- c(na,nh,2)

deriv(P[1, 1:nh,1:2]) <- rT*T[i,j,k] - rP*P[i,j,k] - (eta+age_rate[i])*P[i,j,k]
deriv(P[2:na, 1:nh,1:2]) <- rT*T[i,j,k] - rP*P[i,j,k] - (eta+age_rate[i])*P[i,j,k] +
  age_rate[i-1]*P[i-1,j,k]


# The number of individuals able to acquire clinical malaria
dim(Y) <- c(na,nh,2)
Y[1:na, 1:nh,1:2] <- S[i,j,k]+A[i,j,k]+U[i,j,k]

# The number of new cases at this timestep
dim(clin_inc) <- c(na,nh,2)
clin_inc[1:na, 1:nh,1:2] <- phi[i,j,k]*effective_FOI[i,j,k]*Y[i,j,k]
output(clin_inc)<-clin_inc
output(phi)<-phi
#output(FOI)<-FOI
output(Y)<-Y
# Sum compartments over all age, heterogeneity and intervention categories
Sh <- sum(S[,,])
Th <- sum(T[,,])
Dh <- sum(D[,,])
Ah <- sum(A[,,])
Uh <- sum(U[,,])
Ph <- sum(P[,,])
H <- Sh + Th + Dh + Ah + Uh + Ph
output(pop)<-H
##------------------------------------------------------------------------------
#####################
## IMMUNITY STATES ##
#####################
##------------------------------------------------------------------------------

# See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6

# ICM - Maternally acquired immunity acquired by babies from mothers (assumed to be proportional to the immunity of a 15 to 30 year old woman)
# ICA - Immunity acquired due to exposure to previous infection, increases with age
# IC - Clinical immunity. Upon infection, immunity against clinical case. IC = ICA + ICM
# IB - Infection blocking immunity, chances of preventing infection upon receiving infectious bite
# ID - Detection immunity, when immunity suppresses parasite densities this makes it less likely that diagnostics will detect parasite infection

# fitted immunity parameters:
uCA <- user() # scale parameter (see Supplementary mats. 3.1.2)
dCA <- user() # decay for clinical immunity
dB <- user() # decay for infection blocking immunity
uB <- user() # scale param for IB immunity
dID <- user() # decay for detection immunity
uD <- user() # scale param for ID immunity
x_I[] <- user() # intermediate variable for calculating immunity functions
dim(x_I) <- na
age20l <- user(integer=TRUE) # lower index of age 20 age compartment
age20u <- user(integer=TRUE) # upper index of age 20 age compartment
age_20_factor <- user() # factor calculated in equilibrium solution
PM <- user() # immunity constant

# ICM - maternally acquired immunity
dim(init_ICM_pre) <- c(nh,2)
init_ICM_pre[1:nh,1:2] <- PM*(ICA[age20l,i,j] + age_20_factor*(ICA[age20u,i,j]-ICA[age20l,i,j]))
ICM_age[]<-user()
dim(ICM_age)<-na
dim(ICM) <- c(na,nh,2)
ICM[1:na, 1:nh,1:2]<-ICM_age[i]*init_ICM_pre[j,k]


# ICA - exposure driven immunity
init_ICA[,] <- user()
dim(init_ICA) <- c(na,nh)
initial(ICA[1:na,1:nh,1:2]) <- init_ICA[i,j]
dim(ICA) <- c(na,nh,2)
deriv(ICA[1, 1:nh,1:2]) <- effective_FOI[i,j,k]/(effective_FOI[i,j,k]*uCA + 1) - ICA[i,j,k]/dCA - ICA[i,j,k]/x_I[i]
deriv(ICA[2:na, 1:nh,1:2]) <- effective_FOI[i,j,k]/(effective_FOI[i,j,k]*uCA + 1) - ICA[i,j,k]/dCA - (ICA[i,j,k]-ICA[i-1,j,k])/x_I[i]

# clinical immunity is a combination of maternal and exposure-driven immunity
dim(IC) <- c(na,nh,2)
IC[,,] <- ICM[i,j,k] + ICA[i,j,k]

# phi - probability of clinical disease, dependent on clinical immunity
phi0 <- user()
phi1 <- user() # these parameters characterise the hill function
IC0 <- user() # for probability of clinical disease
kC <- user() # See supplementary materials 1.1.3
dim(phi) <- c(na,nh,2)
phi[1:na,1:nh,1:2] <- phi0*((1-phi1)/(1+(IC[i,j,k]/IC0)^kC) + phi1)

# IB - infection blocking immunity
init_IB[,] <- user()
dim(init_IB) <- c(na,nh)
initial(IB[1:na,1:nh,1:2]) <- init_IB[i,j]
dim(IB) <- c(na,nh,2)

deriv(IB[1, 1:nh,1:2]) <- EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - IB[i,j,k]/x_I[i]
deriv(IB[2:na, 1:nh,1:2]) <- EIR[i,j,k]/(EIR[i,j,k]* uB + 1) - IB[i,j,k]/dB - (IB[i,j,k]-IB[i-1,j,k])/x_I[i]

# b - probability of disease from infectious bite, depends on infection blocking immunity
b0 <- user() # these parameters characterise the hill function for b
b1 <- user() # prob of infection from bite with zero immunity
kB <- user() #
IB0 <- user()
dim(b) <- c(na,nh,2)
b[1:na, 1:nh,1:2] <- b0 * ((1-b1)/(1+(IB[i,j,k]/IB0)^kB)+b1)

# detection immunity
init_ID[,] <- user()
dim(init_ID) <- c(na,nh)
initial(ID[1:na,1:nh,1:2]) <- init_ID[i,j]
dim(ID) <- c(na,nh,2)

deriv(ID[1, 1:nh,1:2]) <- effective_FOI[i,j,k]/(effective_FOI[i,j,k]*uD + 1) - ID[i,j,k]/dID - ID[i,j,k]/x_I[i]
deriv(ID[2:na, 1:nh,1:2]) <- effective_FOI[i,j,k]/(effective_FOI[i,j,k]*uD + 1) - ID[i,j,k]/dID - (ID[i,j,k]-ID[i-1,j,k])/x_I[i]

# p_det - probability of detection by microscopy, immunity decreases chances of
# infection because it pushes parasite density down
aD <- user()
fD0 <- user()
gammaD <- user()
d1 <- user()
ID0 <- user()
kD <- user()
dim(age) <- na
age[] <- user() # vector of age categories supplied by user

dim(fd) <- na
fd[1:(na-1)] <- 1-(1-fD0)/(1+((age[i]+age[i+1])/2/aD)^gammaD)
fd[na]<-1-(1-fD0)/(1+(age[i]/aD)^gammaD)
dim(p_det) <- c(na,nh,2)
p_det[,,] <- d1 + (1-d1)/(1 + fd[i]*(ID[i,j,k]/ID0)^kD)

# # Force of infection, depends on level of infection blocking immunity
# dim(FOI) <- c(na,nh)
# FOI[1:na, 1:nh] <- EIR[i,j] * (if(IB[i,j]==0) b0 else b[i,j])

# Force of infection, depends on level of infection blocking immunity
dim(FOI_lag) <- c(na,nh,2)
FOI_lag[1:na, 1:nh,1:2] <- EIR[i,j,k] * (if(IB[i,j,k]==0) b0 else b[i,j,k])

# Current FOI depends on humans that have been through the latent period
dE <- user() # latent period of human infection.
lag_rates <- user()

FOI_eq[,] <- user()
dim(FOI_eq) <- c(na,nh)
init_FOI[,,] <- FOI_eq[i,j]
dim(init_FOI) <- c(na,nh,lag_rates)
initial(FOI[1:na,1:nh,1:lag_rates,1]) <- init_FOI[i,j,k]
initial(FOI[1:na,1:nh,1:lag_rates,2]) <- init_FOI[i,j,k]
dim(FOI) <- c(na,nh,lag_rates,2)

deriv(FOI[1:na,1:nh,1,1]) <- (lag_rates/dE)*FOI_lag[i,j,1] - (lag_rates/dE)*FOI[i,j,1,1]
deriv(FOI[1:na,1:nh,1,2]) <- (lag_rates/dE)*FOI_lag[i,j,2] - (lag_rates/dE)*FOI[i,j,1,2]
deriv(FOI[1:na,1:nh,2:lag_rates,1:2]) <- (lag_rates/dE)*FOI[i,j,k-1,l] - (lag_rates/dE)*FOI[i,j,k,l]

dim(effective_FOI)<-c(na,nh,2)
effective_FOI[1:na,1:nh,1]<-FOI[i,j,lag_rates,1]
effective_FOI[1:na,1:nh,2]<-FOI[i,j,lag_rates,2]*(1-(1-SMC_prot)*on_vect[i])
# EIR -rate at which each age/het/int group is bitten
# rate for age group * rate for biting category * FOI for age group * prop of
DY<-user()
# init_EIR <- user()
# max_EIR <- user()
# EIR[,] <- exp(log_EIR) * rel_foi[j] * foi_age[i]
# output(EIR_out) <- exp(log_EIR)*DY
output(EIR_out) <- (av * Iv/omega)*DY

av0 <- user()
dim(EIR) <- c(na,nh,2)
EIR[,,] <- rel_foi[j] * foi_age[i] * Iv*av0/omega

#EIR_td<-interpolate(EIR_times, EIR_valsd, "constant")
# EIR_times[]<-user()
# EIR_vals[]<-user()
#EIR_valsd[]<-EIR_vals[i]/DY
#dim(EIR_times)<-user()
# dim(EIR_vals)<-user()
#dim(EIR_valsd)<-length(EIR_vals)
dim(foi_age) <- na
foi_age[] <- user()
dim(rel_foi) <- nh
rel_foi[] <- user()

#####Feed in EIR values; COMMENT OUT IF FEEDING IN BETAA VALS#########
  # DY <- user()
  # dim(EIR) <- c(na,nh)
  # EIR_times[]<-user()
  # EIR_seq[]<-user()
  # EIR_valsd[]<-EIR_seq[i]/DY
  # EIR_td<-interpolate(EIR_times, EIR_valsd, "constant")
  #
  # dim(EIR_times)<-user()
  # dim(EIR_seq)<-user()
  # dim(EIR_valsd)<-length(EIR_seq)
  #
  # EIR[,] <- EIR_td * rel_foi[j] * foi_age[i]
  # output(EIR) <- EIR
  # output(EIR_td) <- EIR_td
  # output(Ivout) <- Iv
  #
  # output(omega) <- omega
  #

  # ##------------------------------------------------------------------------------
  # #####################
  # ## MOSQUITO STATES ##
  # #####################
  # ##------------------------------------------------------------------------------

  # See supplementary materials S1 from http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324#s6
  # fitted entomological parameters:
  mv0 <- user() # initial mosquito density
  mu0 <- user() # baseline mosquito death rate
  tau1 <- user() # duration of host-seeking behaviour
  tau2 <- user() # duration of resting behaviour
  p10 <- user() # prob of surviving 1 feeding cycle
  p2 <- user() #prob of surviving one resting cycle
  fv <- 1/( tau1 + tau2 ) # mosquito feeding rate (zbar from intervention param.)
  mu <- -fv*log(p10*p2) # mosquito death rate
  omega <- user() #normalising constant for biting rates
  Q0 <- user() # proportion of anthropophagy
  av <- fv*Q0
  
  # cA is the infectiousness to mosquitoes of humans in the asmyptomatic compartment broken down
  # by age/het/int category, infectiousness depends on p_det which depends on detection immunity
  cU <- user() # infectiousness U -> mosq
  cD <- user() # infectiousness D -> mosq
  cT <- user() # T -> mosq
  gamma1 <- user() # fitted value of gamma1 characterises cA function
  dim(cA) <- c(na,nh,2)
  cA[,,] <- cU + (cD-cU)*p_det[i,j,k]^gamma1
  # Current hum->mos FOI depends on the number of individuals now producing gametocytes (12 day lag)
  delayGam <- user() # Lag from parasites to infectious gametocytes
  delayMos <- user() # Extrinsic incubation period.

  # Number of mosquitoes that become infected at each time point
  surv <- exp(-mu*delayMos)

  # Force of infection from humans to mosquitoes

  FOIv_eq <- user()
  initial(FOIv[]) <- FOIv_eq*delayGam/lag_rates
  dim(FOIv) <- lag_rates

  FOIvijk[1:na, 1:nh,1:2] <- av0 * (cT*T[i,j,k] + cD*D[i,j,k] + cA[i,j,k]*A[i,j,k] + cU*U[i,j,k]) * rel_foi[j] *foi_age[i]/omega
  dim(FOIvijk) <- c(na,nh,2)
  lag_FOIv <- sum(FOIvijk)
  output(lag_FOIv)<-lag_FOIv
  deriv(FOIv[1]) <- lag_FOIv - (lag_rates/delayGam)*FOIv[1]
  deriv(FOIv[2:lag_rates]) <- (lag_rates/delayGam)*FOIv[i-1] -
    (lag_rates/delayGam)*FOIv[i]

  ince <- FOIv[lag_rates] * lag_rates/delayGam * Sv

  initial(ince_delay[]) <- FOIv_eq*init_Sv*mv0*delayMos/lag_rates
  dim(ince_delay) <- lag_rates

  deriv(ince_delay[1]) <- ince - (lag_rates/delayMos)*ince_delay[1]
  deriv(ince_delay[2:lag_rates]) <- (lag_rates/delayMos)*ince_delay[i-1] -
    (lag_rates/delayMos)*ince_delay[i]

  incv <- ince_delay[lag_rates]*lag_rates/delayMos *surv

  # feed in betaa values from a random walk
####Comment in if reading in betaa values#####
  betaa_times[]<-user()
  dim(betaa_times)<-user()
  betaa_vals[]<-user()
  dim(betaa_vals)<-user()
  # Interpolate constant betaa values between user-defined switch points
  betaa_td <- interpolate(betaa_times, betaa_vals, "constant")
  output(betaa_out) <- betaa_td
  dim(SMC_times)<-user()
  dim(SMC_vals)<-user()
  SMC_times[]<-user()
  SMC_vals[]<-user()
  SMC_prot<-interpolate(SMC_times, SMC_vals, "constant")
#####Comment next line out if reading in betaa values#####
   # betaa_td <- mv * mu

  # Sv - Susceptible mosquitoes
  # Ev - latently infected (exposed) mosquitoes. Number of compartments used to simulate delay in becoming infectious
  # Iv - Infectious mosquitoes

  # initial state values:
  init_Sv <- user()
  init_Ev <- user()
  init_Iv <- user()
  initial(Sv) <- init_Sv * mv0
  initial(Ev) <- init_Ev * mv0
  initial(Iv) <- init_Iv * mv0

  deriv(Sv) <- -ince - mu*Sv + betaa_td
  deriv(Ev) <- ince - incv - mu*Ev
  deriv(Iv) <- incv - mu*Iv

  # Total mosquito population
  mv = Sv+Ev+Iv

  ##------------------------------------------------------------------------------
  ###################
  ## MODEL OUTPUTS ##
  ###################
  ##------------------------------------------------------------------------------
  output(Sout) <- sum(S[,,])
  output(Tout) <- sum(T[,,])
  output(Dout) <- sum(D[,,])
  output(Aout) <- sum(A[,,])
  output(Uout) <- sum(U[,,])
  output(Pout) <- sum(P[,,])

  # Outputs for clinical incidence and prevalence on a given day
  # population densities for each age category
  den[] <- user()
  dim(den) <- na
  output(den)<-den
  # index of the age vector above 59 months
  age59 <- user(integer=TRUE)
  age05 <- user(integer=TRUE)
  # slide positivity in 0 -5 year age bracket
  dim(prev0to59) <- c(age59,nh,2)
  prev0to59[1:age59,,] <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]
  output(prev) <- sum(prev0to59[,,])/sum(den[1:age59])
  output(age59)<-age59
  
  dim(prevall) <- c(na,nh,2)
  prevall[,,] <- T[i,j,k] + D[i,j,k]  + A[i,j,k]*p_det[i,j,k]
  output(prev_all) <- sum(prevall[,,])/sum(den[])
  
  # clinical incidence
  dim(clin_inc0tounder5) <- c(age59,nh,2)
  clin_inc0tounder5[1:age59,,] <- clin_inc[i,j,k]
  output(incunder5) <- sum(clin_inc0tounder5)/sum(den[1:age59])
  
  dim(clin_inc0to5) <- c(age05,nh,2)
  clin_inc0to5[1:age05,,] <- clin_inc[i,j,k]
  output(inc05) <- sum(clin_inc0to5)/sum(den[1:age05])
  output(cov_vect)<-cov_vect
  output(inc) <- sum(clin_inc[,,])
  output(inc1) <- sum(clin_inc[1,,])/den[1]
  output(inc_smc) <- sum(clin_inc[2:age05,,])/sum(den[2:age05])
  # Param checking outputs
  # output(mu) <- mu
  # output(beta_larval) <- beta_larval
  # output(KL) <- KL
   output(mv) <- mv
  # output(K0) <- K0
  # output(agestart) <- agestart
  # output(ageend) <- ageend
