library("odin.dust")
library("odin")
library("patchwork")
library('mcstate')
library(didehpc)
library(pkgdepends)
library(dplyr)
library("coda")
library(binom)
library(ggplot2)
library(bayesplot)

#Required functions
source('nnp/run_pmcmc.R')
source('nnp/run_pmcmc_pg.R')
source('shared/plot_particle_filter.R')
source('shared/addCIs.R')
source('shared/model_parameters.R')
source('shared/equilibrium-init-create-stripped.R')

##Import NNP data (incomplete data sets)
##All
data_raw_ng_asa <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_NG_asa.RDS')
data_raw_ng_ifenorth <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_NG_ifenorth.RDS')
data_raw_ng_ejigbo <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_NG_ejigbo.RDS')
data_raw_ng_moro <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_NG_moro.RDS')

data_raw_bf_banfora <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_BF_banfora.RDS')
data_raw_bf_orodara <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_BF_orodara.RDS')
data_raw_bf_gaoua <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_BF_gaoua.RDS')

data_raw_mz_guro <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_MZ_guro.RDS')
data_raw_mz_chemba <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_MZ_chemba.RDS')
data_raw_mz_changara <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_MZ_changara.RDS')

nnp_list <- list(data_raw_bf_banfora,data_raw_bf_gaoua,data_raw_bf_orodara,
                 data_raw_mz_changara,data_raw_mz_chemba,data_raw_mz_guro,
                 data_raw_ng_asa,data_raw_ng_ejigbo,data_raw_ng_ifenorth,data_raw_ng_moro)

##Primigrav
data_raw_ng_pg_asa <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_ng_pg_asa.RDS')
data_raw_ng_pg_ifenorth <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_ng_pg_ifenorth.RDS')
data_raw_ng_pg_ejigbo <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_ng_pg_ejigbo.RDS')
data_raw_ng_pg_moro <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_ng_pg_moro.RDS')

data_raw_bf_pg_banfora <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_bf_pg_banfora.RDS')
data_raw_bf_pg_orodara <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_bf_pg_orodara.RDS')
data_raw_bf_pg_gaoua <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_bf_pg_gaoua.RDS')

data_raw_mz_pg_guro <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_mz_pg_guro.RDS')
data_raw_mz_pg_chemba <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_mz_pg_chemba.RDS')
data_raw_mz_pg_changara <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/data_raw_mz_pg_changara.RDS')

nnp_pg_list <- list(data_raw_bf_pg_banfora,data_raw_bf_pg_gaoua,data_raw_bf_pg_orodara,
                    data_raw_mz_pg_changara,data_raw_mz_pg_chemba,data_raw_mz_pg_guro,
                    data_raw_ng_pg_asa,data_raw_ng_pg_ejigbo,data_raw_ng_pg_ifenorth,data_raw_ng_pg_moro)

##Test run_pmcmc function##
test_run <- run_pmcmc(data = data_raw_bf_banfora,
                      n_particles = 10,
                      proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                      max_EIR=1000,
                      max_steps = 1e7,
                      atol = 1e-5,
                      rtol = 1e-6,
                      n_steps = 10,
                      n_threads = 2)
plot_particle_filter(test_run$history,true_history=data_raw_bf_banfora,times=data_raw_bf_banfora$t)


##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("nnp/run_pmcmc.R",
             "shared/model_parameters.R",
             "shared/equilibrium-init-create-stripped.R")

ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

config_16 <- didehpc::didehpc_config(cores = 16, parallel = TRUE)
obj_16 <- didehpc::queue_didehpc(ctx,config = config_16)
obj_16$login()

##Fitting all gravidity NNP as <5 yo prevalence
nnp_all_bulk <- obj_16$enqueue_bulk(1:10, function(i,data_site){
  run_pmcmc(data = data_site[[i]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 16)
  },data_site=nnp_list)

nnp_all_result_list <- lapply(1:10, function(id){
  nnp_all_bulk$tasks[[id]]$result()
})


##Fit to primigrav data to <5 yo
nnp_pg_bulk <- obj_16$enqueue_bulk(1:10, function(i,data_site){
  run_pmcmc(data = data_site[[i]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 16)
},data_site=nnp_pg_list)

nnp_pg_bulk$status()

nnp_pg_result_list <- lapply(1:10, function(id){
  nnp_pg_bulk$tasks[[id]]$result()
})

