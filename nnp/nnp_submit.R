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

##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("run_pmcmc.R","run_pmcmc.R",
             "model_parameters.R","equilibrium-init-create-stripped.R")


ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

config_16 <- didehpc::didehpc_config(cores = 16, parallel = TRUE)
obj_16 <- didehpc::queue_didehpc(ctx,config = config_16)
obj_16$login()

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
nnp_bulk <- obj_16$enqueue_bulk(1:10, function(i,data_site){
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

nnp_bulk_allage <- obj_16$task_bundle_get('deistic_grison')
nnp_bulk$status() #'deistic_grison'
nnp_bulk$tasks$`4d6c8776a4d59c636bca57982d462330`$log()
test <- nnp_bulk$tasks$`734c6d65b473efd1f031f1268d65c1c6`$result()
task_ids <- nnp_bulk$ids
bf_banfora <-  nnp_bulk_allage$tasks$`734c6d65b473efd1f031f1268d65c1c6`$result()
bf_gaoua <- nnp_bulk_allage$tasks$`5a97412fcf9e9787ac3b6bdff17b962a`$result()
bf_orodara <- nnp_bulk_allage$tasks$c73ef29830068b971dfcd6ab11048438$result()
mz_changara <- nnp_bulk_allage$tasks$`4d6c8776a4d59c636bca57982d462330`$result()
mz_chemba <- nnp_bulk_allage$tasks$fb6717770a8f503846edbdd5edb0ec3c$result()
mz_guro <- nnp_bulk_allage$tasks$`8a16be32b0444c24847021bef114c851`$result()
ng_asa <- nnp_bulk_allage$tasks$`91f1846e274d751ea13143febeec2e05`$result()
ng_ejigbo <- nnp_bulk_allage$tasks$`56c25737682075594a81364213b60c85`$result()
ng_ifenorth <- nnp_bulk_allage$tasks$fc9daf98af27d12cef1180f327324dc0$result()
ng_moro <- nnp_bulk_allage$tasks$`98d1902129744dbb782b7f4340b3cafc`$result()
nnp_result_list_aa <- list(bf_banfora ,
                        bf_gaoua ,
                        bf_orodara ,
                        mz_changara ,
                        mz_chemba ,
                        mz_guro ,
                        ng_asa ,
                        ng_ejigbo ,
                        ng_ifenorth ,
                        ng_moro)

##Fitting as <5 yo prevalence, bundle='blessed_greathornedowl'
nnp_bulk$status() 
task_ids <- nnp_bulk$ids
id <- task_ids[1]
nnp_result_list <- lapply(1:10, function(id){
  nnp_bulk$tasks[[id]]$result()
})

##Fit to primigrav data to <5 yo
##Set up cluster##
root <- "T:/jth/contexts"
sources <- c("run_pmcmc.R",
             "model_parameters.R","equilibrium-init-create-stripped.R")


ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

##Set up and run scenarios on 64 node core
config_16 <- didehpc::didehpc_config(cores = 16, parallel = TRUE)
obj_16 <- didehpc::queue_didehpc(ctx,config = config_16)
obj_16$login()

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
nnp_pg_bulk$name #"granitic_seriema"
nnp_pg_result_list <- lapply(1:10, function(id){
  nnp_pg_bulk$tasks[[id]]$result()
})

##Test
run_pmcmc(data = data_raw_bf_pg_banfora,
          n_particles = 200,
          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
          max_EIR=1000,
          max_steps = 1e7,
          atol = 1e-5,
          rtol = 1e-6,
          n_steps = 1000,
          n_threads = 16)
source('run_pmcmc_pg.R')
run_pmcmc_pg(data = data_raw_bf_pg_banfora,
          n_particles = 200,
          proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
          max_EIR=1000,
          max_steps = 1e7,
          atol = 1e-5,
          rtol = 1e-6,
          n_steps = 1000,
          n_threads = 16)
pg_corr_sample <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_corr/pg_corr_sample.RDS')
saveRDS(pg_corr_sample,'pg_corr_sample.RDS')
root <- "T:/jth/contexts"
sources <- c("run_pmcmc_pg.R","run_pmcmc.R",
             "model_parameters.R","equilibrium-init-create-stripped.R")


ctx <- context::context_save("T:/jth/contexts", sources = sources,
                             packages = c('statmod','coda','dplyr'),
                             package_sources = conan::conan_sources(c("mrc-ide/odin.dust",'mrc-ide/mcstate')))

##Set up and run scenarios on 64 node core
config_16 <- didehpc::didehpc_config(cores = 16, parallel = TRUE)
obj_16 <- didehpc::queue_didehpc(ctx,config = config_16)

obj_16$cluster_load(TRUE)
obj_16$login()
test_pg <- obj_16$enqueue(run_pmcmc_pg(data = data_raw_bf_pg_banfora,
                                       n_particles = 200,
                                       proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
                                       max_EIR=1000,
                                       max_steps = 1e7,
                                       atol = 1e-5,
                                       rtol = 1e-6,
                                       n_steps = 1000,
                                       n_threads = 16))
test_pg$status()

test_pg$log()
obj_16$unsubmit(test_pg$id)

nnp_pgcorr_bulk <- obj_16$enqueue_bulk(1:10, function(i,data_site){
  run_pmcmc_pg(data = data_site[[i]],
            n_particles = 200,
            proposal_matrix = matrix(c(0.0336,-0.000589,-0.000589,0.049420),nrow=2),
            max_EIR=1000,
            max_steps = 1e7,
            atol = 1e-5,
            rtol = 1e-6,
            n_steps = 1000,
            n_threads = 16)
},data_site=nnp_pg_list)
nnp_pgcorr_bulk$status() #'dendrological_yorkshireterrier'
nnp_pgcorr_result_list <- lapply(1:10, function(id){
  nnp_pgcorr_bulk$tasks[[id]]$result()
})

##Run pg with relaxed prior on volatility
nnp_pg_list <- list(data_raw_bf_pg_banfora,data_raw_bf_pg_gaoua,data_raw_bf_pg_orodara,
                    data_raw_mz_pg_changara,data_raw_mz_pg_chemba,data_raw_mz_pg_guro,
                    data_raw_ng_pg_asa,data_raw_ng_pg_ejigbo,data_raw_ng_pg_ifenorth,data_raw_ng_pg_moro)
nnp_pgvol_bulk <- obj_16$enqueue_bulk(1:10, function(i,data_site){
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
config_16
nnp_pgvol_bulk$status()#'untiring_zenaida'
nnp_pgvol_bulk$tasks$cfd5748a9c445fa9bf3a8c1e8b6b765c$log()
nnp_pgvol_result_list <- lapply(1:10, function(id){
  nnp_pgvol_bulk$tasks[[id]]$result()
})
