#Data prep for pmcmc
library(excel.link)
library(zoo)
library(plyr)
library(dplyr)
library(ggpubr)
library(binom)
library(gridExtra)
library(ggplot2)
library(viridis)
library(lubridate)


##Primigrav data sets
##Nigeria
NG_anc <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/NG_ANC_mother_grouped_sitegrav.rds')

NG_pg_asa <- NG_anc[NG_anc$site=='Asa'&NG_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_pg_asa,'data_raw_NG_pg_asa.RDS')

NG_pg_ejigbo <- NG_anc[NG_anc$site=='Ejigbo'&NG_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_pg_ejigbo,'data_raw_NG_pg_ejigbo.RDS')

NG_pg_ifenorth <- NG_anc[NG_anc$site=='Ife North'&NG_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_pg_ifenorth,'data_raw_NG_pg_ifenorth.RDS')

NG_pg_moro <- NG_anc[NG_anc$site=='Moro'&NG_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_pg_moro,'data_raw_NG_pg_moro.RDS')

##Burkina Faso
BF_anc <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/BF_ANC_mother_grouped_sitegrav.rds')

BF_pg_banfora <- BF_anc[BF_anc$site=='Banfora'&BF_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_pg_banfora,'data_raw_BF_pg_banfora.RDS')


BF_pg_gaoua <- BF_anc[BF_anc$site=='Gaoua'&BF_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_pg_gaoua,'data_raw_BF_pg_gaoua.RDS')

BF_pg_orodara <- BF_anc[BF_anc$site=='Orodara'&BF_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_pg_orodara,'data_raw_BF_pg_orodara.RDS')

##Mozambique
MZ_anc <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/MZ_ANC_mother_grouped_sitegrav.rds')

MZ_pg_changara <- MZ_anc[MZ_anc$site=='Changara'&MZ_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_pg_changara,'data_raw_MZ_pg_changara.RDS')

MZ_pg_chemba <- MZ_anc[MZ_anc$site=='Chemba'&MZ_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_pg_chemba,'data_raw_MZ_pg_chemba.RDS')

MZ_pg_guro <- MZ_anc[MZ_anc$site=='Guro'&MZ_anc$grav_cat=='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_pg_guro,'data_raw_MZ_pg_guro.RDS')

##Multigrav
NG_mg_asa <- NG_anc[NG_anc$site=='Asa'&NG_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_mg_asa,'data_raw_NG_mg_asa.RDS')

NG_mg_ejigbo <- NG_anc[NG_anc$site=='Ejigbo'&NG_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_mg_ejigbo,'data_raw_NG_mg_ejigbo.RDS')

NG_mg_ifenorth <- NG_anc[NG_anc$site=='Ife North'&NG_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_mg_ifenorth,'data_raw_NG_mg_ifenorth.RDS')

NG_mg_moro <- NG_anc[NG_anc$site=='Moro'&NG_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(NG_mg_moro,'data_raw_NG_mg_moro.RDS')

##Burkina Faso
BF_anc <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/BF_ANC_mother_grouped_sitegrav.rds')

BF_mg_banfora <- BF_anc[BF_anc$site=='Banfora'&BF_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_mg_banfora,'data_raw_BF_mg_banfora.RDS')

BF_mg_gaoua <- BF_anc[BF_anc$site=='Gaoua'&BF_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_mg_gaoua,'data_raw_BF_mg_gaoua.RDS')

BF_mg_orodara <- BF_anc[BF_anc$site=='Orodara'&BF_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(BF_mg_orodara,'data_raw_BF_mg_orodara.RDS')

##Mozambique
MZ_anc <- readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/MZ_ANC_mother_grouped_sitegrav.rds')

MZ_mg_changara <- MZ_anc[MZ_anc$site=='Changara'&MZ_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_mg_changara,'data_raw_MZ_mg_changara.RDS')

MZ_mg_chemba <- MZ_anc[MZ_anc$site=='Chemba'&MZ_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_mg_chemba,'data_raw_MZ_mg_chemba.RDS')

MZ_mg_guro <- MZ_anc[MZ_anc$site=='Guro'&MZ_anc$grav_cat!='Gravidities 1',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(total),
            positive=sum(positive))
saveRDS(MZ_mg_guro,'data_raw_MZ_mg_guro.RDS')
