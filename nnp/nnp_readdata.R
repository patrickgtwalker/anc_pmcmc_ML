library (bayesplot)
library(excel.link)
library(zoo)
library(plyr)
library(dplyr)
library(ggpubr)
library(binom)
library(gridExtra)
library(ggplot2)
library(viridis)

source('shared/addCIs.R')
#######################
########NIGERIA########
#######################
#Read in Nigeria data and rename most used variables
NG_ANC_mother <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Nigeria/ANC-based surveillance/data_anc_mother_nigeria.xlsx')%>%
  dplyr::rename(primigrav = q1_preg,
         prev_pregs = q2_preg,
         mal_symp = q1_mal,
         rdt = q2_mal,
         month = date) %>%
  dplyr::mutate(month = as.yearmon(month),
         grav = prev_pregs+1,
         ward = toupper(ward))
#Read in Nigeria data and rename most used variables
NG_CS_child_2020 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2020.xlsx')%>%
  dplyr::rename(rdt = q85c_result,
         mal_symp = q86a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2020 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2020.xlsx')%>%
  rename(individual_id = hh_id,
         date = start) %>%
  mutate(date = as.Date(date))

#Merge child to hh to get date of survey
NG_CS_2020 <- merge(NG_CS_child_2020,NG_CS_hh_2020,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))

##Read in 2021 data
NG_CS_child_2021 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2021.xlsx')%>%
  rename(rdt = q90c_result_6to59months,
         mal_symp = q91_a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2021 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2021.xlsx')%>%
  rename(individual_id = hh_id,
         date = start) %>%
  mutate(date = as.Date(date))

#Merge child to hh to get date of survey
NG_CS_2021 <- merge(NG_CS_child_2021,NG_CS_hh_2021,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))

#Combine 2 years of data
NG_CS_all <- plyr::rbind.fill(NG_CS_2020,NG_CS_2021)
NG_CS_all$district_code <- as.integer(substr(NG_CS_all$cluster_n,1,1))
NG_CS_all$district_n <- unlist(lapply(NG_CS_all$district_code, FUN = function(x) switch(x,'Ejigbo','Ife North','Asa','Moro')))

#By site and gravidity#
NG_CS_all_grouped_site <- NG_CS_all %>%
  rename(site=district_n)%>%
  mutate(month=as.yearmon(date),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(site))%>%
  group_by(site,month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_CS_all_grouped_site <- addCIs(NG_CS_all_grouped_site,NG_CS_all_grouped_site$positive,NG_CS_all_grouped_site$total)
saveRDS(NG_CS_all_grouped_site,'NG_CS_all_grouped_site.rds')
#Group by grav
NG_ANC_mother_grouped_sitegrav <- NG_ANC_mother %>%
  dplyr::rename(site = lga) %>%
  dplyr::mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")))%>%
  filter(!is.na(site))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_ANC_mother_grouped_sitegrav <- addCIs(NG_ANC_mother_grouped_sitegrav,NG_ANC_mother_grouped_sitegrav$positive,NG_ANC_mother_grouped_sitegrav$total)
saveRDS(NG_ANC_mother_grouped_sitegrav,'NG_ANC_mother_grouped_sitegrav.rds')

#######################
######MOZAMBIQUE#######
#######################
#Read in Mozambique data and rename most used variables
MZ_ANC_mother <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Mozambique/ANC-based surveillance/ANC Surveillance_11.2011.xlsx',
                              xl.sheet = 'ANC')%>%
  dplyr::rename(primigrav = q4_primagravidae,
                prev_pregs = q5_num_pregnancy,
                mal_symp = q6_mal_symptoms,
                rdt = q7_rdt_result,
                age = q2_age_n) %>%
  dplyr::mutate(month = as.yearmon(as.Date(date_interview_n)),
                grav = ifelse(primigrav=='Yes',1,ifelse(prev_pregs==0,NA,prev_pregs+1)))
#Read in Mozambique data and rename most used variables
MZ_CS_rdt_base <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Mozambique/Cross-sectional survey/Moz Baseline CSS Datasets 2020.xlsx',
                               xl.sheet = 'RDT')%>%
  mutate(date = as.Date(date))
MZ_CS_rdt_mid <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Mozambique/Cross-sectional survey/Moz Midline CSS Datasets 2021.xlsx',
                              xl.sheet = 'RDT')%>%
  mutate(date = as.Date(date,"%m/%d/%Y"))
MZ_CS_rdt <- rbind(MZ_CS_rdt_base,MZ_CS_rdt_mid)

#Group CS data by site and calculate prevalence and CI by month
MZ_CS_all_grouped_site <- MZ_CS_rdt %>%
  rename(site = district) %>%
  mutate(month=as.yearmon(date),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  filter(site %in% c('Changara','Chemba','Guro'))%>%
  group_by(site,month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_CS_all_grouped_site <- addCIs(MZ_CS_all_grouped_site,MZ_CS_all_grouped_site$positive,MZ_CS_all_grouped_site$total)
saveRDS(MZ_CS_all_grouped_site,'MZ_CS_all_grouped_site_160822.rds')

#Group ANC data by grav and site
MZ_ANC_mother_grouped_sitegrav <- MZ_ANC_mother %>%
  rename(site = district_n) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA))))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age>=12&age<50) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_ANC_mother_grouped_sitegrav <- addCIs(MZ_ANC_mother_grouped_sitegrav,MZ_ANC_mother_grouped_sitegrav$positive,MZ_ANC_mother_grouped_sitegrav$total)
saveRDS(MZ_ANC_mother_grouped_sitegrav,'MZ_ANC_mother_grouped_sitegrav_0822.rds')

#########################
######Burkina Faso#######
#########################
#Read in Burkina Faso data and rename most used variables
BF_ANC_mother <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Burkina Faso/ANC-based surveillance/ANC_BF_5.16.2022.xlsx',
                              xl.sheet = 'ANC_woman')%>%
  dplyr::rename(primigrav = Q04_first_preg,
                prev_pregs = Q05_nbr_prior_preg,
                mal_symp = Q06_mal_symp,
                rdt = RDT,
                age = Age) %>%
  dplyr::mutate(month = as.yearmon(as.Date(Date)),
                grav = ifelse(primigrav=='Oui',1,prev_pregs+1))
#Read in Burkina Faso data and rename most used variables
BF_CS_child_2019 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Burkina Faso/Cross-sectional survey/NNET_CSS_2019.xls',
                                 xl.sheet = 'Child_data')%>%
  dplyr::rename(rdt = Result_RDT) %>%
  dplyr::mutate(month = as.yearmon(as.Date('2019-07-01')))

BF_CS_child_2020 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Burkina Faso/Cross-sectional survey/NNET_CSS_2020.xls',
                                 xl.sheet = 'Child_data')%>%
  dplyr::rename(rdt = Result_RDT) %>%
  dplyr::mutate(month = as.yearmon(as.Date('2020-06-01')))

BF_CS_child_2021 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Burkina Faso/Cross-sectional survey/NNET_CSS_2021.xlsx',
                                 xl.sheet = 'Child_data')%>%
  dplyr::rename(rdt = Result_RDT) %>%
  dplyr::mutate(month = as.yearmon(as.Date('2021-06-01')))

BF_CS_child_all <- rbind.fill(BF_CS_child_2019,BF_CS_child_2020,BF_CS_child_2021)

#Group CS data and calculate prevalence and CI by month
BF_CS_all_grouped_site <- BF_CS_child_all %>%
  rename(site = District) %>%
  mutate(
    rdt=as.numeric(ifelse(rdt=='Négatif',0,1))
  )%>%
  filter(!is.na(rdt)) %>%
  group_by(site,month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_CS_all_grouped_site <- addCIs(BF_CS_all_grouped_site,BF_CS_all_grouped_site$positive,BF_CS_all_grouped_site$total)
saveRDS(BF_CS_all_grouped_site,'BF_CS_all_grouped_site_0822.rds')

#Group ANC data by site
#Group by grav
BF_ANC_mother_grouped_sitegrav <- BF_ANC_mother %>%
  rename(site = District) %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         rdt=as.numeric(ifelse(rdt=='Positif',1,0)))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age!=88) %>%
  filter(!is.na(grav))%>%
  group_by(site,month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_ANC_mother_grouped_sitegrav <- addCIs(BF_ANC_mother_grouped_sitegrav,BF_ANC_mother_grouped_sitegrav$positive,BF_ANC_mother_grouped_sitegrav$total)
saveRDS(BF_ANC_mother_grouped_sitegrav,'BF_ANC_mother_grouped_sitegrav.rds')

########################################
####Create Plots for Each Country#######
########################################
windows(10,8)
#NIGERIA#
#Create data tables to annotate figure
NG_iv <- data.frame(site = c('Ejigbo','Asa','Moro','Ife North'),
                    month = as.yearmon(rep(as.Date('2020-11-01'),4)))
NG_labs <- c('Ejigbo - Standard','Asa - IG2','Moro - RG','Ife North - PBO')
names(NG_labs) <- c('Ejigbo','Asa','Moro','Ife North')
#Grav colored, facetted by site
ggplot(data=NG_ANC_mother_grouped_sitegrav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  annotate("rect", xmin = as.yearmon('Oct 2020'), xmax = as.yearmon('Nov 2020'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.yearmon('Jul 2021'), xmax = as.yearmon('Nov 2021'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(vars(site),nrow=2,ncol=2,labeller = labeller(site = NG_labs))+
  geom_point(data=NG_CS_all_grouped_site,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=NG_CS_all_grouped_site,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")+
  geom_vline(data=NG_iv,aes(xintercept=month),color='darkgrey',size=2,alpha = 0.5)+
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12, face = "bold"))


#MOZAMBIQUE#
#Create data tables to annotate figure
MZ_iv <- data.frame(site = c('Chemba','Guro','Changara'),
                    month = as.yearmon(rep(as.Date('2020-11-01'),3)))
MZ_labs <- c('Chemba - Standard','Guro - IG2','Changara - PBO')
names(MZ_labs) <- c('Chemba','Guro','Changara')

#Grav colored
ggplot(data=MZ_ANC_mother_grouped_sitegrav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  annotate("rect", xmin = as.yearmon('Jan 2021'), xmax = as.yearmon('Jun 2021'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  geom_point(data=MZ_CS_all_grouped_site,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=MZ_CS_all_grouped_site,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")+
  facet_wrap(vars(site),nrow=2,ncol=2,labeller = labeller(site = MZ_labs))+
  geom_vline(data=MZ_iv,aes(xintercept=month),color='darkgrey',size=2,alpha = 0.5)+
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12, face = "bold"))

#BURKINA FASO#
#Create data tables to annotate figure
BF_iv <- data.frame(site = c('Gaoua','Banfora','Orodara'),
                    month = as.yearmon(c(as.Date('2019-08-01'),as.Date('2019-10-01'),as.Date('2019-06-01'))))
BF_labs <- c('Gaoua - Standard','Banfora - IG2','Orodara - PBO')
names(BF_labs) <- c('Gaoua','Banfora','Orodara')

ggplot(data=BF_ANC_mother_grouped_sitegrav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  annotate("rect", xmin = as.yearmon('Jun 2019'), xmax = as.yearmon('Oct 2019'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.yearmon('Jun 2020'), xmax = as.yearmon('Oct 2020'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  annotate("rect", xmin = as.yearmon('Jun 2021'), xmax = as.yearmon('Oct 2021'), ymin = 0, ymax = 1,alpha = .1,fill = "#0072B2")+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  geom_point(data=BF_CS_all_grouped_site,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=BF_CS_all_grouped_site,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")+
  facet_wrap(vars(site),nrow=2,ncol=2,labeller = labeller(site = BF_labs))+
  geom_vline(data=BF_iv,aes(xintercept=month),color='darkgrey',size=2,alpha = 0.5)+
  theme(legend.position="bottom",
        axis.text.x = element_text(size = 10, angle = 45, vjust = 1.2, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 12, face = "bold"))

#######################################
######Process data for pMCMC###########
#######################################
##Primigrav data sets
##Nigeria
NG_anc <- NG_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/NG_ANC_mother_grouped_sitegrav.rds')

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
BF_anc <- BF_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/BF_ANC_mother_grouped_sitegrav.rds')

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
MZ_anc <- MZ_ANC_mother_grouped_sitegrav
#readRDS('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Analysis/nnp_explor/MZ_ANC_mother_grouped_sitegrav.rds')

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

