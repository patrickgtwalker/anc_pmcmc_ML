library(excel.link)
library(zoo)
library(tidyverse)
library(ggpubr)
library(binom)
library(gridExtra)
library(ggplot2)
library(viridis)

##Read in MipMon data
mipmon_raw <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Mozambique/isglobal_cism_data/mipmon_merged.csv')

##Format data
mipmon.all <- mipmon_raw %>%
  mutate(date = as.Date(visdate,tryFormats='%d/%m/%Y'),
         year = format(date, format='%Y'),
         week = week(date),
         month = as.yearmon(date),
         N=1,
         rdt=ifelse(pcrpos=='PCR-',
                    0,
                    ifelse(pcrpos=='PCR+',
                           ifelse(density<100&!is.na(density),
                                  0,
                                  ifelse(!is.na(density),
                                         1,
                                         NA)
                           ),
                           NA)
         ),
         site=recode(posto_code,
                     `Ilha-Josina`='Ilha Josina',
                     `Magude-sede`='Magude Sede',
                     `MOtzae`='Magude Sede',
                     `Manhica-Sede`='Manhica',
                     .default = NA_character_),
         grav=ifelse(gestnum==1,'primi',ifelse(!is.na(gestnum),'multi',NA)))%>%
  filter(visit=='PN'&!is.na(site)&!is.na(gestnum)&!is.na(month)&!is.na(rdt))

data_raw_mipmon_pg_ij <- mipmon.all[mipmon.all$site=='Ilha Josina'&mipmon.all$grav=='primi',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(N),
            positive=sum(rdt,na.rm = TRUE))
saveRDS(data_raw_mipmon_pg_ij,'mipmon/data/data_raw_mipmon_pg_ij.RDS')
data_raw_mipmon_pg_ij <- readRDS('mipmon/data/data_raw_mipmon_pg_ij.RDS')

data_raw_mipmon_mg_ij <- mipmon.all[mipmon.all$site=='Ilha Josina'&mipmon.all$grav=='multi',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(N),
            positive=sum(rdt,na.rm = TRUE))
saveRDS(data_raw_mipmon_mg_ij,'mipmon/data/data_raw_mipmon_mg_ij.RDS')
data_raw_mipmon_mg_ij <- readRDS('mipmon/data/data_raw_mipmon_mg_ij.RDS')

data_raw_mipmon_pg_ms <- mipmon.all[mipmon.all$site=='Magude Sede'&mipmon.all$grav=='primi',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(N),
            positive=sum(rdt,na.rm = TRUE))
saveRDS(data_raw_mipmon_pg_ms,'mipmon/data/data_raw_mipmon_pg_ms.RDS')
data_raw_mipmon_pg_ms <- readRDS('mipmon/data/data_raw_mipmon_pg_ms.RDS')

data_raw_mipmon_mg_ms <- mipmon.all[mipmon.all$site=='Magude Sede'&mipmon.all$grav=='multi',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(N),
            positive=sum(rdt,na.rm = TRUE))
saveRDS(data_raw_mipmon_mg_ms,'mipmon/data/data_raw_mipmon_mg_ms.RDS')

data_raw_mipmon_pg_man <- mipmon.all[mipmon.all$site=='Manhica'&mipmon.all$grav=='primi',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(N),
            positive=sum(rdt,na.rm = TRUE))
saveRDS(data_raw_mipmon_pg_man,'mipmon/data/data_raw_mipmon_pg_man.RDS')

data_raw_mipmon_mg_man <- mipmon.all[mipmon.all$site=='Manhica'&mipmon.all$grav=='multi',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(N),
            positive=sum(rdt,na.rm = TRUE))
saveRDS(data_raw_mipmon_mg_man,'mipmon/data/data_raw_mipmon_mg_man.RDS')

data_raw_mipmon_pg_all <- mipmon.all[mipmon.all$grav=='primi',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(N),
            positive=sum(rdt,na.rm = TRUE))
saveRDS(data_raw_mipmon_pg_all,'mipmon/data/data_raw_mipmon_pg_all.RDS')
data_raw_mipmon_mg_all <- mipmon.all[mipmon.all$grav=='multi',] %>%
  group_by(month)%>%
  summarise(t=cur_group_id()*30,
            tested=sum(N),
            positive=sum(rdt,na.rm = TRUE))
saveRDS(data_raw_mipmon_mg_all,'mipmon/data/data_raw_mipmon_mg_all.RDS')

mipmon_data_pg <- list(`Ilha Josina` = data_raw_mipmon_pg_ij,
                       `Magude Sede` = data_raw_mipmon_pg_ms,
                       `Manhica` = data_raw_mipmon_pg_man,
                       `All Sites` = data_raw_mipmon_pg_all)
mipmon_data_mg <- list(`Ilha Josina` = data_raw_mipmon_mg_ij,
                       `Magude Sede` = data_raw_mipmon_mg_ms,
                       `Manhica` = data_raw_mipmon_mg_man,
                       `All Sites` = data_raw_mipmon_mg_all)
