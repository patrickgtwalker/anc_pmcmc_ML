library(excel.link)
library(zoo)
library(plyr)
library(dplyr)
library(ggpubr)
library(binom)
library(gridExtra)
library(ggplot2)
library(viridis)

tanz_data_all_2017 <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Tanzania/TZ_ANC_data_region_2014_2017.csv')
data<-tanz_data_all_2017
tanz_data_process <- function(data, level = c('District','Region'), remove_before=NULL){
  if(is.null(remove_before)){
    remove_before <- min(data$yearmon)
  }
  level_call <- ifelse(level=='District','Council','Region')
  
  all <- data.frame(month = as.yearmon(data$yearmon),
                    tested = data$tested,
                    positive = data$positive,
                    site = data[,level_call])
  lt20 <- data.frame(month = as.yearmon(data$yearmon),
                       tested = data$test_LT20,
                       positive = data$pos_LT20,
                     site = data[,level_call])
  ge20 <- data.frame(month = as.yearmon(data$yearmon),
                     tested = data$test_GE20,
                     positive = data$pos_GE20,
                     site = data[,level_call])
  
  split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])
  
  all_list <- split_tibble(all[all$month>=remove_before,],'site')
  lt20_list <- split_tibble(lt20[lt20$month>=remove_before,],'site')
  ge20_list <- split_tibble(ge20[ge20$month>=remove_before,],'site')
  
  return(list(all = all_list, lt20 = lt20_list, ge20 = ge20_list))
}

tanz_data_list <- tanz_data_process(tanz_data_all_2017,'Region')
tanz_data_list_15to17 <- tanz_data_process(tanz_data_all_2017,remove_before = as.yearmon('Jan 2015'),level='Region')
rm(tanz_data_list_14to17)
names(tanz_data_list$lt20[1])
tanz_rainfall <- read_csv("./tanz/processed_inputs/ANC_data_Rainfall_ad1.csv")%>%
  mutate(yearmon=as.yearmon(yearmon))
##Lake Malawi Districts
lake_malawi_raw <- read.csv('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/Tanzania/ANC_rainfall_Lake_Malawi_dists.csv')
tanz_lakemal_list_15to17 <- tanz_data_process(lake_malawi_raw,level = 'District', remove_before = as.yearmon('Jan 2015'))
