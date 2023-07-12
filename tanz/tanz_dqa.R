library(foreign)
library(sas7bdat)
library(tidyverse)
library(readxl)
library(devtools)
devtools::install_github("mrc-ide/umbrella")
library(umbrella)
library(ggplot2)
library(lubridate)
library(zoo)
library(terra)
library(stringr)
library(geofacet)

list.files("./Data")
#### read in ANC data ##
library(excel.link)
Full_data_2016<-xl.read.file("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",xl.sheet="2016_2019")
Full_data_2020<-xl.read.file("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",xl.sheet="2020_2021")
Full_data_2022<-xl.read.file("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",xl.sheet="2022")



names(Full_data_2016)[5]<-"Month"
names(Full_data_2016)[6]<-"Before_12wk"
names(Full_data_2016)[7]<-"After_12wk"
names(Full_data_2016)[8]<-"ANC_re_ad"
names(Full_data_2016)[9]<-"Total_re_ad"
names(Full_data_2016)[10]<-"ANC_test"
names(Full_data_2016)[11]<-"ANC_pos"
names(Full_data_2016)[16]<-"Hb_test"
names(Full_data_2016)[17]<-"Anaemia"




names(Full_data_2020)[5]<-"Month"
names(Full_data_2020)[6]<-"Before_12wk"
names(Full_data_2020)[7]<-"After_12wk"
names(Full_data_2020)[8]<-"ANC_re_ad"
names(Full_data_2020)[9]<-"Total_re_ad"
names(Full_data_2020)[10]<-"ANC_test"
names(Full_data_2020)[11]<-"ANC_pos"
names(Full_data_2020)[16]<-"Hb_test"
names(Full_data_2020)[17]<-"Anaemia"





names(Full_data_2022)[5]<-"Month"
names(Full_data_2022)[6]<-"Year"
names(Full_data_2022)[7]<-"Before_12wk"
names(Full_data_2022)[8]<-"After_12wk"
names(Full_data_2022)[9]<-"ANC_re_ad"
names(Full_data_2022)[10]<-"Total_re_ad"
names(Full_data_2022)[11]<-"ANC_test"
names(Full_data_2022)[12]<-"ANC_pos"



names(Full_data_2022)[17]<-"Hb_test"
names(Full_data_2022)[18]<-"Anaemia"





Full_data<-rbind(
  Full_data_2022%>%select(Region,Council,HF,Year,Month,Before_12wk,After_12wk,ANC_re_ad,Total_re_ad,ANC_test,ANC_pos,Hb_test,Anaemia),
  Full_data_2020%>%select(Region,Council,HF,Year,Month,Before_12wk,After_12wk,ANC_re_ad,Total_re_ad,ANC_test,ANC_pos,Hb_test,Anaemia),
  Full_data_2016%>%select(Region,Council,HF,Year,Month,Before_12wk,After_12wk,ANC_re_ad,Total_re_ad,ANC_test,ANC_pos,Hb_test,Anaemia)
)
#### read in original ANC data ##
Full_data<-read.sas7bdat("./tanz/Patrick/Data/anc2014_17.sas7bdat")
Full_data$Region[Full_data$Region=="Songwe Region"]="Mbeya Region"

#### summarise at district level ##
tanz_data_all_2017_dist <-Full_data%>%
  mutate(yearmon=as.yearmon(paste0(Full_data$Year,"-",Full_data$Month),"%Y-%B"))%>%
  group_by(yearmon,Region,Council)%>%
  dplyr::summarise(count=n(),positive=sum(Mal),tested=sum(RDT),total=sum(ANC),
            pos_LT20=sum(MALLT20),test_LT20=sum(RDTLT20),total_LT20=sum(ANC_LT20),
            pos_GE20=sum(MALGE20),test_GE20=sum(RDTGE20),total_GE20=sum(ANC_GE20))%>%
  arrange(Region,Council,yearmon)
write_csv(tanz_data_all_2017_dist,"./tanz/Patrick/processed_inputs/TZ_ANC_data_district_2014_2017.csv")

tanz_data_all_2017_region<-Full_data%>%
  mutate(yearmon=as.yearmon(paste0(Full_data$Year,"-",Full_data$Month),"%Y-%B"))%>%
  group_by(yearmon,Region)%>%
  dplyr::summarise(count=n(),positive=sum(Mal),tested=sum(RDT),total=sum(ANC),
            pos_LT20=sum(MALLT20),test_LT20=sum(RDTLT20),total_LT20=sum(ANC_LT20),
            pos_GE20=sum(MALGE20),test_GE20=sum(RDTGE20),total_GE20=sum(ANC_GE20))%>%
  arrange(Region,yearmon)

write_csv(tanz_data_all_2017_region,"./tanz/Patrick/processed_inputs/TZ_ANC_data_region_2014_2017.csv")

###Updated data March 2023###
#### read in ANC data ##
# Full_data_2016<-readxl::read_excel("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",sheet = '2016_2019')
# Full_data_2020<-readxl::read_excel("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",sheet = '2020_2021')
# Full_data_2022<-readxl::read_excel("./tanz/Patrick/Data/ANC_Data_Compiled.xlsx",sheet = '2022')
# var.names <- c('region','council','hg','periodid','year','month','first_anc_lt12w','first_anc_gt12w','reattend_anc','reattend_total','tested','positive','llin','ipt2','ipt3','ipt4','hb_measured','hb_lt8.5')
# names(Full_data_2016) <- var.names[-4]
# names(Full_data_2020) <- var.names[-4]
# names(Full_data_2022) <- var.names
# ##Fix 2022 months and years
# Full_data_2022 <- Full_data_2022%>%
#   mutate(temp = month,
#          month=year,
#          year=temp)%>%
#   select(-temp)
# tanz_data_all_16to22 <- plyr::rbind.fill(Full_data_2016,Full_data_2020,Full_data_2022)
# 
# tanz_data_all_16to22$region[tanz_data_all_16to22$region=="Songwe Region"]="Mbeya Region"
# 
# ubungo <- tanz_data_all_16to22[tanz_data_all_16to22$council=='Ubungo Municipal Council',]
# 
# ubungo_clinic_count <- ubungo %>%
#   group_by(year,hg)%>%
#   dplyr::summarise(count=n())
  
#### summarise at district level ##
Full_data$Region[Full_data$Region=="Songwe Region"]="Mbeya Region"
tanz_data_all_16to22_dist <-Full_data%>%
  mutate(yearmon=as.yearmon(paste0(Full_data$Year,"-",Full_data$Month),"%Y-%m"))%>%
  # filter(ANC_pos<=ANC_test)%>%
  group_by(yearmon,Region,Council)%>%
  dplyr::summarise(count=n(),positive=sum(ANC_pos, na.rm = TRUE),tested=sum(ANC_test, na.rm=TRUE),total=sum(Total_re_ad, na.rm = TRUE))
write_csv(tanz_data_all_16to22_dist,"./tanz/Patrick/processed_inputs/TZ_ANC_data_district_2016_2022.csv")
tanz_data_all_16to22_dist <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_district_2016_2022.csv")

tanz_data_all_16to22_region <-Full_data%>%
  mutate(yearmon=as.yearmon(paste0(Full_data$Year,"-",Full_data$Month),"%Y-%m"))%>%
    # filter(ANC_pos<=ANC_test)%>%
  group_by(yearmon,Region)%>%
  dplyr::summarise(count=n(),positive=sum(ANC_pos, na.rm = TRUE),tested=sum(ANC_test, na.rm=TRUE),total=sum(Total_re_ad, na.rm = TRUE))
 write_csv(tanz_data_all_16to22_region,"./tanz/Patrick/processed_inputs/TZ_ANC_data_region_2016_2022.csv")

obs_prev_plot <- ggplot(tanz_data_all_16to22_region)+
  geom_point(aes(x=as.Date(as.yearmon(yearmon)),y=positive/tested),pch = 19,color='darkgray')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ Region)+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC Prevalence')

obs_tested_plot <- ggplot(tanz_data_all_16to22_region)+
  geom_point(aes(x=as.Date(as.yearmon(yearmon)),y=tested),pch = 19,color='darkgray')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ Region)+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC Tested')
lindi_mtwara<-Full_data%>%
  mutate(yearmon=as.yearmon(paste0(Full_data$Year,"-",Full_data$Month),"%Y-%m"))%>%
  filter(!is.na(ANC_test)&ANC_re_ad!=0)%>%
  group_by(yearmon,Region)%>%
  dplyr::summarise(count=n(),positive=sum(ANC_pos,na.rm = TRUE),tested=sum(ANC_test,na.rm = TRUE),total=sum(Total_re_ad,na.rm = TRUE))%>%
  filter(Region == 'Lindi Region' | Region == 'Mtwara Region') %>%
  mutate(period = ifelse(yearmon>=as.yearmon('2018-01-01')&yearmon<=as.yearmon('2019-12-31'),'2018-2019',
                         ifelse(yearmon>=as.yearmon('2020-01-01')&yearmon<=as.yearmon('2021-12-31'),'2020-2021',NA)))
windows(5,10)
ylim.prim <- c(0, 23000)   # school prev limits
ylim.sec <- c(0, 500)    # ANC prev limits

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
colors <- c('Number of women tested' = "#FF7F00", 'Number of women attended' = '#FDBF6F',
            'Number of facilities reporting' = "#1F78B4")
obs_tested_plot_dist <- ggplot(lindi_mtwara)+
  geom_line(aes(x=as.Date(as.yearmon(yearmon)),y=tested,color='Number of women tested'),size=1)+
  geom_line(aes(x=as.Date(as.yearmon(yearmon)),y=total,color='Number of women attended'),size=1)+
  geom_line(aes(x=as.Date(as.yearmon(yearmon)),y=count*b+a,color='Number of facilities reporting'),size=1)+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ Region,nrow = 2)+
  scale_color_manual(values=colors)+
  scale_x_date(date_labels = "'%y")+
  scale_y_continuous("# Women Attended and Tested", sec.axis = sec_axis(~ (. - a)/b, name = "# Facilities Reporting"),limits = c(0,NA)) +
  labs(x='Date')



windows(30,15)
obs_tested_plot_dist <- ggplot(tanz_data_all_16to22_dist)+
  geom_point(aes(x=as.Date(as.yearmon(yearmon)),y=tested),pch = 19,color='darkgray')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ Council)+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC Tested')
obs_attend_plot_dist <- ggplot(tanz_data_all_16to22_dist)+
  geom_point(aes(x=as.Date(as.yearmon(yearmon)),y=total),pch = 19,color='darkgray')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ council)+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC Attendance')
obs_clinics_plot_dist <- ggplot(tanz_data_all_16to22_dist)+
  geom_point(aes(x=as.Date(as.yearmon(yearmon)),y=count),pch = 19,color='darkgray')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ council, scales = 'free_y')+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC Clinics Reporting')

windows(20,20)
ggsave('tanz/figures/obs_tested_plot_dist_2.pdf',plot = obs_tested_plot_dist,width = 30,height=15)
ggsave('tanz/figures/obs_attend_plot_dist_2.pdf',plot = obs_attend_plot_dist,width = 30,height=15)
ggsave('tanz/figures/obs_clinics_plot_dist_2.pdf',plot = obs_clinics_plot_dist,width = 30,height=15)


obs_prev_plot_dist <- ggplot(tanz_data_all_16to22_dist)+
  geom_point(aes(x=as.Date(yearmon),y=positive/tested),pch = 19,color='darkgray')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_wrap(~ council)+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC Prevalence')

View(tanz_data_all_2017)
View(tanz_data_all_16to22)
