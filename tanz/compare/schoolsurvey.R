
tanz_data_all_2017_region <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_region_2014_2017.csv")%>%
  select(1:6)%>%
  filter(yearmon <= as.yearmon('Dec 2015')) %>%
  rename(region = Region)
tanz_data_all_16to22_region <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_region_2016_2022.csv")%>%
  rename(region = Region)

tanz_data_all_14to22_region <- bind_rows(tanz_data_all_2017_region,tanz_data_all_16to22_region)
tanz_data_all_14to22_region <- addCIs(tanz_data_all_14to22_region,tanz_data_all_14to22_region$positive,tanz_data_all_14to22_region$tested)
tanz_data_all_14to22_region$Region <- sub(' Region','',tanz_data_all_14to22_region$region)

tanz_data_all_2017_dist <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_district_2014_2017.csv")%>%
  select(1:6)%>%
  filter(yearmon <= as.yearmon('Dec 2015')) %>%
  rename(region = Region,
         council = Council)
tanz_data_all_16to22_dist <- read_csv("./tanz/Patrick/processed_inputs/TZ_ANC_data_district_2016_2022.csv")%>%
  select(-total)%>%
  rename(region = Region,
         council = Council)

tanz_data_all_14to22_dist <- bind_rows(tanz_data_all_2017_dist,tanz_data_all_16to22_dist)%>%
  mutate(yearmon=as.yearmon(yearmon))%>%
  filter(positive<=tested)

##read in school survey dates
survey_dates_2015 <- readxl::read_excel("./tanz/Patrick/Data/Dates of interview for all survey rounds.xlsx",sheet = '2015')
survey_dates_2017 <- readxl::read_excel("./tanz/Patrick/Data/Dates of interview for all survey rounds.xlsx",sheet = '2017')
survey_dates_2019 <- readxl::read_excel("./tanz/Patrick/Data/Dates of interview for all survey rounds.xlsx",sheet = '2019')
survey_dates_2021 <- readxl::read_excel("./tanz/Patrick/Data/Dates of interview for all survey rounds.xlsx",sheet = '2021')

survey_dates_2015_proc <- survey_dates_2015 %>%
  rename(School = Name_primary_school)%>%
  mutate(yearmon = as.yearmon(paste0(survey_dates_2015$year,"-",survey_dates_2015$Month),"%Y-%m"))%>%
  select(c(Region,Council,School,yearmon))
survey_dates_2017_proc <- survey_dates_2017 %>%
  mutate(yearmon = as.yearmon(paste0(survey_dates_2017$Year,"-",survey_dates_2017$Month),"%Y-%m"))%>%
  select(c(Region,Council,School,Tested,Positive,yearmon))
survey_dates_2019_proc <- survey_dates_2019 %>%
  rename(School = Name_primary_school)%>%
  mutate(yearmon = as.yearmon(paste0(survey_dates_2019$year,"-",survey_dates_2019$Month),"%Y-%m"))%>%
  select(c(Region,Council,School,yearmon))
survey_dates_2021_proc <- survey_dates_2021 %>%
  rename(School = Name_primary_school)%>%
  mutate(yearmon = as.yearmon(paste0(survey_dates_2021$year,"-",survey_dates_2021$Month),"%Y-%m"))%>%
  select(c(Region,Council,School,yearmon))

survey_dates_all <- plyr::rbind.fill(survey_dates_2015_proc,survey_dates_2017_proc,survey_dates_2019_proc,survey_dates_2021_proc)
survey_dates_all$Council <- recode(survey_dates_all$Council,
                                   `Arusha MC`='Arusha CC',
                                   `Dodoma CC`='Dodoma MC',
                                   Kondoa='Kondoa DC',
                                   `Bukombe DC` = 'Bukombe',
                                   `Chato DC` = 'Chato',
                                   `Mbogwe DC` = 'Mbogwe',
                                   `Nyang'hwale DC` = "Nyang'hwale",
                                   `Kilolo DC` = "Kilolo",
                                   `Mpanda TC` = "Mpanda MC",
                                   `Kakonko RC` = 'Kakonko',
                                   `Kasulu RC` = 'Kasulu DC',
                                   `Kibondo RC` = 'Kibondo',
                                   `Kigoma Rural` = 'Kigoma DC',
                                   `Kigoma Urban` = 'Kigoma MC',
                                   `Uvinza DC` = 'Uvinza',
                                   `Rombo DC` = 'Rombo',
                                   `Lindi Rural` = 'Lindi DC',
                                   `Lindi Urban` = 'Lindi MC',
                                   `Hanang DC` = 'Hanang',
                                   `Kiteto DC` = 'Kiteto',
                                   `Simanjiro DC` = 'Simanjiro',
                                   Bunda='Bunda DC',
                                   `Mbeya CC`='Mbeya MC',
                                   `Mtwara Rural` = 'Mtwara DC',
                                   `Mtwara urban` = 'Mtwara MC',
                                   Missungwi='Misungwi',
                                   `Nyamagana MC` = 'Nyamagana',
                                   `Makambako TC` = "Makambako",
                                   `Njombe MC` = 'Njombe TC',
                                   `Wang'ing'ombe` = "Wanging'ombe",
                                   Mbinga='Mbinga DC',
                                   `Songea Rural` = 'Songea DC',
                                   `Songea Urban` = 'Songea MC',
                                   `Meatu DC` = 'Meatu',
                                   `Busega DC` = 'Busega',
                                   Nzega='Nzega DC',
                                   `Lushoto DC` = 'Lushoto',
                                   `Muheza DC` = 'Muheza',
                                   `Tanga MC`='Tanga CC',
                                   .default = survey_dates_all$Council)
survey_dates_all$Region <- recode(survey_dates_all$Region,
                                  `Dar es Salaam`='Dar Es Salaam',
                                  Dsm = 'Dar Es Salaam',
                                  .default = survey_dates_all$Region)


school_counts <- survey_dates_all %>%
  group_by(Region,Council,yearmon) %>%
  summarise(n_schools = n(),
            tested = sum(Tested, na.rm = TRUE),
            positive = sum(Positive, na.rm = TRUE))%>%
  arrange(.by_group = TRUE)


##district key
district_school <- unique(survey_dates_all[c('Region','Council')])
write_csv(district_school,"./tanz/Patrick/Data/districts_school.csv")

district_key <- read_csv("./tanz/Patrick/Data/District_key.csv")

school_counts$council_anc <- district_key$`ANC district`[match(as.character(school_counts$Council), as.character(district_key$school))]

##Merge school dates with ANC data to only keep months and councils with school surveys
tanz_data_all_14to22_dist_sch <- left_join(x=school_counts,y=tanz_data_all_14to22_dist,
                                       by = join_by(council_anc == council, yearmon == yearmon),
                                       suffix = c('_sch','_anc'),
                                       multiple = 'warning')
#add zones
zone_key <- read_csv("./tanz/Patrick/Data/zone_key.csv")
tanz_data_all_14to22_dist_sch$zone <- zone_key$Zones[match(as.character(tanz_data_all_14to22_dist_sch$Region), as.character(zone_key$Region))]

source('shared/addCIs.R')
##Unweighted ANC prevalence in months with school surveys
tanz_anc_14to22_regprev_wt <- tanz_data_all_14to22_dist_sch %>%
  mutate(year = factor(ifelse(year(yearmon)==2014,2015,year(yearmon)),levels=c('2015','2017','2019','2021'),labels=c("'15","'17","'19","'21")),
         Region = ifelse(Region=='Songwe','Mbeya',Region),
         positive_wt = positive_anc*n_schools,
         tested_wt = tested_anc*n_schools)%>%
  group_by(Region,year)%>%
  summarise(school_count = sum(n_schools),
            positive = sum(positive_anc, na.rm = TRUE),
            tested = sum(tested_anc, na.rm = TRUE),
            positive_wt = sum(positive_wt, na.rm = TRUE),
            tested_wt = sum(tested_wt, na.rm = TRUE),
            svy_dt = median(yearmon))
tanz_anc_14to22_regprev_wt <- addCIs(tanz_anc_14to22_regprev_wt,tanz_anc_14to22_regprev_wt$positive_wt,tanz_anc_14to22_regprev_wt$tested_wt)%>%
  rename(mean_wt = mean,
         upper_wt = upper,
         lower_wt = lower)
tanz_anc_14to22_regprev_wt <- addCIs(tanz_anc_14to22_regprev_wt,tanz_anc_14to22_regprev_wt$positive,tanz_anc_14to22_regprev_wt$tested)
  
tanz_anc_14to22_zoneprev_wt <- tanz_data_all_14to22_dist_sch %>%
  mutate(year = factor(ifelse(year(yearmon)==2014,2015,year(yearmon)),levels=c('2015','2017','2019','2021'),labels=c("'15","'17","'19","'21")),
         positive_wt = positive_anc*n_schools,
         tested_wt = tested_anc*n_schools)%>%
  group_by(zone,year)%>%
  summarise(school_count = sum(n_schools),
            positive = sum(positive_anc, na.rm = TRUE),
            tested = sum(tested_anc, na.rm = TRUE),
            positive_wt = sum(positive_wt, na.rm = TRUE),
            tested_wt = sum(tested_wt, na.rm = TRUE),
            svy_dt = median(yearmon))
tanz_anc_14to22_allprev_wt <- tanz_data_all_14to22_dist_sch %>%
  mutate(year = factor(ifelse(year(yearmon)==2014,2015,year(yearmon)),levels=c('2015','2017','2019','2021'),labels=c("'15","'17","'19","'21")),
         positive_wt = positive_anc*n_schools,
         tested_wt = tested_anc*n_schools)%>%
  group_by(year)%>%
  summarise(school_count = sum(n_schools),
            positive = sum(positive_anc, na.rm = TRUE),
            tested = sum(tested_anc, na.rm = TRUE),
            positive_wt = sum(positive_wt, na.rm = TRUE),
            tested_wt = sum(tested_wt, na.rm = TRUE),
            svy_dt = median(yearmon))%>%
  mutate(zone='All')
tanz_anc_14to22_zoneprev_wt <- bind_rows(tanz_anc_14to22_zoneprev_wt,tanz_anc_14to22_allprev_wt)
tanz_anc_14to22_zoneprev_wt <- addCIs(tanz_anc_14to22_zoneprev_wt,tanz_anc_14to22_zoneprev_wt$positive_wt,tanz_anc_14to22_zoneprev_wt$tested_wt)%>%
  rename(mean_wt = mean,
         upper_wt = upper,
         lower_wt = lower)
tanz_anc_14to22_zoneprev_wt <- addCIs(tanz_anc_14to22_zoneprev_wt,tanz_anc_14to22_zoneprev_wt$positive,tanz_anc_14to22_zoneprev_wt$tested)

table(tanz_anc_14to22_prev_unwt$Region,tanz_anc_14to22_prev_unwt$year)  
province_grid<-read_csv("./tanz/processed_inputs/Province_grid.csv")
province_grid$name=province_grid$code=province_grid$NAME_1

windows(10,10)
obs_regprev_plot <- ggplot(tanz_anc_14to22_regprev_wt)+
  geom_point(aes(x=as.Date(svy_dt),y=mean),pch = 19,color='darkgray')+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower,ymax=upper),width = 0,color='darkgray')+
  geom_point(aes(x=as.Date(svy_dt),y=mean_wt),pch = 19,color="#414487FF")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_wt,ymax=upper_wt),width = 0,color="#414487FF")+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ Region, grid = province_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC Prevalence')
zone_grid<-read_csv("./tanz/processed_inputs/zone_grid.csv")
zone_grid$name=zone_grid$code=zone_grid$NAME_1

zone_prev <- readxl::read_excel("./tanz/Patrick/Data/zone_school_prev.xlsx")
zone_prev<- pivot_longer(zone_prev,cols = `2015`:`2021`)%>%
  mutate(year = factor(name,levels=c('2015','2017','2019','2021'),labels=c("'15","'17","'19","'21")),
         prev = value/100)%>%
  select(-c(name,value))

tanz_anc_14to22_zoneprev_wt <- left_join(tanz_anc_14to22_zoneprev_wt,zone_prev, by=c('zone','year'))
  
tanz_data_all_14to22_dist <- addCIs(tanz_data_all_14to22_dist,tanz_data_all_14to22_dist$positive,tanz_data_all_14to22_dist$tested)

obs_zoneprev_plot <- ggplot(tanz_anc_14to22_zoneprev_wt)+
  geom_point(aes(x=as.Date(svy_dt),y=mean_wt),pch = 19,color="#414487FF")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_wt,ymax=upper_wt),width = 0,color="#414487FF")+
  geom_point(aes(x=as.Date(svy_dt),y=prev),pch = 19, color = 'darkgrey')+
  geom_line(data=tanz_data_all_14to22_dist,aes(x=as.Date(yearmon),y=mean),color="#414487FF")+
  geom_ribbon(data=tanz_data_all_14to22_dist,aes(x=as.Date(yearmon),ymin=lower,ymax=upper),size=1,fill="#414487FF",alpha=0.2)+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ zone, grid = zone_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y")+
  labs(x='Date',y='ANC (purple)/School (grey) Prevalence')

##Add monthly ANC prev by zone
tanz_data_all_14to22_dist$Region <- sub(' Region','',tanz_data_all_14to22_dist$region)
tanz_data_all_14to22_dist$zone <- zone_key$Zones[match(as.character(tanz_data_all_14to22_dist$Region), as.character(zone_key$Region))]

tanz_data_all_14to22_zone <- tanz_data_all_14to22_dist%>%
  group_by(zone,yearmon)%>%
  summarise(positive = sum(positive, na.rm = TRUE),
            tested = sum(tested, na.rm = TRUE),
            svy_dt = median(yearmon))
tanz_data_all_14to22_all <- tanz_data_all_14to22_dist%>%
  group_by(yearmon)%>%
  summarise(positive = sum(positive, na.rm = TRUE),
            tested = sum(tested, na.rm = TRUE),
            svy_dt = median(yearmon))%>%
  mutate(zone='All')
tanz_data_all_14to22_zone <- bind_rows(tanz_data_all_14to22_zone,tanz_data_all_14to22_all)
tanz_data_all_14to22_zone <- addCIs(tanz_data_all_14to22_zone,tanz_data_all_14to22_zone$positive,tanz_data_all_14to22_zone$tested)
  
ylim.prim <- c(0, 0.6)   # school prev limits
ylim.sec <- c(0, 0.3)    # ANC prev limits

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]

obs_zoneprev_plot <- ggplot(tanz_anc_14to22_zoneprev_wt)+
  geom_line(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),y=a + mean*b),color="#27808EFF",size=1)+
  geom_ribbon(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),ymin=a + lower*b,ymax=a+upper*b),fill="#27808EFF",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_wt*b),pch = 19,size=2.5,color="#27808EFF")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_wt*b,ymax=a+upper_wt*b),width = 0,color="#27808EFF")+
  geom_point(aes(x=as.Date(svy_dt),y=prev),pch = 19, size=2.5,color = '#666666')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ zone, grid = zone_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y")+
  scale_y_continuous("School Survey Prevalence (grey)", sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)")) +
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10))

obs_zoneprev_plot <- ggplot(tanz_anc_14to22_zoneprev_wt)+
  geom_line(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),y=a + mean*b),color="white",size=1)+
  geom_ribbon(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),ymin=a + lower*b,ymax=a+upper*b),fill="white",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_wt*b),pch = 19,size=2.5,color="white")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_wt*b,ymax=a+upper_wt*b),width = 0,color="white")+
  geom_point(aes(x=as.Date(svy_dt),y=prev),pch = 19, size=2.5,color = '#666666')+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ zone, grid = zone_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y")+
  scale_y_continuous("School Survey Prevalence (grey)", sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)")) +
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10))


##Compare DHS and ANC
dhs_zones <- readxl::read_excel("./tanz/Patrick/Data/u5yo_prev_dhs_tanz_2017-2022.xlsx",sheet = 'Zones')
dhs_regions <- readxl::read_excel("./tanz/Patrick/Data/u5yo_prev_dhs_tanz_2017-2022.xlsx",sheet = 'Regions')

dhs_zones <- dhs_zones %>%
  mutate(positive = round(prev*tested/100),
         svy_dt = as.yearmon(ifelse(year==2017,as.yearmon('Nov 2017'),as.yearmon('May 2022'))))
dhs_all <- dhs_zones %>%
  group_by(svy_dt) %>%
  summarise(positive = sum(positive),
            tested = sum(tested)) %>%
  mutate(Zone = 'All')
dhs_zones <- bind_rows(dhs_zones,dhs_all)
dhs_zones <- addCIs(dhs_zones,dhs_zones$positive,dhs_zones$tested)

dhs_regions <- dhs_regions %>%
  mutate(positive = round(prev*tested/100),
         svy_dt = as.yearmon(ifelse(year==2017,as.yearmon('Nov 2017'),as.yearmon('May 2022'))))
dhs_all <- dhs_regions %>%
  group_by(svy_dt) %>%
  summarise(positive = sum(positive),
            tested = sum(tested)) %>%
  mutate(Zone = 'All',
         Regions = 'All')
dhs_zones <- bind_rows(dhs_zones,dhs_all)
dhs_regions <- bind_rows(dhs_regions,dhs_all)

dhs_regions <- addCIs(dhs_regions,dhs_regions$positive,dhs_regions$tested)
dhs_regions$Regions <- ifelse(dhs_regions$Regions=='Dar es Salaam','Dar Es Salaam',dhs_regions$Regions)

tanz_anc_14to22_regprev_dhs <- tanz_data_all_14to22_dist %>%
  filter((yearmon>=as.yearmon('Oct 2017')&yearmon<=as.yearmon('Dec 2017'))|(yearmon>=as.yearmon('Feb 2022')&yearmon<=as.yearmon('Jul 2022')))%>%
  mutate(Region = ifelse(Region=='Songwe','Mbeya',Region),
         svy_dt = as.yearmon(ifelse(year(yearmon)==2017,as.yearmon('Nov 2017'),as.yearmon('May 2022'))))%>%
  group_by(Region,svy_dt)%>%
  summarise(positive = sum(positive, na.rm = TRUE),
            tested = sum(tested, na.rm = TRUE))
tanz_anc_14to22_zoneprev_dhs <- tanz_data_all_14to22_dist %>%
  filter((yearmon>=as.yearmon('Oct 2017')&yearmon<=as.yearmon('Dec 2017'))|(yearmon>=as.yearmon('Feb 2022')&yearmon<=as.yearmon('Jul 2022')))%>%
  mutate(svy_dt = as.yearmon(ifelse(year(yearmon)==2017,as.yearmon('Nov 2017'),as.yearmon('May 2022'))))%>%
  group_by(zone,svy_dt)%>%
  summarise(positive = sum(positive, na.rm = TRUE),
            tested = sum(tested, na.rm = TRUE))
tanz_anc_14to22_allprev_dhs <- tanz_data_all_14to22_dist %>%
  filter((yearmon>=as.yearmon('Oct 2017')&yearmon<=as.yearmon('Dec 2017'))|(yearmon>=as.yearmon('Feb 2022')&yearmon<=as.yearmon('Jul 2022')))%>%
  mutate(svy_dt = as.yearmon(ifelse(year(yearmon)==2017,as.yearmon('Nov 2017'),as.yearmon('May 2022'))))%>%
  group_by(svy_dt)%>%
  summarise(positive = sum(positive, na.rm = TRUE),
            tested = sum(tested, na.rm = TRUE))%>%
  mutate(zone = 'All',
         Region = 'All')
tanz_anc_14to22_zoneprev_dhs <- bind_rows(tanz_anc_14to22_zoneprev_dhs,tanz_anc_14to22_allprev_dhs)
tanz_anc_14to22_zoneprev_dhs <- left_join(tanz_anc_14to22_zoneprev_dhs,dhs_zones,by=c('zone'='Zone','svy_dt'),suffix=c('_anc','_dhs'))
tanz_anc_14to22_zoneprev_dhs <- tanz_anc_14to22_zoneprev_dhs %>%
  rename(mean_dhs = mean,
         upper_dhs = upper,
         lower_dhs = lower)
tanz_anc_14to22_zoneprev_dhs <- addCIs(tanz_anc_14to22_zoneprev_dhs,tanz_anc_14to22_zoneprev_dhs$positive_anc,tanz_anc_14to22_zoneprev_dhs$tested_anc)
tanz_anc_14to22_zoneprev_dhs <- tanz_anc_14to22_zoneprev_dhs %>%
  rename(mean_anc = mean,
         upper_anc = upper,
         lower_anc = lower)

tanz_anc_14to22_regprev_dhs <- bind_rows(tanz_anc_14to22_regprev_dhs,tanz_anc_14to22_allprev_dhs)
tanz_anc_14to22_regprev_dhs <- left_join(tanz_anc_14to22_regprev_dhs,dhs_regions,by=c('Region'='Regions','svy_dt'),suffix=c('_anc','_dhs'))
tanz_anc_14to22_regprev_dhs <- tanz_anc_14to22_regprev_dhs %>%
  rename(mean_dhs = mean,
         upper_dhs = upper,
         lower_dhs = lower)
tanz_anc_14to22_regprev_dhs <- addCIs(tanz_anc_14to22_regprev_dhs,tanz_anc_14to22_regprev_dhs$positive_anc,tanz_anc_14to22_regprev_dhs$tested_anc)
tanz_anc_14to22_regprev_dhs <- tanz_anc_14to22_regprev_dhs %>%
  rename(mean_anc = mean,
         upper_anc = upper,
         lower_anc = lower)

ylim.prim <- c(0, 0.3)   # school prev limits
ylim.sec <- c(0, 0.3)    # ANC prev limits

b <- diff(ylim.prim)/diff(ylim.sec)
a <- ylim.prim[1] - b*ylim.sec[1]
windows(10,10)
obs_zoneprev_plot_dhs <- ggplot(tanz_anc_14to22_zoneprev_dhs)+
  geom_line(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),y=a + mean*b),color="#27808EFF",size=1)+
  geom_ribbon(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),ymin=a + lower*b,ymax=a+upper*b),fill="#27808EFF",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_anc*b),pch = 19,size=2.5,color="#27808EFF")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_anc*b,ymax=a+upper_anc*b),width = 0,color="#27808EFF")+
  geom_point(aes(x=as.Date(svy_dt),y=mean_dhs),pch = 19, size=2.5,color = '#666666')+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs),width = 0,color="#666666")+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ zone, grid = zone_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y",limits = as.Date(c('2016-01-01',Inf)))+
  scale_y_continuous("DHS Prevalence (grey)", sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)")) +
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10))

obs_zoneprev_plot_dhs <- ggplot(tanz_anc_14to22_zoneprev_dhs)+
  geom_line(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),y=a + mean*b),color="white",size=1)+
  geom_ribbon(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),ymin=a + lower*b,ymax=a+upper*b),fill="white",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_anc*b),pch = 19,size=2.5,color="white")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_anc*b,ymax=a+upper_anc*b),width = 0,color="white")+
  geom_point(aes(x=as.Date(svy_dt),y=mean_dhs),pch = 19, size=2.5,color = '#666666')+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs),width = 0,color="#666666")+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ zone, grid = zone_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y",limits = as.Date(c('2016-01-01',Inf)))+
  scale_y_continuous("DHS Prevalence (grey)", sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)")) +
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10))
colors <- c('ANC' = "#27808EFF", 'Under 5 year olds' = '#666666')
windows(4,4)
obs_tzprev_plot <- ggplot(tanz_anc_14to22_zoneprev_dhs[tanz_anc_14to22_zoneprev_dhs$zone=='All',])+
  geom_line(data=tanz_data_all_14to22_zone[tanz_data_all_14to22_zone$zone=='All',],aes(x=as.Date(svy_dt),y=a + mean*b,color='ANC'),size=1)+
  geom_ribbon(data=tanz_data_all_14to22_zone[tanz_data_all_14to22_zone$zone=='All',],aes(x=as.Date(svy_dt),ymin=a + lower*b,ymax=a+upper*b),fill="#27808EFF",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_anc*b,color='ANC'),pch = 19,size=2.5)+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_anc*b,ymax=a+upper_anc*b,color='ANC'),width = 0)+
  geom_point(aes(x=as.Date(svy_dt),y=mean_dhs, color = 'Under 5 year olds'),pch = 19, size=2.5)+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs, color = 'Under 5 year olds'),width = 0)+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  scale_x_date(date_labels = "'%y",limits = as.Date(c('2016-01-01',Inf)))+
  scale_y_continuous("Prevalence",limits = c(0,0.1),breaks = c(0,0.02,0.04,0.06,0.08,0.1))+#, sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)")) +
  scale_color_manual(values=colors)+
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

obs_tzprev_plot <- ggplot(tanz_anc_14to22_zoneprev_dhs[tanz_anc_14to22_zoneprev_dhs$zone=='All',])+
  geom_line(data=tanz_data_all_14to22_zone[tanz_data_all_14to22_zone$zone=='All',],aes(x=as.Date(svy_dt),y=a + mean*b),color='white',size=1)+
  geom_ribbon(data=tanz_data_all_14to22_zone[tanz_data_all_14to22_zone$zone=='All',],aes(x=as.Date(svy_dt),ymin=a + lower*b,ymax=a+upper*b),fill="white",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_anc*b,color='ANC'),pch = 19,size=2.5)+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_anc*b,ymax=a+upper_anc*b,color='ANC'),width = 0)+
  geom_point(aes(x=as.Date(svy_dt),y=mean_dhs, color = 'Under 5 year olds'),pch = 19, size=2.5)+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs, color = 'Under 5 year olds'),width = 0)+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  scale_x_date(date_labels = "'%y",limits = as.Date(c('2016-01-01',Inf)))+
  scale_y_continuous("Prevalence",limits = c(0,0.1),breaks = c(0,0.02,0.04,0.06,0.08,0.1))+#, sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)")) +
  scale_color_manual(values=colors)+
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

obs_tzprev_plot <- ggplot(tanz_anc_14to22_zoneprev_dhs[tanz_anc_14to22_zoneprev_dhs$zone=='All',])+
  geom_line(data=tanz_data_all_14to22_zone[tanz_data_all_14to22_zone$zone=='All',],aes(x=as.Date(svy_dt),y=a + mean*b),color='white',size=1)+
  geom_ribbon(data=tanz_data_all_14to22_zone[tanz_data_all_14to22_zone$zone=='All',],aes(x=as.Date(svy_dt),ymin=a + lower*b,ymax=a+upper*b),fill="white",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_anc*b),color='white',pch = 19,size=2.5)+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_anc*b,ymax=a+upper_anc*b),color='white',width = 0)+
  geom_point(aes(x=as.Date(svy_dt),y=mean_dhs, color = 'Under 5 year olds'),pch = 19, size=2.5)+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs, color = 'Under 5 year olds'),width = 0)+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  scale_x_date(date_labels = "'%y",limits = as.Date(c('2016-01-01',Inf)))+
  scale_y_continuous("Prevalence",limits = c(0,0.1),breaks = c(0,0.02,0.04,0.06,0.08,0.1))+#, sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)")) +
  scale_color_manual(values=colors)+
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

windows(10,10)
obs_regprev_plot_dhs <- ggplot(tanz_anc_14to22_regprev_dhs)+
  geom_line(data=tanz_data_all_14to22_region[!(tanz_data_all_14to22_region$Region=='Mbeya'&tanz_data_all_14to22_region$yearmon=='Oct 2020'),],aes(x=as.Date(as.yearmon(yearmon)),y=a + mean*b),color="#27808EFF",size=0.7)+
  geom_ribbon(data=tanz_data_all_14to22_region[!(tanz_data_all_14to22_region$Region=='Mbeya'&tanz_data_all_14to22_region$yearmon=='Oct 2020'),],aes(x=as.Date(as.yearmon(yearmon)),ymin=a + lower*b,ymax=a+upper*b),fill="#27808EFF",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_anc*b),pch = 19,size=2.5,color="#27808EFF")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_anc*b,ymax=a+upper_anc*b),width = 0,color="#27808EFF")+
  geom_point(aes(x=as.Date(svy_dt),y=mean_dhs),pch = 19, size=2.5,color = '#666666')+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs),width = 0,color="#666666")+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ Region, grid = province_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y",limits = as.Date(c('2016-01-01',Inf)))+
  scale_y_continuous("DHS Prevalence (grey)", sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)"),limits = c(0,0.3)) +
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y.right = element_blank())

obs_regprev_plot_dhs_2017 <- ggplot(tanz_anc_14to22_regprev_dhs)+
  geom_line(data=tanz_data_all_14to22_region[!(tanz_data_all_14to22_region$Region=='Mbeya'&tanz_data_all_14to22_region$yearmon=='Oct 2020'),],aes(x=as.Date(as.yearmon(yearmon)),y=a + mean*b),color="#27808EFF",size=0.7)+
  geom_ribbon(data=tanz_data_all_14to22_region[!(tanz_data_all_14to22_region$Region=='Mbeya'&tanz_data_all_14to22_region$yearmon=='Oct 2020'),],aes(x=as.Date(as.yearmon(yearmon)),ymin=a + lower*b,ymax=a+upper*b),fill="#27808EFF",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_anc*b),pch = 19,size=2.5,color="#27808EFF")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_anc*b,ymax=a+upper_anc*b),width = 0,color="#27808EFF")+
  geom_point(aes(x=as.Date(svy_dt),y=mean_dhs),pch = 19, size=2.5,color = '#666666')+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs),width = 0,color="#666666")+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ Region, grid = province_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "%b",date_breaks='3 months',limits = as.Date(c('2017-01-01','2017-12-31')))+
  scale_y_continuous("DHS Prevalence (grey)", sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)"),limits = c(0,0.3)) +
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8,angle=45,hjust=1),
        axis.text.y.right = element_blank())

obs_regprev_plot_dhs <- ggplot(tanz_anc_14to22_regprev_dhs)+
  geom_line(data=tanz_data_all_14to22_region,aes(x=as.Date(as.yearmon(yearmon)),y=a + mean*b),color="white",size=0.7)+
  geom_ribbon(data=tanz_data_all_14to22_region,aes(x=as.Date(as.yearmon(yearmon)),ymin=a + lower*b,ymax=a+upper*b),fill="white",alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_anc*b),pch = 19,size=2.5,color="#27808EFF")+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_anc*b,ymax=a+upper_anc*b),width = 0,color="#27808EFF")+
  geom_point(aes(x=as.Date(svy_dt),y=mean_dhs),pch = 19, size=2.5,color = '#666666')+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs),width = 0,color="#666666")+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ Region, grid = province_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y",limits = as.Date(c('2016-01-01',Inf)))+
  scale_y_continuous("DHS Prevalence (grey)", sec.axis = sec_axis(~ (. - a)/b, name = "ANC Survey Prevalence (blue)")) +
  labs(x='Date')+
  theme(strip.text.x = element_text(size = 7),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y.right = element_blank())

colors <- c('ANC' = "#27808EFF", 'Under 5 year olds' = '#D95F02', 'School age' = '#666666')
##Combine all
obs_zoneprev_all_plot <- ggplot(tanz_anc_14to22_zoneprev_wt)+
  geom_line(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),y=a + mean*b,color="ANC"),size=1)+
  geom_ribbon(data=tanz_data_all_14to22_zone,aes(x=as.Date(svy_dt),ymin=a + lower*b,ymax=a+upper*b),fill='#27808EFF',alpha=0.2)+
  geom_point(aes(x=as.Date(svy_dt),y=a+mean_wt*b,color="ANC"),pch = 19,size=2.5)+
  geom_errorbar(aes(x=as.Date(svy_dt),ymin=a+lower_wt*b,ymax=a+upper_wt*b,color="ANC"),width = 0)+
  geom_point(aes(x=as.Date(svy_dt),y=prev,color = 'School age'),pch = 19, size=2.5)+
  geom_point(data=tanz_anc_14to22_zoneprev_dhs,aes(x=as.Date(svy_dt),y=mean_dhs,color = 'Under 5 year olds'),pch = 19, size=2.5)+
  geom_errorbar(data=tanz_anc_14to22_zoneprev_dhs,aes(x=as.Date(svy_dt),ymin=lower_dhs,ymax=upper_dhs,color='Under 5 year olds'),width = 0)+
  # scale_color_manual(values=colors)+
  # scale_fill_manual(values=colors)+
  facet_geo(~ zone, grid = zone_grid%>%
              select(row,col,code,name))+
  scale_x_date(date_labels = "'%y")+
  scale_y_continuous("Malaria Prevalence") +
  scale_color_manual(values = colors)+
  labs(x='Date',color='Legend')+
  theme(strip.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
