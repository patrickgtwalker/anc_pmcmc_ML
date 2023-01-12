library(DiagrammeR)
library(greta)
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

addCIs<-function(df,Ys,Ns){
  df$mean<-NA
  df$upper<-NA
  df$lower<-NA
  CIs<-binom.confint(Ys,Ns,method="exact")
  df$mean[Ns>0]<-CIs$mean[Ns>0]
  df$upper[Ns>0]<-CIs$upper[Ns>0]
  df$lower[Ns>0]<-CIs$lower[Ns>0]
  return(df)
}

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
         

#Basic missing check
table(is.na(NG_ANC_mother$age)) #Number missing: 0
table(NG_ANC_mother$age)
table(is.na(NG_ANC_mother$month)) #Number missing: 0
table(NG_ANC_mother$month) #Just month - not exact date
table(is.na(NG_ANC_mother$primigrav)) #Number missing: 0
table(NG_ANC_mother$primigrav)
table(is.na(NG_ANC_mother$prev_pregs)) #Number missing: 0
table(NG_ANC_mother$prev_pregs) #No. of primigrav=1 (1940) does not equal prev_pregs=0 (1947)
table(is.na(NG_ANC_mother$mal_symp)) #Number missing: 0
table(NG_ANC_mother$mal_symp)
table(is.na(NG_ANC_mother$rdt)) #Number missing: 0
table(NG_ANC_mother$rdt)
table(is.na(NG_ANC_mother$date_mens)) #Number missing: 0
table(NG_ANC_mother$date_mens)

#Check noted inconsistency between primigrav and prev_pregs
grav_check <- NG_ANC_mother[NG_ANC_mother$primigrav==0&NG_ANC_mother$prev_pregs==0,]

#Frequency plots
ggplot(NG_ANC_mother,aes(x=age))+
  geom_histogram()
ggplot(NG_ANC_mother,aes(x=age))+
  geom_bar()
ggplot(NG_ANC_mother,aes(x=month))+
  geom_bar()
ggplot(NG_ANC_mother,aes(x=primigrav))+
  geom_bar()
ggplot(NG_ANC_mother,aes(x=grav))+
  geom_bar()
ggplot(NG_ANC_mother,aes(x=month,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(NG_ANC_mother,aes(x=age,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(NG_ANC_mother,aes(x=age,fill=as.factor(rdt)))+
  geom_bar(position='stack')
ggplot(NG_ANC_mother,aes(x=grav,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(NG_ANC_mother,aes(x=grav,fill=as.factor(rdt)))+
  geom_bar(position='stack')
ggplot(NG_ANC_mother,aes(x=date_mens))+
  geom_bar()
ggplot(NG_ANC_mother,aes(x=date_mens,fill=as.factor(rdt)))+
  geom_bar(position='stack')
ggplot(NG_ANC_mother,aes(x=date_mens,fill=as.factor(rdt)))+
  geom_bar(position='fill')

#Missing and values for character variables
table(NG_ANC_mother$hf_name_clean)
table(is.na(NG_ANC_mother$hf_name_clean)) #no missing
table(NG_ANC_mother$lga)
table(is.na(NG_ANC_mother$lga)) #7 missing
table(NG_ANC_mother$ward)
table(is.na(NG_ANC_mother$ward)) #0 missing; many different spelling/spacings

#Group data and calculate prevalence and CI by month and gravidity
NG_ANC_mother_grouped <- NG_ANC_mother %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         age_cat=cut(age,breaks=c(12,15,20,25,30,35,40,Inf),labels=c('12-15 years','16-20 years','21-25 years','26-30 years','31-35 years','36-40 years','Over 40 years'))
  )%>%
  group_by(month,grav_cat,age_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_ANC_mother_grouped <- addCIs(NG_ANC_mother_grouped,NG_ANC_mother_grouped$positive,NG_ANC_mother_grouped$total)

windows(height=40,width=50)
ggplot(data=NG_ANC_mother_grouped,aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~grav_cat,dir="v")+theme(legend.position="none")

#Group by grav
NG_ANC_mother_grouped_grav <- NG_ANC_mother %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")))%>%
  group_by(month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_ANC_mother_grouped_grav <- addCIs(NG_ANC_mother_grouped_grav,NG_ANC_mother_grouped_grav$positive,NG_ANC_mother_grouped_grav$total)

#Grav facetted
ggplot(data=NG_ANC_mother_grouped_grav,aes(x=month,y=mean))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~grav_cat,dir="v")+theme(legend.position="none")

#Grav colored
ggplot(data=NG_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")

#Group by age
NG_ANC_mother_grouped_age <- NG_ANC_mother %>%
  mutate(age_cat=cut(age,breaks=c(12,15,20,25,30,35,40,Inf),labels=c('12-15 years','16-20 years','21-25 years','26-30 years','31-35 years','36-40 years','Over 40 years')))%>%
  group_by(month,age_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_ANC_mother_grouped_age <- addCIs(NG_ANC_mother_grouped_age,NG_ANC_mother_grouped_age$positive,NG_ANC_mother_grouped_age$total)

#Age facetted
ggplot(data=NG_ANC_mother_grouped_age,aes(x=month,y=mean))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~age_cat,dir="v")+theme(legend.position="none")

#Age colored
ggplot(data=NG_ANC_mother_grouped_age[NG_ANC_mother_grouped_age$age_cat!='12-15 years'&NG_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")

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


#Basic missing check
table(is.na(MZ_ANC_mother$age)) #Number missing: 12
table(MZ_ANC_mother$age) #Also age out of reasonable ranges (<12 and >50)
table(is.na(MZ_ANC_mother$month)) #Number missing: 2
table(MZ_ANC_mother$month) #Exact dates
table(is.na(MZ_ANC_mother$primigrav)) #Number missing: 12
table(MZ_ANC_mother$primigrav) #Yes vs No
table(is.na(MZ_ANC_mother$prev_pregs)) #Number missing: 2306 (missing if primigrav=Yes)
table(MZ_ANC_mother$prev_pregs) #Some unreasonable values (>20)
table(is.na(MZ_ANC_mother$mal_symp)) #Number missing: 9
table(MZ_ANC_mother$mal_symp) #3 coded as number 2
table(is.na(MZ_ANC_mother$rdt)) #Number missing: 47
table(MZ_ANC_mother$rdt) #Several other codings besides Negative and positive
table(is.na(MZ_ANC_mother$date_mens)) #Number missing: 331
table(MZ_ANC_mother$date_mens) #Exact date

#Check noted inconsistency between primigrav and prev_pregs
grav_check <- MZ_ANC_mother[MZ_ANC_mother$primigrav=='No'&MZ_ANC_mother$prev_pregs==0,]

#Frequency plots
ggplot(MZ_ANC_mother,aes(x=age))+
  geom_histogram()
ggplot(MZ_ANC_mother,aes(x=age))+
  geom_bar()
ggplot(MZ_ANC_mother,aes(x=month))+
  geom_bar()
ggplot(MZ_ANC_mother,aes(x=primigrav))+
  geom_bar()
ggplot(MZ_ANC_mother,aes(x=grav))+
  geom_bar()
ggplot(MZ_ANC_mother,aes(x=month,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(MZ_ANC_mother,aes(x=age,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(MZ_ANC_mother,aes(x=age,fill=as.factor(rdt)))+
  geom_bar(position='stack')
ggplot(MZ_ANC_mother,aes(x=grav,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(MZ_ANC_mother,aes(x=grav,fill=as.factor(rdt)))+
  geom_bar(position='stack')


#Missing and values for character variables
table(MZ_ANC_mother$HF)
table(is.na(MZ_ANC_mother$HF)) #6 missing
table(MZ_ANC_mother$district_n) # 3 observations have numbers as district name
table(is.na(MZ_ANC_mother$district_n)) #3 missing
table(MZ_ANC_mother$village)
table(is.na(MZ_ANC_mother$village)) #9 missing; many different spelling/spacings

#Group data and calculate prevalence and CI by month and gravidity
MZ_ANC_mother_grouped <- MZ_ANC_mother %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         age_cat=cut(age,breaks=c(12,15,20,25,30,35,40,Inf),labels=c('12-15 years','16-20 years','21-25 years','26-30 years','31-35 years','36-40 years','Over 40 years')),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age>=12&age<50) %>%
  filter(!is.na(age_cat)) %>%
  filter(!is.na(grav))%>%
  group_by(month,grav_cat,age_cat,.drop=FALSE)%>%
  summarise(positive=sum(rdt),total=n())
MZ_ANC_mother_grouped <- addCIs(MZ_ANC_mother_grouped,MZ_ANC_mother_grouped$positive,MZ_ANC_mother_grouped$total)

windows(height=40,width=50)
ggplot(data=MZ_ANC_mother_grouped,aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~grav_cat,dir="v")+theme(legend.position="none")

#Group by grav
MZ_ANC_mother_grouped_grav <- MZ_ANC_mother %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA))))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age>=12&age<50) %>%
  filter(!is.na(grav))%>%
  group_by(month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_ANC_mother_grouped_grav <- addCIs(MZ_ANC_mother_grouped_grav,MZ_ANC_mother_grouped_grav$positive,MZ_ANC_mother_grouped_grav$total)

#Grav facetted
ggplot(data=MZ_ANC_mother_grouped_grav,aes(x=month,y=mean))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~grav_cat,dir="v")+theme(legend.position="none")

#Grav colored
ggplot(data=MZ_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")

#Group by age
MZ_ANC_mother_grouped_age <- MZ_ANC_mother %>%
  mutate(age_cat=cut(age,breaks=c(12,15,20,25,30,35,40,Inf),labels=c('12-15 years','16-20 years','21-25 years','26-30 years','31-35 years','36-40 years','Over 40 years')),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA))))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age>=12&age<50) %>%
  filter(!is.na(age_cat)) %>%
  filter(!is.na(grav))%>%
  group_by(month,age_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_ANC_mother_grouped_age <- addCIs(MZ_ANC_mother_grouped_age,MZ_ANC_mother_grouped_age$positive,MZ_ANC_mother_grouped_age$total)

#Age facetted
ggplot(data=MZ_ANC_mother_grouped_age,aes(x=month,y=mean))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~age_cat,dir="v")+theme(legend.position="none")

#Age colored
ggplot(data=MZ_ANC_mother_grouped_age[MZ_ANC_mother_grouped_age$age_cat!='12-15 years'&MZ_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")

#########################
######Burkina Faso#######
#########################
#Read in Mozambique data and rename most used variables
BF_ANC_mother <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Burkina Faso/ANC-based surveillance/ANC_BF_5.16.2022.xlsx',
                              xl.sheet = 'ANC_woman')%>%
  dplyr::rename(primigrav = Q04_first_preg,
         prev_pregs = Q05_nbr_prior_preg,
         mal_symp = Q06_mal_symp,
         rdt = RDT,
         age = Age) %>%
  dplyr::mutate(month = as.yearmon(as.Date(Date)),
         grav = ifelse(primigrav=='Oui',1,prev_pregs+1))


#Basic missing check
table(is.na(BF_ANC_mother$age)) #Number missing: 0
table(BF_ANC_mother$age) #Also age out of reasonable ranges (age = 88, refused/dk code?)
table(is.na(BF_ANC_mother$Date)) #Number missing: 0
table(BF_ANC_mother$Date) #Exact dates
table(is.na(BF_ANC_mother$primigrav)) #Number missing: 0
table(BF_ANC_mother$primigrav) #Non vs Oui
table(is.na(BF_ANC_mother$prev_pregs)) #Number missing: 3151 (missing if primigrav=Oui)
table(BF_ANC_mother$prev_pregs) 
table(is.na(BF_ANC_mother$mal_symp)) #Number missing: 0
table(BF_ANC_mother$mal_symp) #Non vs Oui
table(is.na(BF_ANC_mother$rdt)) #Number missing: 0
table(BF_ANC_mother$rdt) #Negatif vs Positif
table(is.na(BF_ANC_mother$Q03_last_period)) #Number missing: 0
table(BF_ANC_mother$Q03_last_period) #Date, Q03a_last_period describes whether this date is an estimate or exact

#Check noted inconsistency between primigrav and prev_pregs
grav_check <- BF_ANC_mother[BF_ANC_mother$primigrav=='Non'&BF_ANC_mother$prev_pregs==0,]

#Frequency plots
ggplot(BF_ANC_mother,aes(x=age))+
  geom_histogram()
ggplot(BF_ANC_mother,aes(x=age))+
  geom_bar()
ggplot(BF_ANC_mother,aes(x=month))+
  geom_bar()
ggplot(BF_ANC_mother,aes(x=primigrav))+
  geom_bar()
ggplot(BF_ANC_mother,aes(x=grav))+
  geom_bar()
ggplot(BF_ANC_mother,aes(x=month,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(BF_ANC_mother,aes(x=age,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(BF_ANC_mother,aes(x=age,fill=as.factor(rdt)))+
  geom_bar(position='stack')
ggplot(BF_ANC_mother,aes(x=grav,fill=as.factor(rdt)))+
  geom_bar(position='fill')
ggplot(BF_ANC_mother,aes(x=grav,fill=as.factor(rdt)))+
  geom_bar(position='stack')


#Missing and values for character variables
table(BF_ANC_mother$HF)
table(is.na(BF_ANC_mother$HF)) #0 missing
table(BF_ANC_mother$District)
table(is.na(BF_ANC_mother$District)) #0 missing
table(BF_ANC_mother$Village)
table(is.na(BF_ANC_mother$Village)) #0 missing; many different spelling/spacings

#Group data and calculate prevalence and CI by month and gravidity
BF_ANC_mother_grouped <- BF_ANC_mother %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         age_cat=cut(age,breaks=c(12,15,20,25,30,35,40,Inf),labels=c('12-15 years','16-20 years','21-25 years','26-30 years','31-35 years','36-40 years','Over 40 years')),
         rdt=as.numeric(ifelse(rdt=='Positif',1,0))
  )%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age!=88) %>%
  filter(!is.na(age_cat)) %>%
  filter(!is.na(grav))%>%
  group_by(month,grav_cat,age_cat,.drop=FALSE)%>%
  summarise(positive=sum(rdt),total=n())
BF_ANC_mother_grouped <- addCIs(BF_ANC_mother_grouped,BF_ANC_mother_grouped$positive,BF_ANC_mother_grouped$total)

windows(height=40,width=50)
ggplot(data=BF_ANC_mother_grouped,aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~grav_cat,dir="v")+theme(legend.position="none")

#Group by grav
BF_ANC_mother_grouped_grav <- BF_ANC_mother %>%
  mutate(grav_cat=cut(grav,breaks=c(0,1,3,Inf),labels=c("Gravidities 1","Gravidities 2-3","Gravidities 4+")),
         rdt=as.numeric(ifelse(rdt=='Positif',1,0)))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age!=88) %>%
  filter(!is.na(grav))%>%
  group_by(month,grav_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_ANC_mother_grouped_grav <- addCIs(BF_ANC_mother_grouped_grav,BF_ANC_mother_grouped_grav$positive,BF_ANC_mother_grouped_grav$total)

#Grav facetted
ggplot(data=BF_ANC_mother_grouped_grav,aes(x=month,y=mean))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~grav_cat,dir="v")+theme(legend.position="none")

#Grav colored
ggplot(data=BF_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")

#Group by age
BF_ANC_mother_grouped_age <- BF_ANC_mother %>%
  mutate(age_cat=cut(age,breaks=c(12,15,20,25,30,35,40,Inf),labels=c('12-15 years','16-20 years','21-25 years','26-30 years','31-35 years','36-40 years','Over 40 years')),
         rdt=as.numeric(ifelse(rdt=='Positif',1,0)))%>%
  filter(!is.na(rdt)) %>%
  filter(!is.na(age)) %>%
  filter(age!=88) %>%
  filter(!is.na(age_cat)) %>%
  filter(!is.na(grav))%>%
  group_by(month,age_cat,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_ANC_mother_grouped_age <- addCIs(BF_ANC_mother_grouped_age,BF_ANC_mother_grouped_age$positive,BF_ANC_mother_grouped_age$total)

#Age facetted
ggplot(data=BF_ANC_mother_grouped_age,aes(x=month,y=mean))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~age_cat,dir="v")+theme(legend.position="none")

#Age colored
ggplot(data=BF_ANC_mother_grouped_age[BF_ANC_mother_grouped_age$age_cat!='12-15 years'&BF_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")

#########################################
##Combining Countries into single plot###
#########################################
NG_ANC_mother_grouped_grav$country <- 'Nigeria'
MZ_ANC_mother_grouped_grav$country <- 'Mozambique'
BF_ANC_mother_grouped_grav$country <- 'Burkina Faso'

All_ANC_mother_grouped_grav <- rbind(NG_ANC_mother_grouped_grav,MZ_ANC_mother_grouped_grav,BF_ANC_mother_grouped_grav)

#Country by facet
ggplot(data=All_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  facet_grid(rows=vars(country))

#Country by color, grav by facet
ggplot(data=All_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(country)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(country)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  facet_grid(rows=vars(as.factor(grav_cat)))

#Country by facet, grav by facet
ggplot(data=All_ANC_mother_grouped_grav,aes(x=month,y=mean,))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  facet_grid(rows = vars(country),cols = vars(as.factor(grav_cat)))

##By age
NG_ANC_mother_grouped_age$country <- 'Nigeria'
MZ_ANC_mother_grouped_age$country <- 'Mozambique'
BF_ANC_mother_grouped_age$country <- 'Burkina Faso'

All_ANC_mother_grouped_age <- rbind(NG_ANC_mother_grouped_age,MZ_ANC_mother_grouped_age,BF_ANC_mother_grouped_age)

#Country by facet
ggplot(data=All_ANC_mother_grouped_age[All_ANC_mother_grouped_age$age_cat!='12-15 years'&All_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  facet_grid(rows=vars(country))

#Country by color, age by facet
ggplot(data=All_ANC_mother_grouped_age[All_ANC_mother_grouped_age$age_cat!='12-15 years'&All_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(country)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(country)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  facet_grid(rows=vars(as.factor(age_cat)))

#Country by facet, age by facet
ggplot(data=All_ANC_mother_grouped_age[All_ANC_mother_grouped_age$age_cat!='12-15 years'&All_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  facet_grid(rows = vars(country),cols = vars(as.factor(age_cat)))


###############################################
#######Community Prevalence####################
###############################################
#Read in Nigeria data and rename most used variables
NG_CS_child_2020 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Nigeria/Cross-sectional survey/data_nnp_survey_child_nigeria_2020.xlsx')%>%
  rename(rdt = q85c_result,
         mal_symp = q86a_symptoms) %>%
  mutate(age = ifelse(age==0,age_months/12,age))

NG_CS_hh_2020 <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl300522/Nigeria/Cross-sectional survey/data_nnp_survey_hh_nigeria_2020.xlsx')%>%
  rename(individual_id = hh_id,
         date = start) %>%
  mutate(date = as.Date(date))

table(is.na(NG_CS_hh_2020$village)) #1 missing
table(NG_CS_hh_2020$village) #Lots of different spellings/capitalizations
table(is.na(NG_CS_hh_2020$cluster_n))


#Merge child to hh to get date of survey
NG_CS_2020 <- merge(NG_CS_child_2020,NG_CS_hh_2020,by='submission_id',all.x=TRUE,suffixes=c('.child','.hh'))
table(is.na(NG_CS_2020$cluster_n))

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
table(is.na(NG_CS_2020$cluster_n))

#Combine 2 years of data
NG_CS_all <- plyr::rbind.fill(NG_CS_2020,NG_CS_2021)
NG_CS_all$district_code <- as.integer(substr(NG_CS_all$cluster_n,1,1))
NG_CS_all$district_n <- unlist(lapply(NG_CS_all$district_code, FUN = function(x) switch(x,'Ejigbo','Ife North','Asa','Moro')))

#Basic missing investigation
table(is.na(NG_CS_all$rdt)) #3 missing
table(NG_CS_all$rdt) #2 'invalid'
table(is.na(NG_CS_all$date)) #0 missing
table(NG_CS_all$date)

#Frequency plots
ggplot(NG_CS_all,aes(x=date))+
  geom_bar()
ggplot(NG_CS_all,aes(x=rdt))+
  geom_bar()

#Group data and calculate prevalence and CI by month
NG_CS_all_grouped <- NG_CS_all %>%
  mutate(month=as.yearmon(date),
    rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  group_by(month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
NG_CS_all_grouped <- addCIs(NG_CS_all_grouped,NG_CS_all_grouped$positive,NG_CS_all_grouped$total)

windows(height=40,width=50)
ggplot(data=NG_CS_all_grouped,aes(x=month,y=mean))+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="none")

#Add community prevalence to ANC figures
ggplot(data=NG_ANC_mother_grouped)+
  geom_point(aes(x=month,y=mean,col=as.factor(age_cat)),size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  facet_wrap(~grav_cat,dir="v")+theme(legend.position="none")+
  geom_point(data=NG_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=NG_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")
#Grav colored
ggplot(data=NG_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  geom_point(data=NG_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=NG_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")
#Age colored
ggplot(data=NG_ANC_mother_grouped_age[NG_ANC_mother_grouped_age$age_cat!='12-15 years'&NG_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  geom_point(data=NG_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=NG_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")


#Read in Mozambique data and rename most used variables
MZ_CS_rdt_base <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Mozambique/Cross-sectional survey/Moz Baseline CSS Datasets 2020.xlsx',
                            xl.sheet = 'RDT')%>%
  mutate(date = as.Date(date))
MZ_CS_rdt_mid <- xl.read.file('C:/Users/jthicks/OneDrive - Imperial College London/Imperial_ResearchAssociate/PregnancyModel/PATH/Imperial College (ANC)_data_dl160822/Mozambique/Cross-sectional survey/Moz Midline CSS Datasets 2021.xlsx',
                          xl.sheet = 'RDT')%>%
  mutate(date = as.Date(date,"%m/%d/%Y"))
MZ_CS_rdt <- rbind(MZ_CS_rdt_base,MZ_CS_rdt_mid)


#Basic missing investigation
table(is.na(MZ_CS_rdt$rdt)) #24 missing
table(MZ_CS_rdt$rdt) #5 'Not Done'
table(is.na(MZ_CS_rdt$date)) #0 missing
table(MZ_CS_rdt$date)

#Frequency plots
ggplot(MZ_CS_rdt,aes(x=date))+
  geom_bar()
ggplot(MZ_CS_rdt,aes(x=rdt))+
  geom_bar()

#Group data and calculate prevalence and CI by month
MZ_CS_all_grouped <- MZ_CS_rdt %>%
  mutate(month=as.yearmon(date),
         rdt=as.numeric(ifelse(rdt=='Positive',1,ifelse(rdt=='Negative',0,NA)))
  )%>%
  filter(!is.na(rdt)) %>%
  group_by(month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
MZ_CS_all_grouped <- addCIs(MZ_CS_all_grouped,MZ_CS_all_grouped$positive,MZ_CS_all_grouped$total)

#Add community prevalence to ANC figures
#Grav colored
ggplot(data=MZ_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  geom_point(data=MZ_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=MZ_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")

#Age colored
ggplot(data=MZ_ANC_mother_grouped_age[MZ_ANC_mother_grouped_age$age_cat!='12-15 years'&MZ_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  geom_point(data=MZ_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=MZ_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")

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

#Basic missing investigation
table(is.na(BF_CS_child_all$rdt)) #0 missing
table(BF_CS_child_all$rdt) #Spacing issues with 'Positif

#Frequency plots
ggplot(BF_CS_child_all,aes(x=month))+
  geom_bar()
ggplot(BF_CS_child_all,aes(x=rdt))+
  geom_bar()

#Group data and calculate prevalence and CI by month
BF_CS_all_grouped <- BF_CS_child_all %>%
  mutate(
         rdt=as.numeric(ifelse(rdt=='Négatif',0,1))
  )%>%
  filter(!is.na(rdt)) %>%
  group_by(month,.drop=FALSE)%>%
  dplyr::summarise(positive=sum(rdt),total=n())
BF_CS_all_grouped <- addCIs(BF_CS_all_grouped,BF_CS_all_grouped$positive,BF_CS_all_grouped$total)


#Grav colored
ggplot(data=BF_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  geom_point(data=BF_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=BF_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")

#Age colored
ggplot(data=BF_ANC_mother_grouped_age[BF_ANC_mother_grouped_age$age_cat!='12-15 years'&BF_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  geom_point(data=BF_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=BF_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")

#All three countries
NG_CS_all_grouped$country <- 'Nigeria'
MZ_CS_all_grouped$country <- 'Mozambique'
BF_CS_all_grouped$country <- 'Burkina Faso'

All_CS_all_grouped <- rbind(NG_CS_all_grouped,MZ_CS_all_grouped,BF_CS_all_grouped)

#Country by facet
ggplot(data=All_ANC_mother_grouped_grav,aes(x=month,y=mean,col=as.factor(grav_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(grav_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  geom_point(data=All_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=All_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")+
  facet_grid(rows=vars(country))

#age
ggplot(data=All_ANC_mother_grouped_age[All_ANC_mother_grouped_age$age_cat!='12-15 years'&All_ANC_mother_grouped_age$age_cat!='Over 40 years',],aes(x=month,y=mean,col=as.factor(age_cat)))+
  geom_line()+
  geom_point(size=3)+
  theme_minimal()+
  geom_errorbar(aes(x=month,y=mean,ymax=upper,ymin=lower,width=0,col=as.factor(age_cat)))+
  ylab("Prevalence")+xlab("Month")+
  scale_color_viridis(name="Data source",option = "D",discrete=T,begin=0.2,end=0.9)+
  theme(legend.position="bottom")+
  geom_point(data=All_CS_all_grouped,aes(x=month,y=mean),size=3,color="#D55E00")+
  geom_errorbar(data=All_CS_all_grouped,aes(x=month,y=mean,ymax=upper,ymin=lower,width=0),color="#D55E00")+
  facet_grid(rows=vars(country))


########################################
####Prevalence by site within country###
########################################
#NIGERIA#
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
NG_ANC_mother_grouped_sitegrav <- readRDS('NG_ANC_mother_grouped_sitegrav.rds')
#Create data tables to annotate figure
NG_iv <- data.frame(site = c('Ejigbo','Asa','Moro','Ife North'),
                    month = as.yearmon(rep(as.Date('2020-11-01'),4)))
NG_labs <- c('Ejigbo - Standard','Asa - IG2','Moro - RG','Ife North - PBO')
names(NG_labs) <- c('Ejigbo','Asa','Moro','Ife North')
#Grav colored, facetted by site
windows(height=40,width=50)
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
MZ_CS_all_grouped_site <- readRDS('MZ_CS_all_grouped_site.rds')

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
MZ_ANC_mother_grouped_sitegrav <- readRDS('MZ_ANC_mother_grouped_sitegrav.rds')

table(MZ_CS_all_grouped_site$site)
table(MZ_ANC_mother_grouped_sitegrav$site)

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
BF_CS_all_grouped_site <- readRDS('BF_CS_all_grouped_site.rds')

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
BF_ANC_mother_grouped_sitegrav <- readRDS('BF_ANC_mother_grouped_sitegrav.rds')

table(BF_CS_all_grouped$site)
table(BF_ANC_mother_grouped_sitegrav$site)

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



