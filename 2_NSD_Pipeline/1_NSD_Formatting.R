#set home dir of pipeline
home<-"/Users/kayleigh.chalkowski/Library/CloudStorage/OneDrive-USDA/Projects/NIFA_Analyses/NIFA_Removals_Mvmt/Pipeline"

# Purpose ----------------------------------------------------------------------

# In/Out -----------------------------------------------------------------------

# Setup ------------------------------------------------------------------------

#load libraries

#Set dirs
input=file.path(home,"1_Data","Input",fsep=.Platform$file.sep)
objdir=file.path(home,"1_Data","Objects",fsep=.Platform$file.sep)

#Read geoloc data
geo.tox<-readRDS(file.path(objdir,"geotox.rds"))
geo.aer<-readRDS(file.path(objdir,"geoaer.rds"))
geo.trap<-readRDS(file.path(objdir,"geotrap.rds"))

#source needed functions
func.list=list.files(file.path(home,"2_Scripts","Functions",fsep=.Platform$file.sep),full.names=TRUE)
for(f in 1:length(func.list)){
  source(func.list[f])
}

# * Trim incomplete weeks ------------------------------------------------------

#Aerial
aerchk=geo.aer %>% 
  group_by(animalid,removal.period.akdecalc,week) %>% 
  dplyr::summarise(n_distinct(jDate)) %>%
  as.data.frame()
#remove week 6 and 12, only two days
geo.aer=geo.aer[geo.aer$week!=12,]
geo.aer=geo.aer[geo.aer$week!=6,]

#Trap
trapchk=geo.trap %>% 
  group_by(animalid,removal.period.akdecalc,week) %>% 
  dplyr::summarise(n_distinct(jDate)) %>%
  as.data.frame()
#remove week 25, only 1 days
#geo.trapd2=geo.trapd2[geo.trapd2$week!=25,]
#wks 9 and 26-- all only 1 day
geo.trap=geo.trap[geo.trap$week!=9,]
geo.trap=geo.trap[geo.trap$week!=26,]

#each below combo has <6 days
geo.trap=geo.trap[!(geo.trap$week==18&geo.trap$animalid=="48458_A4_A4"),]
geo.trap=geo.trap[!(geo.trap$week==18&geo.trap$animalid=="48479_T2_T2"),]
geo.trap=geo.trap[!(geo.trap$week==20&geo.trap$animalid=="85411_C3_C3"),]
geo.trap=geo.trap[!(geo.trap$week==24&geo.trap$animalid=="48460_B3_B3"),]
geo.trap=geo.trap[!(geo.trap$week==21&geo.trap$animalid=="85454_1W_1W"),]
geo.trap=geo.trap[!(geo.trap$week==22&geo.trap$animalid=="48468_J6_J6"),]

#Tox
toxchk=geo.tox %>% 
  group_by(animalid,removal.period.akdecalc,week) %>% 
  dplyr::summarise(n_distinct(jDate)) %>%
  as.data.frame()

#very short weeks for 6, 15
geo.tox=geo.tox[geo.tox$week!=6,]
geo.tox=geo.tox[geo.tox$week!=15,]

#short weeks for certain IDs for weeks 14/10/9
geo.tox=geo.tox[!(geo.tox$week==14&geo.tox$animalid=="86063_C6_C6"),]
geo.tox=geo.tox[!(geo.tox$week==10&geo.tox$animalid=="86070_H2_H2"),]
geo.tox=geo.tox[!(geo.tox$week==9&geo.tox$animalid=="86071_Y2_Y2"),]
geo.tox=geo.tox[!(geo.tox$week==9&geo.tox$animalid=="86079_B1_B1"),]
geo.tox=geo.tox[!(geo.tox$week==9&geo.tox$animalid=="85734_T6_T6"),]
geo.tox=geo.tox[!(geo.tox$week==9&geo.tox$animalid=="85735_7K_7K"),]

#View summary
nrow(geo.tox) #320761
nrow(geo.trap) #588735
nrow(geo.aer) #225757

# Calculate NSD ----------------------------------------------------------------

#Make functions to calculate NSD
#pigex is pig/period subset from geo.rem dfs
do.NSDcalcs<-function(pigex){
  pigex=pigex[order(pigex$datetime),]
  
  #lagged loc
  pigex$Xstart=pigex$X[1]
  pigex$Ystart=pigex$Y[1]
  #pigex$datetime2=dplyr::lead(pigex$datetime)
  
  pigex$displacement=sqrt((pigex$X-pigex$Xstart)^2+(pigex$Y-pigex$Ystart)^2)
  pigex$NSD=(pigex$displacement)^2
  return(pigex)
}

#loop through each pig, get NSD, rbind together, output geo.rem with NSD
do.NSDcalcs.georem=function(geo.rem){
  pigID=unique(geo.rem$animalid)
  
  for(i in 1:length(pigID)){
    
    print(paste0("calculating NSD for ",pigID[i]," (pig ",i," of ",length(pigID),")"))
    geo.pig=geo.rem[geo.rem$animalid==pigID[i],]
    geo.pig[order(geo.pig$datetime),]
    
    geo.pig2=do.NSDcalcs(geo.pig)
    
    if(i==1){
      geo.rem.out=geo.pig2
    } else{
      geo.rem.out=rbind(geo.rem.out,geo.pig2)
    }
  }
  
  return(geo.rem.out)
  
}

geo.aerd=do.NSDcalcs.georem(geo.aer)
geo.trapd=do.NSDcalcs.georem(geo.trap)
geo.toxd=do.NSDcalcs.georem(geo.tox)

# Summarize NSD by weekly median -----------------------------------------------

#summarize by week
geo.aerd.wk=geo.aerd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()
geo.trapd.wk=geo.trapd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()
geo.toxd.wk=geo.toxd2 %>% group_by(animalid, Removal.Type, removal.period.akdecalc,sex, week) %>% dplyr::summarise(mNSD=median(NSD),mX=mean(X),mY=mean(Y)) %>% as.data.frame()

# Write out data ---------------------------------------------------------------

saveRDS(geo.aerd.wk,file.path(objdir,"NSDgeoaer.rds",fsep=.Platform$file.sep))
saveRDS(geo.trapd.wk,file.path(objdir,"NSDgeotrap.rds",fsep=.Platform$file.sep))
saveRDS(geo.toxd.wk,file.path(objdir,"NSDgeotox.rds",fsep=.Platform$file.sep))

