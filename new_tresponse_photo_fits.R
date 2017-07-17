
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM"
#path="//ad.uws.edu.au/dfshare/HomesHWK$/90925395/My Documents/Repos/PhotoM"
path<-getwd()

hfe<-read.csv(paste(path,"/Data/hfe_An_Tdata_processed.csv",sep=""))

hfe.mix<-split(hfe,paste(hfe$Species,hfe$Month))


hfe.mix[c(3:5)]<-NULL # E.crebra/January & E.crebra/August E.dunii/August (cannot fit either mix models or standard non-linear; data quality??)

hfe.topts<-get_topts(lapply(hfe.mix,function(x)fit.nlme(x,yvar="Photo",random="Tree")))


#hfe.topts$Names <- names(hfe.mix)

hfe.Names <- strsplit(names(hfe.mix)," ")

hfe.topts$Species<-get_names(hfe.mix)[,1]
hfe.topts$Month<-get_names(hfe.mix)[,2]
#----------------------------
#-------------------------------

#get Topts for few species from ACi data. 

hfe_rest<-read.csv(paste(path,"/Data/hfe_ACidata_processed.csv",sep=""))
hfe_rest.a<-subset(hfe_rest,abs(hfe_rest$CO2R-400)<10)

ecre.feb<-data.frame(do.call(rbind,list(fitquad(subset(hfe_rest.a,hfe_rest.a$Species=="Ecla" & hfe_rest.a$Season=="February")))))
ecre.aug<-data.frame(do.call(rbind,list(fitquad(subset(hfe_rest.a,hfe_rest.a$Species=="Ecla" & hfe_rest.a$Season=="August")))))
edun.aug<-data.frame(do.call(rbind,list(fitquad(subset(hfe_rest.a,hfe_rest.a$Species=="Edun" & hfe_rest.a$Season=="August")))))
#emel.feb<-data.frame(do.call(rbind,list(fitquad(subset(hfe_rest.a,hfe_rest.a$Species=="Emel" & hfe_rest.a$Season=="February")))))

hfe_topt_rest<-rbind(ecre.feb,ecre.aug,edun.aug)
hfe_topt_rest$Species<-c("E.crebra","E.crebra","E.dunnii")
hfe_topt_rest$Month<-c("January","August","August")

hfe.topts<-rbind(hfe.topts,hfe_topt_rest)

#-------------------------------

#get mean gs,VPD and Ci
#hfe.sum<-summaryBy(Cond+VpdL+Ci~Species+ Month,keep.names=TRUE,data=hfe,FUN=c(mean,std.error))

hfe.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci+Tleaf~Species+Month+Tleafnew,data=hfe,FUN=c(mean,std.error),keep.names=TRUE)
hfe.sum<-subset(hfe.sum,hfe.sum$Tleafnew==4)

hfe.topts<-merge(hfe.topts,hfe.sum,by=c("Species","Month"),all=T)

#to change January to February
#to facilitate merging with other dataframes

hfe.topts$Month<-ifelse(hfe.topts$Month=="January","February","August") #just change month to fit for vcmax and jmax dataframe

#add latitude of seed source
hfe.topts$SSloc<-NA
hfe.topts$SSloc[which(hfe.topts$Species=="E.cladocalyx")]<--33.03
hfe.topts$SSloc[which(hfe.topts$Species=="E.crebra")]<--33.61 
hfe.topts$SSloc[which(hfe.topts$Species=="E.dunnii")]<--28.41
hfe.topts$SSloc[which(hfe.topts$Species=="E.melliodora")]<--35.11
hfe.topts$SSloc[which(hfe.topts$Species=="E.saligna")]<--30.37
hfe.topts$SSloc[which(hfe.topts$Species=="E.tereticornis")]<--30.04


hfe.topts$DataSet<-"HFE_CG"


#get g1 and add to the dataframe
hfe<-read.csv(paste(path,"/Data/hfe_An_Tdata_processed.csv",sep=""))
hfe.mix<-split(hfe,paste(hfe$Species,hfe$Month))

hfe_g1<-get_topts(lapply(hfe.mix,FUN=fitStom))
hfe_g1$Species<-get_names(hfe.mix)[,1]
hfe_g1$Month<-get_names(hfe.mix)[,2]
hfe_g1$Month<-ifelse(hfe_g1$Month=="January","February","August") #just change month to fit for vcmax and jmax dataframe

hfe.topts<-merge(hfe.topts,hfe_g1,by=c("Species","Month"))


#Get complete species names

hfe.topts$Species[which(hfe.topts$Species=="E.cladocalyx")]<-"Eucalyptus cladocalyx"
hfe.topts$Species[which(hfe.topts$Species=="E.crebra")]<-"Eucalyptus crebra"
hfe.topts$Species[which(hfe.topts$Species=="E.dunnii")]<-"Eucalyptus dunnii"
hfe.topts$Species[which(hfe.topts$Species=="E.melliodora")]<-"Eucalyptus melliodora" 
hfe.topts$Species[which(hfe.topts$Species=="E.saligna")]<-"Eucalyptus saligna"
hfe.topts$Species[which(hfe.topts$Species=="E.tereticornis")]<-"Eucalyptus tereticornis"

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#wtc3.c<-subset(wtc3.b,wtc3.b$Ci>150)
#wtc3.topts<-get_topts(lapply(wtc3.mix,function(x)fit.nlme(x,yvar="Photo",random="Chamber")))
#wtc3.topts$Temp_Treatment<-get_names(wtc3.mix)[,1]
#wtc3.topts$Month<-get_names(wtc3.mix)[,2]

#E.tereticornis from WTC 3

#wtc3<-read.csv("C:/Dushan/Repos/ACI_ALL/WTC3/WTC3_ACidata_processed.csv")
wtc3<-getwtc3(path=paste0(path,"/Data"))
#wtc3<-getwtc3(path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/Data")

wtc3.a<-wtc3[firstobs(~Plotsite,data=wtc3),] #to get ambient CO2 levels from ACi curves
wtc3.b<-subset(wtc3.a,wtc3.a$CO2S>350) #to remove low CO2R values (first observation < 300 ppm)

wtc3.c<-wtc3.b[!(wtc3.b$MeanTemp=="14-Apr" & wtc3.b$Chamber==5),] #remove 5 datapoints from April ambient (chamber 5)
#very low photosynthetic rates and conductance compared to other chambers
wtc3.mix<-split(wtc3.c,paste(wtc3.c$Treatment,wtc3.c$Season))
wtc3.mix[2]<-NULL #January, Ambient ~ No clear peak
#wtc3.mix[6]<-NULL #May, Elevated cannot fit mix models (data only from one chamber)

#fit non linear mix model
wtc3.topts<-lapply(wtc3.mix,function(x)fit.nlme(x,yvar="Photo",random="Chamber"))


wtc3.topts<-data.frame(do.call(rbind,list(wtc3.topts[[1]],wtc3.topts[[2]],wtc3.topts[[3]],wtc3.topts[[4]],wtc3.topts[[5]]
)))
#wtc3.topts$Names <- names(wtc3.mix)


wtc3.topts$Temp_Treatment<-c(Name[[1]][[1]],Name[[2]][[1]],Name[[3]][[1]],Name[[4]][[1]],
                         Name[[5]][[1]])
wtc3.topts$Temp_Treatment<-as.factor(wtc3.topts$Temp_Treatment)
wtc3.topts$Month<-c(Name[[1]][[2]],Name[[2]][[2]],Name[[3]][[2]],Name[[4]][[2]],
                    Name[[5]][[2]])
wtc3.topts$Month<-as.factor(wtc3.topts$Month)


#get mean gs,VPD and Ci
e.tret.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci+Tleaf~Treatment+TargetTemp+Season,keep.names=TRUE,data=wtc3.b,FUN=c(mean,std.error))
e.tret.sum<-subset(e.tret.sum,e.tret.sum$TargetTemp==25)

#e.tret.sum<-summaryBy(Cond+VpdL+Ci~Treatment+Season,keep.names=TRUE,data=wtc3.b,FUN=c(mean,std.error))
colnames(e.tret.sum)[1]<-"Temp_Treatment"
colnames(e.tret.sum)[3]<-"Month"

wtc3.topts<-merge(wtc3.topts,e.tret.sum,by=c("Temp_Treatment","Month"),all=TRUE)
#colnames(wtc3.topts)[15]<-"Ci_ratio"
#colnames(wtc3.topts)[21]<-"Ci_ratio_se"
#add latitude of seed source
wtc3.topts$SSloc<--33.62


wtc3.topts$Species<-"Eucalyptus tereticornis"
wtc3.topts$DataSet<-"WTC3"

#get g1

wtc3_g1<-get_topts(lapply(wtc3.mix,FUN=fitStom))
wtc3_g1$Temp_Treatment<-get_names(wtc3.mix)[,1]
wtc3_g1$Month<-get_names(wtc3.mix)[,2]

wtc3.topts<-merge(wtc3.topts,wtc3_g1,by=c("Month","Temp_Treatment"))
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#E.globulus from WTC 2

#wtc2<-read.csv("C:/Dushan/Repos/ACI_ALL/WTC2/WTC2_ACidata_processednew.csv")

wtc2<-getwtc2(path=paste0(path,"/Data"))

#get ambient CO2 treatment
wtc2<-subset(wtc2,wtc2$C_treat=="0C")

#get ambient CO2 levels form ACi curvs
wtc2.1<-wtc2[firstobs(~Identity,data=wtc2),]

wtc2.a<-subset(wtc2.1,wtc2.1$CO2R>60)#to remove low CO2R values

#e.glob.a<-subset(e.glob,abs((e.glob$CO2S)-400)<22)

#fit non linear mix model

wtc2.mix<-split(wtc2.a,paste(wtc2.a$Temp.treat,wtc2.a$Season))

wtc2.mix[1]<-NULL #April, measured only at 25C
wtc2.mix[4]<-NULL #April, measured only at 25C

wtc2.topts.1<-lapply(wtc2.mix[1],FUN=fitquad) #model without random effects. Mix model cannot fir for this (August/Ambient)
wtc2.topts.2<-lapply(wtc2.mix[6],FUN=fitquad) #model without random effects. Mix model cannot fir for this (Elevated/February)


#to fit mix models for the rest 

wtc2.mix[1]<-NULL
wtc2.mix[5]<-NULL

wtc2.topts.3<-lapply(wtc2.mix,function(x)fit.nlme(x,yvar="Photo",random="Chamber"))


#to pull out paramters
wtc2.topts<-data.frame(do.call(rbind,list(wtc2.topts.3[[1]],wtc2.topts.3[[2]],wtc2.topts.3[[3]],wtc2.topts.3[[4]]
                                          ,wtc2.topts.1[[1]],wtc2.topts.2[[1]])))

#to add temperature treatment and season to the dataframe

wtc2.topts$Temp_Treatment<-as.factor( c("Ambient","Ambient","Elevated","Elevated","Ambient","Elevated"))
wtc2.topts$Month<-as.factor(wtc2.topts$Month<-c("December","February","September","December","September","February"))
wtc2.topts$Species<-"Eucalyptus globulus"  


#get mean gs,VPD and Ci

wtc2.f<-subset(wtc2.a,wtc2.a$Season !="April")
#e.glob.sum<-summaryBy(Cond+VpdL+Ci~Temp.treat+Season,keep.names=TRUE,data=wtc2.f,FUN=c(mean,std.error))
e.glob.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci+Tleaf~Temp.treat+TargetTemp+Season,keep.names=TRUE,data=wtc2.f,FUN=c(mean,std.error))
e.glob.sum<-subset(e.glob.sum,e.glob.sum$TargetTemp==25)

colnames(e.glob.sum)[1]<-"Temp_Treatment"
colnames(e.glob.sum)[3]<-"Month"
e.glob.sum$Temp_Treatment<-c("Ambient","Ambient","Ambient","Elevated","Elevated","Elevated")

e.glob.sum$Month<-c("September","December","February","September","December","February")

wtc2.topt<-merge(wtc2.topts,e.glob.sum,by=c("Temp_Treatment","Month"))


wtc2.topt$DataSet<-"WTC2"

#get g1
wtc2.mix<-split(wtc2.a,paste(wtc2.a$Temp.treat,wtc2.a$Season))
wtc2.mix[1]<-NULL #April, measured only at 25C
wtc2.mix[4]<-NULL #April, measured only at 25C

wtc2_g1<-get_topts(lapply(wtc2.mix,FUN=fitStom))
wtc2_g1$Temp_Treatment<-capFirst(get_names(wtc2.mix)[,1])
wtc2_g1$Month<-c("September","December","February","September","December","February")

wtc2.topt<-merge(wtc2.topt,wtc2_g1,by=c("Month","Temp_Treatment"))


#wtc2.topt<-subset(wtc2.topt,wtc2.topt$Month!="February")
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------

#WTC1

wtc1<-read.csv(paste(path,"/Data/wtc1_Aci_processed.csv",sep=""))
wtc1$TargetTemp<-cut(wtc1$Tleaf,breaks=c(16,22,28,34,37,40),
                     labels=FALSE)

#to get well watered and ambient CO2 treatments (chamber 1, 3 & 11)
wtc1.a<-subset(wtc1,wtc1$CO2_Treat=="Ambient" & wtc1$Water_treat=="wet")
wtc1.b<-subset(wtc1.a,abs((wtc1.a$CO2S)-400)<20 )
#wtc1.b<-subset(wtc1.b,wtc1.b$Tleaf>21 )

wtc1.mix<-split(wtc1.b,wtc1.b$Season)

wtc1.topts<-lapply(wtc1.mix[1],function(x)fit.nlme(x,yvar="Photo",random="Chamber"))
wtc1.topts.1<-lapply(wtc1.mix[2],FUN=fitquad) #for November, mix models cannot be fitted 

wtc1.topts<-data.frame(do.call(rbind,list(wtc1.topts[[1]],wtc1.topts.1[[1]])))

Treatment<-strsplit(names(wtc1.mix)," ")
wtc1.topts$Month<-c(Treatment[[1]][[1]],Treatment[[2]][[1]])
wtc1.topts$Species<-"Eucalyptus saligna"




#get mean gs,VPD and Ci
#e.sal.sum<-summaryBy(Cond+VpdL+Ci~Season,keep.names=TRUE,data=wtc1.b,FUN=c(mean,std.error))

e.sal.mm<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci+Tleaf~Season+TargetTemp,FUN=c(mean,std.error),data=wtc1.b)
e.sal.mm<-subset(e.sal.mm,e.sal.mm$TargetTemp==2)

colnames(e.sal.mm)[1]<-"Month"

wtc1.topts<-merge(wtc1.topts,e.sal.mm,by=c("Month"))

#add latitude of seed source
wtc1.topts$SSloc<--30.43

wtc1.topts$DataSet<-"WTC1"


#get g1

wtc1_g1<-get_topts(lapply(wtc1.mix,FUN=fitStom))
wtc1_g1$Month<-names(wtc1.mix)

wtc1.topts<-merge(wtc1.topts,wtc1_g1,by=c("Month"))
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------


#E.saligna and E.sideroxylon (S30)

e.ss<-read.csv(paste(path,"/Data/GHS30_Euc2-SxTxCO2xW_GEtempresponse_20090105-20090109_L1.csv",sep=""))
e.ss<-subset(e.ss,e.ss$Cond>0 & e.ss$Ci>0) #remove negative gs and Ci
e.ss<-subset(e.ss,e.ss$Cond>0)

e.ss$Temp_T<-NA
e.ss$Temp_T[which(e.ss$Temp=="Amb")]<-"Ambient"
e.ss$Temp_T[which(e.ss$Temp=="Elv")]<-"Elevated"

e.ss$Spp<-NA
e.ss$Spp[which(e.ss$Species=="E. saligna")]<-"E.saligna"
e.ss$Spp[which(e.ss$Species=="E. sideroxylon")]<-"E.sideroxylon"

#get ambient CO2 treatment
e.ss<-subset(e.ss,e.ss$CO2==400)


#to get average photosynthesis for each seedling
#there are five spot measurements per seedling per Tleaf level

e.ss$Number <- c(1)
Number <- c()
count <- 1  
for (i in 2:length(e.ss$TBlk)){
  ifelse(e.ss$TBlk[i-1] - e.ss$TBlk[i]<1 && e.ss$Potnum[i-1]==e.ss$Potnum[i],count <- count,count <- count + 1)
  Number[i] <- count 
}
e.ss$Number[2:length(e.ss$TBlk)] <- na.omit(Number)


e.ss.1<-summaryBy(.~Number,data=e.ss,FUN=mean,keep.names=TRUE)

#to add factor variables to the dataframe
e.ss.2<-e.ss[,c(23,24,25)] 
e.ss.3 <- merge(e.ss.1,e.ss.2,by="Number")
e.ss.4 <- e.ss.3[!duplicated(e.ss.3[,c("Number")]),]
e.ss.4$TargetTemp<-cut(e.ss.4$TBlk,breaks=c(15,19,21,26,32,36,44),
                       labels=FALSE)

#fit mix models
e.ss.mix<-split(e.ss.4,paste(e.ss.4$Spp,e.ss.4$Temp_T))

e.ss.topts<-lapply(e.ss.mix,function(x)fit.nlme(x,yvar="Photo",random="Potnum"))

e.ss.topts<-data.frame(do.call(rbind,list(e.ss.topts[[1]],e.ss.topts[[2]],e.ss.topts[[3]],e.ss.topts[[4]]
)))


#nam<-names(e.ss.mix)
nam<-strsplit(names(e.ss.mix)," ")

e.ss.topts$Species<-c(nam[[1]][[1]],nam[[2]][[1]],nam[[3]][[1]],nam[[4]][[1]])
e.ss.topts$Temp_Treatment<-c(nam[[1]][[2]],nam[[2]][[2]],nam[[3]][[2]],nam[[4]][[2]])


#get mean gs,VPD and Ci
#e.ss.sum<-summaryBy(Cond+VpdL+Ci~Spp+Temp_T,keep.names=TRUE,data=e.ss.4,FUN=c(mean,std.error))
e.ss.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci+Tleaf~Spp+Temp_T+TargetTemp,data=e.ss.4,FUN=c(mean,std.error),keep.names=TRUE)
e.ss.sum<-subset(e.ss.sum,e.ss.sum$TargetTemp==3)

colnames(e.ss.sum)[1]<-"Species"
colnames(e.ss.sum)[2]<-"Temp_Treatment"
e.ss.topts<-merge(e.ss.topts,e.ss.sum,by=c("Species","Temp_Treatment"))
e.ss.topts$Month<-"January" 
#to add longteram Temperature at seed source

#add latitude of seed source
e.ss.topts$SSloc<-NA
e.ss.topts$SSloc[which(e.ss.topts$Species=="E.saligna")]<--30.56
e.ss.topts$SSloc[which(e.ss.topts$Species=="E.sideroxylon")]<--32.98 



#get g1

e.ss_g1<-get_topts(lapply(e.ss.mix,FUN=fitStom))

e.ss_g1$Species<- get_names(e.ss.mix)[,1]
e.ss_g1$Temp_Treatment<- get_names(e.ss.mix)[,2]
e.ss.topts<-merge(e.ss.topts,e.ss_g1,by=c("Species","Temp_Treatment"))


#get complete species name

e.ss.topts$Species[which(e.ss.topts$Species=="E.saligna")]<-"Eucalyptus saligna"
e.ss.topts$Species[which(e.ss.topts$Species=="E.sideroxylon")]<-"Eucalyptus sideroxylon" 


#e.ss.topts$T_seedS<-NA
#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/ClimOrig/S30/Saligna"
#e.ss.topts$T_seedS[which(e.ss.topts$Species=="E.saligna")]<-met_out(path=path,from="1965-01-01",to="2015-12-31")

# path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/ClimOrig/S30/Sideroxylon"
#e.ss.topts$T_seedS[which(e.ss.topts$Species=="E.sideroxylon")]<-met_out(path=path,from="1965-01-01",to="2015-12-31")

e.ss.topts$DataSet<-"S30"

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#E.tetradonta

e.tr<-read.csv(paste(path,"/Data/SIOP Leaf gas exchange Cernusak et al.csv",sep=""))

e.tr<-subset(e.tr,abs((e.tr$CO2R)-400)<10)
e.tr$TargetTemp<-cut(e.tr$Tleaf,breaks=c(24,26,32,36,42),
                     labels=FALSE)

e.tr.mix<-split(e.tr,e.tr$Canopy_position)

e.tr.topts<-lapply(e.tr.mix,function(x)fit.nlme(x,yvar="Photo",random="Leaf.number"))


e.tr.topts<-data.frame(do.call(rbind,list(e.tr.topts[[1]]))) #extract parameters only for canopy

e.tr.topts$Species<-"Eucalyptus tetrodonta"

#get mean gs,VPD and Ci
#e.tr.sum<-summaryBy(Cond+VpdL+Ci~Canopy_position,keep.names=TRUE,data=e.tr,FUN=c(mean,std.error))


e.tr.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci+Tleaf~Canopy_position+TargetTemp,data=e.tr,FUN=c(mean,std.error),keep.names=TRUE)
e.tr.sum<-subset(e.tr.sum,e.tr.sum$TargetTemp==1)
e.tr.topts<-merge(e.tr.topts,e.tr.sum)
e.tr.topts$Canopy_position<-NULL


#add latitude of seed source
e.tr.topts$SSloc<--14.15

e.tr.topts$DataSet<-"SAVANNA"

e.tr_g1<-get_topts(lapply(e.tr.mix[1],FUN=fitStom))
e.tr.topts<-cbind(e.tr.topts,e.tr_g1)
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#E.globulus (Tasmania;Mike Battaglia)

e.glob2<-read.csv(paste(path,"/Data/GLOB_MB.csv",sep=""))
e.glob2$DATE<-as.factor(e.glob2$DATE)
e.glob3<-split(e.glob2,e.glob2$Month)
e.glob2.fit<-lapply(e.glob3,FUN=fitquad)
e.glob2.fit.para<-data.frame(do.call(rbind,list(e.glob2.fit[[1]],e.glob2.fit[[2]],e.glob2.fit[[3]],e.glob2.fit[[4]],
                                                e.glob2.fit[[5]],e.glob2.fit[[6]])))

e.glob2.fit.para$Month<-names(e.glob3)
e.glob2.fit.para$Species<-"Eucalyptus globulus"

#to drop estimates for month=Dec (this time period is not inclded in MB's paper. Year of measurements is unknown)
e.glob2.fit.para<-subset(e.glob2.fit.para,e.glob2.fit.para$Month!="Dec")

#get mean gs,VPD and Ci
e.glob2.sum<-summaryBy(Photo+Cond+Ci~Month+Tleaf,keep.names=TRUE,data=e.glob2,FUN=c(mean,std.error))
e.glob2.sum<-subset(e.glob2.sum,e.glob2.sum$Tleaf==25)

#colnames(e.glob2.sum)[1]<-"Month"
e.glob2.fit.para<-merge(e.glob2.fit.para,e.glob2.sum,by="Month")

#add latitude of seed source
e.glob2.fit.para$SSloc<--42.81
#to add longteram Temperature at seed source

#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/ClimOrig/Hobart"
#e.glob2.fit.para$T_seedS<-as.numeric(met_out(path=path,from="1965-01-01",to="2015-12-31"))

e.glob2.fit.para$DataSet<-"TASMANIA"


#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#WTC4 

wtc4<-read.csv(paste(path,"/Data/wtc4_photo_temp_all.csv",sep=""))
wtc4$Tempfrac<-cut(wtc4$Tleaf,breaks=c(13,18,22,27,32,37.5,44),labels=FALSE)
wtc4.1<-summaryBy(.~Season+Ttreatment+Chamber+Tempfrac,fun=mean,keep.names=TRUE,data=wtc4)


#fit mix models

wtc4.mix<-split(wtc4.1,paste(wtc4.1$Ttreatment,wtc4.1$Season))
wtc4.topts<-lapply(wtc4.mix,function(x)fit.nlme(x,yvar="Photo",random="Chamber"))

wtc4.topts<-data.frame(do.call(rbind,list(wtc4.topts[[1]],wtc4.topts[[2]],wtc4.topts[[3]],wtc4.topts[[4]]
)))
Treatment<-strsplit(names(wtc4.mix)," ")

wtc4.topts$Temp_Treatment<-c(Treatment[[1]][[1]],Treatment[[2]][[1]],Treatment[[3]][[1]],Treatment[[4]][[1]])
wtc4.topts$Month<-c(Treatment[[1]][[2]],Treatment[[2]][[2]],Treatment[[3]][[2]],Treatment[[4]][[2]])

wtc4.topts$Species<-"Eucalyptus parramattensis"

#get mean gs,VPD and Ci
e.par.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci~Season+Ttreatment+Tempfrac,keep.names=TRUE,data=wtc4.1,FUN=c(mean,std.error))
e.par.sum<-subset(e.par.sum,e.par.sum$Tempfrac==3)

colnames(e.par.sum)[1]<-"Month"
colnames(e.par.sum)[2]<-"Temp_Treatment"

wtc4.topts<-merge(wtc4.topts,e.par.sum,by=c("Temp_Treatment","Month"))

#add latitude of seed source
wtc4.topts$SSloc<--33.62

wtc4.topts$DataSet<-"WTC4"

wtc4_g1<-get_topts(lapply(wtc4.mix,FUN=fitStom))
wtc4_g1$Temp_Treatment<-get_names(wtc4.mix)[,1]
wtc4_g1$Month<-get_names(wtc4.mix)[,2]

wtc4.topts<-merge(wtc4.topts,wtc4_g1,by=c("Temp_Treatment","Month"))

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#Rainforests

dat.rf<-read.csv(paste(path,"/Data/Daintree_ACidata_processed.csv",sep=""))
dat.rf.a<-subset(dat.rf,abs((dat.rf$CO2R)-400)<25)
dat.rf.a<-subset(dat.rf.a,dat.rf.a$Ci>175 )


dat.rf.a$Tempfrac<-cut(dat.rf.a$Tleaf,breaks=c(17,23,27,33,37,42),labels=FALSE)
#dat.rf.a<-summaryBy(.~Spp+Season+Tempfrac,data=dat.rf.a,FUN=mean,keep.names=TRUE)

rf.mix<-split(dat.rf.a,paste(dat.rf.a$Season))
#rf.mix<-split(dat.rf.a,paste(dat.rf.a$Spp,dat.rf.a$Season))

#rf.topts.1<-lapply(rf.mix[2],FUN=fitquad) #cannot fit mix models.
rf.topts<-get_topts(lapply(rf.mix,function(x)fit.nlme(x,yvar="Photo",random="Spp"))) #cannot fit mix models.
#rf.mix[2]<-NULL
#rf.topts.2<-lapply(rf.mix[4],function(x)fit.nlme(x,yvar="Photo",random="Leaf"))


#rf.topts<-data.frame(do.call(rbind,list(rf.topts.2[[1]],rf.topts.1[[1]],rf.topts.2[[2]],rf.topts.2[[3]])))


rf.topts$Month<-get_names(rf.mix)[,2]
rf.topts$Species<-as.factor("RF_AUSTRALIA")
#rf.topts$Month<-c("SUMMER","WINTER","SUMMER","WINTER")
#rf.topts$Month<-as.factor(rf.topts$Month)


#get mean gs,VPD and Ci
rf.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci~Season+Tempfrac,keep.names=TRUE,data=dat.rf.a,FUN=c(mean,std.error))
rf.sum<-subset(rf.sum,rf.sum$Tempfrac==2)

colnames(rf.sum)[1]<-"Month"
#colnames(rf.sum)[1]<-"Species"
#rf.topts<-merge(rf.topts,rf.sum,by="Month")
rf.topts<-cbind(rf.topts,rf.sum[c(3:12)])
rf.topts<-rf.topts[-c(12,13)]
#combine all fitted parameters to a data frame

#add latitude of seed source
rf.topts$SSloc<--16.10

rf.topts$DataSet<-"RF_AUS"

#get g1

rf_g1<-get_topts(lapply(rf.mix,FUN=fitStom))

rf.topts<-cbind(rf.topts,rf_g1)


#get complete species name
#rf.topts$Species<-as.character(rf.topts$Species)
#rf.topts$Species[which(rf.topts$Species=="ACGR")]<-"Syzygium graveolens"
#rf.topts$Species[which(rf.topts$Species=="ARPE")]<-"Argyrodendron peralatum"
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#GWW Species
#- fit all data together

gww<-read.csv(paste(path,"/Data/gwwACidata_processed_new.csv",sep=""))
gwwphoto<-subset(gww,abs(gww$CO2S-400)<=30.5)
gwwphoto<-subset(gwwphoto,!gwwphoto$Curve %in% c(19,23,27,44)) #remove four outliers. Three at highest leaf temperature and one at lower
                                                                # not affect Topt or Aopt significantly but lower standard errors by ~6C
#to get average across leaf temperatures~6
gwwphoto$Tleaffrac<-cut(gwwphoto$Tleaf,breaks=c(17,23,26,32,37),labels=FALSE)
#gwwphoto.1<-summaryBy(.~Spp+Tleaffrac,data=gwwphoto,FUN=mean,keep.names=TRUE)


#fot E.salmonophloia and E.scoparia, mix models was fitted 
gww.topts<-data.frame(do.call(rbind,list(fit.nlme(gwwphoto,yvar="Photo",rand="Tree"))))

#gww.topts.1<-fit.nlme(subset(gwwphoto,gwwphoto$Species=="E.scoparia"),yvar="Photo",rand="Tree") #not an Eucalypt
#gww.topts.2<-fit.nlme(subset(gwwphoto,gwwphoto$Species=="E.salmonophloia"),yvar="Photo",rand="Tree")

#for E.salubris, model was fitted to average across replicate trees. Cannott fit mix models 
#gww.topts.3<-fitquad(subset(gwwphoto.1,gwwphoto.1$Species=="E.salubris"))

#gww.topts<-data.frame(do.call(rbind,list(gww.topts.2,gww.topts.3)))

#gww.topts$Species<-c("E.salmonophloia","E.salubris")
gww.topts$Species<-"ASAW spp"

#get mean gs,VPD and Ci
gww.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci+Tleaf~Tleaffrac,keep.names=TRUE,data=gwwphoto,FUN=c(mean,std.error))
gww.sum<-subset(gww.sum,gww.sum$Tleaffrac==2)

gww.topts<-cbind(gww.topts,gww.sum)
gww.topts$Species<-"ASAW spp"
gww.topts$DataSet<-"GWW"
#gww.topts.4<-fitquad(subset(gwwphoto,gwwphoto$Species=="E.salubris"))
#with(subset(gwwphoto.1,gwwphoto.1$Species=="E.scoparia"),plot(Tleaf,Photo,cex=2,pch=16,col=Tree))
#with(subset(gwwphoto.1,gwwphoto.1$Species=="E.salmonophloia"),plot(Tleaf,Photo,cex=2,pch=16))
#with(subset(gwwphoto.1,gwwphoto.1$Species=="E.salubris"),plot(Tleaf,Photo,cex=2,pch=16))


#add latitude of seed source
gww.topts$SSloc<--30.19
gww.topts$DataSet<-"GWW"


#get_g1
gww_split<-split(gwwphoto.1,gwwphoto.1$Species)
gww_g1<-get_topts(lapply(gww_split,FUN=fitStom))
gww_g1$Species<-names(gww_split)

gww.topts<-merge(gww.topts,gww_g1,by="Species")

#get complete species name

gww.topts$Species[which(gww.topts$Species=="E.salmonophloia")]<-"Eucalyptus salmonophloia"
gww.topts$Species[which(gww.topts$Species=="E.salubris")]<-"Eucalyptus salubris"
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#E.delegatensis (Tumbarumba EC site)

#e.del.1$xterm<-with(e.del.1,Photo/CO2S*sqrt(VpdL))

e.del<-read.csv(paste(path,"/Data/TumbarumbaGasex_Spot_Medlyn.csv",sep=""))
e.del.1<-subset(e.del,e.del$PARi>1200& e.del$Cond>0.09 ) #to get saturated PAR and remove abnormally low datapoints in November

with(e.del.1,plot(Tleaf,Photo,col=factor(Season),pch=16,cex=2))
e.del.topts<-(fit.nlme(e.del.1,yvar="Photo",random="Age.Class")) 

e.del.topt<-data.frame(do.call(rbind,list(e.del.topts))) 
e.del.topt$Species<-as.factor(e.del.topt$Species<-"Eucalyptus delegatensis")
e.del.topt$Facility<-as.factor(e.del.topt$Facility<-"TMB")

#e.del.topt$Month<-"February" #this is not actually data from February. Month make to February to get 
#e.del.topt$Month<-"February"#summer values of Vcmax and Jmax. these Vcmax and Jmax values are similar to 
# the values when fit the response without considering seasons
#e.del.topt$Temp_Treatment<-"Ambient"

#get Tleaf around 25
dat25<-subset(e.del.1,abs(e.del.1$Tleaf-25)<=2)

e.del.topt$Cir<-mean(dat25$Ci/dat25$CO2S)
e.del.topt$Cir.std.error<-std.error(dat25$Ci/dat25$CO2S)


e.del.topt$Photo.mean<-mean(dat25$Photo)
e.del.topt$Photo.std.error<-std.error(dat25$Photo)

e.del.topt$Cond.mean<-mean(dat25$Cond)
e.del.topt$Cond.std.error<-std.error(dat25$Cond)

e.del.topt$VpdL.mean<-mean(dat25$VpdL)
e.del.topt$VpdL.std.error<-std.error(dat25$VpdL)

e.del.topt$Ci.mean<-mean(dat25$Ci)
e.del.topt$Ci.std.error<-std.error(dat25$Ci)

colnames(e.del.topt)[c(14,15)]<-c("Ci/CO2S.mean","Ci/CO2S.std.error")

#add latitude of seed source
e.del.topt$SSloc<--35.65


#e.del.topt$Env<-"Field"

#e.del.topt$Type<-"Eucalypts_Temperate"

e.del.topt$DataSet<-"TMB"

#get g1

tmb_g1<-fitStom(e.del.1)
e.del.topt<-cbind(e.del.topt,tmb_g1)


#fitarh(d4) #fit without consedering seasons
#---------------------------------------------------------------------------------------------------------------------
#e.del.aci<-read.csv(paste(path,"/Tumbarumba_ACidata_processed.csv",sep=""))
#e.del.1.aci<-subset(e.del.aci,abs(e.del.aci$CO2R-350)<=15)
#e.del.1.aci$Leaf.Age<-as.factor( e.del.1.aci$Leaf.Age)



#e.del.topts<-fit.nlme(dat=e.del.1.aci,yvar="Photo",random="Age.Class")
#e.del.topts<-fitquad(e.del.1.aci)

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

#Corymbria (Mike aspinwall)

cor<-read.csv(paste(path,"/Data/mike_aspinwall_corymbia_calophylla.V1.csv",sep=""))

cor.a<-subset(cor,cor$CO2R>415 & cor$CO2R<425) #to get ambient CO2 levels from ACi curves
cor.a<-subset(cor.a,cor.a$Cond<0.7)   #remove 4 points with abnormally high Stomatal Conductance


cor.a$xterm<-with(cor.a,Photo/CO2S*sqrt(VpdL))
with(cor.a,plot(Cond~xterm,col=Provenance))

cor.b<-split(cor.a,paste(cor.a$Growth_temp,cor.a$Provenance))

cor.topts<-get_topts(lapply(cor.b,function(x)fit.nlme(x,yvar="Photo",random="PotID")))
cor.topts$Growth_temp<-get_names(cor.b)[,1]
cor.topts$Provenance<-get_names(cor.b)[,2]


#get mean gs,VPD and Ci
cor.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci~Growth_temp+Provenance+Meas_temp,data=cor.a,FUN=c(mean,std.error))
cor.sum<-subset(cor.sum,cor.sum$Meas_temp==25)
cor.sum<-cor.sum[-c(3)]
cor.topts<-merge(cor.topts,cor.sum,by=c("Growth_temp","Provenance"))


cor.topts$Tavg_30<-NA
cor.topts$Tavg_30[which(cor.topts$Growth_temp==26)]<-19 #average temperature (day and night settings). Source: Aspinwall et al (Tree Phys, submitted)
cor.topts$Tavg_30[which(cor.topts$Growth_temp==32)]<-25
cor.topts$DataSet<-"MIKE_ASPINWALL"


#get g1
cor_g1<-get_topts(lapply(cor.b,FUN=fitStom))
cor_g1$Growth_temp<-get_names(cor.b)[,1]
cor_g1$Provenance<-get_names(cor.b)[,2]

cor.topts<-merge(cor.topts,cor_g1,by=c("Growth_temp","Provenance"))

names(cor.topts)[c(1:2)]<-c("Temp_Treatment","Species")

cor.topts$Species[which(cor.topts$Species=="BOO")]<-"Corymbia calophylla_CW"
cor.topts$Species[which(cor.topts$Species=="CRI")]<-"Corymbia calophylla_CD"
cor.topts$Species[which(cor.topts$Species=="GIN")]<-"Corymbia calophylla_WW"
cor.topts$Species[which(cor.topts$Species=="MOG")]<-"Corymbia calophylla_WD"

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#E. pauciflora (Kirschbaum et al)

epau<-read.csv(paste(path,"/Data/euc_pauciflora_mk.csv",sep=""))
epau.1<-do.call(rbind,by(epau,epau$Curve,function(x) x[4,]))
epau.topts<-fitquad(epau.1)
epau.topts<-data.frame(do.call(rbind,list(epau.topts)))

epau.topts$Species<-"Eucalyptus pauciflora"
epau.topts$DataSet<-"KIRSCHBAUM"

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#GREAT (E. tereticornis)

avt<-read.csv(paste(path,"/Data/great_photo_shortterm.csv",sep=""))
#- fit AvT to estimate Topts at high light, AREA BASED
tofit <- subset(avt, LightFac==4)
tofit$Tempfrac<-cut(tofit$Tleaf,breaks=c(18,23,26,32,37.5,44),labels=FALSE)
with(tofit,plot(Tleaf,Photo,col=Prov,pch=16))


tofit.l <- split(tofit,tofit$location)
gr.topts<-get_topts(lapply(tofit.l,function(x)fit.nlme(x,yvar="Photo",random="Code")))
gr.topts$Species<-names(tofit.l)
gr.topts$DataSet<-"GREAT"
gr.topts$Tavg_30<-21.5

#get mean gs,VPD and Ci
gr.sum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci~Growth_temp+location+Tempfrac,data=tofit,FUN=c(mean,std.error))
gr.sum<-subset(gr.sum,gr.sum$Tempfrac==2)
gr.sum<-gr.sum[-c(2)]
names(gr.sum)[1]<-"Species"
gr.sum<-merge(gr.topts,gr.sum,by=c("Species"))


gr.topts$Species[which(gr.topts$Species=="Central")]<-"Eucalyptus tereticornis_C"
gr.topts$Species[which(gr.topts$Species=="Cold-edge")]<-"Eucalyptus tereticornis_Co"
gr.topts$Species[which(gr.topts$Species=="Warm-edge")]<-"Eucalyptus tereticornis_W"


#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#EucFace spot measurements

eface<-read.csv(paste(path,"/Data/eucFace_spot.csv",sep=""))
eface.a<-subset(eface,eface$CO2_Treat1=="amb")

eface.topts<-data.frame(do.call(rbind,list(fit.nlme(eface.a,yvar="Photo",random="Tree")))) 
eface.topts$Species<-"Eucalyptus tereticornis"
eface.topts$DataSet<-"EucFace"
#eface.topts$Tavg_30<-
  

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------




#this script calculates average temperature across different time periods before photosynthesisi measured
#30 days//60 days//90 days before
#get tmin//tmax and tavg~(tmin+tmax/2)

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#hfe common garden

#to add growth temperature 30 days before measurements 
#January measurements
hfj<-wtcMet2(path=paste0(path,"/MET_DATA/WTC1"),
             from="2010-01-04",fname="HFE_TEMP_CM_WTCMET_L1_v2.csv")

#August measurements
hfa<-wtcMet2(path=paste0(path,"/MET_DATA/WTC1"),
             from="2009-08-27",fname="HFE_TEMP_CM_WTCMET_L1_v2.csv")



#Average daily temperature
hfe.topts$Tavg_30<-NA
hfe.topts$Tavg_30[which(hfe.topts$Month=="February")]<-hfj[[1]]
hfe.topts$Tavg_30[which(hfe.topts$Month=="August")]<-hfa[[1]]

#Average daily minimum
hfe.topts$Tmin_30<-NA
hfe.topts$Tmin_30[which(hfe.topts$Month=="February")]<-hfj[[2]]
hfe.topts$Tmin_30[which(hfe.topts$Month=="August")]<-hfa[[2]]

#Average daily maximum
hfe.topts$Tmax_30<-NA
hfe.topts$Tmax_30[which(hfe.topts$Month=="February")]<-hfj[[3]]
hfe.topts$Tmax_30[which(hfe.topts$Month=="August")]<-hfa[[3]]


#to add longteram Temperature at seed source

ecre<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="hfecre")
edun<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="hfedun")
emel<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="hfemel")
esal<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="hfesal")
eter<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="hfeter")
ecla<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="hfecla")

#TavgSS: Average temperature across entire time duration considered (1975-2015)
hfe.topts$TavgSS<-NA
hfe.topts$TavgSS[which(hfe.topts$Species=="Eucalyptus crebra")]<-ecre[[1]]
hfe.topts$TavgSS[which(hfe.topts$Species=="Eucalyptus dunnii")]<-edun[[1]]
hfe.topts$TavgSS[which(hfe.topts$Species=="Eucalyptus melliodora")]<-emel[[1]]
hfe.topts$TavgSS[which(hfe.topts$Species=="Eucalyptus saligna")]<-esal[[1]]
hfe.topts$TavgSS[which(hfe.topts$Species=="Eucalyptus tereticornis")]<-eter[[1]]
hfe.topts$TavgSS[which(hfe.topts$Species=="Eucalyptus cladocalyx")]<-ecla[[1]]



#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#WTC3

#to add met data

apr<-wtcMet(path=paste0(path,"/MET_DATA/WTC3"),
            fname="WTC_TEMP_CM_WTCMET_L1_v2.csv",from="2014-04-22")

sep<-wtcMet(path=paste0(path,"/MET_DATA/WTC3"),
            fname="WTC_TEMP_CM_WTCMET_L1_v2.csv",from="2013-09-17")

jan<-wtcMet(path=paste0(path,"/MET_DATA/WTC3"),
            fname="WTC_TEMP_CM_WTCMET_L1_v2.csv",from="2014-01-13")


#before 30 day mean temperature
wtc3.topts$Tavg_30<-NA
wtc3.topts$Tavg_30[which(wtc3.topts$Temp_Treatment=="Ambient" & wtc3.topts$Month=="April")]<-apr[[4]][[1]]
wtc3.topts$Tavg_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="April")]<-apr[[4]][[2]]
wtc3.topts$Tavg_30[which(wtc3.topts$Temp_Treatment=="Ambient" & wtc3.topts$Month=="September")]<-sep[[4]][[1]]
wtc3.topts$Tavg_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="September")]<-sep[[4]][[2]]
wtc3.topts$Tavg_30[which(wtc3.topts$Temp_Treatment=="Ambient" & wtc3.topts$Month=="January")]<-jan[[4]][[1]]
wtc3.topts$Tavg_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="January")]<-jan[[4]][[2]]


#before 30 day min temperature
wtc3.topts$Tmin_30<-NA
wtc3.topts$Tmin_30[which(wtc3.topts$Temp_Treatment=="Ambient" & wtc3.topts$Month=="April")]<-apr[[2]][[1]]
wtc3.topts$Tmin_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="April")]<-apr[[2]][[2]]
wtc3.topts$Tmin_30[which(wtc3.topts$Temp_Treatment=="Ambient" & wtc3.topts$Month=="September")]<-sep[[2]][[1]]
wtc3.topts$Tmin_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="September")]<-sep[[2]][[2]]
wtc3.topts$Tmin_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="January")]<-jan[[2]][[2]]

#before 30 day max temperature
wtc3.topts$Tmax_30<-NA
wtc3.topts$Tmax_30[which(wtc3.topts$Temp_Treatment=="Ambient" & wtc3.topts$Month=="April")]<-apr[[3]][[1]]
wtc3.topts$Tmax_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="April")]<-apr[[3]][[2]]
wtc3.topts$Tmax_30[which(wtc3.topts$Temp_Treatment=="Ambient" & wtc3.topts$Month=="September")]<-sep[[3]][[1]]
wtc3.topts$Tmax_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="September")]<-sep[[3]][[2]]
wtc3.topts$Tmax_30[which(wtc3.topts$Temp_Treatment=="Elevated" & wtc3.topts$Month=="January")]<-jan[[3]][[2]]
wtc3.topts$Tmax_30[which(wtc3.topts$Temp_Treatment=="Ambient" & wtc3.topts$Month=="January")]<-jan[[3]][[2]]


#to add longteram Temperature at seed source

rich<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="wtc3")

wtc3.topts$TavgSS<-rich[[1]]

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#WTC2

#to add met data

dec1<-wtcMet(path=paste0(path,"/MET_DATA/WTC2"),
             fname="WTC2_TEMP_CM_WTCMET_L1_v2.csv",from="2010-12-03")

sep1<-wtcMet(path=paste0(path,"/MET_DATA/WTC2"),
             fname="WTC2_TEMP_CM_WTCMET_L1_v2.csv",from="2011-09-01")

feb1<-wtcMet(path=paste0(path,"/MET_DATA/WTC2"),
             fname="WTC2_TEMP_CM_WTCMET_L1_v2.csv",from="2011-02-09")

#before 30 day mean temperature
wtc2.topt$Tavg_30<-NA
wtc2.topt$Tavg_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="December")]<-dec1[[4]][[1]]
wtc2.topt$Tavg_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="December")]<-dec1[[4]][[2]]
wtc2.topt$Tavg_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="September")]<-sep1[[4]][[1]]
wtc2.topt$Tavg_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="September")]<-sep1[[4]][[2]]
wtc2.topt$Tavg_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="February")]<-feb1[[4]][[1]]
wtc2.topt$Tavg_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="February")]<-feb1[[4]][[2]]


#before 30 day min temperature
wtc2.topt$Tmin_30<-NA
wtc2.topt$Tmin_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="December")]<-dec1[[2]][[1]]
wtc2.topt$Tmin_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="December")]<-dec1[[2]][[2]]
wtc2.topt$Tmin_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="September")]<-sep1[[2]][[1]]
wtc2.topt$Tmin_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="September")]<-sep1[[2]][[2]]
wtc2.topt$Tmin_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="February")]<-feb1[[2]][[1]]
wtc2.topt$Tmin_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="February")]<-feb1[[2]][[2]]

#before 30 day max temperature
wtc2.topt$Tmax_30<-NA
wtc2.topt$Tmax_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="December")]<-dec1[[3]][[1]]
wtc2.topt$Tmax_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="December")]<-dec1[[3]][[2]]
wtc2.topt$Tmax_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="September")]<-sep1[[3]][[1]]
wtc2.topt$Tmax_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="September")]<-sep1[[3]][[2]]
wtc2.topt$Tmax_30[which(wtc2.topt$Temp_Treatment=="Ambient" & wtc2.topt$Month=="February")]<-feb1[[3]][[1]]
wtc2.topt$Tmax_30[which(wtc2.topt$Temp_Treatment=="Elevated" & wtc2.topt$Month=="February")]<-feb1[[3]][[2]]

#-----

#to add longteram Temperature at seed source

capeO<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="wtc2")

wtc2.topt$TavgSS<-capeO[[1]]

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#WTC1

jan3<-wtcMet2(path=paste0(path,"/MET_DATA/WTC1"),
              from="2009-01-19",fname="HFE_TEMP_CM_WTCMET_L1_v2.csv")

nov3<-wtcMet2(path=paste0(path,"/MET_DATA/WTC1"),
              from="2008-11-04",fname="HFE_TEMP_CM_WTCMET_L1_v2.csv")

#before 30 day mean temperature
wtc1.topts$Tavg_30<-NA
wtc1.topts$Tavg_30[which(wtc1.topts$Month=="January")]<-jan3[[1]]
wtc1.topts$Tavg_30[which(wtc1.topts$Month=="November")]<-nov3[[1]]

#before 30 day mean minimum temperature
wtc1.topts$Tmin_30<-NA
wtc1.topts$Tmin_30[which(wtc1.topts$Month=="January")]<-jan3[[2]]
wtc1.topts$Tmin_30[which(wtc1.topts$Month=="November")]<-nov3[[2]]

#before 30 day mean maximum temperature
wtc1.topts$Tmax_30<-NA
wtc1.topts$Tmax_30[which(wtc1.topts$Month=="January")]<-jan3[[3]]
wtc1.topts$Tmax_30[which(wtc1.topts$Month=="November")]<-nov3[[3]]

#---

#to add longteram Temperature at seed source

arm<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="wtc1")

wtc1.topts$TavgSS<-arm[[1]]

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#WTC4

win<-wtcMet(path=paste0(path,"/MET_DATA/WTC4"),
            fname="WTC_TEMP_CM_WTCMET_L1_v2.csv",from="2016-07-07")

spr<-wtcMet(path=paste0(path,"/MET_DATA/WTC4"),
            fname="WTC_TEMP_CM_WTCMET_L1_v2.csv",from="2016-10-05")

#before 30 day mean temperature
wtc4.topts$Tavg_30<-NA
wtc4.topts$Tavg_30[which(wtc4.topts$Temp_Treatment=="ambient" & wtc4.topts$Month=="Winter")]<-win[[4]][[1]]
wtc4.topts$Tavg_30[which(wtc4.topts$Temp_Treatment=="elevated" & wtc4.topts$Month=="Winter")]<-win[[4]][[2]]
wtc4.topts$Tavg_30[which(wtc4.topts$Temp_Treatment=="ambient" & wtc4.topts$Month=="Spring")]<-spr[[4]][[1]]
wtc4.topts$Tavg_30[which(wtc4.topts$Temp_Treatment=="elevated" & wtc4.topts$Month=="Spring")]<-spr[[4]][[2]]

#before 30 day mean min temperature
wtc4.topts$Tmin_30<-NA
wtc4.topts$Tmin_30[which(wtc4.topts$Temp_Treatment=="ambient" & wtc4.topts$Month=="Winter")]<-win[[2]][[1]]
wtc4.topts$Tmin_30[which(wtc4.topts$Temp_Treatment=="elevated" & wtc4.topts$Month=="Winter")]<-win[[2]][[2]]
wtc4.topts$Tmin_30[which(wtc4.topts$Temp_Treatment=="ambient" & wtc4.topts$Month=="Spring")]<-spr[[2]][[1]]
wtc4.topts$Tmin_30[which(wtc4.topts$Temp_Treatment=="elevated" & wtc4.topts$Month=="Spring")]<-spr[[2]][[2]]

#before 30 day mean max temperature
wtc4.topts$Tmax_30<-NA
wtc4.topts$Tmax_30[which(wtc4.topts$Temp_Treatment=="ambient" & wtc4.topts$Month=="Winter")]<-win[[3]][[1]]
wtc4.topts$Tmax_30[which(wtc4.topts$Temp_Treatment=="elevated" & wtc4.topts$Month=="Winter")]<-win[[3]][[2]]
wtc4.topts$Tmax_30[which(wtc4.topts$Temp_Treatment=="ambient" & wtc4.topts$Month=="Spring")]<-spr[[3]][[1]]
wtc4.topts$Tmax_30[which(wtc4.topts$Temp_Treatment=="elevated" & wtc4.topts$Month=="Spring")]<-spr[[3]][[2]]

#-----


#to add longteram Temperature at seed source

rich<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="local")

wtc4.topts$TavgSS<-rich[[1]]

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#Savanna: E.tetradonta
#newpath<-"C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/Met from Mingkai/Processed"

sav<-wtcMet2(path=paste0(path,"/MET_DATA/Met from Mingkai/Processed"),
             from="2009-09-08",fname="DalyUncleared_met_forcing.csv")

e.tr.topts$Tavg_30<-sav[[1]]
e.tr.topts$Tmin_30<-sav[[2]]
e.tr.topts$Tmax_30<-sav[[3]]

#to add longteram Temperature at seed source


sav.1<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="daly")


e.tr.topts$TavgSS<-sav.1[[1]]
#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#GWW Species

#newpath<-"C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/Met from Mingkai/Processed"

gww<-wtcMet2(path=paste0(path,"/MET_DATA/Met from Mingkai/Processed"),
             from="2013-04-04",fname="greatwesternwoodlands_met_forcing.csv")

gww.topts$Tavg_30<-gww[[1]]
gww.topts$Tmin_30<-gww[[2]]
gww.topts$Tmax_30<-gww[[3]]

#to add longteram Temperature at seed source


gww.1<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="gww")

gww.topts$TavgSS<-gww.1[[1]]

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#E.delegatensis (Tumbarumba EC site)

#newpath<-"C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/Met from Mingkai/Processed"

tmb<-wtcMet2(path=paste0(path,"/MET_DATA/Met from Mingkai/Processed"),
             from="2002-02-11",fname="Tumbarumba_met_forcing.csv")

e.del.topt$Tavg_30<-tmb[[1]]
e.del.topt$Tmin_30<-tmb[[2]]
e.del.topt$Tmax_30<-tmb[[3]]

#to add longterm Temperature at seed source


tmb.1<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="tmb")

e.del.topt$TavgSS<-tmb.1[[1]]

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#RF species

sum1<-wtcMet2(path=paste0(path,"/MET_DATA/ClimOrig/RF"),
              from="2011-04-07",fname="IDCJAC0011_031037_1800_all.csv")

win1<-wtcMet2(path=paste0(path,"/MET_DATA/ClimOrig/RF"),
              from="2010-07-22",fname="IDCJAC0011_031037_1800_all.csv")

#before 30 day mean temperature
rf.topts$Month<-as.factor(rf.topts$Month)

rf.topts$Tavg_30<-NA
rf.topts$Tavg_30[which(rf.topts$Month=="SUMMER")]<-sum1[[1]]
rf.topts$Tavg_30[which(rf.topts$Month=="WINTER")]<-win1[[1]]

rf.topts$Tmin_30<-NA
rf.topts$Tmin_30[which(rf.topts$Month=="SUMMER")]<-sum1[[2]]
rf.topts$Tmin_30[which(rf.topts$Month=="WINTER")]<-win1[[2]]

rf.topts$Tmax_30<-NA
rf.topts$Tmax_30[which(rf.topts$Month=="SUMMER")]<-sum1[[3]]
rf.topts$Tmax_30[which(rf.topts$Month=="WINTER")]<-win1[[3]]


#to add longterm Temperature at seed source
#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/ClimOrig/RF"
#tmb.1<-met_out(path=path,from="1975-01-01",to="2015-12-31")

rf.1<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="rf")

rf.topts$TavgSS<-rf.1[[1]]

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#S30 GS species

#add Tgrowth. 
#this is from Oula et al, 2010 PCE paper. Need raw data to get other averages
e.ss.topts$Tavg_30<-NA
e.ss.topts$Tavg_30[which(e.ss.topts$Temp_Treatment=="Ambient")]<-22
e.ss.topts$Tavg_30[which(e.ss.topts$Temp_Treatment=="Elevated")]<-26

e.ss.topts$Tmax_30<-NA
e.ss.topts$Tmax_30[which(e.ss.topts$Temp_Treatment=="Ambient")]<-26
e.ss.topts$Tmax_30[which(e.ss.topts$Temp_Treatment=="Elevated")]<-30


s30sa<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="S30ESa")
s30sy<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="S30ESy")

e.ss.topts$TavgSS<-NA
e.ss.topts$TavgSS[which(e.ss.topts$Species=="Eucalyptus saligna")]<-s30sa[[1]]
e.ss.topts$TavgSS[which(e.ss.topts$Species=="Eucalyptus sideroxylon")]<-s30sy[[1]]


#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------
#E. globulus: Tasmania

#one month before
e.glob2.fit.para$Tavg_30<-NA
e.glob2.fit.para$Tavg_30[which(e.glob2.fit.para$Month=="April")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1994-04-01")[[1]]
e.glob2.fit.para$Tavg_30[which(e.glob2.fit.para$Month=="December")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1994-01-01")[[1]]
e.glob2.fit.para$Tavg_30[which(e.glob2.fit.para$Month=="February")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1994-02-01")[[1]]
e.glob2.fit.para$Tavg_30[which(e.glob2.fit.para$Month=="November")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1993-12-01")[[1]]
e.glob2.fit.para$Tavg_30[which(e.glob2.fit.para$Month=="October")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1993-11-01")[[1]]

e.glob2.fit.para$Tmax_30<-NA
e.glob2.fit.para$Tmax_30[which(e.glob2.fit.para$Month=="April")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1994-04-01")[[3]]
e.glob2.fit.para$Tmax_30[which(e.glob2.fit.para$Month=="December")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1994-01-01")[[3]]
e.glob2.fit.para$Tmax_30[which(e.glob2.fit.para$Month=="February")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1994-02-01")[[3]]
e.glob2.fit.para$Tmax_30[which(e.glob2.fit.para$Month=="November")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1993-12-01")[[3]]
e.glob2.fit.para$Tmax_30[which(e.glob2.fit.para$Month=="October")]<-met_eMAST.2(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania",from="1993-11-01")[[3]]

#to add longterm Temperature at seed source

tas<-met_eMAST(path=paste0(path,"/MET_DATA/eMAST"),fname="tasmania")
e.glob2.fit.para$TavgSS<-tas[[1]]

#-----------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#species outside Australia

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#Martijn Slot's data from Panama

panama<-read.csv(paste(path,"/Data/data_out_auz/SlotWinter2017_NewPhyt.Data.csv",sep=""))
panama.1<-subset(panama,panama$Form=="Tree")
with(panama.1,plot(Tleaf,Photo,col=Site))
#panama.2<-split(panama.1,panama.1$Site)

panama.topts<-data.frame(do.call(rbind,list(fit.nlme(panama.1,yvar="Photo",random="Species")))) 
#panama.topts$Site<-names(panama.2)
panama.topts$DataSet<-"Martijn_Slot"
#panama.topts$Species<-c("Fort Sherman Spp", "PNM Spp")
panama.topts$Species<-"Panama Rainforest Spp"
#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM"
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#Amazone

amz.1<-read.csv(paste(path,"/Data/data_out_auz/AmazonACIdata_f.csv",sep=""))

#to get ambient Ci values (every 4th observation in each curve)
amz.2<-subset(amz.1,amz.1$Ci<270 & amz.1$Ci>175)
with(amz.2,plot(Tleaf,Photo,col=Season,cex=2,pch=16))

amz.3<-split(amz.2,amz.2$Season)

amz.topts.1<-get_topts(lapply(amz.3[1],function(x)fit.nlme(x,yvar="Photo",random="Species"))) 
amz.topts.2<-get_topts(lapply(amz.3[2],FUN=fitquad))#mix model didin't converged

amz.topts<-rbind(amz.topts.1,amz.topts.2)
amz.topts$Season<-names(amz.3)
amz.topts$Species<-"Amazon spp"

#amz.topts$Type<-"Tropical_RF"
amz.topts$DataSet<-"RF_AMZ"

amz.topts$Tavg_30<-NA
amz.topts$Tavg_30[which(amz.topts$Season=="dry")]<-27.8 #average temperature during measurement months as data collected over several months
amz.topts$Tavg_30[which(amz.topts$Season=="wet")]<-25.8

amz.topts$Tmax_30<-NA
amz.topts$Tmax_30[which(amz.topts$Season=="dry")]<-32.14 #average temperature during measurement months as data collected over several months
amz.topts$Tmax_30[which(amz.topts$Season=="wet")]<-29.26


#No measurements at 25C Tmin_30
path<-"\\\\ad.uws.edu.au/dfshare/HomesHWK$/90931217/My Documents/Documents/Dushan/Repos/PhotoM/data_with_gs/Amazon/csv"
amz.gs.dat<-read_and_convert_licor(path=path)
amz.gs.dat$Tempfrac<-cut(amz.gs.dat$Tleaf,breaks=c(25,28,32,36,40,45),labels=FALSE)
amz.topts$Cond.mean<-with(subset(amz.gs.dat,amz.gs.dat$Tempfrac==1),mean(Cond))
amz.topts$Photo.mean<-with(subset(amz.gs.dat,amz.gs.dat$Tempfrac==1),mean(Photo))
path<-getwd()

#amz.gs.dat$xterm<-with(amz.gs.dat,Photo/CO2S*sqrt(VpdL))
#with(amz.gs.dat,plot(Cond~xterm))

#get g1
amz_g1<-fitStom(amz.gs.dat)

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#Tropical species in Rwanda (AV)

avr<-read.csv(paste(path,"/Data/data_out_auz/Angelica_Varhammar_tropical_species_with_gs.csv",sep=""))

#to get ambient Ci values (every 4th observation in each curve)
avr.1<-do.call(rbind,by(avr,avr$Curve,function(x) x[5,]))

avr.1<-subset(avr.1,avr.1$Ci<420) #remove one Ci value ~700

avr.2<-split(avr.1,avr.1$Species)
avr.topts<-get_topts(lapply(avr.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))

avr.topts$Species<-names(avr.2)
#avr.topts$Season<-"dry"
#avr.topts$Temp_Treatment<-"Ambient"
#avr.topts$Facility<-"Field"
avr.topts$SSloc<--2.6

avr.topts$Tavg_30<-19.38 # No met data available for 2011. However, as climate is stable, I use preveous year data (average T for July and August)
avr.topts$Tmax_30<-24.28

#avr.topts$Type<-"Tropical_Montane_RF"
avr.topts$DataSet<-"RF_MON"
#with(avr.1,plot(Tleaf,Photo,col=Species,cex=2,pch=16))


#get mean photosynthesis at 25C
avr.1$TargetTemp<-cut(avr.1$Tleaf,breaks=c(19,23.5,26.5,30,40),labels=FALSE)

avr.sum<-summaryBy(Photo+Cond+VpdL~Species+TargetTemp,FUN=c(mean,std.error),data=avr.1)
avr.sum<-subset(avr.sum,avr.sum$TargetTemp==2)
avr.sum[2]<-NULL
avr.topts<-merge(avr.topts,avr.sum,by=c("Species"))

#add gs at 25C
#Source Varhammar et al, 2015 New Phytologist
#avr.topts$Cond.mean<-NA
#avr.topts$Cond.mean[which(avr.topts$Species=="C_g")]<-0.09
#avr.topts$Cond.mean[which(avr.topts$Species=="C_s")]<-0.15
#avr.topts$Cond.mean[which(avr.topts$Species=="E_e")]<-0.04
#avr.topts$Cond.mean[which(avr.topts$Species=="E_i")]<-0.41
#avr.topts$Cond.mean[which(avr.topts$Species=="E_m")]<-0.66
#avr.topts$Cond.mean[which(avr.topts$Species=="H_a")]<-0.36



#get complete species name
avr.topts$Species[which(avr.topts$Species=="C.g.")]<-"Carapa procera"
avr.topts$Species[which(avr.topts$Species=="C.s.")]<-"Toona sinensis"
avr.topts$Species[which(avr.topts$Species=="E.e.")]<-"Entandrophragma excelsum"
avr.topts$Species[which(avr.topts$Species=="E.i.")]<-"Eucalyptus globulus" 
avr.topts$Species[which(avr.topts$Species=="E.m.")]<-"Eucalyptus microcorys"
avr.topts$Species[which(avr.topts$Species=="H.a.")]<-"Hagenia abyssinica"

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
with(wang.1,plot(Tleaf,Photo,col=Species,cex=2,pch=16,main="Finland data"))
with(subset(wang.1,wang.1$Species=="Pinus sylvestris"),plot(Tleaf,Photo,col="black"))
#wang et al: two species (Ci not sure)

wang<-read.csv(paste(path,"/Data/data_out_auz/Betula_pendula_Pinus_sylvestris_wang_etal.csv",sep=""))
wang.1<-do.call(rbind,by(wang,wang$Curve,function(x) x[5,]))

wang.1<-subset(wang.1,wang.1$Ci<400) #to remove 4 points which are not belongs to ambient CO2 levels

wang.2<-split(wang.1,wang.1$Species)

wang.topts<-get_topts(lapply(wang.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
wang.topts$Species<-as.factor(names(wang.2))

wang.topts$DataSet<-"WANG_ET_AL"

#add Photosynthesis at 25C
wang.1$TargetTemp<-cut(wang.1$Tleaf,breaks=c(5,22.5,27,40),labels=FALSE)
wang.sum<-summaryBy(Photo+VPD~Species+TargetTemp,FUN=c(mean,std.error),data=wang.1)
wang.sum<-subset(wang.sum,wang.sum$TargetTemp==2)
wang.sum[2]<-NULL;names(wang.sum)[c(3,5)]<-c("VpdL.mean","VpdL.std.error")

wang.topts<-merge(wang.topts,wang.sum,by="Species",all=TRUE)

wang.topts$Tavg_30<-14

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# Chamaecyparis obtusa

cha<-read.csv(paste(path,"/Data/data_out_auz/Chamaecyparis obtusa_han_etal.csv",sep=""))
colnames(cha)[7]<-"VpdL"

#to get ambient Ci values (every 4th observation)
cha.1<-do.call(rbind,by(cha,cha$Curve,function(x) x[6,]))

cha.2<-split(cha.1,cha.1$Season)
cha.topts<-get_topts(lapply(cha.2,FUN=fitquad)) #No data to fit random intercept


cha.topts$Species<-"Chamaecyparis obtusa"
cha.topts$Season<-as.factor(names(cha.2))

cha.1$Tleaffrac<-cut(cha.1$Tleaf,breaks=c(14,19,22,26,30),labels=FALSE)

chasum<-summaryBy(Ci/CO2S+Photo+Cond+VpdL+Ci+Tleaf~Season+Tleaffrac,data=cha.1,FUN=c(mean))
chasum<-subset(chasum,chasum$Tleaffrac==3)

cha.topts<-merge(cha.topts,chasum,by="Season")

#cha.topts$Type<-"Evergreen_Gymnosperm"
#cha.topts$Type<-"Temperate"
#colnames(cha.topts)[c(1)]<-"Month"
cha.topts$DataSet<-"HAN_ET_AL_A"
cha.topts$Tavg_30<-NA 
cha.topts$Tavg_30[which(cha.topts$Season=="Apr")]<-12.0  #warmest 17.2
cha.topts$Tavg_30[which(cha.topts$Season=="Oct")]<-15.2  


#get g1

cha.1$xterm<-with(cha.1,Photo/CO2S*sqrt(VpdL))
with(cha.1,plot(Cond~xterm,col=Season,cex=2))

cha_g1<-get_topts(lapply(cha.2,FUN=fitStom))
cha_g1$Season<-names(cha.2)

cha.topts<-merge(cha.topts,cha_g1,by="Season")

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#with(dil.1,plot(Tleaf,Photo,cex=2,col=Species,pch=c(16,15,1,5)[Season]))
#for(i in 1:length(dil.2)){
 # toplot<-dil.2[[i]]
 # name<-names(dil.2)[i]
 # plot(Photo~Tleaf,data=toplot,main=name,pch=16,cex=2,col=Leaf)
#}
#with(subset(dil.1,dil.1$Species=="Betulapapyrifera"),plot(Tleaf,Photo,cex=2,col=Season,pch=16))
#with(subset(dil.1,dil.1$Species=="Populustremuloides"),plot(Tleaf,Photo,cex=2,col=Season,pch=16))
#with(subset(dil.1,dil.1$Species=="Liquidambarstyraciflua"),plot(Tleaf,Photo,cex=2,col=Season,pch=16))
#with(subset(dil.1,dil.1$Species=="Populusdeltoides"),plot(Tleaf,Photo,cex=2,col=Season,pch=16))


# Dillaway et al 

dil<-read.csv(paste(path,"/Data/data_out_auz/dillaway_etal.csv",sep=""))
dil.1<- do.call(rbind,by(dil,dil$Curve,function(x) x[4,]))
dil.1<-subset(dil.1,!Curve %in% c(184,170,134,121:126,162)) #remove few outliers. Higher photo compared to ther data 

dil.1$Species<-factor(str_replace_all(dil.1$Species, fixed(" "), ""))
dil.1$Season<-factor(str_replace_all(dil.1$Season, fixed(" "), ""))

dil.2<-split(dil.1,paste(dil.1$Species,dil.1$Season))

#list1<-dil.2[3]
#dil.topts.1<-get_topts(lapply(list1,FUN=fitquad)) #cannot fit mix models


#list2<-dil.2[c(1:2,4:12)]
dil.topts<-get_topts(lapply(dil.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))


#dil.topts.1$Species<-get_names(list1)[,2]
dil.topts$Species<-get_names(dil.2)[,1]

#dil.topts.1$Season<-get_names(dil.2)[,4]
dil.topts$Season<-get_names(dil.2)[,2]

#dil.topts<-rbind(dil.topts.1,dil.topts.2)

#get mean Cond
dil.1$Tleaffrac<-cut(dil.1$Tleaf,breaks=c(20,24,27,30),labels=FALSE)

dilsum<-summaryBy(Photo+Cond~Season+Species+Tleaffrac,keep.names=TRUE,data=dil.1,FUN=c(mean,std.error))
dilsum<-subset(dilsum,dilsum$Tleaffrac==2)

dil.topts<-merge(dil.topts,dilsum,by=c("Season","Species"))

#add Tgrowth
#data source: Dillaway and Kruger PCE 2010
#averages are for the entire study period. not fro preciding 30 days

dil.topts$Tavg_30<-NA
dil.topts$Tavg_30[which(dil.topts$Season=="Illinois")]<-30.6
dil.topts$Tavg_30[which(dil.topts$Season=="northernWisconsin")]<-18.7
dil.topts$Tavg_30[which(dil.topts$Season=="southernWisconsin")]<-22.5
#fu<-subset(dil.1,dil.1$Species=="Populusdeltoides" & dil.1$Season=="northernWisconsin")
#fu<-subset(dil.1,dil.1$Species=="Populustremuloides" & dil.1$Season=="northernWisconsin")
#with(fu,plot(Tleaf,Photo,cex=2))


dil.topts$Species[which(dil.topts$Species=="Betulapapyrifera")]<-"Betula papyrifera"
dil.topts$Species[which(dil.topts$Species=="Liquidambarstyraciflua")]<-"Liquidambar styraciflua"
dil.topts$Species[which(dil.topts$Species=="Populusdeltoides")]<-"Populus deltoides"
dil.topts$Species[which(dil.topts$Species=="Populustremuloides")]<-"Populus tremuloides"

dil.topts$DataSet<-"DILLAWAY"
#dil.topts$DataSet[which(dil.topts$Season=="Illinois")]<-"DILLAWAY_1"
#dil.topts$DataSet[which(dil.topts$Season=="northernWisconsin")]<-"DILLAWAY_2"
#dil.topts$DataSet[which(dil.topts$Season=="southernWisconsin")]<-"DILLAWAY_3"

dil.topts<-subset(dil.topts,dil.topts$b>0)

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# Dreyer et al 
#for(i in 1:length(dre.2)){
 # toplot<-dre.2[[i]]
  #name<-names(dre.2)[i]
#  plot(Photo~Tleaf,data=toplot,main=name)
#}

dre<-read.csv(paste(path,"/Data/data_out_auz/DreyerSevenSpp_final.csv",sep=""))
#dre$Species<-str_replace_all(dre$Species, fixed("_"), " ") 

dre.1<- subset(dre,dre$Ci>=275 & dre$Ci<=325) #get ambient Ci levels

dre.2<-split(dre.1,paste(dre.1$Species))

dre.topts<-get_topts(lapply(dre.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
dre.topts$Species<-as.factor(names(dre.2))

#No Cond/VPD data 
dre.1$Tleaffrac<-cut(dre.1$Tleaf,breaks=c(20,24,26,30),labels=FALSE)

dresum<-summaryBy(Photo+Ci~Species+Tleaffrac,keep.names=TRUE,data=dre.1,FUN=c(mean,std.error))
dresum<-subset(dresum,dresum$Tleaffrac==2)

dre.topts<-merge(dre.topts,dresum,by="Species")
#dre.topts$Season<-"Summer"

dre.topts$Tavg_30<-17.3 #source: Table 4, Dreyer et al 2001, Tree Physiology
dre.topts$Tmax_30<-22.4
#dre.topts$Type<-"Deciduous_Angiosperm"
#dre.topts$Type<-"Temperate"
#dre.1$Species<-as.factor(dre.1$Species)
#subset(dre.1,dre.1$Species=="Acer pseudoplatanus")
#with(subset(dre.1,dre.1$Species=="Quercus robur ED"),plot(Tleaf,Photo,col=Species,pch=16))

dre.topts$DataSet<-"DREYER"

dre.topts$Species<-factor(str_replace_all(dre.topts$Species, fixed("_"), " "))
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#for(i in 1:length(med.2)){
  #toplot<-med.2[[i]]
  #name<-names(med.2)[i]
  #plot(Photo~Tleaf,data=toplot,main=name)
#}

# Medlyn et al
med<-read.csv(paste(path,"/Data/data_out_auz/Medlyn_etal_all.csv",sep=""))
med.1<-subset(med,med$Ci>=175 & med$Ci<=350)
with(med.1,plot(Tleaf,Photo))
med.1<-subset(med.1,med.1$Curve!=62) #remove one outlier

med.2<-split(med.1,paste(med.1$Species,med.1$Season))
med.2[c(4,10,11)]<-NULL #poor fits. remove from further analysis

med.topts<-get_topts(lapply(med.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
med.topts$Species<-get_names(med.2)[,1]
med.topts$Season<-get_names(med.2)[,2]


#med.topts$Type<-"Evergreen_Gymnosperm"
#med.topts$Type<-"Temperate"
med.topts$DataSet<-"MEDLYN"
#path="//ad.uws.edu.au/dfshare/HomesHWK$/90925395/My Documents/Repos/PhotoM"

#
#get g1

#add met data

med_met<-read.csv(paste0(path,"/MET_DATA/MET_DATA/medlyn_pinus_pinaster_met.csv"))
med.topts<-merge(med.topts,med_met,by="Season")

#add photosynthesis at 25C

#med.1$Tleaffrac<-cut(med.1$Tleaf,breaks=c(20,24,27,30),labels=FALSE)
#med.sum<-summaryBy(Photo~Species+Season+Tleaffrac,FUN=c(mean,std.error),data=med.1)
#med.sum<-subset(med.sum,med.sum$Tleaffrac==2)
#med.topts<-merge(med.topts,med.sum,by=c("Species","Season"))


med.topts$Species[which(med.topts$Species=="Landes")]<-"Pinus pinaster_L"
med.topts$Species[which(med.topts$Species=="Tamjoute")]<-"Pinus pinaster_T"

#add gs estimates
gs_dat<-read.csv(paste0(path,"/Data/data_out_auz/Medlyn_gs_data_processed.csv"))
gs_dat<-subset(gs_dat,!is.na(gs_dat$Season))
gs_dat$VpdL<-with(gs_dat,VpdL/10)
gs_dat$Cond<-with(gs_dat,Cond/10^3)

gs_dat<-subset(gs_dat,gs_dat$CO2S>300 & gs_dat$CO2S<400) #get ambient CO2 levels

gs_dat$Tempfrac<-cut(gs_dat$Tleaf,breaks=c(20,25,27.5,36,40,45),labels=FALSE)

gs_sum<-summaryBy(Photo+Cond~Season+Species+Tempfrac,keep.names=TRUE,data=gs_dat,FUN=c(mean,std.error))

gs_sum<-subset(gs_sum,gs_sum$Tempfrac==2)

gs_sum<-gs_sum[-c(2,3)]
gs_sum$Species<-factor(c("Pinus pinaster_L","Pinus pinaster_T","Pinus pinaster_L","Pinus pinaster_T","Pinus pinaster_L","Pinus pinaster_T"))

med.topts<-merge(med.topts,gs_sum,by=c("Species","Season"),all=T)
med.topts<-subset(med.topts,!is.na(med.topts$topt))
#med.topts$Cond.mean=with(med.topts,Cond.mean/10^3)
#med.topts<-subset(med.topts, ! Season %in% c("Jul")) #no sufficient data points.

med.split<-split(gs_dat,paste(gs_dat$Species,gs_dat$Season))
med_g1<-get_topts(lapply(med.split,FUN=fitStom))

#gs_dat$xterm<-with(gs_dat,Photo/CO2S*sqrt(VpdL))
#with(gs_dat,plot(Cond~xterm))

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# Picea_mariana_way_etal

picea<-read.csv(paste(path,"/Data/data_out_auz/Picea_mariana_way_etal.csv",sep=""))
picea.1<-picea[firstobs(~Curve,data=picea),]
picea.2<-split(picea.1,picea.1$Season)
with(picea.1,plot(Tleaf,Photo,pch=c(15,16)[Season]))


picea.topts<-get_topts(lapply(picea.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
picea.topts$Temp_Treatment<-names(picea.2)
picea.topts$Species<-"Picea mariana"

#picea.topts$Type<-"Evergreen_Gymnosperm"
#picea.topts$Type<-"Boreal"

picea.topts$Tavg_30<-NA
picea.topts$Tavg_30[which(picea.topts$Temp_Treatment=="Cool")]<-18 #source: Way et al 2008, PCE materials and methods
picea.topts$Tavg_30[which(picea.topts$Temp_Treatment=="warm")]<-26

picea.topts$Tmax_30<-NA
picea.topts$Tmax_30[which(picea.topts$Temp_Treatment=="Cool")]<-22 #source: Way et al 2008, PCE materials and methods
picea.topts$Tmax_30[which(picea.topts$Temp_Treatment=="warm")]<-30


picea.topts$DataSet<-"WAY_ET_AL"
#with(subset(picea.1,picea.1$Season=="warm"),plot(Tleaf,Photo,col=Leaf,cex=2,pch=16))
#No data avilable around 25C
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

#Onoda et al (poor fit, drop from furtehr analysis???)
fagus<-read.csv(paste(path,"/Data/data_out_auz/Onoda_etal.csv",sep=""))
fagus.1<-subset(fagus,fagus$GrowthCa=="amb")
fagus.2<-do.call(rbind,by(fagus.1,fagus.1$Curve,function(x) x[3,]))

fagus.2<-subset(fagus.2,fagus.2$Photo<8.2) #Remove some outliers

with(fagus.2,plot(Tleaf,Photo,col=Season,pch=16))
fagus.3<-split(fagus.2,paste(fagus.2$Season))
#fagus.topts<-get_topts(lapply(fagus.3,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
fagus.topts<-get_topts(lapply(fagus.3,FUN=fitquad))

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# Picea_abies_tarvainen_etal
piceaal<-read.csv(paste(path,"/Data/data_out_auz/Pieca_abies_tarvainen_etal.csv",sep=""))
piceaal.1<-subset(piceaal,abs(piceaal$CO2R-400)<5)
with(piceaal.1,plot(Tleaf,Photo,col=Season,pch=16,cex=2,main="Sweden data"))

piceaal.2<-split(piceaal.1,piceaal.1$Season)
piceaal.topts<-get_topts(lapply(piceaal.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
piceaal.topts$Species<-"Picea abies"
piceaal.topts$Season<-names(piceaal.2)

picsum<-summaryBy(Photo+Cond+Ci+Ci/CO2S+VpdL~Season+TargetT,keep.names=TRUE,data=piceaal.1,FUN=c(mean,std.error))
picsum<-subset(picsum,picsum$TargetT==25)

piceaal.topts<-merge(piceaal.topts,picsum,by="Season")
#piceaal.topts$Type<-"Evergreen_Gymnosperm"
#piceaal.topts$Type<-"Boreal"

piceaal.topts$DataSet<-"TARVAINEN"
piceaal.topts<-piceaal.topts[-c(12)]

piceaal.topts$Tavg_30<-6.9 #source: Tarvainen et al 2013 Oecologia
piceaal.topts$Tmax_30<-15.1 #this is monthly average for May for the closest met station: Source https://www.yr.no/place/Sweden/V%C3%A4stra_G%C3%B6taland/V%C3%A4nersborg/statistics.html

piceaal.topts<-subset(piceaal.topts,piceaal.topts$Season=="Year3")
piceaal.topts<-piceaal.topts[-1]
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

# Pinus_densiflora_han_etal (high standard errors in fits for NOV)

pden<-read.csv(paste(path,"/Data/data_out_auz/Pinus_densiflora_han_etal.csv",sep=""))
colnames(pden)[7]<-"VpdL"

pden.1<-subset(pden,abs(pden$CO2S-340)<10)


with(pden.1,plot(Tleaf,Photo,col=Season,pch=16,cex=2))
pden.2<-split(pden.1,pden.1$Season)
pden.topts<-get_topts(lapply(pden.2[c(2,3)],FUN=fitquad)) #Data not sufficient to fit mix models
pden.topts$Season<-names(pden.2[c(2,3)])
pden.topts$Species<-"Pinus densiflora"

pden.1$Tleaffrac<-cut(pden.1$Tleaf,breaks=c(20,24,26,30),labels=FALSE)
pdensum<-summaryBy(Photo+Cond+Ci+Ci/CO2S+VpdL~Season+Tleaffrac,keep.names=TRUE,data=pden.1,FUN=c(mean,std.error))
pdensum<-subset(pdensum,pdensum$Tleaffrac==2 & pdensum$Season!="July")

pden.topts<-merge(pden.topts,pdensum,by="Season")

#pden.topts$Type<-"Evergreen_Gymnosperm"
#pden.topts$Type<-"Temperate"

pden.topts$DataSet<-"HAN_2"
#pden.topts$Tavg_30<-15 #From Yan Shih Lin's data

pden.topts$Tavg_30<-NA
pden.topts$Tavg_30[which(pden.topts$Season=="May")]<-11.9 #source: calculated based on 1 degree resolution data for the site
pden.topts$Tavg_30[which(pden.topts$Season=="Nov")]<-14.68


#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------

# Pinus_radiata_Walcroft
prad<-read.csv(paste(path,"/Data/data_out_auz/Pinus_radiata_Walcroft.csv",sep=""))
colnames(prad)[10]<-"VpdL"

prad.1<-subset(prad,prad$CO2S>330 & prad$CO2S<400)
with(prad.1,plot(Tleaf,Photo,col=Season,pch=16,cex=2))

prad.2<-split(prad.1,prad.1$Season)
prad.topts<-get_topts(lapply(prad.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
prad.topts$Species<-"Pinus radiata"


prad.1$Tleaffrac<-cut(prad.1$Tleaf,breaks=c(20,24,26,30),labels=FALSE)
pradnsum<-summaryBy(Photo+Cond+Ci+Ci/CO2S+VpdL~Tleaffrac,keep.names=TRUE,data=prad.1,FUN=c(mean,std.error))
pradnsum<-subset(pradnsum,pradnsum$Tleaffrac==2)

prad.topts<-merge(prad.topts,pradnsum)
#prad.topts$Season<-"No Season"

#prad.topts$Type<-"Evergreen_Gymnosperm"
#prad.topts$Type<-"Temperate"


prad.topts$DataSet<-"WALCROFT"

prad.topts$Tavg_30<-23.7 #walcroft et al, 1997 PCE: this is the average temperature in glasshouse
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# Pinus_taeda_Ellsworth


ptad<-read.csv(paste(path,"/Data/data_out_auz/Pteada_Ellsworth.csv",sep=""))
colnames(ptad)[10]<-"VpdL"

ptad.1<-subset(ptad,ptad$CO2S>340 & ptad$CO2S<370)
with(ptad.1,plot(Tleaf,Photo,col=Season,pch=16,cex=2))
ptad.2<-split(ptad.1,ptad.1$Season)
#ptad.2<-split(ptad.1,ptad.1$Season)
ptad.topts<-get_topts(lapply(ptad.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
ptad.topts$Season<-names(ptad.2)
ptad.topts$Species<-"Pinus taeda"


ptad.1$Tleaffrac<-cut(ptad.1$Tleaf,breaks=c(20,26,28,30),labels=FALSE)
ptadsum<-summaryBy(Photo+Cond+Ci+Ci/CO2S+VpdL~Season+Tleaffrac,keep.names=TRUE,data=ptad.1,FUN=c(mean,std.error))
ptadsum<-subset(ptadsum,ptadsum$Tleaffrac==2)

ptad.topts<-merge(ptad.topts,ptadsum,by="Season")

#ptad.topts$Type<-"Evergreen_Gymnosperm"
#ptad.topts$Type<-"Temperate"

#Add met data paste0(path,"/MET_DATA/DUKE")
ptad.topts$Tavg_30<-NA
ptad.topts$Tavg_30[which(ptad.topts$Season=="summer")]<-getavgduke(path=paste0(path,"/MET_DATA/DUKE"),fname="duke_temp_data.csv",from="1999-8-1")[[1]]
ptad.topts$Tavg_30[which(ptad.topts$Season=="winter")]<-getavgduke(path=paste0(path,"/MET_DATA/DUKE"),fname="duke_temp_data.csv",from="1999-12-1")[[1]]

ptad.topts$Tmax_30<-NA
ptad.topts$Tmax_30[which(ptad.topts$Season=="summer")]<-getavgduke(path=paste0(path,"/MET_DATA/DUKE"),fname="duke_temp_data.csv",from="1999-8-1")[[2]]
ptad.topts$Tmax_30[which(ptad.topts$Season=="winter")]<-getavgduke(path=paste0(path,"/MET_DATA/DUKE"),fname="duke_temp_data.csv",from="1999-12-1")[[2]]

ptad.topts<-ptad.topts[-c(12)]
ptad.topts$DataSet<-"ELLSWORTH"


#get g1
ptad.1$xterm<-with(ptad.1,Photo/CO2S*sqrt(VpdL))
with(ptad.1,plot(Cond~xterm,col=Season,cex=2,pch=16))


ells_g1<-get_topts(lapply(ptad.2,FUN=fitStom))
ells_g1$Season<-names(ptad.2)

ptad.topts<-merge(ptad.topts,ells_g1,by="Season")

#An<-nls(Photo~Aopt-(b*(Tleaf-Topt)^2),data=subset(ptad.1,ptad.1$Photo>2.4 ),start=list(Aopt=3.5,Topt=25,b=0.05))
#plot_nls(An)
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#Strassmayer et al
pstras<-read.csv(paste(path,"/Data/data_out_auz/Strassemeyer.csv",sep=""))
pstras.1<-do.call(rbind,by(pstras,pstras$Curve,function(x) x[3,]))
pstras.2<-split(pstras.1,pstras.1$Species)

pstras.topts<-get_topts(lapply(pstras.2,function(x)fit.nlme(x,yvar="Photo",random="Leaf")))
pstras.topts$Species<-names(pstras.2)

pstras.1$Tleaffrac<-cut(pstras.1$Tleaf,breaks=c(20,24,27.5,30),labels=FALSE)
pstrassum<-summaryBy(Photo+Cond+Ci+Ci/CO2S+VpdL~Species+Tleaffrac,keep.names=TRUE,data=pstras.1,FUN=c(mean,std.error))
pstrassum<-subset(pstrassum,pstrassum$Tleaffrac==2)
pstras.topts<-merge(pstras.topts,pstrassum,by="Species")

pstras.topts$DataSet<-"STRASSEMEYER_ET_AL"

pstras.topts$Tavg_30<-20
with(pstras.1,plot(Tleaf,Photo,col=Species))
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#Arctic species

#---------------------------------------------------------------------------------------------
#to add curve number
#dat_arc<-read.csv(paste(path,"/Data/data_out_auz/Arctic_A-Ci_curves_2012-2015.csv",sep=""))
#dat_arc$Curve<-c(1)
#Curve <- c()
#count <- 1  

#for (i in 2:length(dat_arc$Sample_ID)){
  
  #ifelse(dat_arc$Sample_ID[i-1]==dat_arc$Sample_ID[i],count<-count,count <- count + 1)
  #Curve[i] <- count 
  
  #}
#dat_arc$Curve[2:length(dat_arc$Sample_ID)] <- na.omit(Curve)

#dat_arc$CO2frac<-cut(dat_arc$CO2R,breaks=c(35,46,60,80,105,130,160,210,310,410,610,910,1310,1510,1750,2100),
                # labels=FALSE)

#dat_arc_summary<-summaryBy(.~Curve+USDA_Species_Code+CO2frac,data=dat_arc,FUN=c(mean),keep.names=T)
#write.csv(dat_arc_summary,paste0(path,"/Data/data_out_auz/Arctic_A-Ci_curves_2012-2015_V2.csv"),row.names=FALSE,sep=",")

#--------------------------------------------------------------------------------------------

dat_arc<-read.csv(paste(path,"/Data/data_out_auz/Arctic_A-Ci_curves_2012-2015_V2.csv",sep=""))
dat_arc.a<-subset(dat_arc,dat_arc$CO2S>370 & dat_arc$CO2S<410)  #get abmient CO2 levels

with(dat_arc, table(USDA_Species_Code, Year))

#fit mix models using year as a random variable

arc.topts<-data.frame(do.call(rbind,list(fit.nlme(dat_arc.a,yvar="Photo",rand="Year"))))
arc.topts$DataSet<-"ARCTIC"
arc.topts$Species<-"Arctic spp"


#dat_arc.b<-split(dat_arc.a,dat_arc.a$USDA_Species_Code)
#get_topts(lapply(dat_arc.b,FUN=fitquad))

#windows()
#par(mfrow=c(4,2))

#for(i in 1:length(dat_arc.b)){
#toplot<-dat_arc.b[[i]]
#name<-names(dat_arc.b)[i]
#plot(Photo~Tleaf,data=toplot,main=name,pch=16)
#}

#fit.nlme(dat_arc.a,yvar="Photo",random="Year")

#with(subset(dat_arc.a,dat_arc.a$Year==2012),plot(Tleaf,Photo,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data 2012"))
#with(subset(dat_arc.a,dat_arc.a$Year==2013),plot(Tleaf,Photo,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data 2013"))
#with(subset(dat_arc.a,dat_arc.a$Year==2014),plot(Tleaf,Photo,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data 2014"))
#with(subset(dat_arc.a,dat_arc.a$Year==2015),plot(Tleaf,Photo,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data 2015"))


#with(dat_arc.a,plot(Tleaf,Cond,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data"))
#with(dat_arc.a,plot(VpdL,Cond,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data"))
#with(dat_arc.a,plot(Tleaf,Ci,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data"))


#with(dat_arc.a,plot(Sample_Date,Tleaf,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data"))

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#Spruce Anna et al 2015

spr<-read.csv(paste(path,"/Data/data_out_auz/SPRUCE_3_cohort_ACi_data.csv",sep=""))
spr.a<-subset(spr,spr$CO2R>395 & spr$CO2R<410 & spr$Month %in% c(4,8,10))  #get abmient CO2 levels and remove 4 outliers

spr.b<-summaryBy(.~Year+Month+Cohort_age+Curve,data=spr.a,FUN=mean,keep.names=T) #get averages across multiple measurements

spr.c<-split(spr.b,paste(spr.b$Month))

spr.topts<-get_topts(lapply(spr.c,function(x)fit.nlme(x,yvar="Photo",random="Cohort_age")))
spr.topts$Season<-names(spr.c)
spr.topts$DataSet<-"ANNA"
spr.topts$Species<-"Picea mariana"


with(spr.b,plot(Tleaf,Photo,col=Month,cex=2,pch=c(1,16)[Cohort_age]))
with(spr.b,plot(Tleaf,Cond,col=Month))

with(subset(spr.a,spr.a$Month==4),plot(Tleaf,Photo))
with(subset(spr.a,spr.a$Month==8),plot(Tleaf,Photo,col=Cohort_age,pch=16,cex=2))

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------


#Walcroft et al, Peach

pea<-read.csv(paste(path,"/Data/data_out_auz/walcroft peach.csv",sep=""))
pea.1<-pea[firstobs(~Curve,data=pea),]
pea.topts<-fitquad(pea.1)


#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------


#combine all dataframes together
tfits<-data.frame(rbind.fill(hfe.topts,wtc3.topts,wtc2.topt,wtc1.topts,wtc4.topts,e.ss.topts,e.tr.topts,rf.topts,e.del.topt,gww.topts,e.glob2.fit.para,cor.topts,epau.topts,gr.topts,eface.topts))
names(tfits)[1:9] <- c("Species","Season","Aopt","Topt","b","Aopt.se","Topt.se","b.se","R2")
tfits<-tfits[, !(colnames(tfits) %in% c("W","V9","Tleafnew","Tleaf.mean","Tleaf.std.error","TargetTemp","Facility" ,"Cir","Tleaffrac","Tleaf","Tempfrac","Cir.std.error" ))] 
#names(tfits)[22]<-"Temp_Treatment"
tfits["Season"][is.na(tfits["Season"])] <-"No_Seasonal_Measurements"
tfits["Temp_Treatment"][is.na(tfits["Temp_Treatment"])] <-"No_Temperature_Treatments"

#get seasons 
tfits$Season_New<-NA
tfits$Season_New[which(tfits$Season=="April")]<-"Autumn"
tfits$Season_New[which(tfits$Season=="August")]<-"Winter"
tfits$Season_New[which(tfits$Season=="December")]<-"Summer"
tfits$Season_New[which(tfits$Season=="February")]<-"Summer"
tfits$Season_New[which(tfits$Season=="January")]<-"Summer"
tfits$Season_New[which(tfits$Season=="May")]<-"Autumn"
tfits$Season_New[which(tfits$Season=="November")]<-"Spring"
tfits$Season_New[which(tfits$Season=="October")]<-"Spring"
tfits$Season_New[which(tfits$Season=="September")]<-"Winter"
tfits$Season_New[which(tfits$Season=="Spring")]<-"Spring"
tfits$Season_New[which(tfits$Season=="SUMMER")]<-"Summer"
tfits$Season_New[which(tfits$Season=="Winter")]<-"Winter"
tfits$Season_New[which(tfits$Season=="WINTER")]<-"Winter"
tfits$Season_New[which(tfits$Season=="No_Seasonal_Measurements")]<-"No_Seasonal_Measurements"
tfits$Season_New<-as.factor(tfits$Season_New)



tfits_oaz<-data.frame(rbind.fill(spr.topts,arc.topts,amz.topts,avr.topts,wang.topts,cha.topts,dil.topts,dre.topts,med.topts,picea.topts,piceaal.topts,pden.topts,prad.topts,ptad.topts,pstras.topts,panama.topts))
names(tfits_oaz)[1:9] <- c("Aopt","Topt","b","Aopt.se","Topt.se","b.se","R2","S","pvalue")
#tfits_oaz<-subset(tfits_oaz,tfits_oaz$Topt.se<10)
tfits_oaz["Season"][is.na(tfits_oaz["Season"])] <-"No_Seasonal_Measurements"
tfits_oaz["Temp_Treatment"][is.na(tfits_oaz["Temp_Treatment"])] <-"No_Temperature_Treatments"

tfits_oaz$Season_New<-NA
tfits_oaz$Season_New[which(tfits_oaz$Season=="Autumn")]<-"Autumn"
tfits_oaz$Season_New[which(tfits_oaz$Season=="winter")]<-"Winter"
tfits_oaz$Season_New[which(tfits_oaz$Season=="summer")]<-"Summer"
tfits_oaz$Season_New[which(tfits_oaz$Season=="Spring")]<-"Spring"

tfits_oaz$Season_New[which(tfits_oaz$Season=="Jan")]<-"Winter"
tfits_oaz$Season_New[which(tfits_oaz$Season=="Apr")]<-"Spring"
tfits_oaz$Season_New[which(tfits_oaz$Season=="Jun")]<-"Summer"
tfits_oaz$Season_New[which(tfits_oaz$Season=="Jul")]<-"Summer"
tfits_oaz$Season_New[which(tfits_oaz$Season=="AugSept")]<-"Summer"
tfits_oaz$Season_New[which(tfits_oaz$Season=="Oct")]<-"Autumn"
tfits_oaz$Season_New[which(tfits_oaz$Season=="Nov")]<-"Autumn"
tfits_oaz$Season_New[which(tfits_oaz$Season=="Mar")]<-"Spring"
tfits_oaz$Season_New[which(tfits_oaz$Season=="dry")]<-"Dry"
tfits_oaz$Season_New[which(tfits_oaz$Season=="wet")]<-"Wet"
tfits_oaz$Season_New[which(tfits_oaz$Season=="May")]<-"Spring"
tfits_oaz$Season_New[which(tfits_oaz$Season==10)]<-"Autumn"
tfits_oaz$Season_New[which(tfits_oaz$Season==4)]<-"Spring"
tfits_oaz$Season_New[which(tfits_oaz$Season==8)]<-"Summer"
tfits_oaz$Season_New[which(tfits_oaz$Season=="No_Seasonal_Measurements")]<-"No_Seasonal_Measurements"


tfits_oaz$Season_New<-as.factor(tfits_oaz$Season_New)

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#to get a dataframe with estimated parameters for temperature response of net photosynthesis
#this dataframe contains fitte parameters of quadratic model, average gs,VPD,ci,Ci/Ca at 25C and their standard errors
#and mean growth temperature of preceding 30 days 

tfits_photo_met<-rbind.fill(tfits_oaz,tfits)
#tfits_photo_met<-tfits_photo_met[-c()] #to remove unwanted columns

#tfits_photo_met["Temp_Treatment"][is.na(tfits_photo_met["Temp_Treatment"])] <-"Ambient"

#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#Remove estimates with poor fits
#tfits_photo_met<-subset(tfits_photo_met,tfits_photo_met$b>0) #to remove negative b (two points)
#tfits_photo_met<-subset(tfits_photo_met,tfits_photo_met$Topt.se<10) #to remove Topt with se=80 (one point)
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
dat_info<-read.csv("dataset_seedsource_info.csv")
#tfits_photo_met<-merge(tfits_photo_met,dat_info,by=c("DataSet","Species"),all=TRUE)
#tfits_photo_met<-subset(tfits_photo_met,!is.na(tfits_photo_met$Topt))
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
#load("spp_clim_env,RData")
#add temperature data for species's climate envelop
#names(spp_clim_env)[1]<-"Species"
#tfits_photo_met<-merge(tfits_photo_met,spp_clim_env,by=c("Species"),all=TRUE)

#add temperature data for species's seed source location

#topath<-"//ad.uws.edu.au/dfshare/HomesHWK$/90931217/My Documents/bioclim"
#metdat_ss_new<-get_bioclim_seedsource(data=dat_info,return="summary",topath=topath)

#save(metdat_ss_new,file="metdat_ss_new,RData")
#topath<-"//ad.uws.edu.au/dfshare/HomesHWK$/90931217/My Documents/Documents/Dushan/WorldClim"
#worldclim <- get_worldclim_rasters_new(topath)
#metdat_ss<-get_worldclim_temp_seedsource(data=dat_info,return="summary",topath=topath)

#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
#dat_info_new<-read.csv("dataset_seedsource_info_additions.csv")
#topath<-"//ad.uws.edu.au/dfshare/HomesHWK$/90931217/My Documents/bioclim"
#metdat_ss_add<-get_bioclim_seedsource(data=dat_info_new,return="summary",topath=topath)


#-------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
load("metdat_ss_new,RData")
tfits_photo_met<-merge(tfits_photo_met,metdat_ss_new,by=c("Species","DataSet"),all=TRUE)
tfits_photo_met<-subset(tfits_photo_met,!is.na(tfits_photo_met$DataSet))
write.csv(tfits_photo_met,paste0(path,"/Tables/netphoto_parameters_with_met_data.csv"),row.names=FALSE,sep=",")
#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
t_growth_photo<-tfits_photo_met[c("Species","DataSet","Season" ,"Temp_Treatment","Tavg_30" )]

with(tfits_photo_met,plot(BIO7/10,Topt))
