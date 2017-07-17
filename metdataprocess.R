#install.packages("dismo")
#install.packages('chron')
#install.packages('sp')
#install_bitbucket("remkoduursma/HIEv")

#library(chron)
#library(sp)
#library(devtools)
#library(HIEv)
#library(dismo)

#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

#-- Script contains functions to download,read and process the met data required for chapter 1
#-- 
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#function to download WTC met data
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

metD<-function(path,pattern){
  setToPath(path)
  metsearch<-searchHIEv(pattern,quiet=TRUE)
  downloadHIEv(metsearch)
}


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#function to download tmin and tmax from eMAST database
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------


#set path to the directory where downloaded csv to be saved
#this function download CSV of monthly Tmin and Tmax for the period of 1970 to 2012

#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/MetNew"


getmet<-function(path,var1,var2,lat,lon,site){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  
  base <-"http://dapds00.nci.org.au/thredds/ncss/eMAST_TERN/Climate/eMAST/ANUClimate/0_01deg/v1m0_aus/mon/land/"
  base.1 <- paste(base,sprintf("%s/e_01/1970_2012/eMAST_ANUClimate_mon_%s_v1m0_1970_2012_agg.ncml?var=air_temperature&latitude=%f&longitude=%f",var1,var1,lat,lon),sep = "")
  url.1<-paste(base.1,"&time_start=1970-01-01T00%3A00%3A00Z&time_end=2012-12-01T00%3A00%3A00Z&accept=csv_file")
  
  
  
  base.2 <- paste(base,sprintf("%s/e_01/1970_2012/eMAST_ANUClimate_mon_%s_v1m0_1970_2012_agg.ncml?var=air_temperature&latitude=%f&longitude=%f",var2,var2,lat,lon),sep = "")
  url.2<-paste(base.2,"&time_start=1970-01-01T00%3A00%3A00Z&time_end=2012-12-01T00%3A00%3A00Z&accept=csv_file")
  
  
  desf.1<-paste0(path,sprintf("/eMAST.%s.%s.csv",site,var1))
  desf.2<-paste0(path,sprintf("/eMAST.%s.%s.csv",site,var2))
  
  download.file(url=url.1, destfile=desf.1,cacheOK = TRUE)
  download.file(url=url.2, destfile=desf.2,cacheOK = TRUE)
  
  #names(fl1)[1:4]<-c("DateTime","Latitude", "Longitude","var1")
  #names(fl2)[1:4]<-c("DateTime","Latitude", "Longitude","var2")
  
  #write.csv(fl1,paste0(path,sprintf("/%s.%s.csv",site,var1),collapse = " "))
  #write.csv(fl2,paste0(path,sprintf("/%s.%s.csv",site,var2),collapse = " "))
  
}

#getmet(path=path,var="tmin",lat=-30.4335,lon=152.0421,site="StyxRiver")

#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

#function to get running averages for any defined time (wtc2, wtc3 and wtc4)

wtcMet<- function(path,fname,from){
  
  #set WD to the data directory
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read all csv files and bind to one file
  
  wtcmet <- read.csv(paste0(path,"/",fname) )
  #wtcmet<-read.csv("WTC_TEMP_CM_WTCMET_L1_v2.csv")
  
  wtcmet$DateTime <- as.POSIXct(wtcmet$DateTime,format="%Y-%m-%d")
  
  from=as.Date(from)
  
  end30<-as.Date(from-30)
  end60<-as.Date(from-15)
  end90<-as.Date(from-3)
  
  #get a subset for required period
  dat30<-wtcmet[wtcmet$DateTime >= end30 & wtcmet$DateTime <= from,]
  dat60<-wtcmet[wtcmet$DateTime >= end60 & wtcmet$DateTime <= from,]
  dat90<-wtcmet[wtcmet$DateTime >= end90 & wtcmet$DateTime <= from,]
  

  
  #get average across given time period (set only chamber temperature)
  
  t30 <- dplyr::summarize(group_by(dat30,Treat),
                          #Tair_SP=mean(Tair_SP,na.rm=T),
                          #RH_al=mean(RH_al,na.rm=T),
                          #DP_al=mean(DP_al,na.rm=T),
                          Tmin=mean(Tmin,na.rm=T),
                          Tmax=mean(Tmax,na.rm=T),
                          Tmean=mean(Tmean,na.rm=T))
  
  
  t60 <- dplyr::summarize(group_by(dat60,Treat),
                          #Tair_SP=mean(Tair_SP,na.rm=T),
                          #RH_al=mean(RH_al,na.rm=T),
                          #DP_al=mean(DP_al,na.rm=T),
                          Tmin=mean(Tmin,na.rm=T),
                          Tmax=mean(Tmax,na.rm=T),
                          Tmean=mean(Tmean,na.rm=T))
  
  t90 <- dplyr::summarize(group_by(dat90,Treat),
                          #Tair_SP=mean(Tair_SP,na.rm=T),
                          #RH_al=mean(RH_al,na.rm=T),
                          #DP_al=mean(DP_al,na.rm=T),
                          Tmin=mean(Tmin,na.rm=T),
                          Tmax=mean(Tmax,na.rm=T),
                          Tmean=mean(Tmean,na.rm=T))
  
  
  taverage <- data.frame(t30,t60,t90)
  
  
  
  return(taverage)
}

#wtcMet(path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/WTC4",
        #fname="WTC_TEMP_CM_WTCMET_L1_v2.csv",from="2016-04-01")
        
#------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------
#to get monthly averages for WTC1 HFE and OZ flux sites

wtcMet2 <- function(path,from,fname){
  
  #set WD to the data directory
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read all csv files and bind to one file
  
  wtcmet <- read.csv(paste0(path,"/",fname) )
  #wtcmet <- read.csv("WTC_TEMP_CM_WTCMET_L1_v2.csv")
  
  wtcmet$DateTime <- as.POSIXct(wtcmet$DateTime,format="%Y-%m-%d")
  
  #get a subset for required period
  from=as.Date(from)
  
  end30<-as.Date(from-30)
  end60<-as.Date(from-15)
  end90<-as.Date(from-3)
  
  #get a subset for required period
  dat30<-wtcmet[wtcmet$DateTime >= end30 & wtcmet$DateTime <= from,]
  dat60<-wtcmet[wtcmet$DateTime >= end60 & wtcmet$DateTime <= from,]
  dat90<-wtcmet[wtcmet$DateTime >= end90 & wtcmet$DateTime <= from,]

  
  
  #get average/min/max temperature across given time period (set only chamber temperature)
  #30 day back
  t30 <- mean(dat30$Tmean,na.rm=T)
  t30min<-mean(dat30$Tmin,na.rm=T)
  t30max<-mean(dat30$Tmax,na.rm=T)
  

  #60 day back
  t60 <- mean(dat60$Tmean,na.rm=T)
  t60min<-mean(dat60$Tmin,na.rm=T)
  t60max<-mean(dat60$Tmax,na.rm=T)
  

  
  #90 day back
  t90 <- mean(dat90$Tmean,na.rm=T)
  t90min<-mean(dat90$Tmin,na.rm=T)
  t90max<-mean(dat90$Tmax,na.rm=T)
  
  
  
  taverage<-data.frame(t30,t30min,t30max,t60,t60min,t60max,t90,t90min,t90max)
  #taverage<-data.frame(t30d,t60d,t90d)
  #taverage$Index<-c("Tmean","Tmin","Tmax")
  return(taverage)

}

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#function to read and process met data from ABO stations
#path="C:/Users/90931217/Documents/MetData/MetNew"
#met_out(path=path,from("1965-01-01"),to=("2000-01-01"))

met_out<-function(path,from,to){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read Tmax
  maxfiles <- list.files(pattern="\\Max.csv")
  maxf<-list()
 
  for (i in seq_along(maxfiles)){
    dat.max <- suppressWarnings(read.csv(maxfiles[i]))
    dat.max <- as.data.frame(dat.max)
    maxf[[i]] <- dat.max
  }
  
  #read Tmin
  minfiles <- list.files(pattern="\\Min.csv")
  minf<-list()
  
  for (i in seq_along(minfiles)){
    dat.min <- suppressWarnings(read.csv(minfiles[i]))
    dat.min <- as.data.frame(dat.min)
    minf[[i]] <- dat.min
  }
  
  tmax<-do.call(rbind.data.frame,maxf)
  tmin<-do.call(rbind.data.frame,minf)
  
  tmet<-merge(tmax,tmin,by=c("Bureau.of.Meteorology.station.number","Year","Month","Day"))
  
  #to add Date Time column
  
  tmet<-within(tmet, DateTime <- sprintf("%d-%02d-%02d", Year, Month,Day))
  
 tmet$DateTime <- as.POSIXct(tmet$DateTime,format="%Y-%m-%d")
  
 
 #to get required time frame to be get the averages
 dat<-tmet[tmet$DateTime >= as.Date(from) & tmet$Date <= as.Date(to),]
 
 
 #get monthly averages from daily data
 
 tmonthly <- dplyr::summarize(group_by(dat,Year,Month),
                         #Tair_SP=mean(Tair_SP,na.rm=T),
                         #RH_al=mean(RH_al,na.rm=T),
                         #DP_al=mean(DP_al,na.rm=T),
                         #Tmin=mean(,na.rm=T),
                         Tmax=mean(Maximum.temperature..Degree.C.,na.rm=T),
                         Tmin=mean(Minimum.temperature..Degree.C.,na.rm=T))
 
 
 tmonthly$Season<-NA
 tmonthly$Season[which(tmonthly$Month %in% c(12,1,2))]<-"Summer"
 tmonthly$Season[which(tmonthly$Month %in% c(3,4,5))]<-"Autumn"
 tmonthly$Season[which(tmonthly$Month %in% c(6,7,8))]<-"Winter"
 tmonthly$Season[which(tmonthly$Month %in% c(9,10,11))]<-"Spring"
 tmonthly$Season<-as.factor(tmonthly$Season)
 
 
 #to get monthly average temperature
 tmonthly$Tmean<-(tmonthly$Tmax+tmonthly$Tmin)/2
 
 
 
 #to add precipitation: this is need to get bioclim variables. I added random Precip values until get real RF data
 #Cannot use RF related bioclim variables
 
 #tmonthly$Precip<-tmonthly$Tmean+5
 
 dat.1<-subset(tmonthly,tmonthly$Season=="Summer")
 dat.2<-subset(tmonthly,tmonthly$Season=="Winter")
  
 #get statistics
 Tmean<-mean(tmonthly$Tmean,na.rm=T) #average temperature over the total period
 
 SummerTmean <- mean(dat.1$Tmean,na.rm=T) #average temperature during summer over the total period
 SummerTmax <- mean(dat.1$Tmax,na.rm=T) #average maximum temperature during summer over the total period 
 WinterTmean <- mean(dat.2$Tmean,na.rm=T) #average temperature during winter over the total period
 WinterTmin <- mean(dat.2$Tmin,na.rm=T) #average maximum temperature during winter over the total period 
 
 taverage <- c(Tmean,SummerTmean,SummerTmax,WinterTmean,WinterTmin)
 names(taverage)[1:5]<-c("Tmean","SummerTmean","SummerTmax","WinterTmean","WinterTmin")
 
  return(taverage)
  
}

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#function to read,process and get averagees from eMAST temperature data
#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/Met_Data/eMAST"
#met_eMAST(path=path,fname="wtc3")

met_eMAST<-function(path,fname){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read Tmax

  maxfiles <- list.files(pattern=sprintf("\\.%s.tmax.csv",fname))
  maxf<-list()
  
  for (i in seq_along(maxfiles)){
    dat.max <- suppressWarnings(read.csv(maxfiles[i]))
    dat.max <- as.data.frame(dat.max)
    maxf[[i]] <- dat.max
  }
  
  #read Tmin
  minfiles <- list.files(pattern=sprintf("\\.%s.tmin.csv",fname))
  minf<-list()
  
  for (i in seq_along(minfiles)){
    dat.min <- suppressWarnings(read.csv(minfiles[i]))
    dat.min <- as.data.frame(dat.min)
    minf[[i]] <- dat.min
  }
  
  #to remove unnecessary columns and add column names
  tmax<-do.call(rbind.data.frame,maxf)
  tmax<-tmax[-c(2,3)]
  names(tmax)[1:2]<-c("DateTime","Tmax")
  
  #to remove unnecessary columns and add column names
  tmin<-do.call(rbind.data.frame,minf)
  tmin<-tmin[-c(2,3)]
  names(tmin)[1:2]<-c("DateTime","Tmin")
  
  
  #merge tmax and tmin together
  
  tmet<-merge(tmax,tmin,by=c("DateTime"))
  
 
  tmet$DateTime <- as.POSIXct(tmet$DateTime,format="%Y-%m-%d")
  
  #tmet$Month<-months(tmet$DateTime, abbreviate = FALSE)
  
  tmet$Month<-month(tmet$DateTime)
  
  #to get month of the year to a seperate column
  
  #to get required time frame to be get the averages
  #dat<-tmet[tmet$DateTime >= as.Date(from) & tmet$Date <= as.Date(to),]
  
  #to define seasons
  tmet$Season<-NA
  tmet$Season[which(tmet$Month %in% c(12,1,2))]<-"Summer"
  tmet$Season[which(tmet$Month %in% c(3,4,5))]<-"Autumn"
  tmet$Season[which(tmet$Month %in% c(6,7,8))]<-"Winter"
  tmet$Season[which(tmet$Month %in% c(9,10,11))]<-"Spring"
  tmet$Season<-as.factor(tmet$Season)
  
  
  #to add monthly mean temperature
  
  tmet$Tmean<-(tmet$Tmax+tmet$Tmin)/2
 
  
  #tmonthly$Precip<-tmonthly$Tmean+5
  
  dat.1<-subset(tmet,tmet$Season=="Summer")
  dat.2<-subset(tmet,tmet$Season=="Winter")
  
  #get statistics
  Tmean<-mean(tmet$Tmean,na.rm=T) #average temperature over the total period
  
  SummerTmean <- mean(dat.1$Tmean,na.rm=T) #average temperature during summer over the total period
  SummerTmax <- mean(dat.1$Tmax,na.rm=T) #average maximum temperature during summer over the total period 
  WinterTmean <- mean(dat.2$Tmean,na.rm=T) #average temperature during winter over the total period
  WinterTmin <- mean(dat.2$Tmin,na.rm=T) #average maximum temperature during winter over the total period 
  
  taverage <- c(Tmean,SummerTmean,SummerTmax,WinterTmean,WinterTmin)
  names(taverage)[1:5]<-c("Tmean","SummerTmean","SummerTmax","WinterTmean","WinterTmin")
  
  return(taverage)
  
}


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------


#function to compile all monthly met files together, process and get daily averages for WTC3 and WTC4
#this save a csv with daily average temperature in given directory (path)
#set path to the directory where the individual csv files located

#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/WTC4"
wtcMet_process <- function(path){
  
  #set WD to the data directory
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read all csv files and bind to one file
  metfiles <- list.files(pattern="\\.csv")
  allfiles<-list()
  
  for (i in seq_along(metfiles)){
    dat <- suppressWarnings(read.csv(metfiles[i]))
    dat <- as.data.frame(dat)
    allfiles[[i]] <- dat
  }
  
  wtcmet<-do.call(rbind.data.frame,allfiles)
  
  wtcmet.1<-subset(wtcmet,wtcmet$Tair_al>-5 & wtcmet$Tair_al <55 ) # to remove unusual temperatures
  
  #get date (yyyy-mm-dd)
  wtcmet.1$DateTime <- as.POSIXct(wtcmet.1$DateTime,format="%Y-%m-%d")
  
  wtcmet.2<-dplyr::summarize(group_by(wtcmet.1,chamber,DateTime),
                             #Tair_SP=mean(Tair_SP,na.rm=T),
                             #RH_al=mean(RH_al,na.rm=T),
                             #DP_al=mean(DP_al,na.rm=T),
                             Tmin=min(Tair_al,na.rm=T),
                             Tmax=max(Tair_al,na.rm=T),
                             Tmean=mean(Tair_al,na.rm=T))
  
  
  wtcmet.2$chamber<-as.numeric(wtcmet.2$chamber)
  
  #assign temperature treatments
  is.even <- function(x) x %% 2 == 0
  is.odd <- function(x) x %% 2 != 0
  
  wtcmet.2$Treat<-NA
  wtcmet.2$Treat[which(is.even(wtcmet.2$chamber))]<-"Elevated"
  wtcmet.2$Treat[which(is.odd(wtcmet.2$chamber))]<-"Ambient"
  
  write.csv(wtcmet.2,paste0(path,"/WTC_TEMP_CM_WTCMET_L1_v2.csv",collapse = " "))
  
  return(wtcmet.2)
}

#wtcMet_process(path=path)
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------


#read and process wtc2 met data

readwtcmet<-function(path){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  wtcmet<-read.csv(paste0(path,"/Raw/WTC_TEMP_CM_WTCMET_201007-201111_L1_v1.csv",sep=""))
  
  wtcmet$chamber<-as.numeric(str_sub(wtcmet$Chamber, start = -2, end = -1))
  wtcmet$DateTime<-as.POSIXlt(as.character(wtcmet$DateTime),format="%d/%m/%Y")
  
  wtcmet$Treat<-NA
  wtcmet$Treat[which(wtcmet$chamber %in% c(1,3,4,6,8,11))]<-"Elevated"
  wtcmet$Treat[which(wtcmet$chamber %in% c(7,5,9,10,2,12))]<-"Ambient"
  
  colnames(wtcmet)<-c("DateTime","Chamber","Tmean","Tmin","Tmax","AvgRH_al","MinRH_al",	"MaxRH_al","chamber","Treat")
  
  write.csv(wtcmet,paste0(path,"/WTC_TEMP_CM_WTCMET_L1_v2.csv",collapse = " "))
  
  return(wtcmet)
}
#wtc2mettest<-readwtcmet(path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/WTC2")
#readwtcmet(path=path)


#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#function to read and process WTC1 met data

wtc1metprocess <- function(path){
  
  #set WD to the data directory
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read all csv files and bind to one file
  metfiles <- list.files(pattern="\\.csv")
  allfiles<-list()
  
  for (i in seq_along(metfiles)){
    dat <- suppressWarnings(read.csv(metfiles[i]))
    dat <- as.data.frame(dat)
    allfiles[[i]] <- dat
  }
  
  wtcmet<-do.call(rbind.fill,allfiles)
  
  wtcmet$DateTime <- as.POSIXct(wtcmet$DateTime,format="%Y-%m-%d")
  wtcmet.1<-subset(wtcmet,wtcmet$TairMet>-5 & wtcmet$TairMet <55 ) # to remove unusual temperatures
  
  wtcmet.2<-dplyr::summarize(group_by(wtcmet.1,DateTime),
                             #Tair_SP=mean(Tair_SP,na.rm=T),
                             #RH_al=mean(RH_al,na.rm=T),
                             #DP_al=mean(DP_al,na.rm=T),
                             Tmin=min(TairMet,na.rm=T),
                             Tmax=max(TairMet,na.rm=T),
                             Tmean=mean(TairMet,na.rm=T))
  
  
  write.csv(wtcmet.2,paste0(path,"/WTC_TEMP_CM_WTCMET_L1_v2.csv",collapse = " "))
  
  return(wtcmet.2)
}
#wtc1metprocess(path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/Rwanda")

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

#process hfe met data 

#function to read and process hfe met data for 2010. 

#download required data file from HIEv

#methfe<-metD(path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/HFE",pattern="CG_StomatalSensitivity_MetData_2010-2011")

hfemet<- function(path){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read excel file
  dat <- suppressWarnings(read_excel("hfemetdata20102011.xls",sheet=1))
  
  #change col names to fit for other met data files
  
  colnames(dat)<-c("Year","Day","Month","Time","Rsol","TairMet","RHmet","winddirection",
                   "SoilHeatFlux","VPDmet","Accumilated Rain (mm)",
                   "windmet")
  
  
  #get daily average/min/max temperattures
  
  dat.2<-dplyr::summarize(group_by(dat,Year,Month,Day),
                          #Tair_SP=mean(Tair_SP,na.rm=T),
                          #RH_al=mean(RH_al,na.rm=T),
                          #DP_al=mean(DP_al,na.rm=T),
                          Tmin=min(TairMet,na.rm=T),
                          Tmax=max(TairMet,na.rm=T),
                          Tmean=mean(TairMet,na.rm=T))
  
  #convert month to numeric
  
  dat.2$month_new<-NA
  dat.2$month_new[which(dat.2$Month=="Jan")]<-1
  dat.2$month_new[which(dat.2$Month=="Feb")]<-2
  dat.2$month_new[which(dat.2$Month=="Mar")]<-3
  dat.2$month_new[which(dat.2$Month=="Apr")]<-4
  dat.2$month_new[which(dat.2$Month=="May")]<-5
  dat.2$month_new[which(dat.2$Month=="June")]<-6
  dat.2$month_new[which(dat.2$Month=="July")]<-7
  dat.2$month_new[which(dat.2$Month=="Aug")]<-8
  dat.2$month_new[which(dat.2$Month=="Sept")]<-9
  dat.2$month_new[which(dat.2$Month=="Oct")]<-10
  dat.2$month_new[which(dat.2$Month=="Nov")]<-11
  dat.2$month_new[which(dat.2$Month=="Dec")]<-12
  dat.2$month_new<-as.numeric(dat.2$month_new)
  
  
  
  dat.3<-dat.2[ order(dat.2$Year,dat.2$month_new), ]
  dat.3$DateTime <-seq(from=as.Date("2009-01-01"),to=as.Date("2011-12-31"),by="day")
  
  #to merge with 2007-2009 data
  
  metrest<-read.csv(paste(path,"/WTC_TEMP_CM_WTCMET_L1_v2.csv",sep=""))
  
  dat.4<-rbind.fill(metrest,dat.3)
  dat.4<-dat.4[ ,-c(1,6,7,8,9) ]
  dat.4$DateTime<-as.POSIXct(dat.4$DateTime,format="%Y-%m-%d")
  
  write.csv(dat.4,paste0(path,"/HFE_TEMP_CM_WTCMET_L1_v2.csv",collapse = " "))
  
}

#hfemet(path=path)

#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/WTC1"


#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

#function to process met data from flux sites

#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/Met from Mingkai"

readfluxmet<-function(path,filename,startdate,enddate){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  flxmet<-read.csv(paste0(path,"/",filename)) 
  
  flxmet.1<-dplyr::summarize(group_by(flxmet,year,doy),
                             
                             Tmin=min(tair,na.rm=T),
                             Tmax=max(tair,na.rm=T),
                             Tmean=mean(tair,na.rm=T))
  
  
  flxmet.1$DateTime<-seq(as.POSIXct(startdate), as.POSIXct(enddate), by="day")
  
  newpath<-"C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/Met from Mingkai/Processed"
  write.csv(flxmet.1,paste0(newpath,"/",filename,collapse = " "))
  
  #dat<-flxmet.1[flxmet.1$DateTime >= as.Date(from) & flxmet.1$Date <= as.Date(to),]
  #dat$yvar <- dat[[yvar]]
  
  #taverage <- mean(dat$yvar,na.rm=T)
  #tmin<-min(dat$yvar,na.rm=T)
  #tmax<-max(dat$yvar,na.rm=T)
  
  #taverage <- c(taverage,tmin,tmax)
  #names(taverage)[1:3]<-c("Tmean","Tmin","Tmax")
  return(flxmet.1)
  
}

#readfluxmet(path=path,
#filename="DalyPasture_met_forcing.csv",startdate=as.Date("2008-01-01"),enddate=as.Date("2012-12-31"))


#readfluxmet(path=path,
#filename="Tumbarumba_met_forcing.csv",startdate=as.Date("2002-01-01"),enddate=as.Date("2014-12-31"))

#readfluxmet(path=path,
#filename="DalyUncleared_met_forcing.csv",startdate=as.Date("2008-01-01"),enddate=as.Date("2014-12-31"))

#readfluxmet(path=path,
#filename="greatwesternwoodlands_met_forcing.csv",startdate=as.Date("2013-01-01"),enddate=as.Date("2014-12-31"))


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------

#function to process met data from ABM sites (this is only for use in RF site as daily met data during 
#gas exchange measurements is not available)

#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/ClimOrig/RF"

readabm<-function(path,filename,startdate,enddate){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  flxmet.3<-read.csv(paste0(path,"/",filename)) 
  
  flxmet.3$DateTime<-seq(as.POSIXct(startdate), as.POSIXct(enddate), by="day")
  
  
  flxmet.3<-flxmet.3[,-c(1,2)]
  colnames(flxmet.3)<-c("year","Month","day","Tmin","Tmax","DateTime")
  
  flxmet.3$Tmean<-(flxmet.3$Tmin+flxmet.3$Tmax)/2
  
  
  write.csv(flxmet.3,paste0(path,"/",filename,collapse = " "))
  
  #dat<-flxmet.1[flxmet.1$DateTime >= as.Date(from) & flxmet.1$Date <= as.Date(to),]
  #dat$yvar <- dat[[yvar]]
  
  #taverage <- mean(dat$yvar,na.rm=T)
  #tmin<-min(dat$yvar,na.rm=T)
  #tmax<-max(dat$yvar,na.rm=T)
  
  #taverage <- c(taverage,tmin,tmax)
  #names(taverage)[1:3]<-c("Tmean","Tmin","Tmax")
  return(flxmet.3)
  
}

#readabm(path=path,filename="IDCJAC0011_031037_1800_all.csv",startdate="1967-01-01",enddate="2016-11-02")






#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#process HFE met station, 2009 data
#path="C:/Users/90931217/Documents/MetData/HFE"

#hfemet2009<-read.csv("Copyof2009MetDataHFE.csv")
#hfemet2009.1<-summaryBy(.~Year+Month+Day,data=hfemet2009,FUN=mean,keep.names=TRUE)

#hfemet2009$DateTime<-seq(as.POSIXct("2009-01-01"), as.POSIXct("2009-12-31"), by="day")
#colnames(hfemet2009)<-c("year","day","month","Time","solrad","TairMet","RHmet","winddirection","SoilHeatFlux","Battery.Voltage",
#                        "VPDmet","Rain_Accu","windmet","DateTime")

#write.csv(hfemet2009,"hfecommongardenmet_2009.csv",row.names=FALSE,sep=",")
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#wtcMet1 <- function(path,from,to){
  
  #o <- getwd()
  #on.exit(setwd(o))
  #setwd(path)
  
  #wtc1metdata<-read.csv(paste(path,"/WTC_CO2DROUGHT_CM_WTCFLUX_20080220-20090316_L0_V1.csv",sep=""))
  
  
  #wtc1metdata$Year <- as.factor(format(wtc1metdata$DateTime, format="%y"))
  #wtc1metdata$month <- as.factor(format(wtc1metdata$DateTime, format="%m"))
  
  #wtc1metdata.1<-subset(wtc1metdata,wtc1metdata$TAref>-5 & wtc1metdata$TAref <55 ) # to remove unusual temperatures
  #get date (yyyy-mm-dd)
  #wtc1metdata.1$Date <- as.POSIXct(wtc1metdata.1$Date,format="%Y-%m-%d")
  
  
  #get a subset for required period
  #dat<-wtc1metdata.1[wtc1metdata.1$Date >= from & wtc1metdata.1$Date <= to,]
  #dat$chamber<-as.numeric(dat$chamber)
  
  #dat$yvar <- dat[[yvar]]
  
  #get average across given time period (set only chamber temperature)
  
  #taverage <- mean(dat$TAref,na.rm=T)
  
  
  #taverage <- as.data.frame(taverage)
  
  
  
  #return(taverage)
#}
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#function to get wtc2 met data (this is not in use)

#wtc2Met <- function(path,from,to,yvar){
  
  #set WD to the data directory
  #o <- getwd()
  #on.exit(setwd(o))
  #setwd(path)
  
  #read all csv files and bind to one file
  #metfiles <- list.files(pattern="\\.csv")
  #allfiles<-list()
  
  #for (i in seq_along(metfiles)){
    #dat <- suppressWarnings(read.csv(metfiles[i]))
    #dat <- as.data.frame(dat)
    #allfiles[[i]] <- dat
  #}
  
  #wtcmet<-do.call(rbind.data.frame,allfiles)
  #wtcmet.1<-subset(wtcmet,wtcmet$Tair_al>-5 & wtcmet$Tair_al <55 ) # to remove unusual temperatures
  
  #get date (yyyy-mm-dd)
  
  
  #wtcmet.1$Date <- as.POSIXct(wtcmet.1$DateTime,format="%Y-%m-%d")
  
  
  #get a subset for required period
  #dat<-wtcmet.1[wtcmet.1$Date >= from & wtcmet.1$Date <= to,]
  #dat$chamber<-as.numeric(dat$chamber)
  
  #dat$yvar <- dat[[yvar]]
  
  
  #assign temperature treatments
  
  #dat$Treat<-NA
  #dat$Treat[which(dat$chamber %in% c(1,3,4,6,8,11))]<-"Elevated"
  #dat$Treat[which(dat$chamber %in% c(7,5,9,10,2,12))]<-"Ambient"
  
  
  
  #get average across given time period (set only chamber temperature)
  
  #taverage <- dplyr::summarize(group_by(dat,Treat),
                               #Tair_SP=mean(Tair_SP,na.rm=T),
                               #RH_al=mean(RH_al,na.rm=T),
                               #DP_al=mean(DP_al,na.rm=T),
                               #Tmin=min(yvar,na.rm=T),
                               #Tmax=max(yvar,na.rm=T),
                               #Tmean=mean(yvar,na.rm=T))
  
  
  #taverage <- as.data.frame(taverage)
  
  
  
  #return(taverage)
#}

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#wtc.test<-wtc2Met(path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/MET_DATA/WTC2",
                  #from=as.Date("2010-12-01"),to=as.Date("2011-02-28"),yvar="Tair_al")   

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#function to get average temperatures for E. globulus/Tasmania data.
#This function get average Tmin/Tmax for the closest past 1, 2 and 3 months (not by dates)

met_eMAST.2<-function(path,fname,from){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  #read Tmax
  
  maxfiles <- list.files(pattern=sprintf("\\.%s.tmax.csv",fname))
  maxf<-list()
  
  for (i in seq_along(maxfiles)){
    dat.max <- suppressWarnings(read.csv(maxfiles[i]))
    dat.max <- as.data.frame(dat.max)
    maxf[[i]] <- dat.max
  }
  
  #read Tmin
  minfiles <- list.files(pattern=sprintf("\\.%s.tmin.csv",fname))
  minf<-list()
  
  for (i in seq_along(minfiles)){
    dat.min <- suppressWarnings(read.csv(minfiles[i]))
    dat.min <- as.data.frame(dat.min)
    minf[[i]] <- dat.min
  }
  
  #to remove unnecessary columns and add column names
  tmax<-do.call(rbind.data.frame,maxf)
  tmax<-tmax[-c(2,3)]
  names(tmax)[1:2]<-c("DateTime","Tmax")
  
  #to remove unnecessary columns and add column names
  tmin<-do.call(rbind.data.frame,minf)
  tmin<-tmin[-c(2,3)]
  names(tmin)[1:2]<-c("DateTime","Tmin")
  
  
  #merge tmax and tmin together
  
  tmet<-merge(tmax,tmin,by=c("DateTime"))
  
  
  tmet$DateTime <- as.POSIXct(tmet$DateTime,format="%Y-%m-%d")
  
  #tmet$Month<-months(tmet$DateTime, abbreviate = FALSE)
  
  tmet$Month<-month(tmet$DateTime)
  
  #to get month of the year to a seperate column
  
  #to add monthly mean temperature
  
  tmet$Tmean<-(tmet$Tmax+tmet$Tmin)/2
  
  
  from=as.Date(from)
  end30<-as.Date(from-30)
  end60<-as.Date(from-60)
  end90<-as.Date(from-90)
  
  dat30<-tmet[tmet$DateTime >= end30 & tmet$DateTime <= from,]
  dat60<-tmet[tmet$DateTime >= end60 & tmet$DateTime <= from,]
  dat90<-tmet[tmet$DateTime >= end90 & tmet$DateTime <= from,]
  
  #to get required time frame to be get the averages
  #dat<-tmet[tmet$DateTime >= as.Date(from) & tmet$Date <= as.Date(to),]
  
  #to define seasons
  tmet$Season<-NA
  tmet$Season[which(tmet$Month %in% c(12,1,2))]<-"Summer"
  tmet$Season[which(tmet$Month %in% c(3,4,5))]<-"Autumn"
  tmet$Season[which(tmet$Month %in% c(6,7,8))]<-"Winter"
  tmet$Season[which(tmet$Month %in% c(9,10,11))]<-"Spring"
  tmet$Season<-as.factor(tmet$Season)
  
  

  
  #get average/min/max temperature across given time period (set only chamber temperature)
  #30 day back
  t30 <- mean(dat30$Tmean,na.rm=T)
  t30min<-mean(dat30$Tmin,na.rm=T)
  t30max<-mean(dat30$Tmax,na.rm=T)
  
  
  #60 day back
  t60 <- mean(dat60$Tmean,na.rm=T)
  t60min<-mean(dat60$Tmin,na.rm=T)
  t60max<-mean(dat60$Tmax,na.rm=T)
  
  
  
  #90 day back
  t90 <- mean(dat90$Tmean,na.rm=T)
  t90min<-mean(dat90$Tmin,na.rm=T)
  t90max<-mean(dat90$Tmax,na.rm=T)
  
  
  
  taverage<-data.frame(t30,t30min,t30max,t60,t60min,t60max,t90,t90min,t90max)
  #taverage<-data.frame(t30d,t60d,t90d)
  #taverage$Index<-c("Tmean","Tmin","Tmax")
  return(taverage)
  
  

  
}

#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM/Met_Data/eMAST"
#met_eMAST.2(path=path,fname="tasmania",from="2000-01-01")

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#function to download Mean annual average/min/max temperature for a given location
#data should be a dtaframe with latitude and longitude of a/or different species
#this is just an extention of Remko's functions to get met data for species occurence locations

get_worldclim_temp<-function(data, topath=topath, return=c("all","summary"), worldclim=NULL){
  
  return <- match.arg(return)
  
  if(is.null(worldclim)){
    worldclim <- get_worldclim_rasters_new(topath)
  }
  tmean_raster <- worldclim$tmean_raster
  tmax_raster <- worldclim$tmax_raster
  tmin_raster <- worldclim$tmin_raster
  
  #extract worldclim data; extract the gridCell ID for observations
  tmeanm <- tmaxm <- tminm <- matrix(ncol=12, nrow=nrow(data))
  for(i in 1:12){
    tmeanm[,i] <- 0.1 * extract(tmean_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
    tmaxm[,i] <- 0.1 * extract(tmax_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
    tminm[,i] <- 0.1 * extract(tmin_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
  }
  colnames(tmeanm) <- paste0("tmean_",1:12)
  colnames(tmaxm) <- paste0("tmax_",1:12)
  colnames(tminm) <- paste0("tmin_",1:12)
  
  pxy <- cbind(data, as.data.frame(tmeanm), as.data.frame(tmaxm),as.data.frame(tminm))
  #names(pxy)[2:3] <- c("longitude","latitude")
  
  pxy$MAT <- apply(pxy[,grep("tmean_",names(pxy))],1,mean)
  pxy$MATmax <- apply(pxy[,grep("tmax_",names(pxy))],1,mean)
  pxy$MATmin <- apply(pxy[,grep("tmin_",names(pxy))],1,mean)
  
  # 
  if(return == "all")return(pxy)
  
  if(return == "summary"){
    
    #dfr<-summaryBy(MAT+MATmax+MATmin~species,FUN=mean,data=pxy,na.rm=TRUE)
    
    dfr <- suppressWarnings(with(pxy, data.frame(species=unique(data$species),
    n=nrow(data),
    lat_mean=mean(latitude,na.rm=TRUE),
    long_mean=mean(longitude,na.rm=TRUE),
    MAT_mean=mean(MAT,na.rm=TRUE),
    MAT_q05=quantile(MAT,0.05,na.rm=TRUE),
    MAT_q95=quantile(MAT,0.95,na.rm=TRUE),
    
    MATmax_mean=mean(MATmax,na.rm=TRUE),
    MATmax_q05=quantile(MATmax,0.05,na.rm=TRUE),
    MATmax_q95=quantile(MATmax,0.95,na.rm=TRUE),
    
    MATmin_mean=mean(MATmin,na.rm=TRUE),
    MATmin_q05=quantile(MATmin,0.05,na.rm=TRUE),
    MATmin_q95=quantile(MATmin,0.95,na.rm=TRUE))))
    rownames(dfr) <- NULL
    #colnames(dfr)<-c("species","MAT","MATmax","MATmin")
    return(dfr)
  }
}

#-----------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#function to download worldclim rasters
#this download Tmin,Tmax and Tmean
#resolution was set to 10 min (as Remko did in original functions)

get_worldclim_rasters_new<- function(topath, clean=FALSE){
  
  download_worldclim <- function(basen, topath){
    
    wc_fn_full <- file.path(topath, basen)
    
    if(!file.exists(wc_fn_full)){
      message("Downloading WorldClim 10min layers ... ", appendLF=FALSE)
      download.file(file.path("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur",basen),
                    wc_fn_full, mode="wb")
      message("done.")
    }
    
    u <- unzip(wc_fn_full, exdir=topath)
    
    return(u)
  }
  
  download_worldclim("tmean_30s_esri.zip", topath)
  download_worldclim("tmax_30s_esri.zip", topath) 
  download_worldclim("tmin_30s_esri.zip", topath) 
  
  if(clean){
    unlink(c(wc_fn_full,dir(file.path(topath,"tmean"),recursive=TRUE)))
    unlink(c(wc_fn_full,dir(file.path(topath,"tmax"),recursive=TRUE)))
    unlink(c(wc_fn_full,dir(file.path(topath,"tmin"),recursive=TRUE)))
  }
  
  # Read the rasters into a list
  tmean_raster <- list()
  tmax_raster <- list()
  tmin_raster <- list()
  
  #message("Reading Worldclim rasters ... ", appendLF = FALSE)
  for(i in 1:12){
    tmean_raster[[i]] <- raster(file.path(topath, sprintf("tmean/tmean_%s", i)))
    tmax_raster[[i]] <- raster(file.path(topath, sprintf("tmax/tmax_%s", i)))
    tmin_raster[[i]] <- raster(file.path(topath, sprintf("tmin/tmin_%s", i)))
    
  }
  #message("done.")
  
  return(list(tmean_raster=tmean_raster, tmax_raster=tmax_raster,tmin_raster=tmin_raster))
}

#topath<-"C:/Users/90931217/Documents/Dushan/WorldClim"

#data<-read.csv(paste0(topath,"/locations_species_oaz.csv"))

#get_worldclim_temp(topath=topath,data=data,return="summary")            


#get_worldclim_rasters_new(topath)
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#function to download daily Tmin,Tmax and Taverage for any location 
#resolution of the grided dataset is 1degree
#data source<-NASA Prediction of Worldwide Energy Resource (POWER)/Climatology Resource for Agroclimatology/Global coverage on a 1° latitude by 1° longitude grid 
#https://power.larc.nasa.gov

#path<-set directory to a place where downloaded files to be saved
#lat and lon<- in degree decimals
#currently, the function set to download daily data from 1990-01-01 to 2015-12-31


getmetworld<-function(path,lat,lon,site){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  
  base <-"https://power.larc.nasa.gov/cgi-bin/"
  url.1 <- paste(base,sprintf("agro.cgi?email=&step=1&lat=%f&lon=%f&ms=1&ds=1&ys=1990&me=31&de=12&ye=2015&p=T2M&p=T2MN&p=T2MX&submit=Submit",lat,lon),sep = "")
  
  desf.1<-paste0(path,sprintf("/CRA.%s.txt",site))
  download.file(url=url.1, destfile=desf.1,cacheOK = TRUE)
  
  data<-read.table(sprintf("CRA.%s.txt",site),header = FALSE,skip=15,sep="")
  names(data)[1:5]<-c("Year","DOY", "Tmean","Tmin","Tmax")
  data$DateTime <-seq(from=as.Date("1990-01-01"),to=as.Date("2015-12-31"),by="day")
  
  
  write.csv(data,paste0(path,sprintf("/%s.csv",site),collapse = " "))
  
  #write.csv(fl2,paste0(path,sprintf("/%s.%s.csv",site,var2),collapse = " "))
  
}



#function to get average temperature for DUKE FACE site (Pinus taeda)
getavgduke<-function(path,fname,from){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  
  wtcmet <- read.csv(paste0(path,"/",fname) )
  #wtcmet <- read.csv("WTC_TEMP_CM_WTCMET_L1_v2.csv")
  
  wtcmet$DateTime <- as.POSIXct(wtcmet$DateTime,format="%Y-%m-%d")
  
  #get a subset for required period
  from=as.Date(from)
  
  end30<-as.Date(from-30)
  end60<-as.Date(from-15)
  end90<-as.Date(from-3)
  
  #get a subset for required period
  dat30<-wtcmet[wtcmet$DateTime >= end30 & wtcmet$DateTime <= from,]
  dat60<-wtcmet[wtcmet$DateTime >= end60 & wtcmet$DateTime <= from,]
  dat90<-wtcmet[wtcmet$DateTime >= end90 & wtcmet$DateTime <= from,]
  
  
  #get average/min/max temperature across given time period (set only chamber temperature)
  #30 day back
  t30 <- mean(dat30$Tmean,na.rm=T)
  tmax30 <- max(dat30$Tmean,na.rm=T)
  
  #60 day back
  t60 <- mean(dat60$Tmean,na.rm=T)
 
  
  #90 day back
  t90 <- mean(dat90$Tmean,na.rm=T)
   
  
  
  taverage<-data.frame(t30,tmax30,t60,t90)
  return(taverage)
  
  
  }
#getavgduke(path=path,fname="duke_temp_data.csv",from="1999-9-30")

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#function to get long-term average temperatures at seed source

get_worldclim_temp_seedsource<-function(data, topath=topath, return=c("all","summary"), worldclim=NULL){
  
  return <- match.arg(return)
  
  if(is.null(worldclim)){
    worldclim <- get_worldclim_rasters_new(topath)
  }
  tmean_raster <- worldclim$tmean_raster
  tmax_raster <- worldclim$tmax_raster
  tmin_raster <- worldclim$tmin_raster
  
  #extract worldclim data; extract the gridCell ID for observations
  tmeanm <- tmaxm <- tminm <- matrix(ncol=12, nrow=nrow(data))
  for(i in 1:12){
    tmeanm[,i] <- 0.1 * extract(tmean_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
    tmaxm[,i] <- 0.1 * extract(tmax_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
    tminm[,i] <- 0.1 * extract(tmin_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
  }
  colnames(tmeanm) <- paste0("tmean_",1:12)
  colnames(tmaxm) <- paste0("tmax_",1:12)
  colnames(tminm) <- paste0("tmin_",1:12)
  
  pxy <- cbind(data, as.data.frame(tmeanm), as.data.frame(tmaxm),as.data.frame(tminm))
  #names(pxy)[2:3] <- c("longitude","latitude")
  
  pxy$MAT <- apply(pxy[,grep("tmean_",names(pxy))],1,mean)
  pxy$MATmax <- apply(pxy[,grep("tmax_",names(pxy))],1,mean)
  pxy$MATmin <- apply(pxy[,grep("tmin_",names(pxy))],1,mean)
  
  # 
  if(return == "all")return(pxy)
  
  if(return == "summary"){
    
    dfr<-summaryBy(MAT+MATmax+MATmin~species,FUN=mean,data=pxy,na.rm=TRUE)
    
        rownames(dfr) <- NULL
    #colnames(dfr)<-c("species","MAT","MATmax","MATmin")
    return(dfr)
  }
}


#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#functions to download BIOCLIM data
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
topath<-"//ad.uws.edu.au/dfshare/HomesHWK$/90931217/My Documents/bioclim"
#function to download bioclim rasters
get_bioclim<-function(topath, clean=FALSE){
  
  download_bioclim <- function(basen, topath){
    
    wc_fn_full <- file.path(topath, basen)
    
    if(!file.exists(wc_fn_full)){
      message("Downloading Bioclim 0.5min layers ... ", appendLF=FALSE)
      download.file(file.path("http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur",basen),
                    wc_fn_full, mode="wb")
      message("done.")
    }
    
    u <- unzip(wc_fn_full, exdir=topath)
    
    return(u)
  }
  
  download_bioclim("bio_30s_esri.zip", topath)
  
  if(clean){
    unlink(c(wc_fn_full,dir(file.path(topath,"bio"),recursive=TRUE)))
    
  }
  
  # Read the rasters into a list
  bio_raster <- list()
  
  
  #message("Reading Worldclim rasters ... ", appendLF = FALSE)
  for(i in 1:19){
    bio_raster[[i]] <- raster(file.path(topath, sprintf("bio/bio_%s", i)))
    
  }
  #message("done.")
  
  return(list(bio_raster=bio_raster))
}



#function to extract bioclim variables

get_bioclim_seedsource<-function(data, topath=topath, return=c("all","summary"), worldclim=NULL){
  
  return <- match.arg(return)
  
  if(is.null(worldclim)){
    worldclim <- get_bioclim(topath)
  }
  bio_raster <- worldclim$bio_raster
  
  
  #extract worldclim data; extract the gridCell ID for observations
  bio<- matrix(ncol=19, nrow=nrow(data))
  for(i in 1:19){
    bio[,i] <-  extract(bio_raster[[i]], cbind(data$longitude,data$latitude), method='simple')
  }
  
  colnames(bio) <- paste0("BIO",1:19)
  
  pxy <- cbind(data, as.data.frame(bio))
  
  
  if(return == "all")return(pxy)
  
  if(return == "summary"){
    
    dfr<-pxy[c(1:20)]
    
    rownames(dfr) <- NULL
    #colnames(dfr)<-c("species","MAT","MATmax","MATmin")
    return(dfr)
  }
}

#get_bioclim_seedsource(data=data,return="summary",topath=topath)


