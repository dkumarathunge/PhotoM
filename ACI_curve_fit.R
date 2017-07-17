

##----------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------
#FIT ACi Curves

#Note: Change the path to the derectory containing data
##----------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------------------------------
#fit ACi curves
#----------------------------------------------------------------------------------------------------------------------
#path="C:/Users/90931217/Documents/Dushan/Repos/PhotoM"
path=getwd()
#----------------------------------------------------------------------------------------------------------------------

source("T response functions.R")
source("metdataprocess.R")
source("rsq_mm.R")
source("new_tresponse_photo_fits.R")
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Daintree 
fa <- makecurves(path,"/Data/Daintree_ACidata_processed_summer.csv")
da <- makedata(path,"/Data/Daintree_ACidata_processed_summer.csv",fa)
write.csv(da,paste0(path,"/DATA_JIM/Daintree.csv"),row.names=FALSE,sep=",")

#HFECommon garden
fc <- makecurves(path,"/Data/hfe_ACidata_processed.csv")
dc <- makedata.1(path,"/Data/hfe_ACidata_processed.csv",fc)
write.csv(dc,paste0(path,"/DATA_JIM/HFE_CG.csv"),row.names=FALSE,sep=",")


#GWW Henrique's
fd <- makecurves(path,"/Data/gwwACidata_processed_new.csv")
dd <- makedata(path,"/Data/gwwACidata_processed_new.csv",fd)
write.csv(dd,paste0(path,"/DATA_JIM/Western_woodlands.csv"),row.names=FALSE,sep=",")


#WTC1
f1<-makecurves(path,"/Data/wtc1_Aci_processed.csv")
d1 <- makedata.1(path,"/Data/wtc1_Aci_processed.csv",f1)
write.csv(d1,paste0(path,"/DATA_JIM/wtc1.csv"),row.names=FALSE,sep=",")


#WTC2
f2 <- makecurves(path,"/Data/WTC2_ACidata_processed.csv")
d2 <- makedata.1(path,"/Data/WTC2_ACidata_processed.csv",f2)
write.csv(d2,paste0(path,"/DATA_JIM/wtc2.csv"),row.names=FALSE,sep=",")


#WTC3
f3 <- makecurves(path,"/Data/WTC3_ACidata_processed.csv")
d3 <- makedata(path,"/Data/WTC3_ACidata_processed.csv",f3)
write.csv(d3,paste0(path,"/DATA_JIM/wtc3.csv"),row.names=FALSE,sep=",")



#Tumbarumba site
f4 <- makecurves(path,"/Data/Tumbarumba_ACidata_processed.csv")
d4 <- makedata(path,"/Data/Tumbarumba_ACidata_processed.csv",f4)
write.csv(d4,paste0(path,"/DATA_JIM/tumbarumba.csv"),row.names=FALSE,sep=",")



#WTC4
f5 <- makecurves(path,"/Data/acitwtc4_cleaned.csv")
d5 <- makedata(path,"/Data/acitwtc4_cleaned.csv",f5)
write.csv(d5,paste0(path,"/DATA_JIM/wtc4.csv"),row.names=FALSE,sep=",")



#SAVANA
f6<- makecurves(path,"/Data/SIOP Leaf gas exchange Cernusak et al.csv")
d6 <- makedata.1(path,"/Data/SIOP Leaf gas exchange Cernusak et al.csv",f6)
write.csv(d6,paste0(path,"/DATA_JIM/savanna.csv"),row.names=FALSE,sep=",")


#Mike Aspinwall Data

f7 <- makecurves(path,"/Data/mike_aspinwall_corymbia_calophylla.V1.csv")
d7 <- makedata(path,"/Data/mike_aspinwall_corymbia_calophylla.V1.csv",f7)
write.csv(d7,paste0(path,"/DATA_JIM/Mike_aspinwall.csv"),row.names=FALSE,sep=",")


#GREAT

fg<-makecurves(path,"/Data/Great_Aci_data_Processed.csv")
dg <- makedata.1(path,"/Data/Great_Aci_data_Processed.csv",fg)
write.csv(dg,paste0(path,"/DATA_JIM/Great.csv"),row.names=FALSE,sep=",")

#with(dg,plot(Ts,Vcmax,col=c("black","red","green")[factor(Room)],pch=c(16,15,2)[factor(Povanance)]))
#with(dg,plot(Ts,Jmax,col=c("black","red","green")[factor(Room)],pch=c(16,15,2)[factor(Povanance)]))

#-----------------------------------------------------------------------------------------------------------------
#species outside Australia

#Amazone
amz<- makecurves(path,"/Data/data_out_auz/AmazonACIdata_f.csv")
damz <- makedata(path,"Data/data_out_auz/AmazonACIdata_f.csv",amz)
write.csv(damz,paste0(path,"/DATA_JIM/amazon.csv"),row.names=FALSE,sep=",")


#Anjelica:tropical species from Rwanda

rwa<- makecurves(path,"/Data/data_out_auz/Angelica_Varhammar_tropical_species.csv")
drwa <- makedata(path,"Data/data_out_auz/Angelica_Varhammar_tropical_species.csv",rwa)
write.csv(drwa,paste0(path,"/DATA_JIM/RF_MONTANE.csv"),row.names=FALSE,sep=",")


#wang et al: two species

wang<- makecurves(path,"/Data/data_out_auz/Betula_pendula_Pinus_sylvestris_wang_etal.csv")
dwang <- makedata(path,"Data/data_out_auz/Betula_pendula_Pinus_sylvestris_wang_etal.csv",wang)
write.csv(dwang,paste0(path,"/DATA_JIM/wang et al.csv"),row.names=FALSE,sep=",")



# Chamaecyparis obtusa

cob<- makecurves(path,"/Data/data_out_auz/Chamaecyparis obtusa_han_etal.csv")
dcob <- makedata.1(path,"Data/data_out_auz/Chamaecyparis obtusa_han_etal.csv",cob)
write.csv(dcob,paste0(path,"/DATA_JIM/cham.obtusa.csv"),row.names=FALSE,sep=",")


# Dillaway et al 

dil<- makecurves(path,"/Data/data_out_auz/dillaway_etal.csv")
ddil <- makedata(path,"Data/data_out_auz/dillaway_etal.csv",dil)
write.csv(ddil,paste0(path,"/DATA_JIM/dillaway.csv"),row.names=FALSE,sep=",")


# Dreyer et al 

dre<- makecurves(path,"/Data/data_out_auz/DreyerSevenSpp_final.csv")
ddre <- makedata(path,"Data/data_out_auz/DreyerSevenSpp_final.csv",dre)
write.csv(ddre,paste0(path,"/DATA_JIM/dreyer.csv"),row.names=FALSE,sep=",")


# Haley et al (Glycine_max)

hal<- makecurves(path,"/Data/data_out_auz/Glycine_max_Harley_f2.csv")
dhal <- makedata(path,"Data/data_out_auz/Glycine_max_Harley_f2.csv",hal)
write.csv(dhal,paste0(path,"/DATA_JIM/Haley et al.csv"),row.names=FALSE,sep=",")


# Medlyn et al

med<- makecurves(path,"/Data/data_out_auz/Medlyn_etal_all.csv")
dmed <- makedata(path,"Data/data_out_auz/Medlyn_etal_all.csv",med)

write.csv(dmed,paste0(path,"/DATA_JIM/medlyn.csv"),row.names=FALSE,sep=",")


# Onoda et al

on<- makecurves(path,"/Data/data_out_auz/Onoda_etal.csv")
don <- makedata.1(path,"Data/data_out_auz/Onoda_etal.csv",on)
write.csv(don,paste0(path,"/DATA_JIM/onoda.csv"),row.names=FALSE,sep=",")


# Picea_mariana_way_etal

way<- makecurves(path,"/Data/data_out_auz/Picea_mariana_way_etal.csv")
dway <- makedata.1(path,"Data/data_out_auz/Picea_mariana_way_etal.csv",way)
write.csv(dway,paste0(path,"/DATA_JIM/way.csv"),row.names=FALSE,sep=",")


# Picea_abies_tarvainen_etal

tar<- makecurves(path,"/Data/data_out_auz/Pieca_abies_tarvainen_etal.csv")
dtar <- makedata(path,"Data/data_out_auz/Pieca_abies_tarvainen_etal.csv",tar)
write.csv(dtar,paste0(path,"/DATA_JIM/Picea_abies_tarvainen_etal.csv"),row.names=FALSE,sep=",")


# Pinus_densiflora_han_etal

han<- makecurves(path,"/Data/data_out_auz/Pinus_densiflora_han_etal.csv")
dhan <- makedata.1(path,"Data/data_out_auz/Pinus_densiflora_han_etal.csv",han)
write.csv(dhan,paste0(path,"/DATA_JIM/Pinus_densiflora_han_etal.csv"),row.names=FALSE,sep=",")


# Pinus_radiata_Walcroft

wal<- makecurves(path,"/Data/data_out_auz/Pinus_radiata_Walcroft.csv")
dwal <- makedata(path,"Data/data_out_auz/Pinus_radiata_Walcroft.csv",wal)
write.csv(dwal,paste0(path,"/DATA_JIM/Pinus_radiata_Walcroft.csv"),row.names=FALSE,sep=",")


# Pinus_teada_Ellsworth

elw<- makecurves(path,"/Data/data_out_auz/Pteada_Ellsworth.csv")
delw <- makedata.1(path,"Data/data_out_auz/Pteada_Ellsworth.csv",elw)
write.csv(delw,paste0(path,"/DATA_JIM/Pinus_teada_Ellsworth.csv"),row.names=FALSE,sep=",")


# Strassmayer et al

str<-makecurves(path,"/Data/data_out_auz/Strassemeyer.csv")
dstr <- makedata.1(path,"Data/data_out_auz/Strassemeyer.csv",str)
write.csv(dstr,paste0(path,"/DATA_JIM/Strassmayer et al.csv"),row.names=FALSE,sep=",")


# Walcroft et al_2

wpea<- makecurves(path,"/Data/data_out_auz/walcroft peach.csv")
dpea <- makedata(path,"Data/data_out_auz/walcroft peach.csv",wpea)
write.csv(dpea,paste0(path,"/DATA_JIM/Walcroft et al.csv"),row.names=FALSE,sep=",")


#E. pauciflora 

wep<- makecurves(path,"/Data/euc_pauciflora_mk.csv")
dwep <- makedata(path,"Data/euc_pauciflora_mk.csv",wep)
write.csv(dwep,paste0(path,"/DATA_JIM/pauciflora.csv"),row.names=FALSE,sep=",")


#Arctic species

arc<-makecurves(path,"/Data/data_out_auz/Arctic_A-Ci_curves_2012-2015_V2.csv")
darc <- makedata.1(path,"Data/data_out_auz/Arctic_A-Ci_curves_2012-2015_V2.csv",arc)

with(darc,plot(Ts,Rd,col=USDA_Species_Code))
with(darc,plot(Ts,Jmax,col=USDA_Species_Code))



#Spruce Anna et al

spru<-makecurves(path,"/Data/data_out_auz/SPRUCE_3_cohort_ACi_data.csv")
dspru <- makedata.1(path,"Data/data_out_auz/SPRUCE_3_cohort_ACi_data.csv",spru)

with(dspru,plot(Ts,lowAs,col=Year))
with(subset(dspru,dspru$Month==4),plot(Ts,Jmax,col=Year))

#source("Topt_Photo_300.R")
#to save fitted parametes
#save(da,dc,dd,d1,d2,d3,d4,d5,d6,file="acifits,RData")

#save("acifits,RData")
#load("acifits,RData")
no_curves<-c(23,123,39,66,156,130,37,125,20,148,71,138,28,8,247,164,177,78,47,31,26,21,56,12)

#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

get_all_parameters <- function(pathtodata){
  
  o <- getwd()
  on.exit(setwd(o))
  setwd(pathtodata)
  
  files.summer <- list.files(pattern="\\.csv")
  allxlfiles.s<-list()
  
  for (i in seq_along(files.summer)){
    
    
    dat <- suppressWarnings(read.csv(files.summer[i]))
   
    dat <- as.data.frame(dat)
    
    #dat$Obs<-as.numeric(dat$Obs)
    
    
    
    allxlfiles.s[[i]] <- dat
    
  }
  
  summer<-data.frame(rbind.fill(allxlfiles.s))
  
  
  return(summer)
  
}

path<-getwd()
pathtodata<-paste0(path,"/DATA_JIM")

aci_parameters<-get_all_parameters(pathtodata)
aci_parameters<-aci_parameters[c("Vcmax","Jmax","Rd","Ts","TsK","JVr","Treatment","Season","T_Treatment","Povanance","season","Ttreatment")]
#fit Tresponse

euc_bio.b$TsK<-with(euc_bio.b,Ts+273.15)

vcm<-fitpeaked(aci_parameters,return="Peak")
jma<-fitjmax_mm(subset(euc_bio.b,!is.na(Jmax)),random="Species",return="Peak")



euc_bio.b$VcmaxFit<-with(euc_bio.b, vcm[[1]] * exp((vcm[[2]] * (TsK - 298.15))/(298.15 * 0.008314 * 
                                                                                  TsK)) * (1 + exp((298.15 * vcm[[3]] - 200)/(298.15 * 0.008314)))/(1 + 
                                                                                                                                                      exp((TsK * vcm[[3]] - 200)/(TsK * 0.008314))))

euc_bio.b$JmaxFit<-with(euc_bio.b, jma[[1]] * exp((jma[[2]] * (TsK - 298.15))/(298.15 * 0.008314 * 
                                                                                 TsK)) * (1 + exp((298.15 * jma[[3]] - 200)/(298.15 * 0.008314)))/(1 + 
                                                                                                                                                     exp((TsK * jma[[3]] - 200)/(TsK * 0.008314))))

with(aci_parameters,plot(JVr~Ts))
