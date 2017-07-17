


fit.nlme<-function(dat,random,yvar){
  
  dat$rand <- dat[[random]]
  dat$yvar <- dat[[yvar]]
  
  testfit2<-nlme(yvar~Aopt-(b*(Tleaf-Topt)^2),fixed=list(Aopt + Topt + b ~ 1),random = Aopt+Topt ~ 1 | rand,
                 start=list(fixed=c(Aopt=max(dat$yvar),Topt=25,b=0.05)),data=dat)
  
  test<-testfit2$tTable
  
  aopt<-summary(testfit2)$tTable[[1]]
  topt<-summary(testfit2)$tTable[[2]]
  b<-summary(testfit2)$tTable[[3]]
  aopt.se<-summary(testfit2)$tTable[[4]]
  topt.se<-summary(testfit2)$tTable[[5]]
  b.se<-summary(testfit2)$tTable[[6]]
  
  #to get R2 between fitted and observed photosynthesis
  r<-cor(fitted(testfit2),dat$yvar)
  r2<-r*r
  
  #test for normality of residuals
  rest<-residuals(testfit2)
  norm<-shapiro.test(rest)
  s<-norm$statistic
  pvalue<-norm$p.value
  
  param<-cbind(aopt,topt,b,aopt.se,topt.se,b.se)
  
  names(param)[1:6] <- c("aopt","topt","b","aopt.se","topt.se","b.se")
  
  #return(param)
  
  
  return(c(param,r2,s,pvalue))
  
}

#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------


#fit quadratic model to estimate Topt
#------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------
#WTC 1
d1.1<-subset(d1,d1$CO2_Treat=="Ambient" & d1$Water_treat=="wet")
d1.2<-split(d1.1,d1.1$Season)


#fit quadratic model to estimate Topt

#e.sal.300<-lapply(d1.2[1],function(x)fit.nlme(x,yvar="lowAs",rand="Chamber"))
wtc1.300<-get_topts(lapply(d1.2,FUN=fitquad300)) #cannot fit mix model. Igrore, negative b
wtc1.300$Month<-names(d1.2)
wtc1.300$Species<-as.factor("E.saligna")
wtc1.300<-subset(wtc1.300,wtc1.300$b>0) #remove poor fit with negative b
wtc1.300$DataSet<-"WTC1"


wtc1.800<-get_topts(lapply(d1.2,FUN=fitquad800))
wtc1.800$Month<-names(d1.2)
wtc1.800$Species<-as.factor("E.saligna")
wtc1.800$DataSet<-"WTC1"
wtc1.800<-subset(wtc1.800,wtc1.800$b>0) #remove poor fit with negative b

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#WTC2

d2.1<-subset(d2,d2$C_treat=="0C")
d2.2<-split(d2.1,paste(d2.1$Temp.treat,d2.1$Season))
list1<-d2.2[3]
wtc2.300.1<-get_topts(lapply(list1,FUN=fitquad300)) #model without random effects (Ambient/September) # few data points. cannot fit a mix model
list2<-d2.2[c(1:2,4:6)]
wtc2.300.2<-get_topts(lapply(list2,function(x)fit.nlme(x,yvar="lowAs",rand="Chamber")))
wtc2.300<-data.frame((rbind(wtc2.300.1,wtc2.300.2)))

#wtc2.300$Ttreatment<-get_names(d2.2)[,1]
wtc2.300$Month<-get_names(d2.2)[,2]
wtc2.300$Species<-"E.globulus"  
wtc2.300$DataSet<-"WTC2"
wtc2.300$Treatment<-capFirst(get_names(d2.2)[,1])

wtc2.800<-get_topts(lapply(d2.2,function(x)fit.nlme(x,yvar="highAs",rand="Chamber")))
wtc2.800$Treatment<-capFirst(get_names(d2.2)[,1])
wtc2.800$Month<-get_names(d2.2)[,2]
wtc2.800$Species<-"E.globulus"  
wtc2.800$DataSet<-"WTC2"


#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------

#WTC3
d3.2<-split(d3,paste(d3$Treatment,d3$Season))

list1<-d3.2[1]
wtc3.300.1<-get_topts(lapply(list1,FUN=fitquad300)) #model without random effects. Mix model didn't converged


list2<-d3.2[c(2:6)]
wtc3.300.2<-get_topts(lapply(list2,function(x)fit.nlme(x,yvar="lowAs",rand="Chamber")))

wtc3.300<-data.frame(rbind(wtc3.300.1,wtc3.300.2))

wtc3.300$Treatment<-get_names(d3.2)[,1]
wtc3.300$Month<-get_names(d3.2)[,2]
wtc3.300$Species<-"E.tereticornis"
wtc3.300$DataSet<-"WTC3"


list3<-d3.2[4]
wtc3.800.1<-get_topts(lapply(list3,FUN=fitquad800)) #model without random effects. Mix model didn't converged
wtc3.800.1$Treatment<-get_names(list3)[,2]
wtc3.800.1$Month<-get_names(list3)[,4]


list4<-d3.2[c(1:3,5,6)]
wtc3.800.2<-get_topts(lapply(list4,function(x)fit.nlme(x,yvar="highAs",rand="Chamber")))
wtc3.800.2$Treatment<-get_names(list4)[,1]
wtc3.800.2$Month<-get_names(list4)[,2]

wtc3.800<-data.frame(rbind(wtc3.800.1,wtc3.800.2))
wtc3.800$Species<-"E.tereticornis"
wtc3.800$DataSet<-"WTC3"


#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#WTC4

d5.1<-split(d5,paste(d5$Ttreatment,d5$Season))


wtc4.300<-get_topts(lapply(d5.1,function(x)fit.nlme(x,yvar="lowAs",rand="Chamber")))

wtc4.300$Treatment<-get_names(d5.1)[,1]
wtc4.300$Month<-get_names(d5.1)[,2]


wtc4.300$Species<-"E.parramattensis"
wtc4.300$DataSet<-"WTC4"


wtc4.800<-suppressWarnings(get_topts(lapply(d5.1,function(x)fit.nlme(x,yvar="highAs",rand="Chamber"))))

wtc4.800$Treatment<-get_names(d5.1)[,1]
wtc4.800$Month<-get_names(d5.1)[,2]


wtc4.800$Species<-"E.parramattensis"
wtc4.800$DataSet<-"WTC4"


#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#for(i in 1:length(dc.2)){
 # toplot<-dc.2[[i]]
#  name<-names(dc.2)[i]
 # plot(lowAs~Ts,data=toplot,main=name,col=Species,pch=16,cex=2)
#}

#HFE

dc.2<-split(dc,paste(dc$Spp,dc$Season))
list1<-dc.2[-c(14)]

hfe.300<-get_topts(lapply(list1,FUN=fitquad300))

hfe.300$Species<-get_names(list1)[,1]
hfe.300$Month<-get_names(list1)[,2]

hfe.300<-subset(hfe.300,hfe.300$topt.se<100) #remove a datapoint wit SE=287 (poor fit for Topt)

#hfe.300<-subset(hfe.300,hfe.300$topt.se<10) #remove poor fits with very high SEs (SE>15)
#to remove poor fits (negative b)
#hfe.300<-subset(hfe.300,hfe.300$b>0)
hfe.300$DataSet<-"HFE_CG"

#lapply(dc.2[15],FUN=fitquad300)
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

#Tumbarumba site

tmb.300<-fit.nlme(d4,yvar="lowAs",random="Leaf.Age")
tmb.300<-data.frame(do.call(rbind,list(tmb.300)))
tmb.300$Species<-as.factor(tmb.300$Species<-"E.delegatensis")


tmb.300$Month<-"November" #Just to merge with November Vcmax and Jmax data



tmb.300$DataSet<-"TMB"

#No clear peak for highAs
#with(d4.1,plot(Ts,lowAs,col=Month,pch=16,cex=2))
#legend("bottomright",levels(d4$Month),pch=16,col=COL)

#d4.1<-subset(d4,d4$Curve!=17)
#d4.2<-split(d4.1,d4.1$Month)
#lapply(d4.2[2],FUN=fitquad)

#lapply(d4.2[1],function(x)fit.nlme(x,yvar="lowAs",rando="Leaf.Age"))

with(d4,plot(Ts,lowAs,main="Tumbarumba",pch=16,cex=2))
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

#SAVANA

e.tr.boot<-split(d6,d6$Canopy_position)

e.tr.300<-lapply(e.tr.boot,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf.number")) #fit only for canopy

e.tr.300<-data.frame(do.call(rbind,list(e.tr.300[[1]]))) #extract parameters only for canopy

e.tr.300$Species<-"E.tetradonta"
e.tr.300$Month<-"September"
e.tr.300$DataSet<-"SAVANNA"

e.tr.800<-lapply(e.tr.boot,function(x)fit.nlme(x,yvar="highAs",rand="Leaf.number")) #fit only for canopy
e.tr.800<-data.frame(do.call(rbind,list(e.tr.800[[1]]))) #extract parameters only for canopy

e.tr.800$Species<-"E.tetradonta"
e.tr.800$Month<-"September"
e.tr.800$DataSet<-"SAVANNA"

with(subset(d6,d6$Canopy_position=="Canopy"),plot(Ts,lowAs,main="Savanna",pch=16,cex=2))
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

#GWW
dd.1<-subset(dd,!dd$Curve %in% c(19,23,27,44)) #remove four outliers. Three at highest leaf temperature and one at lower
                                                          # not affect Topt or Aopt significantly but lower standard error 


#E.salmonophloia<-subset(dd,dd$Species=="E.salmonophloia"& Curve !=7)#remove abnormally high lowAs value at Tleaf=26C
#E.salubris<-subset(dd,Species=="E.salubris"& Curve !=17)#remove abnormally low lowAs value at Tleaf=26C

#E.scoparia.fit<-fit.nlme(E.scoparia,yvar="lowAs",rand="Tree")
#E.salmonophloia.fit<-fitquad300(E.salmonophloia)
#E.salubris.fit<-fitquad300(E.salubris)


gww.300<-data.frame(do.call(rbind,list(fitquad300(dd.1))))

gww.300$Species<-"ASAW spp"
gww.300$DataSet<-"GWW"


#E.scoparia.800<-fit.nlme(E.scoparia,yvar="highAs",rand="Tree")
E.salmonophloia.800<-fitquad800(E.salmonophloia) #topt out of the range
E.salubris.800<-fitquad800(E.salubris)

gww.800<-data.frame(do.call(rbind,list(fitquad800(dd.1))))
gww.800$Species<-"ASAW spp"
gww.800$DataSet<-"GWW"


with(dd.1,plot(Ts,lowAs,main="GWW",pch=16,cex=2))
#with(E.scoparia,plot(Ts,highAs,col=Tree,pch=16,cex=2))
#with(E.salmonophloia,plot(Ts,highAs,col=Tree,pch=16,cex=2))
#with(E.salubris,plot(Ts,highAs,col=Tree,pch=16,cex=2))
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------
#
#Daintree

#rf.300<-split(da,paste(da$Spp,da$Season))
rf.300.list<-split(da,paste(da$Spp))
#for(i in 1:length(rf.300)){
  
  #dat.plot<-rf.300[[i]]
  #plot(lowAs~Ts,data=dat.plot,pch=19,cex=2,col=Leaf,main=names(rf.300[i]))
#}

#rf.300[2]<-NULL # Data not sufficient to fit models

#rf.300<-get_topts(lapply(rf.300.list,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf")))
rf.300<-data.frame(do.call(rbind,list(fit.nlme(da,yvar="lowAs",rand="Spp"))))

rf.300$Species<-"RF_AUSTRALIA"

rf.300$DataSet<-"RF_AUS"
rf.300$Month<-"SUMMER"

rf.800<-data.frame(do.call(rbind,list(fit.nlme(da,yvar="highAs",rand="Spp"))))
rf.800$Species<-"RF_AUSTRALIA"
rf.800$Month<-"SUMMER"

rf.800$DataSet<-"RF_AUS"

with(da,plot(Ts,lowAs,main="Daintree",pch=16,cex=2))
#----------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------------

#Corymbia provanances: Mike Aspinwall

d7.2<-split(d7,paste(d7$Provenance,d7$Growth_temp))
cor300<-get_topts(lapply(d7.2,function(x)fit.nlme(x,yvar="lowAs",rand="PotID")))
cor300$Treatment<-get_names(d7.2)[,2]
cor300$Species<-get_names(d7.2)[,1]
cor300$DataSet<-"MIKE_ASPINWALL"


d7.2[6]<-NULL
cor800<-get_topts(lapply(d7.2,function(x)fit.nlme(x,yvar="highAs",rand="PotID")))
cor800$Treatment<-get_names(d7.2)[,2]
cor800$Species<-get_names(d7.2)[,1]
cor800$DataSet<-"MIKE_ASPINWALL"


#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#E. pauciflora (Kirschbaum et al)

ep.fit<-fitquad300(dwep)
ep.300<-data.frame(do.call(rbind,list(ep.fit)))
ep.300$Species<-"Eucalyptus pauciflora"
ep.300$DataSet<-"KIRSCHBAUM"


ep.fit.800<-fitquad800(dwep)
ep.800<-data.frame(do.call(rbind,list(ep.fit.800)))
ep.800$Species<-"Eucalyptus pauciflora"
ep.800$DataSet<-"KIRSCHBAUM"

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

#species outside Australia

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
with(damz,plot(Ts,lowAs,cex=2,pch=16,main="Amazon"))
#Amazone

#Photosynthesis at Ci=275 ppm
#remove one outlier
#damz<-subset(damz,damz$Curve!=77)

damz.1<-split(damz,paste(damz$Season))

amz300.1<-get_topts(lapply(damz.1[1],function(x)fitquad300(x)))
amz300.2<-get_topts(lapply(damz.1[2],function(x)fit.nlme(x,yvar="lowAs",rand="Species")))

amz300<-rbind(amz300.1,amz300.2)

amz300$Season<-c("dry","wet")
amz300$Species<-"Amazon spp"
amz300$DataSet<-"RF_AMZ"

#Photosynthesis at Ci=800 ppm
amz800.2<-get_topts(lapply(damz.1[2],function(x)fit.nlme(x,yvar="highAs",rand="Species")))
amz800.1<-get_topts(lapply(damz.1[1],function(x)fitquad800(x)))

amz800<-rbind(amz800.1,amz800.2)


amz800$Season<-names(damz.1)
amz800$Species<-"Amazon spp"
amz800$DataSet<-"RF_AMZ"

with(damz,plot(Ts,lowAs,main="Amazon"))
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

#Anjelica:tropical species from Rwanda

drwa.1<-split(drwa,paste(drwa$Species))

#Photosynthesis at Ci=275 ppm
list.1<-drwa.1[c(1,4)] #for these species, mix models cannot be fitted. So fit standard quadratic form
rwa300.1<-get_topts(lapply(list.1,FUN=fitquad300))
rwa300.1$Species<-names(list.1)

list.2<-drwa.1[c(2,3,5,6)]
rwa300.2<-get_topts(lapply(list.2,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf")))
rwa300.2$Species<-names(list.2)

rwa300<-rbind(rwa300.1,rwa300.2)
#rwa300$Season<-"dry"
rwa300$DataSet<-"RF_MON"

#Photosynthesis at Ci=800 ppm
list.3<-drwa.1[c(2,4:6)]
rwa800.1<-get_topts(lapply(list.3,function(x)fit.nlme(x,yvar="highAs",rand="Leaf")))
rwa800.1$Species<-names(list.3)

list.4<-drwa.1[c(1,3)]
rwa800.2<-get_topts(lapply(list.4,FUN=fitquad300))
rwa800.2$Species<-names(list.4)

rwa800<-rbind(rwa800.1,rwa800.2)
#rwa800$Season<-"dry"
rwa800$DataSet<-"RF_MON"

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
with(subset(dwang,dwang$Species=="Pinus sylvestris"),plot(Ts,lowAs,col="black",main="Finland sp",cex=2,pch=16))

#wang et al: two species

#Photosynthesis at Ci=275 ppm
dwang.1<-split(dwang,paste(dwang$Species))
wang300<-get_topts(lapply(dwang.1,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf")))
wang300$Species<-names(dwang.1)
#wang300$Season<-"Summer"
wang300$DataSet<-"WANG_ET_AL"


#Photosynthesis at Ci=800 ppm
wang800<-get_topts(lapply(dwang.1,FUN=fitquad800))
wang800$Species<-names(dwang.1)
#wang800$Season<-"Summer"
wang800$DataSet<-"WANG_ET_AL"


#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#with(dcob,plot(Ts,lowAs,col=Season))

# Chamaecyparis obtusa

#Photosynthesis at Ci=275 ppm
dcob.1<-split(dcob,dcob$Season)
cob300<-get_topts(lapply(dcob.1,FUN=fitquad300))
cob300$Season<-names(dcob.1)
cob300$Species<-"Chamaecyparis obtusa"
cob300$DataSet<-"HAN_ET_AL_A"


#Photosynthesis at Ci=800 ppm
cob800<-get_topts(lapply(dcob.1,FUN=fitquad800))
cob800$Season<-names(dcob.1)
cob800$Species<-"Chamaecyparis obtusa"
cob800$DataSet<-"HAN_ET_AL_A"
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

# Dillaway et al 
for(i in 1:length(ddil.1)){
toplot<-ddil.1[[i]]
name<-names(ddil.1)[i]
plot(lowAs~Ts,data=toplot,main=name,col=Leaf,pch=16,cex=2)
}
#to remove space between species names and season names

#Photosynthesis at Ci=275 ppm

ddil$Species<-str_replace_all(ddil$Species, fixed(" "), "") 
ddil$Season<-str_replace_all(ddil$Season, fixed(" "), "") 

ddil.1<-split(ddil,paste(ddil$Species,ddil$Season))

listd.1<-ddil.1[c(3,4,6,9,10)]
dil300.1<-get_topts(lapply(listd.1,FUN=fitquad300))
dil300.1$Species<-get_names(listd.1)[,1]
dil300.1$Season<-get_names(listd.1)[,2]


listd.2<-ddil.1[c(1,2,5,7,8,11,12)]
dil300.2<-get_topts(lapply(listd.2,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf")))
dil300.2$Species<-get_names(listd.2)[,1]
dil300.2$Season<-get_names(listd.2)[,2]

dil300<-rbind(dil300.1,dil300.2)
#dil300<-subset(dil300,dil300$b>0) #remove Topts with negative b
#dil300<-subset(dil300,dil300$topt.se<10) #remove Topt with very high SE (SE=98, Topt=38)
dil300$DataSet<-"DILLAWAY"
dil300<-dil300[-c(8,11),]
#Photosynthesis at Ci=800 ppm

dil800<-get_topts(lapply(ddil.1,FUN=fitquad800)) #mix model didin't converged for most cases
dil800$Species<-get_names(ddil.1)[,1]
dil800$Season<-get_names(ddil.1)[,2]
dil800$DataSet<-"DILLAWAY"

#toplot<-subset(ddil,ddil$Species=="Liquidambarstyraciflua" & ddil$Season=="Illinois")
#with(toplot,plot(Ts,lowAs,cex=2))
#fitquad300(toplot)
#get_topts(lapply(listd.2,FUN=fitquad300))

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
# Dreyer et al
#Photosynthesis at Ci=275 ppm
ddre.1<-split(ddre,ddre$Species)
dre300<-get_topts(lapply(ddre.1,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf")))
dre300$Species<-names(ddre.1)
#dre300$Season<-"Summer"
dre300$DataSet<-"DREYER"


#Photosynthesis at Ci=800 ppm
dre800<-get_topts(lapply(ddre.1,function(x)fit.nlme(x,yvar="highAs",rand="Leaf")))
dre800$Species<-names(ddre.1)
#dre800$Season<-"Summer"
dre800$DataSet<-"DREYER"
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

# Medlyn et al

#Photosynthesis at Ci=275 ppm
dmed.1<-split(dmed,paste(dmed$Species,dmed$Season))

med300.1<-get_topts(lapply(dmed.1[10],FUN=fitquad300))
med300.1$Species<-get_names(dmed.1[10])[,2]
med300.1$Season<-get_names(dmed.1[10])[,4]



dmed.1[c(3,10)]<-NULL

med300.2<-get_topts(lapply(dmed.1,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf")))
med300.2$Species<-get_names(dmed.1)[,1]
med300.2$Season<-get_names(dmed.1)[,2]


#Photosynthesis at Ci=800 ppm

dmed.1<-split(dmed,paste(dmed$Species,dmed$Season))

listm.1<-dmed.1[c(3,7,10)]
med800.1<-get_topts(lapply(listm.1,FUN=fitquad800))
med800.1$Species<-get_names(listm.1)[,1]
med800.1$Season<-get_names(listm.1)[,2]

listm.2<-dmed.1[c(1:2,4:6,8:9,11:14)]
med800.2<-suppressWarnings(get_topts(lapply(listm.2,function(x)fit.nlme(x,yvar="highAs",rand="Leaf"))))
med800.2$Species<-get_names(listm.2)[,1]
med800.2$Season<-get_names(listm.2)[,2]


med300<-rbind(med300.1,med300.2)
med300$DataSet<-"MEDLYN"
med300$Species[which(med300$Species=="Landes")]<-"Pinus pinaster_L"
med300$Species[which(med300$Species=="Tamjoute")]<-"Pinus pinaster_T"
#med300<-subset(med300, Season %in% c("AugSept","Jan")) #get summer and winter 



med800<-rbind(med800.1,med800.2)
med800$DataSet<-"MEDLYN"
med800$Species[which(med800$Species=="Landes")]<-"Pinus pinaster_L"
med800$Species[which(med800$Species=="Tamjoute")]<-"Pinus pinaster_T"
#med800<-subset(med800, Season %in% c("AugSept","Jan")) #get summer and winter 

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#with(don.2,plot(Ts,lowAs,col=Season,pch=16))
# Onoda et al (poor fit/drop from further analysis??): No do not drop, they are important

don.1<-subset(don,don$GrowthCa=="amb")
don.2<-subset(don.1,don.1$lowAs<9.5) #remove few outliers


don.3<-split(don.2,paste(don.2$Season))

#Photosynthesis at Ci=300 ppm
don300<-suppressWarnings(get_topts(lapply(don.3,FUN=fitquad300)))
don300$Season<-names(don.3)
don300$Species<-"Fagus crenata"
don300$DataSet<-"OND"

#Photosynthesis at Ci=800 ppm (No clear peak, increasing trend)
#don800<-get_topts(lapply(don.2,fitquad800)) #mix model didn't converged
#don800$Season<-names(don.2)
#don800$Species<-"Fagus crenata"
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
# Picea_mariana_way_etal

dway.1<-split(dway,dway$Season)
dway300<-suppressWarnings(get_topts(lapply(dway.1,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf"))))
dway300$Treatment<-names(dway.1)
dway300$Species<-"Picea mariana"
dway300$DataSet<-"WAY_ET_AL"


dway800<-suppressWarnings(get_topts(lapply(dway.1,function(x)fit.nlme(x,yvar="highAs",rand="Leaf"))))
dway800$Treatment<-names(dway.1)
dway800$Species<-"Picea mariana"
dway800$DataSet<-"WAY_ET_AL"
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

# Picea_abies_tarvainen_etal
# fit only for year 1 data

dtar.1<-split(dtar,dtar$Season)
dtar300<-suppressWarnings(get_topts(lapply(dtar.1[1],function(x)fit.nlme(x,yvar="lowAs",rand="Leaf"))))
#dtar300$Season<-names(dtar.1)[1]
dtar300$Species<-"Picea abies"
dtar300$DataSet<-"TARVAINEN"


dtar800<-suppressWarnings(get_topts(lapply(dtar.1[1],function(x)fit.nlme(x,yvar="highAs",rand="Leaf"))))
#dtar800$Season<-names(dtar.1)[1]
dtar800$Species<-"Picea abies"
dtar800$DataSet<-"TARVAINEN"

with(subset(dtar,dtar$Season=="Year1"),plot(Ts,lowAs,col="black",main="Sweden",cex=2,pch=16))
with(dtar,plot(Ts,lowAs,col="black",main="Sweden"))
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
# Pinus_radiata_Walcroft
dwal.1<-split(dwal,dwal$Species)

dwal300<-get_topts(lapply(dwal.1,function(x)fit.nlme(x,yvar="lowAs",rand="Leaf"))) #topt is 9C out from the highest temperatture 
dwal300$Species<-"Pinus radiata"

dwal300$DataSet<-"WALCROFT"

dwal800<-get_topts(lapply(dwal.1,function(x)fit.nlme(x,yvar="highAs",rand="Leaf")))#negative b parameter
dwal800$Species<-"Pinus radiata"

dwal800$DataSet<-"WALCROFT"
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
# Pinus_teada_Ellsworth

delw.1<-split(delw,delw$Season)

delw300<-get_topts(lapply(delw.1,FUN=fitquad300)) #poor fit for winter. dropped from further analysis
delw300$Season<-names(delw.1)
delw300$Species<-"Pinus taeda"
delw300$DataSet<-"ELLSWORTH"

#Photo @ Ci=800 linearly increase with Tleaf. Quadratic model gave very high Topt
#dropped from further analysis

#delw800<-get_topts(lapply(delw.1[1],FUN=fitquad800)) 
#delw800$Season<-names(delw.1[1])
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#Strassmayer et al
dstr.1<-split(dstr,dstr$Species)
dstr.300<-get_topts(lapply(dstr.1,FUN=fitquad300))
dstr.300$Species<-names(dstr.1)
dstr.300$DataSet<-"STRASSEMEYER_ET_AL"

dstr.800<-get_topts(lapply(dstr.1,FUN=fitquad800))
dstr.800$Species<-names(dstr.1)
dstr.800$DataSet<-"STRASSEMEYER_ET_AL"

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#Han et al: P. densiflora


han.300<-data.frame(do.call(rbind,list(fitquad300(dhan))))
han.300$DataSet<-"HAN_2"
han.300$Species<-"Pinus densiflora"

han.800<-data.frame(do.call(rbind,list(fitquad800(dhan))))
han.800$DataSet<-"HAN_2"
han.800$Species<-"Pinus densiflora"


with(dhan,plot(Ts,lowAs,col="black",main="Japan (Han et al)",cex=2,pch=16))

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------
#Arctic Species
darc.n<-subset(darc,darc$maxCi>400) #remove few bad ACi curves

with(darc.n,plot(Ts,lowAs,col=USDA_Species_Code,pch=16,cex=2,main="Arctic All"))
with(subset(darc.n,darc.n$Year==2012),plot(Ts,lowAs,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data 2012"))
with(subset(darc.n,darc.n$Year==2013),plot(Ts,lowAs,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data 2013"))
with(subset(darc.n,darc.n$Year==2014),plot(Ts,lowAs,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data 2014"))
with(subset(darc.n,darc.n$Year==2015),plot(Ts,lowAs,col=USDA_Species_Code,pch=16,cex=2,main="Arctic Data 2015"))


fitquad300(subset(darc.n,darc.n$Year==2012))
fitquad300(subset(darc.n,darc.n$Year==2013))
fitquad300(subset(darc.n,darc.n$Year==2014))


arc.300<-data.frame(do.call(rbind,list(fit.nlme(darc.n,yvar="lowAs",rand="Year"))))
arc.300$DataSet<-"ARCTIC"
arc.300$Species<-"Arctic spp"
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

#-Spruce Anna et al 2015

spruce.300<-data.frame(do.call(rbind,list(fit.nlme(spru.b,yvar="lowAs",rand="Cohort_age"))))

#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

p_300<-data.frame(rbind.fill(arc.300,amz300,rwa300,wang300,cob300,dil300,dre300,med300,dway300,dtar300,dwal300,delw300,dstr.300,han.300))
names(p_300)[1:9] <- c("Aopt_275","Topt_275","b_275","Aopt_275.se","Topt_275.se","b.275.se","R2_275","S_275","pvalue_275")
p_300["Season"][is.na(p_300["Season"])] <-"No_Seasonal_Measurements"
p_300["Treatment"][is.na(p_300["Treatment"])] <-"No_Temperature_Treatments"


p_800<-data.frame(rbind.fill(amz800,rwa800,wang800,dil800,dre800,med800,dway800,dtar800,dwal800,dstr.800,han.800))
names(p_800)[1:9] <- c("Aopt_800","Topt_800","b_800","Aopt_800.se","Topt_800.se","b.800.se","R2_800","S_800","pvalue_800")
p_800["Treatment"][is.na(p_800["Treatment"])] <-"No_Temperature_Treatments"
p_800["Season"][is.na(p_800["Season"])] <-"No_Seasonal_Measurements"

tfits_275_oaz<-merge(p_300,p_800,by=c("DataSet","Species","Season","Treatment"),all=TRUE)

#to combine both Auz and out Auz datasets

#tfits_275_all<-rbind.fill(tfits_275_oaz,tfits_275)
#tfits_275_all["Temp_Treatment"][is.na(tfits_275_all["Temp_Treatment"])] <-"Ambient"
#-----------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------

tfits_275<-data.frame(rbind.fill(wtc1.300,wtc2.300,wtc3.300,wtc4.300,hfe.300,tmb.300,e.tr.300,rf.300,gww.300,cor300,ep.300))
names(tfits_275)[1:9] <- c("Aopt_275","Topt_275","b_275","Aopt_275.se","Topt_275.se","b.275.se","R2_275","S_275","pvalue_275")
tfits_275["Treatment"][is.na(tfits_275["Treatment"])] <-"No_Temperature_Treatments"
tfits_275["Month"][is.na(tfits_275["Month"])] <-"No_Seasonal_Measurements"


tfits_800<-data.frame(rbind.fill(wtc1.800,wtc2.800,wtc3.800,wtc4.800,e.tr.800,rf.800,gww.800,cor800,ep.800))
names(tfits_800)[1:9] <- c("Aopt_800","Topt_800","b_800","Aopt_800.se","Topt_800.se","b.800.se","R2_800","S_800","pvalue_800")
tfits_800["Treatment"][is.na(tfits_800["Treatment"])] <-"No_Temperature_Treatments"
tfits_800["Month"][is.na(tfits_800["Month"])] <-"No_Seasonal_Measurements"



#tfits_275_800<-tfits_275_800[-14]
tfits_275_Aus<-merge(tfits_275,tfits_800,by=c("Treatment","Month","Species","DataSet"),all=TRUE)

