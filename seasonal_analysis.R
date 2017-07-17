#seasonal responces

#function to fit fitted line in a plot 

plot_fit_line<-function(data,yvar,xvar,byvar=NULL,fitoneline=FALSE,linecol){
  
  data$yvar<-data[[yvar]]
  data$xvar<-data[[xvar]]
  
  
  
  if(!fitoneline){
    data$byvar<-data[[byvar]]
    
    lmlis1 <- lmList(yvar ~ xvar|byvar, data=data)
    liscoef <- coef(lmlis1)
    dfr<-split(data,factor(data[[byvar]]))
    
    for (i in 1:length(dfr)) {
      
      xmin <- min(dfr[[i]][[xvar]])
      xmax <- max(dfr[[i]][[xvar]])
      
      ablineclip(liscoef[i, 1], liscoef[i, 2], x1 = xmin, x2 = xmax,col=linecol,lwd=1)
    }
  }
  
  
  else{
    
    
    lmfit<-lm(yvar~xvar,data=data)
    
    xmin <- min(data[[xvar]])
    xmax <- max(data[[xvar]])
    
    ablineclip(a=coef(lmfit)[[1]],b=coef(lmfit)[[2]],x1=xmin,x2=xmax,col=linecol,lwd=1)
    
  }
  
  
}


#get a subset of data with seasonal measurements
path<-getwd()
tfits<-read.csv(paste0(path,"/Tables/netphoto_parameters_with_met_data.csv"))

tfits$MATSS<-tfits$BIO1/10
tfits$MAXSS<-tfits$BIO5/10

tfits$Tdiff<-tfits$MATSS-tfits$Tavg_30

tfits$PFT<-NA
tfits$PFT[tfits$Leafspan=="Deciduous" & tfits$Type=="Angiosperm"]<-"BDE_TE"
tfits$PFT[tfits$Leafspan=="Evergreen" & tfits$Type=="Angiosperm"]<-"BET_TE"
tfits$PFT[tfits$Tregion=="Tropical" & tfits$Type=="Angiosperm"]<-"TET"
tfits$PFT[tfits$Leafspan=="Evergreen" & tfits$Type=="Gymnosperm"]<-"NET"
tfits$PFT<-factor(tfits$PFT)

tfits$Tdiff<-with(tfits,MATSS-Tavg_30)

#---------------------
#---------------------
#non_seasonal<-subset(tfits, ! DataSet %in% c("S30",)) #not seasonal data

#tfits.1<-tfits[c("Species", "DataSet","Season_New","Temp_Treatment","Aopt","Topt","b","Aopt.se","Topt.se","Tavg_30",          
               #"b.se","PFT","BIO1","BIO5","MATSS","Cond.mean","Tdiff")]

tfits.1<-tfits[c("Species", "DataSet","Season_New","Temp_Treatment","Aopt","Topt","b","Aopt.se","Topt.se","Tavg_30",          
                 "b.se","PFT","BIO1","BIO5","MATSS")]

seasons<-subset(tfits.1,tfits.1$Season_New %in% c("Summer","Autumn","Spring","Winter","Dry","Wet"))
seasons.1<-subset(seasons,seasons$Temp_Treatment %in% c("ambient","Ambient","No_Temperature_Treatments")) #only ambient treatments
seasons.1$DataSet<-factor(seasons.1$DataSet)

seasons.2<-subset(seasons.1, ! DataSet %in% c("HAN_2","S30")) #not seasonal data

seasons.2<-subset(seasons.2,!is.na(Topt))
seasons.2$DataSet<-factor(seasons.2$DataSet)
seasons.2$Species<-factor(seasons.2$Species)

#seasons.2$A_std<-with(seasons.2,Aopt/Photo.mean)
#########################################################################################################################################
#Read biochemical data

tfits_VJ<-read.csv(paste0(path,"/Tables/biochemical_parameters_with_met_data.csv"))
#tfits_VJ<-tfits_VJ[c("Species", "DataSet","Season_New","Temp_Treatment","Vcmax25","Topt","b","Aopt.se","Topt.se","Tavg_30",          
#"b.se","PFT","BIO1","BIO5","BIO7")]


tfits_VJ$TPUFrac<-with(tfits_VJ,nTPU/nVcmax)
tfits_VJ$MATSS<-tfits_VJ$BIO1/10
tfits_VJ$MAXSS<-tfits_VJ$BIO5/10

tfits_VJ$Tdiff<-tfits_VJ$MATSS-tfits_VJ$Tavg_30

tfits_VJ$PFT<-NA
tfits_VJ$PFT[tfits_VJ$Leafspan=="Deciduous" & tfits_VJ$Type=="Angiosperm"]<-"BDE_TE"
tfits_VJ$PFT[tfits_VJ$Leafspan=="Evergreen" & tfits_VJ$Type=="Angiosperm"]<-"BET_TE"
tfits_VJ$PFT[tfits_VJ$Tregion=="Tropical" & tfits_VJ$Type=="Angiosperm"]<-"TET"
tfits_VJ$PFT[tfits_VJ$Leafspan=="Evergreen" & tfits_VJ$Type=="Gymnosperm"]<-"NET"
tfits_VJ$PFT<-factor(tfits_VJ$PFT)

#Get seasonal data seperately

seasons_vj<-subset(tfits_VJ,tfits_VJ$Season_New %in% c("Summer","Autumn","Spring","Winter","Dry","Wet"))
seasons_vj<-subset(seasons_vj,seasons_vj$Temp_Treatment %in% c("ambient","Ambient","No_Temperature_Treatments"))

seasons_vj.1<-subset(seasons_vj, ! DataSet %in% c("RF_AUS","SAVANNA","OND")) #not seasonal data
seasons_vj.1$DataSet<-factor(seasons_vj.1$DataSet)

#get subset of dataframe without NAs and remove TMB (no Topt275 data seasonally)

seasons_vj.2<-seasons_vj.1[c("Species","DataSet","Season_New","Tavg_30","Topt_275","TPUFrac")]
seasons_vj.2<-subset(seasons_vj.2,!is.na(Topt_275)&seasons_vj.2$DataSet!="TMB")
seasons_vj.2$DataSet<-factor(seasons_vj.2$DataSet)

###########################################################################################################################################

#Figure 1
#plot species measured at their home climates

t_home<-subset(tfits,  DataSet %in% c("GWW","RF_AMZ","SAVANNA","TARVAINEN","TMB",
                                             "WANG_ET_AL","RF_AUS","HAN_2","EucFace","Martijn_Slot","ARCTIC")) #these datasets are from
t_home.1<-subset(t_home,!t_home$Species %in% c("Eucalyptus globulus","Eucalyptus microcorys","Betula pendula","Toona sinensis") & !t_home$Season %in% c("WINTER","dry","summer","Nov"))
t_home.1$DataSet<-factor(t_home.1$DataSet)
library(RColorBrewer)
palette(rev(brewer.pal(12,"Paired")))
COL<-palette()[1:12]


#plot topt for photosynthesis
windows(50,50);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))

#with(t_home.1,plot(MAXSS,Topt,col=COL[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(10,35)))

with(t_home.1,plot(MATSS,Topt,col=COL[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(-20,35),ylim=c(10,35)))

ae<-subset(t_home.1,t_home.1$Topt.se<8) #remove one Se=9.4 to improve the clarity
adderrorbars(x=t_home.1$MATSS,y=t_home.1$Topt,SE=ae$Topt.se ,direction="updown")

#adderrorbars(x=t_home.1$MAXSS,y=t_home.1$Topt,SE=ae$Topt.se ,direction="updown")

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)

#with(t_home.1,plot(MATSS,Topt,col=COL[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,40)))
#adderrorbars(x=t_home.1$MATSS,y=t_home.1$Topt,SE=t_home.1$Topt.se ,direction="updown")
#magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T)
#legend(8,11,expression(italic(Eucalyptus)*~spp^1,Amazon~spp^2,Tropical~Rainforest~spp^3,Tropical~montane~rain~forest~spp,
                           # italic(E.~tetrodonta),italic(P.~abies),italic(E.~delegatensis),italic(P.~radiata),italic(P.~sylvestris)
                            #),col=COL,pch=16,ncol=2,bg="white",cex=0.8,title="Species",bty="n")


#title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T)
title(ylab=expression(Topt[A]~(degree*C)),
      xpd=NA,cex.lab=2)
title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T)

#add EucFace Topt (data from Jim)
#Topt estimated across seasons. THome similar to WTC3 

#points(17.1,21.88982,col=COL[9],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(10,35))
#adderrorbars(x=17.1,y=21.88982,SE=1.067941 ,direction="updown")

#add Tropical Panama Topts (data from Martijn Slot)
#Topt estimated across sites. 

points(3.6,20.4,col="black",pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(10,35))
adderrorbars(x=3.6,y=20.4,SE=2.6 ,direction="updown")


#add data from Austria. Source:Ann.For.Science (67, Wieser et al)
#Pinus cembra from central Austrian Alps
#Topt averaged over two seasons Fall and Summer (16.3 & 17.1 C respectively)
#T growth is ~8C > Tseedsource)

points(2.5,16.7,col=COL[12],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(10,35))
adderrorbars(x=2.5,y=16.7,SE=0.4 ,direction="updown")


legend("topleft",c("Arctic spp","Forest Red Gum","Australian semi-arid woodlands","Japanese Red Pine","Panamanian rainforest","Brazilian rainforest","Australian rainforest","Australian savanna","Spruce (Sweden)",
                "Alpine Ash","Scots pine (Finland)","Austrian stone pine"),
       col=COL,pch=16,ncol=2,bg="white",cex=.9,title="Dataset",bty="n")

#----------------------------------------------------------------------------
#par(new=T)
#smoothplot(Tavg_30,Topt,fittype="gam",pointcols=c(alpha("black",0.0)),
          # polycolor=c(alpha("lightgrey",0.7))
           #,cex=1,main="",
          # xlim=c(0,40),ylim=c(0,40),xlab="",ylab="",
          # data=t_home.1, kgam=6,axes=F)

#-----------------------------------------------------------------------------------------------------------------------------------

#plot relationships
#Seasonal Topts
palette(rev(brewer.pal(12,"Paired")))
COL<-palette()[1:12]


#plot topt for photosynthesis
windows(40,40);
par(cex.lab=1.5,mar=c(.5,.5,.5,.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))
#windows(40,40);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))


with(seasons.2,plot(Tavg_30,Topt,col=COL[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(10,40)))



#adderrorbars(x=seasons.2$Tavg_30,y=seasons.2$Topt,SE=seasons.2$Topt.se ,direction="updown")

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,majorn=5)
#axis.break(2, 20, style = "zigzag")

#add fittedline s to HFE species
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="ELLSWORTH"),yvar="Topt",xvar="Tavg_30",linecol=COL[1],fitoneline=T)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="HAN_ET_AL_A"),yvar="Topt",xvar="Tavg_30",linecol=COL[2],fitoneline=T)

plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="HFE_CG"),yvar="Topt",xvar="Tavg_30",byvar="Species",linecol=COL[3],fitoneline=F)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="MEDLYN"),yvar="Topt",xvar="Tavg_30",byvar="Species",linecol=COL[4],fitoneline=F)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="RF_AMZ"),yvar="Topt",xvar="Tavg_30",linecol=COL[5],fitoneline=T)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="RF_AUS"),yvar="Topt",byvar="Species",xvar="Tavg_30",linecol=COL[6],fitoneline=T)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="TASMANIA"),yvar="Topt",xvar="Tavg_30",linecol=COL[7],fitoneline=T)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="WTC1"),yvar="Topt",xvar="Tavg_30",linecol=COL[8],fitoneline=T)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="WTC2"),yvar="Topt",xvar="Tavg_30",linecol=COL[9],fitoneline=T)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="WTC3"),yvar="Topt",xvar="Tavg_30",linecol=COL[10],fitoneline=T)
plot_fit_line(data=subset(seasons.2,seasons.2$DataSet=="WTC4"),yvar="Topt",xvar="Tavg_30",linecol=COL[11],fitoneline=T)


#add HAN_2 dataset (due to changes in color, add seperately)

with(subset(seasons,seasons$DataSet=="HAN_2"),points(Tavg_30,Topt,pch=16,col=COL[12],cex=1.5))
plot_fit_line(data=subset(seasons,seasons$DataSet=="HAN_2"),yvar="Topt",xvar="Tavg_30",linecol=COL[12],fitoneline=T)


legend(0,41,c("Duke pine","Hinoki cypress","CG-Au","Maritime pine","Brazilian rainforest",
              "Australian rainforest","Tasmanian Blue Gum (A)","Sydney Blue Gum","Tasmanian Blue Gum (B)","Forest Red Gum","Parramatta Red Gum","Japanese Red Pine") ,col=COL,pch=16,ncol=2,bg="white",
       cex=.8,title="Dataset",bty="n") 



#legend(0,40,c("Duke forest NC, USA","Japan (b)","CG-Au","France (a)","Brazilian rainforest",
                 #  "Australian rainforest","Tasmania","WTC-V1","WTC-V2","WTC-V3","WTC-V4","Japan (a)") ,col=COL,pch=16,ncol=3,bg="white",
      # cex=.8,title="Dataset",bty="n") 
                                                                          

#legend("topright",expression(bold(a)),cex=1.2,bty="n")

#title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T)
title(ylab=expression(Topt[A]~(degree*C)),
      xpd=NA,cex.lab=2)
title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T)


abline(0,1,lty=2)
#add Topts from selected species which Topt doesn't exixts

# WTC3 summer

points(22.63,18,pch=2,cex=1.5,col=COL[10],lwd=3)
#lmx<-lm(c(wtc3.topts[[4]][[1]],18)~c(wtc3.topts[[28]][[1]],wtc3.topts[[28]][[2]]))

#ablineclip(lmx[1], lmx[2], x1 = wtc3.topts[[28]][[1]], x2 = wtc3.topts[[28]][[2]],col=COL[10],lwd=3,lty=2)


#problem in Medlyn et al fit. manually add fittedline for Tamjout
#data<-subset(seasons.2,seasons.2$DataSet=="MEDLYN")
#lmlis1 <- lmList(Topt ~ Tavg_30|Species, data=data)
#liscoef <- coef(lmlis1)
#xmin <- min(data$Tavg_30)
#xmax <- max(data$Tavg_30)
#ablineclip(liscoef[2, 1], liscoef[2, 2], x1 = xmin, x2 = xmax,col=COL[4],lwd=3)


aeb.1<-subset(seasons.2,seasons.2$Topt.se<8) #remove two SEs to improve clarity
adderrorbars(x=seasons.2$Tavg_30,y=seasons.2$Topt,SE=aeb.1$Topt.se ,direction="updown")

#add best fitted linear mix model
lmm_topt<-lmer(Topt~Tavg_30+(1|Species),data=subset(seasons.2,seasons.2$DataSet != "HAN_ET_AL_A")) #random intercept model. Dataset as random factor
lmm_topt.1<-lmer(Topt~Tavg_30+(Tavg_30|Species),data=subset(seasons.2,seasons.2$DataSet != "HAN_ET_AL_A")) #model with random slope and intercept  
                                                                                                                #HAN ET AL A dataset not included in model fitting
                                                                                                                # due to uncertainty in Tgrowth 
Anova(lmm_topt,test="F")
Anova(lmm_topt.1,test="F")

anova(lmm_topt,lmm_topt.1) #random slope no need
rsquared.glmm(lmm_topt)

ablineclip(round(summary(lmm_topt)$coefficients[1],2), round(summary(lmm_topt)$coefficients[2],2), x1 = min(seasons.2$Tavg_30-2), x2 = max(seasons.2$Tavg_30+2),col="black",lwd=3,lty=2)
#text(34,13,expression(hat(beta)==0.38~(degree*C~degree*C^-1)))
#text(37,11,expression(italic(R^2)==0.36))

text(34,13,expression(T[optA]==19.1+0.34~T[Growth])) #95%CI
text(32,12,expression(R^2==0.39)) 
############################################################################################################################################
############################################################################################################################################

#------------------------------------------------

#plot Topt for photosynthesis at a common Ci vs Tgrowth

COL.1<-palette()[c(1:5,8:11)]

windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(seasons_vj.2,plot(Tavg_30,Topt_275,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(10,40)))

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
#axis.break(2, 20, style = "zigzag")
#add best fitted linear mix model



plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="ELLSWORTH"),yvar="Topt_275",xvar="Tavg_30",linecol=COL.1[1],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="HAN_ET_AL_A"),yvar="Topt_275",xvar="Tavg_30",linecol=COL.1[2],fitoneline=T)

#add HFE (some issue, function not work)
hfe<-subset(seasons_vj.2,seasons_vj.2$DataSet=="HFE_CG" & seasons_vj.2$Specie!="Eucalyptus saligna")
dfr_hfe<-split(hfe,factor(hfe[["Species"]]))

lmlis1_hfe <- lmList(Topt_275 ~ Tavg_30|Species, data=hfe)

liscoef_hfe <- coef(lmlis1_hfe)
xmin <- min(hfe$Tavg_30)
xmax <- max(hfe$Tavg_30)


for (i in 1:length(dfr_hfe)) {
  
  xmin <- min(dfr_hfe[[i]]$Tavg_30)
  xmax <- max(dfr_hfe[[i]]$Tavg_30)
  
  ablineclip(liscoef_hfe[i, 1], liscoef_hfe[i, 2], x1 = xmin, x2 = xmax,col=COL.1[3],lwd=1)
}


#add Medlyn et al (some issue, function not work)
med<-subset(seasons_vj.2,seasons_vj.2$DataSet=="MEDLYN")
dfr_med<-split(med,factor(med[["Species"]]))

lmlis1_med <- lmList(Topt_275 ~ Tavg_30|Species, data=med)

liscoef_med <- coef(lmlis1_med)
xmin <- min(med$Tavg_30)
xmax <- max(med$Tavg_30)


for (i in 1:length(dfr_med)) {
  
  xmin <- min(dfr_med[[i]]$Tavg_30)
  xmax <- max(dfr_med[[i]]$Tavg_30)
  
  ablineclip(liscoef_med[i, 1], liscoef_med[i, 2], x1 = xmin, x2 = xmax,col=COL.1[4],lwd=1)
}
#plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="HFE_CG"),yvar="Topt_275",xvar="Tavg_30",byvar="Species",linecol=COL[3],fitoneline=F)
#plot_fit_line(data=subset(seasons_vj.1,seasons_vj.1$DataSet=="MEDLYN"),yvar="Topt",xvar="Tavg_30",byvar="Species",linecol=COL[4],fitoneline=F)
plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="RF_AMZ"),yvar="Topt_275",xvar="Tavg_30",linecol=COL.1[5],fitoneline=T)
#plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="TMB"),yvar="Topt_275",xvar="Tavg_30",linecol=COL[7],fitoneline=T)

plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="WTC1"),yvar="Topt_275",xvar="Tavg_30",linecol=COL.1[6],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="WTC2"),yvar="Topt_275",xvar="Tavg_30",linecol=COL.1[7],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="WTC3"),yvar="Topt_275",xvar="Tavg_30",linecol=COL.1[8],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.2,seasons_vj.2$DataSet=="WTC4"),yvar="Topt_275",xvar="Tavg_30",linecol=COL.1[9],fitoneline=T)

#legend("bottomright",expression(italic(P.~taeda),italic(C.~obtusa),italic(Eucalyptus~spp^1),italic(P.~pinaster^2),Amazon~spp,
                            #italic(E.~saligna),italic(E.~globulus),italic(E.~tereticornis),italic(E.~parramattensis)),col=COL.1,pch=16,ncol=2,bg="white",cex=.8,title="Species",bty="n")

legend(12,16,c("Duke pine","Hinoki cypress","CG-Au","Maritime pine","Brazilian rainforest",
               "Sydney Blue Gum","Tasmanian Blue Gum (B)","Forest Red Gum","Parramatta Red Gum") ,col=COL.1,pch=16,ncol=2,bg="white",
       cex=.8,title="Dataset",bty="n") 

#legend(0,41,c("Duke pine","Hinoki cypress","CG-Au","Maritime pine","Brazilian rainforest",
              #"Australian rainforest","Tasmanian Blue Gum (A)","Sydney Blue Gum","Tasmanian Blue Gum (B)","Forest Red Gum","Parramatta Red Gum","Japanese Red Pine") ,col=COL,pch=16,ncol=2,bg="white",
      # cex=.8,title="Dataset",bty="n") 



#legend("topright",expression(bold(b)),cex=1.2,bty="n")
title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T)

#title(ylab=expression(Topt[A275]~(degree*C)),
      #xpd=NA,cex.lab=2)
title(ylab=expression(Topt[A275]~(degree*C)),
      xpd=NA,cex.lab=2)

#par(new=T)
#add errorbars (remove very high error bars to improve clarity)
ebd<-subset(seasons_vj.1,seasons_vj.1$Topt_275.se<10)
adderrorbars(x=ebd$Tavg_30,y=ebd$Topt_275,SE=ebd$Topt_275.se ,direction="updown")

abline(0,1,lty=2)
#add best fitted linear mix model)
lmm_275<-lmer(Topt_275~Tavg_30+(1|Species),data=subset(seasons_vj.2,seasons_vj.2$DataSet != "HAN_ET_AL_A")) #random intercept model. Dataset as random factor
lmm_275.1<-lmer(Topt_275~Tavg_30+(Tavg_30|Species),data=seasons_vj.2) #model with random slope and intercept 

Anova(lmm_275,test="F")
Anova(lmm_275.1,test="F")

anova(lmm_275,lmm_275.1) #random slope no need
rsquared.glmm(lmm_275)

#add best fitted linear mix model
ablineclip(22.1976 , 0.2517, x1 = min(seasons_vj.2$Tavg_30-2), x2 = max(seasons_vj.2$Tavg_30+2),col="black",lwd=3,lty=2)
#text(6,40,expression(hat(beta)==0.27~(degree*C~degree*C^-1)))
#text(6,38,expression(italic(R^2)==0.18))

text(8,40,expression(Topt[A275]==22.1+0.25~T[Growth])) #95%CI
text(6,38,expression(R^2==0.18)) 
#-------------------------------------------------------------------------------------------------------
#to get fitted GAM plots fot Topt vs Tgrowth
windows(5,5)
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))
#add fitted GAM function to the same plot
#par(new=T)
smoothplot(Tavg_30,Topt,pointcols=c(alpha("black",0.0)),fitmethod="lm",
           linecol=c("black"),polycolor=c(alpha("lightgrey",0.5))
           ,cex=1,main="",
           xlim=c(0,40),ylim=c(10,40),xlab="",ylab="",
           data=seasons.2, kgam=3,axes=F)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
title(ylab=expression(Topt[A]~(degree*C)),
      xpd=NA,cex.lab=2)
title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.9,line=3)

#plot_gam(Tavg_30,Topt,data=seasons.2,randommethod="gamm", R="Species",k=6,fittype="gam",linecols="black",pointcols = c(alpha("black",0.0)),cex=2,band=T,axes=F,fitoneline = T)
#par(new=T)
#to get fitted GAM plot Topt275 vs Tgrowth
windows(5,5)
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))
smoothplot(Tavg_30,Topt_275,pointcols=c(alpha("black",0.0)),fitmethod="lm",
           linecol=c("black"),polycolor=c(alpha("lightgrey",0.5))
           ,cex=1,main="",
           xlim=c(0,40),ylim=c(10,40),xlab="",ylab="",
           data=seasons_vj.2, kgam=3,axes=F)
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
title(ylab=expression(Topt[A275]~(degree*C)),
      xpd=NA,cex.lab=2)
title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.9,line=3)

############################################################################################################################################
############################################################################################################################################

#Plot data from commongardens
#test for adaptation of Topt of Photosynthesis
#adapt<-subset(tfits,  DataSet %in% c("HFE_CG","DILLAWAY","S30","MEDLYN","DREYER","RF_MON",
                                     #"WANG_ET_AL","GWW","MIKE_ASPINWALL","STRASSEMEYER_ET_AL")) #these datasets are from


adapt<-subset(tfits,  DataSet %in% c("HFE_CG","DILLAWAY","S30","MEDLYN","DREYER","RF_MON",
                                     "MIKE_ASPINWALL","GREAT")) #these datasets are from commongardens or species grown in comon Tgrowth
adapt<-adapt[order(adapt$DataSet),]
adapt$DataSet<-as.character(adapt$DataSet)

adapt$DataSet[1:12]<-rep(c("DILLAWAY_1","DILLAWAY_2","DILLAWAY_3"),4)

#get only summer data is seasonally measured and ambient temperature treatments from warming experiments
# either commongardens or contrasting spp grown in a common temperatures.
adapt.1<-subset(adapt,!Season_New %in% c("Autumn","Winter","Spring") & !Temp_Treatment %in% c(32,"Elevated"))
adapt.1<-subset(adapt.1,adapt.1$Species!="Toona sinensis") #Anjelica's dataset. Climate of origine not sure
adapt.1$DataSet<-factor(adapt.1$DataSet)

#adapt.1$MATSS<-with(adapt.1,BIO1/10)
#adapt.1$MAXWQ<-with(adapt.1,BIO5/10)
COL.2<-palette()[c(1,2,5,6,7,3,4,8,9,10,11,12)]

windows(100,100);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

#windows(40,40);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))


with(adapt.1,plot(MAXSS,Topt,col=COL.2[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(15,40),ylim=c(10,40)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)


plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="DILLAWAY_1"),yvar="Topt",xvar="MAXSS",linecol=COL.2[1],fitoneline=T)
plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="DILLAWAY_2"),yvar="Topt",xvar="MAXSS",linecol=COL.2[2],fitoneline=T)
plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="DILLAWAY_3"),yvar="Topt",xvar="MAXSS",linecol=COL.2[3],fitoneline=T)
plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="DREYER"),yvar="Topt",xvar="MAXSS",linecol=COL.2[4],fitoneline=T)
plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="GREAT"),yvar="Topt",xvar="MAXSS",linecol=COL.2[5],fitoneline=T)

plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="HFE_CG"),yvar="Topt",xvar="MAXSS",linecol=COL.2[6],fitoneline=T)
plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="MEDLYN"),yvar="Topt",xvar="MAXSS",linecol=COL.2[7],fitoneline=T)
plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="MIKE_ASPINWALL"),yvar="Topt",xvar="MAXSS",linecol=COL.2[8],fitoneline=T)
plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="RF_MON"),yvar="Topt",xvar="MAXSS",linecol=COL.2[9],fitoneline=T)
plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="S30"),yvar="Topt",xvar="MAXSS",linecol=COL.2[10],fitoneline=T)

#plot_fit_line(data=subset(adapt.1,adapt.1$DataSet=="WANG_ET_AL"),yvar="Topt",xvar="MATSS",linecol=COL.2[12],fitoneline=T)


#legend("bottomright",levels(adapt.1$DataSet),col=COL.2,pch=16,ncol=2,bg="white",cex=1,title="DataSet",bty="n")

title(ylab=expression(Topt[A]~(degree*C)),
      xpd=NA,cex.lab=2)
title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T)


#title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T)

ebd.1<-subset(adapt.1,adapt.1$Topt.se<7)
adderrorbars(x=ebd.1$MAXSS,y=ebd.1$Topt,SE=ebd.1$Topt.se ,direction="updown")

abline(0,1,lty=2)

#add other data from different home climate
#adapt_rest<-subset(tfits,  DataSet %in% c("RF_AMZ","TMB","HAN_2","SAVANNA","TARVAINEN","WALCROFT")) #these datasets are from
#adapt_rest$DataSet<-factor(adapt_rest$DataSet)
#Add Amazone (No seasonal difference in Topt)

#with(adapt_rest,points(MATSS,Topt,col="black",pch=c(0,1,2,5,7,9)[DataSet],cex=1.5,ylab="",xlab=""))
#adderrorbars(x=adapt_rest$MATSS,y=adapt_rest$Topt,
            #SE=adapt_rest$Topt.se,direction="updown")

#legend("bottomright",expression(Common~garden(D^1),Common~garden(D^2),Common~garden(D^3),Deciduous^1~spp,
                            #italic(Eucalyptus~spp^2),italic(P.~pinaster^3),italic(Corymbia)*~provenances^4,Montane~rain~forest~spp,
                            #italic(Eucalyptus~spp^5)),col=COL.2,pch=16,ncol=2,bg="white",cex=.8,title="Species",bty="n")

legend(29,21,c("CG-Il (USA)","CG-NW (USA)","CG-SW (USA)","France (B)","Australian Eucalyptus (A)","CG-Au","Maritime Pine","Australian Corymbia","Tropical montane forest", 
               "Australian Eucalyptus (B)") ,col=COL.2,pch=16,ncol=1,bg="white",
       cex=.8,title="Dataset",bty="n") 

#par(new=T)
#smoothplot(MATSS,Topt,pointcols=c(alpha("black",0.0)),
          # linecol=c("black"),polycolor=c(alpha("lightgrey",0.5))
          # ,cex=1,main="",xlim=c(0,30),ylim=c(10,40),
          # xlab="",ylab="",
           #data=adapt.1, kgam=2,axes=F)
#fit mix model

lmm_adap<-lmer(Topt~MAXSS+(1|DataSet),data=adapt.1) #random intercept model. Dataset as random factor
summary(lmm_adap)
Anova(lmm_adap,test="F")
rsquared.glmm(lmm_adap)

ablineclip(24.4745, 0.1232, x1 = min(adapt.1$MAXSS-2), x2 = max(adapt.1$MAXSS+2),col="black",lwd=3,lty=2)
#text(18,40,expression(hat(beta)==0.09~(-0.16-0.34)~degree*C^-1))

text(20,40,expression(T[optA]==24.47+0.12~T[Home])) #95%CI
text(18,38,expression(R^2==0.03^NS)) 
#text(18,38,expression(italic(R^2)==0.02))


############################################################################################################################################
#############################################################################################################

# Plot Biochemical parameters

adapt_bio<-subset(tfits_VJ,  DataSet %in% c("HFE_CG","DILLAWAY","MEDLYN","DREYER","RF_MON",
                                     "MIKE_ASPINWALL")) #these datasets are from
adapt_bio<-adapt_bio[order(adapt_bio$DataSet,adapt_bio$Season),]
adapt_bio$DataSet<-as.character(adapt_bio$DataSet)

adapt_bio$DataSet[1:12]<-c(rep("DILLAWAY_1",4),rep("DILLAWAY_2",4),rep("DILLAWAY_3",4))

adapt_bio.1<-subset(adapt_bio,!Season_New %in% c("Autumn","Winter","Spring") & !Temp_Treatment %in% c(32,"Elevated"))
adapt_bio.1<-subset(adapt_bio.1,!Season %in% c("Jul","Jun") & !is.na(Topt_275)) #remove Jun n Jul data in Medlyn et al

adapt_bio.1$DataSet<-factor(adapt_bio.1$DataSet)
#adapt_bio.1$JVr<-with(adapt_bio.1,Jmax25/Vcmax25)
adapt_bio.1<-subset(adapt_bio.1,adapt_bio.1$Species!="Toona sinensis") #Anjelica's dataset. Climate of origine not sure
#windows(40,40);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))

windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

COL.2<-palette()[c(1,2,5,6,7,3,4,8,9,11,12,10)]
with(adapt_bio.1,plot(MAXSS,JVr,col=COL.2[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(10,35),ylim=c(0,3)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)



plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_1"),yvar="JVr",xvar="MAXSS",linecol=COL.2[1],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_2"),yvar="JVr",xvar="MAXSS",linecol=COL.2[2],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_3"),yvar="JVr",xvar="MAXSS",linecol=COL.2[3],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="DREYER"),yvar="JVr",xvar="MAXSS",linecol=COL.2[4],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="HFE_CG"),yvar="JVr",xvar="MAXSS",linecol=COL.2[6],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="MEDLYN"),yvar="JVr",xvar="MAXSS",linecol=COL.2[7],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="MIKE_ASPINWALL"),yvar="JVr",xvar="MAXSS",linecol=COL.2[8],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="RF_MON"),yvar="JVr",xvar="MAXSS",linecol=COL.2[9],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="WANG_ET_AL"),yvar="JVr",xvar="MAXSS",linecol=COL.2[11],fitoneline=T)

legend("topright",expression(bold(b)),cex=1.2,bty="n")

#title(ylab=expression(Topt[A275]~(degree*C)),
#xpd=NA,cex.lab=2)
title(ylab=expression(Topt[A275]~(degree*C)),
      xpd=NA,line=-32,srt=90,cex.lab=2)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.9,line=3)



#title(ylab=expression(Topt[A275]~(degree*C)),
      #xpd=NA,cex.lab=2)


#title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T)

ebd.2<-subset(adapt_bio.1,adapt_bio.1$Topt_275.se<7)
adderrorbars(x=ebd.2$MATSS,y=ebd.2$Topt_275,SE=ebd.2$Topt_275.se ,direction="updown")

abline(0,1,lty=2)

#add Topts which were above the measured temperature range

#HFE E. saligna in summer
points(12.7,32,pch=1,cex=1.5,col=COL.2[6],lwd=2)

#Dillaway et al Populus tremuloides Nothern Wisconsin commongarden

points(7.9,32,pch=1,cex=1.5,col=COL.2[2],lwd=2)


#Dillaway et al Liquidambar styraciflua Nothern Wisconsin commongarden

points(12.2,32,pch=1,cex=1.5,col=COL.2[2],lwd=2)

legend("bottomright",expression(Common~garden(D^1),Common~garden(D^2),Common~garden(D^3),Deciduous^1~spp,Woodlands~spp,
                                italic(Eucalyptus~spp^2),italic(P.~pinaster^3),italic(Corymbia)*~provenances^4,Montane~rain~forest~spp,
                                italic(Eucalyptus~spp^5),Deciduous^6~spp,italic(Betula)*~and~italic(Pinus~spp)),col=COL.2,pch=16,ncol=2,bg="white",cex=.8,title="Species",bty="n")


#add other data from different home climate
adapt_bio_rest<-subset(tfits_VJ,  DataSet %in% c("RF_AMZ","HAN_2","SAVANNA","TARVAINEN","WALCROFT")) #these datasets are from

adapt_bio_rest$DataSet<-factor(adapt_bio_rest$DataSet)
#Add Amazone (No seasonal difference in Topt)

with(adapt_bio_rest,points(MATSS,Topt_275,col="black",pch=c(0,1,2,5,9)[DataSet],cex=1.5,ylab="",xlab=""))
adderrorbars(x=adapt_bio_rest$MATSS,y=adapt_bio_rest$Topt_275,
             SE=adapt_bio_rest$Topt_275.se,direction="updown")

c(0,1,2,5,7,9)
#############################################################################################################
#############################################################################################################

#Plot JV ratio (seasonal)

seasons_vj.3<-seasons_vj.1[c("Species","DataSet","Season_New","Topt_275","Topt_275.se","Tavg_30","Vcmax25","Jmax25","EaV","EaV.se","EaJ","EaJ.se","Rday","delsV",
                             "delsV.se","ToptV","delsJ","delsJ.se","ToptJ","JVr","JVr_SE","MATSS","MAXSS","Vcmax25.se","Jmax25.se","TPUFrac" )]
seasons_vj.3$JVr.1<-with(seasons_vj.3,Jmax25/Vcmax25)
#seasons_vj.3<-subset(seasons_vj.3,!is.na(JVr))
seasons_vj.3$DataSet<-factor(seasons_vj.3$DataSet)

seasons_vj.3$JVr<-ifelse(is.na(seasons_vj.3$JVr),seasons_vj.3$JVr.1,seasons_vj.3$JVr)
#------------------------------------------------

#plot JV Ratio vs Tgrowth

COL.1<-palette()[c(1:4,6,5,8:11)]

windows(100,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,5),cex.axis=1.5,las=1,mfrow=c(1,2))

with(seasons_vj.3,plot(Tavg_30,JVr,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,4)))

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
#axis.break(2, 20, style = "zigzag")
adderrorbars(x=seasons_vj.3$Tavg_30,y=seasons_vj.3$JVr,
SE=seasons_vj.3$JVr_SE,direction="updown")

plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="ELLSWORTH"),yvar="JVr",xvar="Tavg_30",linecol=COL.1[1],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="HAN_ET_AL_A"),yvar="JVr",xvar="Tavg_30",linecol=COL.1[2],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="HFE_CG"),yvar="JVr",xvar="Tavg_30",byvar="Species",linecol=COL.1[3],fitoneline=F)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="MEDLYN"),yvar="JVr",xvar="Tavg_30",byvar="Species",linecol=COL.1[4],fitoneline=F)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="TMB"),yvar="JVr",xvar="Tavg_30",linecol=COL.1[6],fitoneline=T)

plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="RF_AMZ"),yvar="JVr",xvar="Tavg_30",linecol=COL.1[5],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC1"),yvar="JVr",xvar="Tavg_30",linecol=COL.1[7],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC2"),yvar="JVr",xvar="Tavg_30",linecol=COL.1[8],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC3"),yvar="JVr",xvar="Tavg_30",linecol=COL.1[9],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC4"),yvar="JVr",xvar="Tavg_30",linecol=COL.1[10],fitoneline=T)



legend("topleft",c("Duke pine","Hinoki cypress","CG-Au","Maritime pine","Brazilian rainforest","Australian alpine",
              "Sydney Blue Gum","Tasmanian Blue Gum (B)","Forest Red Gum","Parramatta Red Gum") ,col=COL.1,pch=16,ncol=2,bg="white",
       cex=1,title="Dataset",bty="n") 


#legend("topleft",c("Duke forest NC, USA","Japan (b)","CG-Au","France (a)","Brazilian rainforest","Australian alpine",
              # "WTC-V1","WTC-V2","WTC-V3","WTC-V4") ,col=COL.1,pch=16,ncol=2,bg="white",
       #cex=1,title="Dataset",bty="n") 

title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.2,line=3)
title(ylab=expression(JV[r]),
      xpd=NA,cex.lab=2)

#lm5<-lmer(JVr~Tavg_30+(1|Species),data=seasons_vj.3) 
#Anova(lm5,test="F")
#qqPlot(residuals(lm5))
#plot(resid(lm5) ~ fitted(lm5))
#abline(h=0)
#Rsq_marginal<-rsquared.glmm(lm5)[[4]]


lm6<-lmer(JVr~Tavg_30+(Tavg_30|Species),data=seasons_vj.3) 
Anova(lm6,test="F")
#qqPlot(residuals(lm6))
#plot(resid(lm5) ~ fitted(lm6))
#abline(h=0)
Rsq_marginal<-rsquared.glmm(lm6)[[4]]

ablineclip(2.3207, -0.0338, x1 = min(seasons_vj.3$Tavg_30-2), x2 = max(seasons_vj.3$Tavg_30+2),col="black",lwd=3,lty=2)
#text(30,0.5,expression(hat(beta)==-0.035~(-0.02-0.05)~degree*C^-1)) #95%CI
#text(36,0.3,expression(italic(R^2)==0.31))

text(30,0.5,expression(JV[r]==2.32-0.034~T[Growth])) #95%CI
text(30,0.3,expression(R^2==0.30)) 

legend("topright",expression(bold(a)),bty="n",cex=2)
#plot JVr for P. taeda in winter (from fitted Vcmax25 and Jmax25 values)
#points(11.294165,55.97584/30.95969,pch=16,cex=1.5,col=COL.1[1])
#par(new=T)
#smoothplot(Tavg_30,JVr,pointcols=c(alpha("black",0.0)),
           #linecol=c("black"),polycolor=c(alpha("lightgrey",0.5))
           #,cex=1,main="",
           #xlim=c(0,40),ylim=c(0,5),xlab="",ylab="",
           #data=seasons_vj.3, kgam=3,axes=F)
#windows(20,20);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

#with(seasons_vj.3,plot(ToptJ~JVr,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,5),ylim=c(20,50)))
#magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
#magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)

#adderrorbars(x=seasons_vj.3$Vcmax25,y=seasons_vj.3$Jmax25,
             #SE=seasons_vj.3$Jmax25.se ,direction="updown")

#abline(0,1,lty=2)

#lm_vj<-lmer(Jmax25~Vcmax25+(1|DataSet),data=seasons_vj.3)
#rsquared.glmm(lm_vj)[[4]]

#ablineclip(22.2254, 1.3705, x1 = min(seasons_vj.3$Vcmax25-10), x2 = max(seasons_vj.3$Vcmax25+10),col="black",lwd=3,lty=2)
#text(10,300,expression(hat(beta)==0.27~(degree*C~degree*C^-1)))
#text(10,290,expression(italic(R^2)==0.18))
#plot relationship between Topt_275 vs JVr

#-------------------------------------------------------------------------------------------------------
#plot JVr vs Topt275

#windows(40,40);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(seasons_vj.3,plot(Topt_275~JVr,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0.5,3),ylim=c(10,40)))

magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T,at=1:13)

aeb.1<-subset(seasons_vj.3,seasons_vj.3$Topt_275.se<6)
adderrorbars(x=seasons_vj.3$JVr,y=seasons_vj.3$Topt_275,
             SE=aeb.1$Topt_275.se,direction="updown")


adderrorbars(x=seasons_vj.3$JVr,y=seasons_vj.3$Topt_275,
             SE=seasons_vj.3$JVr_SE,direction="leftright")

#title(ylab=expression(Topt[A275]~(degree*C)),
      #xpd=NA,cex.lab=2)

title(ylab=expression(Topt[A275]~(degree*C)),
      xpd=NA,cex.lab=2,srt=-10,line=-27)

title(xlab=expression(JV[r]),
      xpd=NA,cex.lab=2,adj=0.5,line=3)

lmm_tev<-lmer(Topt_275~JVr+(1|Species),data=seasons_vj.3) #random intercept model. Dataset as random factor
lmm_tev.1<-lmer(Topt_275~JVr+(Tavg_30|Species),data=seasons_vj.3) #model with random slope and intercept 

Anova(lmm_tev,test="F")
Anova(lmm_tev.1,test="F")
rsquared.glmm(lmm_tev)

anova(lmm_tev,lmm_tev.1) #random slope no need
#qqPlot(residuals(lmm_tev)) #nice
#plot(resid(lmm_tev) ~ fitted(lmm_tev)) #not too bad
#abline(h=0)

ablineclip(37.315, -6.387, x1 = min(seasons_vj.3$JVr,na.rm=T), x2 = max(seasons_vj.3$JVr,na.rm=T),col="black",lwd=3,lty=2)
#text(2.3,40,expression(hat(beta)==-6.39~(-4.04-8.95)~degree*C)) #95%CI
#text(2.5,38,expression(italic(R^2)==0.40))

text(2.4,13,expression(Topt[A275]==37.3-6.4~JV[r])) #95%CI
text(2.4,12,expression(R^2==0.40))
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="ELLSWORTH"),yvar="Topt_275",xvar="JVr",linecol=COL.1[1],fitoneline=T)
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="HAN_ET_AL_A"),yvar="Topt_275",xvar="JVr",linecol=COL.1[2],fitoneline=T)
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="HFE_CG"),yvar="Topt_275",xvar="JVr",byvar="Species",linecol=COL.1[3],fitoneline=F)
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="MEDLYN"),yvar="Topt_275",xvar="JVr",byvar="Species",linecol=COL.1[4],fitoneline=F)
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="TMB"),yvar="Topt_275",xvar="JVr",linecol=COL.1[6],fitoneline=T)

#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="RF_AMZ"),yvar="Topt_275",xvar="JVr",linecol=COL.1[5],fitoneline=T)
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC1"),yvar="Topt_275",xvar="JVr",linecol=COL.1[7],fitoneline=T)
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC2"),yvar="Topt_275",xvar="JVr",linecol=COL.1[8],fitoneline=T)
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC3"),yvar="Topt_275",xvar="JVr",linecol=COL.1[9],fitoneline=T)
#plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC4"),yvar="Topt_275",xvar="JVr",linecol=COL.1[10],fitoneline=T)

legend("topright",expression(bold(b)),bty="n",cex=2)

############################################################################################################
#############################################################################################################
#plot delsV

windows(100,100);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,5),cex.axis=1.5,las=1,mfrow=c(2,2))


with(seasons_vj.3,plot(Tavg_30,delsV*1000,col=COL.1[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(600,670)))
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T,at=1:13)
adderrorbars(x=seasons_vj.3$Tavg_30,y=seasons_vj.3$delsV*1000,
             SE=seasons_vj.3$delsV.se*1000,direction="updown")

title(ylab=expression(Delta*S[Vcmax]~(J~mol^-1~degree*C^-1)),
      xpd=NA,cex.lab=1.5)
legend("topright",expression(bold(a)),cex=2,bty="n")


#plot delsJ

with(seasons_vj.3,plot(Tavg_30,delsJ*1000,col=COL.1[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(600,670)))

magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T,at=1:13)
adderrorbars(x=seasons_vj.3$Tavg_30,y=seasons_vj.3$delsJ*1000,
             SE=seasons_vj.3$delsJ.se*1000,direction="updown")

title(ylab=expression(Delta*S[Jmax]~(J~mol^-1~degree*C^-1)),
      xpd=NA,cex.lab=1.5,srt=-660,line=-19)
legend("topright",expression(bold(b)),cex=2,bty="n")

#plot ToptV
#windows(100,100);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))

with(seasons_vj.3,plot(Tavg_30,ToptV,col=COL.1[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(20,50)))

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)

title(ylab=expression(Topt[Vcmax]~(degree*C)),
      xpd=NA,cex.lab=1.5)
legend("topright",expression(bold(c)),cex=2,bty="n")
#windows(100,100);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))

#plot ToptV

with(seasons_vj.3,plot(Tavg_30,ToptJ,col=COL.1[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(20,50)))

magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T,at=1:13)

title(ylab=expression(Topt[Jmax]~(degree*C)),
      xpd=NA,cex.lab=1.5,srt=-660,line=-19)

title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.2,line=3)

title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.9,line=3)
legend("topright",expression(bold(d)),cex=2,bty="n")

#------------------------------------------------------------------------------------------------------------


#plot Tgrowth vs EaV

#windows(100,100);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))
windows(100,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,5),cex.axis=1.5,las=1,mfrow=c(1,2))

with(subset(seasons_vj.3,seasons_vj.3$EaV>20),plot(Tavg_30,EaV,col=COL.1[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,30),ylim=c(0,100)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)

title(ylab=expression(Ea[Vcmax]~(kJ~mol^-1)),
      xpd=NA,cex.lab=2)



lmm_ea<-lmer(EaV~Tavg_30+(1|Species),data=subset(seasons_vj.3,seasons_vj.3$EaV>20)) #random intercept model. Dataset as random factor
lmm_ea.1<-lmer(EaV~Tavg_30+(Tavg_30|Species),data=subset(seasons_vj.3,seasons_vj.3$EaV>20)) #model with random slope and intercept 

Anova(lmm_ea,test="F")
Anova(lmm_ea.1,test="F")
rsquared.glmm(lmm_ea)

ablineclip(45.606, 0.719, x1 = min(seasons_vj.3$Tavg_30,na.rm=T)-2, x2 = max(seasons_vj.3$Tavg_30,na.rm=T)+2,col="black",lwd=3,lty=2)
text(7,100,expression(Ea[Vcmax]==45.6-0.7~T[growth])) #95%CI
text(7,96,expression(italic(R^2)==0.11))


adderrorbars(x=seasons_vj.3$Tavg_30,y=seasons_vj.3$EaV,
             SE=aeb.1$EaV.se,direction="updown")
legend("topright",expression(bold(a)),cex=2,bty="n")

#--------------------------------------------------------------------------------------------------------------------------------------

#############################################################################################################
#with(seasons_vj.3,plot(Topt_275~Ear,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,1),ylim=c(15,40)))
#magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
#lmm_ear<-lmer(Topt_275~Ear+(1|Species),data=subset(seasons_vj.3,seasons_vj.3$EaV>20)) #random intercept model. Dataset as random factor
#rsquared.glmm(lmm_ear)
#Anova(lmm_ear,test="F")
#--------------------------------------------------------------------------------------------------------------------------------------

#windows(40,40);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(seasons_vj.3,plot(Topt_275~EaV,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(20,100),ylim=c(10,40)))

magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T,at=1:13)

aeb.1<-subset(seasons_vj.3,seasons_vj.3$Topt_275.se<6)
adderrorbars(x=seasons_vj.3$EaV,y=seasons_vj.3$Topt_275,
             SE=aeb.1$Topt_275.se,direction="updown")


adderrorbars(x=seasons_vj.3$EaV,y=seasons_vj.3$Topt_275,
             SE=seasons_vj.3$EaV.se,direction="leftright")


title(ylab=expression(Topt[A275]~(degree*C)),
      xpd=NA,cex.lab=2,line=-27)



title(xlab=expression(Ea[Vcmax]~(kJ~mol^-1)),
      xpd=NA,cex.lab=2)

lmm_tae<-lmer(Topt_275~EaV+(1|Species),data=seasons_vj.3) #random intercept model. Dataset as random factor
lmm_tae.1<-lmer(Topt_275~EaV+(Tavg_30|Species),data=seasons_vj.3) #model with random slope and intercept 

Anova(lmm_tae,test="F")
Anova(lmm_tae.1,test="F")

anova(lmm_tae,lmm_tae.1) #random slope no need
rsquared.glmm(lmm_tae)


#qqPlot(residuals(lmm_tae)) #nice
#plot(resid(lmm_tae) ~ fitted(lmm_tae)) #not too bad
#abline(h=0)


ablineclip(16.8587, 0.1666, x1 = min(seasons_vj.3$EaV,na.rm=T), x2 = max(seasons_vj.3$EaV,na.rm=T),col="black",lwd=3,lty=2)
text(40,40,expression(Topt[A275]==16.8-0.17~Ea[Vcmax])) #95%CI
text(40,38,expression(italic(R^2)==0.25))

title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.2,line=3)


legend("topright",expression(bold(b)),cex=2,bty="n")

legend("bottomright",c("Duke pine","Hinoki cypress","CG-Au","Maritime pine","Brazilian rainforest","Australian alpine",
                   "Sydney Blue Gum","Tasmanian Blue Gum (B)","Forest Red Gum","Parramatta Red Gum") ,col=COL.1,pch=16,ncol=2,bg="white",
       cex=1,title="Dataset",bty="n") 

#--------------------------------------------------------------------------------------------------------------------------------------
seasons_vj.3$Ear<-with(seasons_vj.3,EaJ/EaV)
windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(seasons_vj.3,plot(Topt_275~EaJ,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,70),ylim=c(10,40)))

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)

lmm_the<-lmer(Topt_275~EaJ+(1|Species),data=subset(seasons_vj.3,seasons_vj.3$EaJ<90)) #random intercept model. Dataset as random factor
lmm_the.1<-lmer(Topt_275~EaJ+(Tavg_30|Species),data=subset(seasons_vj.3,seasons_vj.3$EaJ<90)) #model with random slope and intercept 

Anova(lmm_tae,test="F")
Anova(lmm_tae.1,test="F")

anova(lmm_tae,lmm_tae.1) #random slope no need
rsquared.glmm(lmm_the)

ablineclip(16.5417, 0.1725, x1 = min(seasons_vj.3$EaV,na.rm=T), x2 = max(seasons_vj.3$EaV,na.rm=T),col="black",lwd=3,lty=2)
text(75,40,expression(hat(beta)==0.17~(0.09-0.26)~degree*C~kJ^-1~mol^-2)) #95%CI
text(80,38,expression(italic(R^2)==0.26))

#############################################################################################################
#Plot JV ratio (adaptation)

COL.3<-palette()[c(1,2,5,6,7,3,4,8,9,10,12,11)]

windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))


with(adapt_bio.1,plot(MAXSS,JVr,col=COL.3[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,50),ylim=c(0,5)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)



plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_1"),yvar="JVr",xvar="MAXSS",linecol=COL.2[1],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_2"),yvar="JVr",xvar="MAXSS",linecol=COL.2[2],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_3"),yvar="JVr",xvar="MAXSS",linecol=COL.2[3],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="DREYER"),yvar="JVr",xvar="MAXSS",linecol=COL.2[4],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="HFE_CG"),yvar="JVr",xvar="MAXSS",linecol=COL.2[6],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="MEDLYN"),yvar="JVr",xvar="MAXSS",linecol=COL.2[7],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="MIKE_ASPINWALL"),yvar="JVr",xvar="MAXSS",linecol=COL.2[8],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="RF_MON"),yvar="JVr",xvar="MAXSS",linecol=COL.2[9],fitoneline=T)
plot_fit_line(data=subset(adapt_bio.1,adapt_bio.1$DataSet=="WANG_ET_AL"),yvar="JVr",xvar="MAXSS",linecol=COL.2[12],fitoneline=T)



legend("topleft",expression(Common~garden(D^1),Common~garden(D^2),Common~garden(D^3),Deciduous^1~spp,Woodlands~spp,
                                italic(Eucalyptus~spp^2),italic(P.~pinaster^3),italic(Corymbia)*~provenances^4,Montane~rain~forest~spp,
                                Deciduous^2~spp,italic(Betula)*~and~italic(Pinus~spp)),col=COL.3,pch=16,ncol=2,bg="white",cex=.8,title="Species",bty="n")





title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T)
title(ylab=expression(J[max25]/V[cmax25]~(degree*C)),
      xpd=NA,cex.lab=2)



##########################################################################################################################################################

#plot TPU ftraction vs Tgrowth

COL.1<-palette()[c(1:4,6,5,8:11)]

windows(80,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(seasons_vj.3,plot(Tavg_30,TPUFrac,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,1)))

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
#axis.break(2, 20, style = "zigzag")

lmm_the<-lmer(TPUFrac~Tavg_30+(1|Species),data=subset(seasons_vj.3,seasons_vj.3$EaJ<90)) #random intercept model. Dataset as random factor
lmm_the.1<-lmer(TPUFrac~Tavg_30+(Tavg_30|Species),data=subset(seasons_vj.3,seasons_vj.3$EaJ<90)) #model with random slope and intercept 

Anova(lmm_the,test="F")
Anova(lmm_the.1,test="F")

anova(lmm_tae,lmm_tae.1) #random slope no need
rsquared.glmm(lmm_the)


##########################################################################################################################################################




##########################################################################################################################################################
##########################################################################################################################################################

#analyse with warming experiments

##########################################################################################################################################################
##########################################################################################################################################################


warm<-subset(tfits.1,tfits.1$Temp_Treatment %in% c(26,32, "ambient", "Ambient", "Cool", "elevated", "Elevated", "warm")) 
warm.1<-subset(warm,!warm$Season_New %in% c("Autumn","Winter"))

warm.1$Treatment<-NA
warm.1$Treatment[warm.1$Temp_Treatment==26]<-"Ambient"
warm.1$Treatment[warm.1$Temp_Treatment==32]<-"Warmed"
warm.1$Treatment[warm.1$Temp_Treatment=="ambient"]<-"Ambient"
warm.1$Treatment[warm.1$Temp_Treatment=="Ambient"]<-"Ambient"
warm.1$Treatment[warm.1$Temp_Treatment=="Cool"]<-"Ambient"
warm.1$Treatment[warm.1$Temp_Treatment=="elevated"]<-"Warmed"
warm.1$Treatment[warm.1$Temp_Treatment=="Elevated"]<-"Warmed"
warm.1$Treatment[warm.1$Temp_Treatment=="warm"]<-"Warmed"
warm.1$Treatment<-factor(warm.1$Treatment)

warm.1$DataSet<-factor(warm.1$DataSet)

#plot relationships
#Seasonal Topts
palette(rev(brewer.pal(12,"Paired")))
COL<-palette()[6:12]


#plot topt for photosynthesis
windows(40,40);
par(cex.lab=1.5,mar=c(.5,.5,.5,.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))
#windows(40,40);
#par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))


with(warm.1,plot(Tavg_30,Topt,col=COL[factor(DataSet)],pch=c(16,17)[Treatment],cex=1.5,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(10,40)))

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)


plot_fit_line(data=subset(warm.1,warm.1$DataSet=="MIKE_ASPINWALL"),yvar="Topt",xvar="Tavg_30",byvar="Species",linecol=COL[1],fitoneline=F)
plot_fit_line(data=subset(warm.1,warm.1$DataSet=="S30"),yvar="Topt",xvar="Tavg_30",byvar="Species",linecol=COL[2],fitoneline=F)
plot_fit_line(data=subset(warm.1,warm.1$DataSet=="WAY_ET_AL"),yvar="Topt",xvar="Tavg_30",linecol=COL[3],fitoneline=T)
plot_fit_line(data=subset(warm.1,warm.1$DataSet=="WTC2"),yvar="Topt",xvar="Tavg_30",linecol=COL[4],fitoneline=T)
#plot_fit_line(data=subset(warm.1,warm.1$DataSet=="WTC3"),yvar="Topt",xvar="Tavg_30",linecol=COL[5],fitoneline=T)
plot_fit_line(data=subset(warm.1,warm.1$DataSet=="WTC4"),yvar="Topt",xvar="Tavg_30",linecol=COL[6],fitoneline=T)

points(22.63,18,pch=1,cex=1.5,col=COL[5],lwd=3)


adderrorbars(x=warm.1$Tavg_30,y=warm.1$Topt,
             SE=warm.1$Topt.se,direction="updown")

title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T)

title(ylab=expression(Topt[A]~(degree*C)),
      xpd=NA,cex.lab=2)


legend(10,40,c("Mike Aspinwall","S30_2spp","Dani_Way","WTC-V2","WTC-V3","WTC-V4") ,col=COL,pch=16,ncol=3,bg="white",
       cex=.8,title="Dataset",bty="n") 

legend(10,34,c("Ambient","Warmed") ,col="red",pch=c(16,17),ncol=1,bg="white",
       cex=.8,title="Temperature treatment",bty="n") 

##########################################################################################################################################################
##########################################################################################################################################################

#biochemistry

warm.vj<-subset(tfits_VJ,tfits_VJ$Temp_Treatment %in% c(26,32, "ambient", "Ambient", "Cool", "elevated", "Elevated", "warm")) #only ambient treatments
warm.vj.1<-subset(warm.vj,!warm.vj$Season_New %in% c("Autumn","Winter"))
warm.vj.1$DataSet<-factor(warm.vj.1$DataSet)

warm.vj.1$Treatment<-NA
warm.vj.1$Treatment[warm.vj.1$Temp_Treatment==26]<-"Ambient"
warm.vj.1$Treatment[warm.vj.1$Temp_Treatment==32]<-"Warmed"
warm.vj.1$Treatment[warm.vj.1$Temp_Treatment=="ambient"]<-"Ambient"
warm.vj.1$Treatment[warm.vj.1$Temp_Treatment=="Ambient"]<-"Ambient"
warm.vj.1$Treatment[warm.vj.1$Temp_Treatment=="Cool"]<-"Ambient"
warm.vj.1$Treatment[warm.vj.1$Temp_Treatment=="elevated"]<-"Warmed"
warm.vj.1$Treatment[warm.vj.1$Temp_Treatment=="Elevated"]<-"Warmed"
warm.vj.1$Treatment[warm.vj.1$Temp_Treatment=="warm"]<-"Warmed"
warm.vj.1$Treatment<-factor(warm.vj.1$Treatment)

warm.vj.1$DataSet<-factor(warm.vj.1$DataSet)



#Vcmax and Jmax

COL.1<-palette()[c(6,8,9,10,11)]
windows(100,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,8,1,8),cex.axis=1.5,las=1,mfrow=c(1,2))

with(warm.vj.1,plot(Tavg_30,Vcmax25,col=COL.1[factor(DataSet)],pch=c(16,17)[Treatment],cex=2,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(0,250)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
adderrorbars(x=warm.vj.1$Tavg_30,y=warm.vj.1$Vcmax25,SE=warm.vj.1$Vcmax25.se ,direction="updown")

title(ylab=expression(V[cmax]^25~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=2)
title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.2,line=3)

legend(10,250,c("Mike Aspinwall","Dani_Way","WTC-V2","WTC-V3","WTC-V4") ,col=COL.1,pch=16,ncol=1,bg="white",
       cex=.8,title="Dataset",bty="n") 

legend(20,250,c("Ambient","Warmed") ,col="red",pch=c(16,17),ncol=1,bg="white",
       cex=.8,title="Temperature treatment",bty="n") 

#-----------

with(warm.vj.1,plot(Tavg_30,Jmax25,col=COL.1[factor(DataSet)],pch=c(16,17)[Treatment],cex=2,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(0,250)))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T)
adderrorbars(x=warm.vj.1$Tavg_30,y=warm.vj.1$Jmax25,SE=warm.vj.1$Jmax25.se ,direction="updown")

title(ylab=expression(J[max]^25~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=2,line=-27)
title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.9,line=3)


#-----------
#JV ratio
windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(warm.vj.1,plot(Tavg_30,JVr,col=COL.1[factor(DataSet)],pch=c(16,17)[Treatment],cex=2,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(1,3)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
adderrorbars(x=warm.vj.1$Tavg_30,y=warm.vj.1$JVr,SE=warm.vj.1$JVr_SE ,direction="updown")

#add Dani Way data. No measurements at 25. So get fitted Jmax:Vcmax 
with(subset(warm.vj.1,warm.vj.1$DataSet=="WAY_ET_AL"),points(Tavg_30,Jmax25/Vcmax25,col=COL.1[2],cex=2,pch=c(16,17)[Treatment]))


legend(10,3,c("Mike Aspinwall","Dani_Way","WTC-V2","WTC-V3","WTC-V4") ,col=COL.1,pch=16,ncol=3,bg="white",
       cex=1,title="Dataset",bty="n") 

legend(10,2.6,c("Ambient","Warmed") ,col="red",pch=c(16,17),ncol=1,bg="white",
       cex=1,title="Temperature treatment",bty="n") 


title(ylab=expression(JV[r]),
      xpd=NA,cex.lab=2)
title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.6,line=3)

##########################################################################################################################################################
##########################################################################################################################################################
windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))
with(warm.vj.1,plot(Topt_275~Tavg_30,col=COL.1[factor(DataSet)],pch=c(16,17)[Treatment],cex=2,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(10,30)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)


title(ylab=expression(Topt[A275]~(degree*C)),
      xpd=NA,cex.lab=2)

title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T,line=3)

legend(10,15,c("Mike Aspinwall","Dani_Way","WTC-V2","WTC-V3","WTC-V4") ,col=COL.1,pch=16,ncol=3,bg="white",
       cex=1,title="Dataset",bty="n") 

legend(10,12,c("Ambient","Warmed") ,col="red",pch=c(16,17),ncol=1,bg="white",
       cex=.8,title="Temperature treatment",bty="n") 

abline(a=0,b=1,lty=2)

lm_jvr<-lmer(Topt_275~Tavg_30+(1|Species),data=warm.vj.1)
lm_jvr.1<-lmer(Jmax25~Vcmax25+(Vcmax25|Species),data=warm.vj.1)
anova(lm_jvr,lm_jvr.1)

abline(a=1.90769,b=1.64528,lty=1)
rsquared.glmm(lm_jvr)

##########################################################################################################################################################
##########################################################################################################################################################

windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,5),cex.axis=1.5,las=1,mfrow=c(1,1))
with(warm.vj.1,plot(EaV~Tavg_30,col=COL.1[factor(DataSet)],pch=c(16,17)[Treatment],cex=2,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(10,100)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)



windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,5),cex.axis=1.5,las=1,mfrow=c(1,1))
with(warm.vj.1,plot(EaJ~Tavg_30,col=COL.1[factor(DataSet)],pch=c(16,17)[Treatment],cex=2,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(10,100)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)



windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,5),cex.axis=1.5,las=1,mfrow=c(1,1))
with(warm.vj.1,plot(ToptV~Tavg_30,col=COL.1[factor(DataSet)],pch=c(16,17)[Treatment],cex=2,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(30,45)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)


windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,5),cex.axis=1.5,las=1,mfrow=c(1,1))
with(warm.vj.1,plot(ToptJ~Tavg_30,col=COL.1[factor(DataSet)],pch=c(16,17)[Treatment],cex=2,ylab="",xlab="",axes=F,xlim=c(10,30),ylim=c(30,45)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)

##########################################################################################################################################################
##########################################################################################################################################################

#getting back to figure 1

# define the temperatures to model, assume a constant dewpoint (so VPD increases with temperature)
Temps <- seq(15,40,by=0.1)
Tdew=15
VPD <- DewtoVPD(Tdew=Tdew,TdegC=Temps)

#get Vcmax

vcmax<-subset(tfits_VJ,  DataSet %in% c("GWW","RF_AMZ","SAVANNA","TARVAINEN","TMB",
                                      "WANG_ET_AL","RF_AUS") ) 

vcmax.1<-vcmax[c("Species", "DataSet","Season_New","Temp_Treatment","Vcmax25","BIO1","Tavg_30","Topt_275","Topt_275.se","EaV","EaJ","delsV","delsJ","MATSS","Rday")]

vcmax.1<-subset(vcmax.1,!is.na(Vcmax25) & !vcmax.1$Season %in% c("Autumn","Dry","Spring")& !vcmax.1$Species %in% "Betula pendula")
vcmax.1$DataSet<-factor(vcmax.1$DataSet)
vcmax.1$Species<-factor(vcmax.1$Species)


#- model photosynthesis. 

#-. 1 Vcmax, Jmax and their temperature response parameters specified for each dataset

model_photo<-function(data){
  
  Temps <- seq(15,40,by=0.1)
  Tdew=15
  VPD <- DewtoVPD(Tdew=Tdew,TdegC=Temps)
  
 
  #Tgr<-data$BIO1/10
  Tg<-data$Tavg_30
  
  if(is.na(Tg))
    {Tgr<-with(data,BIO1/10)}
 
   else{Tgr<-Tg}
  
  #Vcmax25<-data$Vcmax25
  #Jmax25<-data$Vcmax25
  #EaV<-data$EaV*1000
  #EaJ<-data$EaJ*1000
  delsV<-data$delsV*1000
  delsJ<-data$delsJ*1000
  #Rday<-data$Rday
  Vcmax25<-68.77453 #from fits to all data
  #EaV=66398
  #EaJ=31857
  #Jmax25<-123.9384
  Jmax25<-Vcmax25*(2.35-0.034*Tgr)

  out <- Photosyn(Tleaf=Temps,VPD=VPD,Ci=275,
                       Vcmax=Vcmax25,EaV=66398,EdVC=2e5,delsC=delsV,
                       Jmax = Jmax25,EaJ=31857,EdVJ=2e5,delsJ=delsJ)
  
  return(out)
}

dat_mod<-split(vcmax.1,vcmax.1$Species)
mod_photo<-lapply(dat_mod,FUN=model_photo)
topt_mod<-lapply(mod_photo,FUN=topt_from_photosyn)
topts<-get_topts(topt_mod)
topts$Species<-names(dat_mod)

topts_cpm<-merge(vcmax.1,topts,by="Species")
topts_cpm$Topt_diff<-with(topts_cpm,topts-Topt_275)






palette(rev(brewer.pal(12,"Paired")))
COL<-palette()[3:12]
windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))


with(topts_cpm,plot(Topt_275,topts,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(15,35),ylim=c(15,35)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
abline(a=0,b=1,lty=2)

title(ylab=expression(Topt[275]^Modelled~(degree*C)),
      xpd=NA,cex.lab=1.5)


title(xlab=expression(Topt[275]~(degree*C)),
      xpd=NA,cex.lab=1.5)


lm_fit<-lm(topts~Topt_275,data=topts_cpm)
#summary(lm_fit)
ablineclip(a=coef(lm_fit)[[1]],b=coef(lm_fit)[[2]],x1=16,x2=35,col="black",lwd=2)
legend("bottomright",expression(R^2==0.60),bty="n",cex=1.5)
#abline(a=15.15875,b=0.43019)
legend("top",c("Australian semi-arid woodlands","Brazilian rainforest","Australian rainforest","Australian savanna","Sweden","Australian alpine","Finland"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")



windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))

with(t_home_bio,plot(MATSS,TPUFrac,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,30),ylim=c(0,1)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)

title(ylab=expression(TPU[ratio]),
      xpd=NA,cex.lab=1.5)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T)

legend("topleft",c("Australian semi-arid woodlands","Brazilian rainforest","Australian rainforest","Australian savanna","Sweden","Australian alpine","Finland"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")




with(topts_cpm,plot(BIO1/10,topts,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(15,35)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)













topt_from_photosyn<-function(dat){
  #dat$ALEAF<-round(dat$ALEAF, digits = 3)
  topt<-c()
  for(i in 1:length(dat$ALEAF)){
    
    topt[i]<-ifelse(abs(dat$ALEAF[i]-dat$ALEAF[i+1])<0.001,dat$Tleaf[i],0)
  }
  
  topt<-data.frame(topt)
  return(max(topt$topt,na.rm=TRUE))
}

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

tfits_VJ$Jmax25_mod<-with(tfits_VJ,Vcmax25*(2.35-0.034*Tavg_30))

with(tfits_VJ,plot(Jmax25,Jmax25_mod))
abline(a=0,b=1)
summary(lm(Jmax25~Jmax25_mod,data=tfits_VJ))

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

dat_mod<-split(vcmax.1,vcmax.1$Species)
mod_photo<-lapply(dat_mod,FUN=model_photo)
topt_mod<-lapply(mod_photo,FUN=topt_from_photosyn)
topts<-get_topts(topt_mod)
topts$Species<-names(dat_mod)

topts_cpm<-merge(vcmax.1,topts,by="Species")
topts_cpm$Topt_diff<-with(topts_cpm,topts-Topt_275)

palette(rev(brewer.pal(12,"Paired")))
COL<-palette()[3:12]


with(topts_cpm,plot(Topt_275,topts,pch=16,cex=2,col=COL[DataSet]))
abline(a=0,b=1)
summary(lm(topts~Topt_275,topts_cpm))
abline(a=12.84391,b=0.45864,lty=2)


windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,5),cex.axis=1.5,las=1,mfrow=c(1,1))

mod_photo<-lapply(dat_mod[4],FUN=model_photo)
with(mod_photo[[1]],plot(Tleaf,ALEAF,ylim=c(0,25),ylab="",xlab="",axes=F))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)

mod_photo<-lapply(dat_mod[1],FUN=model_photo)
with(mod_photo[[1]],points(Tleaf,ALEAF,col="red"))

with(mod_photo[3][[1]],points(Tleaf,ALEAF,col="green"))
with(mod_photo[7][[1]],points(Tleaf,ALEAF,col="pink"))
with(mod_photo[8][[1]],points(Tleaf,ALEAF,col="brown"))

##########################################################################################################################################################
##########################################################################################################################################################

palette(rev(brewer.pal(12,"Paired")))
COL<-palette()[3:12]

windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(vcmax.1,plot(BIO1/10,Topt_275,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(0,35)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)

title(ylab=expression(Topt[A275]~(degree*C)),
      xpd=NA,cex.lab=2)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T,line=3)
legend(2,10,c("Australian semi-arid woodlands","Brazilian rainforest","Australian rainforest","Australian savannah","Sweden","Australian alpine","Finland"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")

vcmax.1$Thome<-with(vcmax.1,BIO1/10)
adderrorbars(x=vcmax.1$BIO1/10,y=vcmax.1$Topt_275,SE=vcmax.1$Topt_275.se ,direction="updown")
lm_fit<-lm(Topt_275~Thome,data=vcmax.1)

abline(a=coef(lm_fit)[1],b=coef(lm_fit)[2],lty=3)
legend("bottomright",expression(R^2==0.73),bty="n",cex=2)


#--------
windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(topts_cpm,plot(Topt_275,topts,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(15,35),ylim=c(15,35)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)

title(ylab=expression(Topt[A275]~(Predicted)~(degree*C)),
      xpd=NA,cex.lab=2)

title(xlab=expression(Topt[A275]~(Observed)~(degree*C)),
      xpd=NA,cex.lab=2)

adderrorbars(x=topts_cpm$Topt_275,y=topts_cpm$topts,SE=topts_cpm$Topt_275.se ,direction="leftright")

lm_fit<-lm(topts~Topt_275,data=topts_cpm)
abline(a=coef(lm_fit)[1],b=coef(lm_fit)[2],lty=3)
abline(a=0,b=1,lty=1)


legend(15,35,c("Australian semi-arid woodlands","Brazilian rainforest","Australian rainforest","Australian savannah","Sweden","Australian alpine","Finland"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")
legend("bottomright",expression(R^2==0.75),bty="n",cex=2)

#--------
windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))
with(tfits_VJ,plot(Jmax25,Jmax25_mod,pch=16,col=alpha("black",.4),cex=2,ylab="",xlab="",axes=F,xlim=c(0,260),ylim=c(0,260)))

title(ylab=expression(J[max]~(Predicted)~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=2)

title(xlab=expression(J[max]~(Observed)~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=2)

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
abline(a=0,b=1,lwd=3)
lm_jm<-lm(Jmax25_mod~Jmax25,data=tfits_VJ)
legend("bottomright",expression(R^2==0.76),bty="n",cex=2)
abline(a=coef(lm_jm)[1],b=coef(lm_jm)[2],lwd=3,lty=3,col="red")
legend("topleft",c("Fitted line","1:1 line"),bty="n",cex=2,lty=c(3,1),lwd=3,col=c("red","black"))

############################################################################################################################################
############################################################################################################################################



























Temps <- seq(15,40,by=0.1)
Tdew=15
VPD <- DewtoVPD(Tdew=Tdew,TdegC=Temps)

am<-Photosyn(Jmax=63,Tleaf=Temps,VPD=VPD)
el<-Photosyn(Jmax=101,Tleaf=Temps,VPD=VPD)

with(am,plot(Tleaf,ALEAF,col="red",ylim=c(0,30)))
with(el,points(Tleaf,ALEAF,col="green"))


lm_jvr<-lmer(Jmax25~Tavg_30+(1|Species),data=warm.vj.1)

lm_jvr.1<-lmer(Jmax25~Tavg_30+(Tavg_30|Species),data=warm.vj.1)

anova(lm_jvr,lm_jvr.1) #no need of random slope

Anova(lm_jvr)



seasons_vj.3<-subset(seasons_vj.3,seasons_vj.3$Rday>0)
windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,5,5),cex.axis=1.5,las=1,mfrow=c(1,1))

with(seasons_vj.3,plot(Tavg_30,Rday,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(-2,5)))

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)
#axis.break(2, 20, style = "zigzag")

plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="ELLSWORTH"),yvar="Rday",xvar="Tavg_30",linecol=COL.1[1],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="HAN_ET_AL_A"),yvar="Rday",xvar="Tavg_30",linecol=COL.1[2],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="HFE_CG"),yvar="Rday",xvar="Tavg_30",byvar="Species",linecol=COL.1[3],fitoneline=F)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="MEDLYN"),yvar="Rday",xvar="Tavg_30",byvar="Species",linecol=COL.1[4],fitoneline=F)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="TMB"),yvar="Rday",xvar="Tavg_30",linecol=COL.1[6],fitoneline=T)

plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="RF_AMZ"),yvar="Rday",xvar="Tavg_30",linecol=COL.1[5],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC1"),yvar="Rday",xvar="Tavg_30",linecol=COL.1[7],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC2"),yvar="Rday",xvar="Tavg_30",linecol=COL.1[8],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC3"),yvar="Rday",xvar="Tavg_30",linecol=COL.1[9],fitoneline=T)
plot_fit_line(data=subset(seasons_vj.3,seasons_vj.3$DataSet=="WTC4"),yvar="Rday",xvar="Tavg_30",linecol=COL.1[10],fitoneline=T)


title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=2,outer=T)
title(ylab=expression(R[day]~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=2)



legend("bottomright",expression(italic(P.~taeda),italic(C.~obtusa),italic(Eucalyptus~spp^1),italic(P.~pinaster^2),Amazon~spp,italic(E.~delegatensis),
                            italic(E.~saligna),italic(E.~globulus),italic(E.~tereticornis),italic(E.~parramattensis)),col=COL.1,pch=16,ncol=2,bg="white",cex=.8,title="Species",bty="n")


##########################################################################################################################################################
##########################################################################################################################################################
seasons_vj.3<-seasons_vj.1[c("Species","DataSet","Season_New","Tavg_30","Vcmax25","EaV","delsV")]
seasons_vj.3<-subset(seasons_vj.2,!is.na(Topt_275)&seasons_vj.2$DataSet!="TMB")
seasons_vj.3$DataSet<-factor(seasons_vj.2$DataSet)

#Biochemical parameters

#Vcmax and Jmax at 25C

windows(60,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,8,1,8),cex.axis=1.5,las=1,mfrow=c(2,2))

with(seasons_vj.1,plot(Tavg_30,Vcmax25,col=COL[factor(DataSet)],pch=16,cex=1.2,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,300)))
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T)
adderrorbars(x=tfits_VJ$MATSS,y=tfits_VJ$Vcmax25,SE=tfits_VJ$Vcmax25.se ,direction="updown")

#fit for NET
#lmm.1<-lm(Vcmax25~MATSS,data=subset(tfits_VJ_sel.1,tfits_VJ_sel.1$PFT=="NET"))
#ablineclip(a=coef(lmm.1)[[1]],b=coef(lmm.1)[[2]],x1=min(tfits_VJ$MATSS,na.rm=T),x2=max(tfits_VJ$MATSS,na.rm=T),y1=min(tfits_VJ$Vcmax25,na.rm=T),
#y2=max(tfits_VJ$Vcmax25,na.rm=T),lwd=3,col=COL[3])

lmm.2<-lm(Vcmax25~MATSS,data=subset(tfits_VJ_sel.1,tfits_VJ_sel.1$PFT=="BET_TE"))


ablineclip(a=coef(lmm.2)[[1]],b=coef(lmm.2)[[2]],x1=min(tfits_VJ$MATSS,na.rm=T),x2=max(tfits_VJ$MATSS,na.rm=T),y1=min(tfits_VJ$Vcmax25,na.rm=T),
           y2=max(tfits_VJ$Vcmax25,na.rm=T),lwd=3,col=COL[2])


title(ylab=expression(V[cmax]~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=1.5)

legend("topleft",levels(tfits_VJ$PFT),col=COL,pch=16,ncol=2,bg="white",cex=1.3,title="PFT",bty="n")

#-------------

with(tfits_VJ,plot(Tavg_30,Vcmax25,col=COL[factor(PFT)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,300)))
magaxis(side=c(1,2,4),labels=c(0,0,1),frame.plot=T)
adderrorbars(x=tfits_VJ$Tavg_30,y=tfits_VJ$Vcmax25,SE=tfits_VJ$Vcmax25.se ,direction="updown")

title(ylab=expression(V[cmax]~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=2,srt=90,line=-100)

#fit for NET
tfits_VJ_sel.1<-subset(tfits_VJ, ! DataSet %in% c("HFE_CG","DILLAWAY"))

#fit for NET
lmm.1<-lm(Vcmax25~Tavg_30,data=subset(tfits_VJ_sel.1,tfits_VJ_sel.1$PFT=="NET"))
ablineclip(a=coef(lmm.1)[[1]],b=coef(lmm.1)[[2]],x1=min(tfits_VJ$Tavg_30,na.rm=T),x2=max(tfits_VJ$Tavg_30,na.rm=T),y1=min(tfits_VJ$Vcmax25,na.rm=T),
           y2=max(tfits_VJ$Vcmax25,na.rm=T),lwd=3,col=COL[3])

lmm.2<-lm(Vcmax25~Tavg_30,data=subset(tfits_VJ_sel.1,tfits_VJ_sel.1$PFT=="BET_TE"))
ablineclip(a=coef(lmm.2)[[1]],b=coef(lmm.2)[[2]],x1=min(tfits_VJ$Tavg_30,na.rm=T),x2=max(tfits_VJ$Tavg_30,na.rm=T),y1=min(tfits_VJ$Vcmax25,na.rm=T),
           y2=max(tfits_VJ$Vcmax25,na.rm=T),lwd=3,col=COL[2])


#------------------

##Jmax

with(tfits_VJ,plot(MATSS,Jmax25,col=COL[factor(PFT)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,300)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
adderrorbars(x=tfits_VJ$MATSS,y=tfits_VJ$Jmax25,SE=tfits_VJ$Jmax25.se ,direction="updown")

lmm<-lm(Jmax25~MATSS,data=tfits_VJ_sel)
ablineclip(a=coef(lmm)[[1]],b=coef(lmm)[[2]],x1=min(tfits$MATSS),x2=max(tfits$MATSS),y1=min(tfits$Topt),
           y2=max(tfits$Topt),lwd=3)
title(ylab=expression(J[max]~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=1.5)

#-------------

with(tfits_VJ,plot(Tavg_30,Jmax25,col=COL[factor(PFT)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,300)))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T)
adderrorbars(x=tfits_VJ$Tavg_30,y=tfits_VJ$Jmax25,SE=tfits_VJ$Jmax25.se ,direction="updown")

#fit for NET
tfits_VJ_sel.1<-subset(tfits_VJ, ! DataSet %in% c("HFE_CG","DILLAWAY"))

#fit for NET
lmm.1<-lm(Jmax25~Tavg_30,data=subset(tfits_VJ_sel.1,tfits_VJ_sel.1$PFT=="NET")) #no HFE_CG and Dillaway
ablineclip(a=coef(lmm.1)[[1]],b=coef(lmm.1)[[2]],x1=min(tfits_VJ$Tavg_30,na.rm=T),x2=max(tfits_VJ$Tavg_30,na.rm=T),y1=min(tfits_VJ$Jmax25,na.rm=T),
           y2=max(tfits_VJ$Jmax25,na.rm=T),lwd=3,col=COL[3])

lmm.2<-lm(Jmax25~Tavg_30,data=subset(tfits_VJ_sel.1,tfits_VJ_sel.1$PFT=="BET_TE"))
ablineclip(a=coef(lmm.2)[[1]],b=coef(lmm.2)[[2]],x1=min(tfits_VJ$Tavg_30,na.rm=T),x2=max(tfits_VJ$Tavg_30,na.rm=T),y1=min(tfits_VJ$Jmax25,na.rm=T),
           y2=max(tfits_VJ$Jmax25,na.rm=T),lwd=3,col=COL[2])


title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.2,line=3)

title(xlab=expression(Growth~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.9,line=3)
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------

COL.1<-palette()[c(1:5,8:11)]

windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))

with(seasons_vj.2,plot(Tavg_30,TPUFrac,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,1.5)))

magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)


title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T)

title(ylab=expression(TPU[ratio]),
      xpd=NA,cex.lab=1.5)

legend("top",c("Duke forest NC, USA","Japan (b)","CG-Au","France (a)","Brazilian rainforest",
               "WTC-V1","WTC-V2","WTC-V3","WTC-V4") ,col=COL,pch=16,ncol=3,bg="white",
       cex=.8,title="Dataset",bty="n") 









with(subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_1"),plot(MATSS,JVr,col=Species,pch=16,cex=2))

with(subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_2"),points(MATSS,JVr,col=Species,pch=1,cex=2))

with(subset(adapt_bio.1,adapt_bio.1$DataSet=="DILLAWAY_3"),points(MATSS,JVr,col=Species,pch=2,cex=2))


lm1<-lmer(Topt~Season_New*MATSS+(1|DataSet),data=seasons.1)
Anova(lm1,test="F")

lm2<-lmer(Aopt~Season_New*MATSS+(1|DataSet),data=seasons.1)
Anova(lm2,test="F")

lm3<-lmer(b~Season_New*MATSS+(1|DataSet),data=seasons.1)
Anova(lm3,test="F")

lm4<-lmer(Topt~Season_New*MATSS+(1|DataSet),data=seasons.1)
Anova(lm4,test="F")

lm5<-lmer(Cond.mean~Season_New*MATSS+(1|DataSet),data=subset(seasons.1,seasons.1$Cond.mean<2)) #remove one outlier
Anova(lm5,test="F")


#############################################################################################################

seasons_vj.4<-subset(seasons_vj.3, seasons_vj.3$Species %in% c("Pinus pinaster_L","Pinus pinaster_T"))

lm5<-lmer(JVr~Tavg_30+(1|Species),data=seasons_vj.3) 
Anova(lm5,test="F")
qqPlot(residuals(lm5))
plot(resid(lm5) ~ fitted(lm5))
abline(h=0)
Rsq_marginal<-rsquared.glmm(lm5)[[4]]


lm6<-lmer(JVr~Tavg_30+(Tavg_30|Species),data=seasons_vj.3) #remove one outlier
Anova(lm6,test="F")
qqPlot(residuals(lm6))
plot(resid(lm5) ~ fitted(lm6))
abline(h=0)
Rsq_marginal<-rsquared.glmm(lm6)[[4]]


windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))
with(subset(seasons_vj.3,seasons_vj.3$EaV>20),plot(Tavg_30,EaV,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,100)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)

lm.eav<-lmer(EaV~Tavg_30+(1|Species),data=subset(seasons_vj.3,seasons_vj.3$EaV>20)) 
Anova(lm.eav,test="F")
rsquared.glmm(lm.eav)

windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,1,1),cex.axis=1.5,las=1,mfrow=c(1,1))
with(subset(seasons_vj.3,seasons_vj.3$EaJ<80),plot(Tavg_30,EaJ,col=COL.1[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(0,40),ylim=c(0,100)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)

lm.eaj<-lmer(EaJ/EaV~Tavg_30+(1|Species),data=subset(seasons_vj.3,seasons_vj.3$EaJ<80)) 
Anova(lm.eaj,test="F")
rsquared.glmm(lm.eaj)
