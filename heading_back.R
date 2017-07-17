#Figure 1
#plot species measured at their home climates

t_home_bio<-subset(tfits_VJ,  DataSet %in% c("ARCTIC","GWW","RF_AMZ","SAVANNA","TARVAINEN","TMB",
                                      "WANG_ET_AL","RF_AUS","HAN_2")) #these datasets are from

t_home_bio<-subset(t_home_bio,!t_home_bio$Species %in% c("Eucalyptus globulus","Eucalyptus microcorys","Betula pendula","Toona sinensis") 
                   & !t_home_bio$Season %in% c("WINTER","dry","May","February","summer"))

t_home_bio$DataSet<-factor(t_home_bio$DataSet)


#-----------------------------------------------------------------

library(RColorBrewer)
palette(rev(brewer.pal(12,"Paired")))
COL<-palette()[3:12]


#plot topt for photosynthesis
windows(100,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,5),cex.axis=1.5,las=1,mfrow=c(1,2))

with(t_home_bio,plot(MATSS,Vcmax25,col=COL[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(10,250)))
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$Vcmax25,SE=t_home_bio$Vcmax25.se ,direction="updown")
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
title(ylab=expression(V[cmax]~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=1.5)


legend("topleft",c("Arctic spp","Australian semi-arid woodlands","Japanese Red Pine","Brazilian rainforest","Australian rainforest","Australian savanna","Spruce (Sweden)","Alpine Ash","Scots pine (Finland)"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")


with(t_home_bio,plot(MATSS,Jmax25,col=COL[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(10,250)))
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$Jmax25,SE=t_home_bio$Jmax25.se ,direction="updown")
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T)
title(ylab=expression(J[max]~(mu*mol~m^-2~s^-1)),
      xpd=NA,cex.lab=1.5,line=-27)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.2,line=3)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.9,line=3)



with(t_home_bio,plot(MATSS,EaV,col=COL[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(10,150)))
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$EaV,SE=t_home_bio$EaV.se ,direction="updown")
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
title(ylab=expression(E[a]~(kJ~mol^-1)),
      xpd=NA,cex.lab=1.5)


with(t_home_bio,plot(MATSS,EaJ,col=COL[factor(DataSet)],pch=16,cex=1.5,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(10,100)))
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$EaJ,SE=t_home_bio$EaJ.se ,direction="updown")
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T)
title(ylab=expression(H[a]~(kJ~mol^-1)),
      xpd=NA,cex.lab=1.5,line=-19)
title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.2,line=3)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.9,line=3)


#--------------------------------------------
summary(lmfit)$r.squared

#plot topt for vcmax and jmax, dels for vcmax and jmax
t_home_bio$delsv.new<-t_home_bio$delsV*1000
t_home_bio$delsj.new<-t_home_bio$delsJ*1000

windows(100,100);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,5),cex.axis=1.5,las=1,mfrow=c(2,2))

#- dels Vcmax
with(t_home_bio,plot(MATSS,delsV*1000,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(610,680)))
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$delsV*1000,SE=t_home_bio$delsV.se*1000 ,direction="updown")
magaxis(side=c(1,2,4),labels=c(0,1,0),frame.plot=T)
title(ylab=expression(Delta*S[Vcmax]~(J~mol^-1~K^-1)),
      xpd=NA,cex.lab=1.5)

#add Tumbarumba  delsv
#points(8.5,tmb_sum[[3]]*1000,pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(620,680),col=COL[6])
#adderrorbars(x=8.5,y=tmb_sum[[3]]*1000,SE=tmb_sum[[6]]*1000 ,direction="updown")

plot_fit_line(t_home_bio,yvar="delsv.new",xvar="MATSS",linecol=COL[1],fitoneline=T)


#- dels Jmax
with(t_home_bio,plot(MATSS,delsJ*1000,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(620,680)))
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$delsJ*1000,SE=t_home_bio$delsJ.se*1000 ,direction="updown")
magaxis(side=c(1,2,4),labels=c(0, 0,1),frame.plot=T)
title(ylab=expression(Delta*S[Jmax]~(J~mol^-1~K^-1)),
      xpd=NA,cex.lab=1.5,line=-19)

plot_fit_line(t_home_bio,yvar="delsj.new",xvar="MATSS",linecol=COL[1],fitoneline=T)

#add Tumbarumba  delsJ
#points(8.5,tmb_sum.j[[3]]*1000,pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(620,680),col=COL[6])
#adderrorbars(x=8.5,y=tmb_sum.j[[3]]*1000,SE=tmb_sum.j[[6]]*1000 ,direction="updown")


#- Topt for Vcmax
with(t_home_bio,plot(MATSS,ToptV,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(20,45)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
title(ylab=expression(Topt[Vcmax]~(degree*C)),
      xpd=NA,cex.lab=1.5)

plot_fit_line(t_home_bio,yvar="ToptV",xvar="MATSS",linecol=COL[1],fitoneline=T)

#add Tumbarumba  Toptv
#points(8.5,tmb_sum[[7]],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(20,45),col=COL[6])


#- Topt for Jmax
with(t_home_bio,plot(MATSS,ToptJ,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(20,45)))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T)
title(ylab=expression(Topt[Jmax]~(degree*C)),
      xpd=NA,cex.lab=1.5,line=-19)
#points(8.5,tmb_sum.j[[7]],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(20,45),col=COL[6])
plot_fit_line(t_home_bio,yvar="ToptJ",xvar="MATSS",linecol=COL[1],fitoneline=T)



title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.2,line=3)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T,adj=0.9,line=3)



#------------------------------------------------------------------------------------------------------------------------


t_home_bio$JVr_Es<-with(t_home_bio,Jmax25/Vcmax25)
t_home_bio$JVr<-ifelse(is.na(t_home_bio$JVr),t_home_bio$JVr_Es,t_home_bio$JVr)

#-plot  

windows(100,60);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,5),cex.axis=1.5,las=1,mfrow=c(1,2))


with(t_home_bio,plot(BIO1/10,Topt_275,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(0,40)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T,at=1:13)

title(ylab=expression(Topt[A275]~(degree*C)),
      xpd=NA,cex.lab=2,adj=0.5)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T,line=3,adj=0.25)

legend("bottomleft",c("Arctic spp","Australian semi-arid woodlands","Japanese Red Pine","Brazilian rainforest","Australian rainforest","Australian savanna","Spruce (Sweden)","Alpine Ash","Scots pine (Finland)"),
       col=COL,pch=16,ncol=2,bg="white",cex=.8,title="Dataset",bty="n")

t_home_bio$Thome<-with(t_home_bio,BIO1/10)
adderrorbars(x=t_home_bio$BIO1/10,y=t_home_bio$Topt_275,SE=t_home_bio$Topt_275.se ,direction="updown")

lm_fit<-lm(Topt_275~Thome,data=t_home_bio)
abline(a=coef(lm_fit)[1],b=coef(lm_fit)[2],lty=3,lwd=2)
legend("topright",expression(R^2==0.66),bty="n",cex=1.5)

----------------------------
with(t_home_bio,plot(MATSS,JVr,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(0,3)))

#with(t_home_bio,plot(MATSS,JVr_Es,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(0,3)))
magaxis(side=c(1,2,4),labels=c(1,0,1),frame.plot=T)
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$JVr,SE=t_home_bio$JVr_SE ,direction="updown")
title(ylab=expression(JV[ratio~at~25~degree*C]),
      xpd=NA,cex.lab=2,line=-27)

#with(subset(t_home_bio,t_home_bio$Species %in% c("Eucalyptus salmonophloia")),points(MATSS,Jmax25/Vcmax25
                   # ,col=COL[c(1,4)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(0,3)))

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T,adj=0.9,line=3)

lm_fit_JVr<-lm(JVr~Thome,data=subset(t_home_bio,t_home_bio$DataSet !="WANG_ET_AL"))
summary(lm_fit_JVr)

abline(a=coef(lm_fit_JVr)[1],b=coef(lm_fit_JVr)[2],lty=3)
abline(a=2.32,b=-0.034,lty=1,lwd=2) # acclimation curve


legend("topright",expression(R^2~(Long~term)==0.76),bty="n",cex=1.5)
abline(a=coef(lm_fit)[1],b=coef(lm_fit)[2],lty=3,lwd=2)
legend("bottomleft",c("Long term T_response","Short term T_response"),lty=c(3,1),lwd=3,bty="n")

#Add Anna's data (just to see)

points(3.6,2.093414,col=COL[12],pch=1,cex=2,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(0,3),lwd=3)
adderrorbars(x=3.6,y=2.093414,SE=0.06441516 ,direction="updown")


#EucFace
points(17.1,1.56,col=COL[11],pch=1,cex=2,ylab="",xlab="",axes=F,xlim=c(-15,35),ylim=c(0,3),lwd=3)
adderrorbars(x=17.1,y=1.56,SE=0.05 ,direction="updown")

#---------------

#plot photosynthesis at common Ci
windows(60,60);

par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))
with(t_home_bio,plot(MATSS,Photo_800,col=COL[factor(DataSet)],pch=17,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(0,35)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$Photo_800,SE=t_home_bio$Photo_800_SE ,direction="updown")


title(ylab=expression(A[common~Ci]),
      xpd=NA,cex.lab=2)


title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T)


with(t_home_bio,points(MATSS,Photo_275,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(0,35)))
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$Photo_275,SE=t_home_bio$Photo_275_SE ,direction="updown")


legend("topright",c("Australian semi-arid woodlands","Brazilian rainforest","Australian rainforest","Australian savanna","Sweden","Australian alpine","Finland"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")

legend("bottom",c(expression(C[i]==275~ppm),expression(C[i]==800~ppm)),pch=c(16,17),bty="n",cex=1,bg="white")

#------------------------------------------------------------


windows(60,60);

par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))
with(t_home_bio,plot(MATSS,A_ratio,col=COL[factor(DataSet)],pch=17,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(0,3)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$A_ratio,SE=t_home_bio$A_ratio_SE ,direction="updown")


legend("topright",c("Australian semi-arid woodlands","Brazilian rainforest","Australian rainforest","Australian savanna","Sweden","Australian alpine","Finland"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")
title(ylab=expression(A[800]/A[275]),
      xpd=NA,cex.lab=2)


title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T)

#-----------------------------------------------------------

#fitted and measured Vcmax25 values
windows(60,60);

par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))
with(t_home_bio,plot(MATSS,Rday,,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,35),ylim=c(0,5)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)
adderrorbars(x=t_home_bio$MATSS,y=t_home_bio$Rday,SE=t_home_bio$Rday_SE ,direction="updown")

title(ylab=expression(R[day]~25~degree*C),
      xpd=NA,cex.lab=2)
title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=2,outer=T)


legend(6,5,c("Australian semi-arid woodlands","Brazilian rainforest","Australian rainforest","Australian savanna","Sweden","Australian alpine","Finland"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")


#----------------------------------------------------------

#- plot TPUfraction

windows(40,40);
par(cex.lab=1.5,mar=c(0.5,0.5,0.5,0.5),oma=c(5,5,0,0),cex.axis=1.5,las=1,mfrow=c(1,1))

with(t_home_bio,plot(MATSS,TPUFrac,col=COL[factor(DataSet)],pch=16,cex=2,ylab="",xlab="",axes=F,xlim=c(0,30),ylim=c(0,1)))
magaxis(side=c(1,2,4),labels=c(1,1,0),frame.plot=T)

title(ylab=expression(TPU[ratio]),
      xpd=NA,cex.lab=1.5)

title(xlab=expression(Home~T[air]~(degree*C)),cex.lab=1.5,outer=T)

legend("topleft",c("Australian semi-arid woodlands","Brazilian rainforest","Australian rainforest","Australian savanna","Sweden","Australian alpine","Finland"),
       col=COL,pch=16,ncol=2,bg="white",cex=1,title="Dataset",bty="n")


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------


Temps <- seq(15,40,by=0.1)
Tdew=15
VPD <- DewtoVPD(Tdew=Tdew,TdegC=Temps)

Vcmax25<-84.14628
Jmax25<-150.9790
EaV<-59.22934*1000
EaJ<-34.75250*1000
delsV<-0.6389840*1000
delsJ<-0.6402895*1000
Rday<-0.9669966

out <- Photosyn(Tleaf=Temps,Ci=275,
                Vcmax=Vcmax25 ,EaV=,EdVC=2e5,delsC=delsV,
                Jmax = Jmax25,EaJ=EaJ,EdVJ=2e5,delsJ=delsJ)

with(out,plot(Tleaf,ALEAF,xlim=c(20,40)))

topt_from_photosyn(out)
abline(v=27.2)
An<-nls(ALEAF~Aopt-(b*(Tleaf-Topt)^2),data=out,start=list(Aopt=max(out$ALEAF),Topt=25,b=0.05))
summary(An)
