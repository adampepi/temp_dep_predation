library(deSolve)
library(rootSolve)
library(reshape2)
library(ggplot2)
library(ggstance)
library(cowplot)
#time vector
Tmax = 1000 # time horizon  
TimeStep = 1 # integration time step
Time <- seq(0, Tmax, by = TimeStep)
#initial
A0 = 10
P0 = 11
J0 = 100
Y0 <- c(A0,P0,J0)
#parameters
a =.1 
c =.1 #
m =.1 #predators live longer
b = 2 #must be greater than 1
g =.4 
n =.3 
o =.3 
K=100
parms <- c(a,c,m,b,g,n,o,K )


#Empirical parameters  --- on a t= 1 day scale
a =.013333 ## From field experiment, predation rate per ant per day
c =.1 # Use base value, from generalised trophic conversion
m =0.0027 # Predator lives for one year, 1/365
b = 2.666 #From normal clutch size/adult life span (80/30)
g =0.02  ## 50 days til invulnerable size at Topt 
n =0.033 ### From life span of adults, 1/30
o =0.003 ## From lab rearing -- death rate per day at Topt
K=1000 ##much more realistic K
parms2 <- c(a,c,m,b,g,n,o,K )


predprey_equations_DD_JD  <- function(t, y, parms) {
  a=parms[1]; c=parms[2] ; m= parms[3]; b=parms[4]; g=parms[5]; n=parms[6];o=parms[7]; K=parms[8]
  A0 <- y[1]; P0 <- y[2]; J0 <- y[3]
  dAdt <- g*J0 - n*A0            #adult ode
  dPdt <- (c*a*J0*P0) - m*P0     #predator ode
  dJdt <- b*A0*(1-(J0+A0)/K) - g*J0 - a*J0*P0 -o*J0  #juvenile ode
  return(list(c(dAdt,dPdt,dJdt)));
}

LV.out <- lsoda(Y0, Time,predprey_equations_DD,parms2)
par(las=1,bty="l")
matplot(LV.out[,1],LV.out[,-1], type='l', xlab="time", ylab="density") 


##Temperature response function
taylor<-function(temp,Rm=10,Topt=20,Tsd=1)Rm*exp(-0.5*((temp-Topt)/Tsd)^2)
ts<-seq(20,40,0.1)

tempresponse2<-function(ts, ##sequence of temperatures
                       ToptJ,  ##topt for juveniles
                       ToptP,  ##topt for predators
                       Tsd=5,  ##thermal niche breadth
                       RmP=.1, ##max predation rate
                       Rmn=.1, ##max (?) adult death rate
                       Rmm=.1, ##max (?) predator death rate
                       Rmo=.1, ##max (?) juvenile rate
                       RmJ=.4, ##max growth rate
                       parms=parms,
                       ODE=predprey_equations ##choose which model version to run
)
{
  gs<-taylor(ts,Topt = ToptJ,Rm=RmJ,Tsd=Tsd)  ##values of G
  as<-taylor(ts,Topt = ToptP,Rm=RmP,Tsd=Tsd)
  ns<-taylor(ts,Topt = ToptJ,Rm=1-Rmn,Tsd=Tsd)
  ms<-taylor(ts,Topt = ToptP,Rm=1-Rmm,Tsd=Tsd)
  os<-taylor(ts,Topt = ToptJ,Rm=1-Rmn,Tsd=Tsd)
  ns<-1-ns
  ms<-1-ms
  os<-1-os
  
  js<-numeric(length(as))
  As<-numeric(length(as))
  ps<-numeric(length(as))
  
  for(i in 1:length(as)){
    
    parms1 <- c(a=as[i],c=parms[2],m=ms[i],b=parms[4],g=gs[i],n=ns[i],o=os[i],K=parms[8] )
    rp<-as.data.frame(steady(y=Y0,time=c(0,1e5),func=ODE,parms=parms1,method='runsteady')) #steady function from rootsolve package gets steady state
    js[i]<-rp[1,1]
    As[i]<-rp[2,1]
    ps[i]<-rp[3,1]
  }
  return(cbind(js,As,ps))
}

###version using empirical values
##Separate Tsds and mortality functions

tempresponse3<-function(ts, ##sequence of temperatures
                        ToptJg=23.42,  ##topt for juvenile
                        ToptJo=17.66,  ##topt for juveniles
                        ToptA=23.42,  ##topt for juveniles
                        ToptP=23.83,  ##topt for predators
                        TsdP=9.14,  ##thermal niche breadth
                        Tsdg=6.7, 
                        Tsdo=7.9,
                        TsdA=7.9,
                        RmP=0.0133, ##max predation rate
                        Rmn=0.005, ##min adult death rate
                        Rmm=0.0027, ##min predator death rate
                        Rmo=0.003, ##min  juvenile death rate
                        RmJ=0.02, ##max growth rate
                        parms=parms2,
                        ODE=predprey_equations_DD_JD ##choose which model version to run
)
{
  gs<-taylor(ts,Topt = ToptJg,Rm=RmJ,Tsd=Tsdg)  ##values of G
  as<-taylor(ts,Topt = ToptP,Rm=RmP,Tsd=TsdP)
  ns<-taylor(ts,Topt = ToptA,Rm=1-Rmn,Tsd=TsdA)
  ms<-taylor(ts,Topt = ToptP,Rm=1-Rmm,Tsd=TsdP)
  os<-taylor(ts,Topt = ToptJo,Rm=1-Rmo,Tsd=Tsdo)
  ns<-1-ns
  ms<-1-ms
  os<-1-os
  
  js<-numeric(length(as))
  As<-numeric(length(as))
  ps<-numeric(length(as))
  
  for(i in 1:length(as)){
    
    parms1 <- c(a=as[i],c=parms[2],m=ms[i],b=parms[4],g=gs[i],n=ns[i],o=os[i],K=parms[8] )
    rp<-as.data.frame(steady(y=Y0,time=c(0,1e5),func=ODE,parms=parms1,method='runsteady')) #steady function from rootsolve package gets steady state
    As[i]<-rp[1,1]
    ps[i]<-rp[2,1]
    js[i]<-rp[3,1]
  }
  return(cbind(As,ps,js))
}

ts<-seq(15,35,0.1)
##Change pred and prey topt here
same<-tempresponse2(ts,ToptJ=25,ToptP=25,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms,ODE=predprey_equations_DD_JD)
same
preym1<-tempresponse2(ts,ToptJ=25,ToptP=23,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)

preyp1<-tempresponse2(ts,ToptJ=25,ToptP=27,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
###Base model plots

preym11<-as.data.frame(preym1)
preym11$Temperature<-ts
preym12<-melt(preym11,id.vars=c("Temperature"))
str(preym12)
preym1plot<-ggplot(preym12,aes(x=Temperature,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density ')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+xlim(18,32)+ylim(-.1,42)+geom_vline(xintercept=25,lty=2)+geom_vline(xintercept=23,lty=2)+
  annotate('text',x=24,y=41,label=expression(paste(T[opt],' ',Predators)),size=2.5)+
  annotate('text',x=25.75,y=41,label=expression(paste(T[opt],' ',Prey)),size=2.5)
preym1plot

preyp11<-as.data.frame(preyp1)
preyp11$Temperature<-ts
preyp12<-melt(preyp11,id.vars=c("Temperature"))
str(preyp12)
preyp1plot<-ggplot(preyp12,aes(x=Temperature,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density ')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+xlim(18,32)+xlim(18,32)+ylim(-.1,42)+geom_vline(xintercept=25,lty=2)+geom_vline(xintercept=27,lty=2)+
  annotate('text',x=28,y=41,label=expression(paste(T[opt],' ',Predators)),size=2.5)+
  annotate('text',x=25.75,y=41,label=expression(paste(T[opt],' ',Prey)),size=2.5)
preyp1plot

same1<-as.data.frame(same)
same1$Temperature<-ts
same2<-melt(same1,id.vars=c("Temperature"))
str(same2)
sameplot<-ggplot(same2,aes(x=Temperature,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density ')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+xlim(18,32)+xlim(18,32)+ylim(-.1,42)+ylim(-.1,42)+geom_vline(xintercept=25,lty=2)+
  annotate('text',x=26.5,y=41,label=expression(paste(T[opt],' ',Predators,' & ' ,Prey)),size=2.5)
  
sameplot

plot_grid(sameplot,preym1plot,preyp1plot,nrow=3,labels=c('A','B','C'))

##Parameterized versions

ts<-seq(15,30,0.1)
avirginalis<-tempresponse3(ts=ts)
avirginalis1<-as.data.frame(avirginalis)
avirginalis1$Temperature<-ts

avirginalis2<-melt(avirginalis1,id.vars=c("Temperature"))
str(avirginalis2)
arrow<-data.frame(x =20,y = 200,xend =22,yend =200,variable="As")
str(arrow)
ggplot(avirginalis2,aes(x=Temperature,y=value+1,color=variable))+geom_line()+scale_y_log10()+theme_classic()+ylab('Density + 1')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_vline(xintercept = 20,lty=2)+geom_vline(xintercept = 22,lty=3)+geom_segment(data=arrow,aes(x=arrow$x,y = arrow$y,xend=arrow$xend,yend=arrow$yend),arrow=arrow(),color="Black")


Y0
Y0[2]<-0
nopred_reponse<-tempresponse2(ts,ToptJ=25,ToptP=25,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms,ODE=predprey_equations_DD_JD)

Y0
Y0[2]<-11

diff1<-preym1-nopred_reponse
matplot(diff1,type="l",x=ts)
diff2<-preyp1-nopred_reponse
matplot(diff2,type="l",x=ts)

###Show differences by asymmetry level
predp.0<-tempresponse2(ts,ToptJ=25,ToptP=25,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp.25<-tempresponse2(ts,ToptJ=25,ToptP=25.25,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp.5<-tempresponse2(ts,ToptJ=25,ToptP=25.5,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp1<-tempresponse2(ts,ToptJ=25,ToptP=26,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp2<-tempresponse2(ts,ToptJ=25,ToptP=27,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp3<-tempresponse2(ts,ToptJ=25,ToptP=28,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp4<-tempresponse2(ts,ToptJ=25,ToptP=29,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp5<-tempresponse2(ts,ToptJ=25,ToptP=30,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)

diff0<-predp.0-predp.0
matplot(diff0,type="l",x=ts)
diff.25<-predp.25-predp.0
matplot(diff.25,type="l",x=ts)
diff.5<-predp.5-predp.0
matplot(diff.5,type="l",x=ts)
diff1<-predp1-predp.0
matplot(diff1,type="l",x=ts)
diff2<-predp2-predp.0
matplot(diff2,type="l",x=ts)
diff3<-predp3-predp.0
matplot(diff3,type="l",x=ts)
diff4<-predp4-predp.0
matplot(diff4,type="l",x=ts)
diff5<-predp5-predp.0
matplot(diff5,type="l",x=ts)



diff0<-as.data.frame(diff0)
diff0$Temperature<-ts
diff0$Asymmetry<-"0°C"
str(diff0)
diff0m<-melt(diff0,id.vars=c("Temperature","Asymmetry"))
str(diff0m)

diff.25<-as.data.frame(diff.25)
diff.25$Temperature<-ts
diff.25$Asymmetry<-"0.25°C"
str(diff.25)
diff.25m<-melt(diff.25,id.vars=c("Temperature","Asymmetry"))
str(diff.25m)

diff.5<-as.data.frame(diff.5)
diff.5$Temperature<-ts
diff.5$Asymmetry<-".5°C"
str(diff.5)
diff.5m<-melt(diff.5,id.vars=c("Temperature","Asymmetry"))
str(diff.5m)

diff1<-as.data.frame(diff1)
diff1$Temperature<-ts
diff1$Asymmetry<-"1°C"
str(diff1)
diff1m<-melt(diff1,id.vars=c("Temperature","Asymmetry"))
str(diff1m)

diff2<-as.data.frame(diff2)
diff2$Temperature<-ts
diff2$Asymmetry<-"2°C"
str(diff2)
diff2m<-melt(diff2,id.vars=c("Temperature","Asymmetry"))
str(diff2m)

diff3<-as.data.frame(diff3)
diff3$Temperature<-ts
diff3$Asymmetry<-"3°C"
str(diff3)
diff3m<-melt(diff3,id.vars=c("Temperature","Asymmetry"))
str(diff3m)

diff4<-as.data.frame(diff4)
diff4$Temperature<-ts
diff4$Asymmetry<-"4°C"
str(diff4)
diff4m<-melt(diff4,id.vars=c("Temperature","Asymmetry"))
str(diff4m)

diff5<-as.data.frame(diff5)
diff5$Temperature<-ts
diff5$Asymmetry<-"5°C"
str(diff5)
diff5m<-melt(diff5,id.vars=c("Temperature","Asymmetry"))
str(diff5m)

diff<-rbind(diff0m, diff1m,diff2m,diff3m,diff4m,diff5m)
str(diff)
diff$Asymmetry<-as.factor(diff$Asymmetry)

diff$variable<-factor(diff$variable,levels=c('js','As','ps'),labels=c("Adults","Predators","Juveniles"))

diff$variable<-relevel(diff$variable,"Juveniles")
diff$variable<-relevel(diff$variable,"Adults")
diff$variable
ggplot(diff,aes(x=Temperature,y=value,color=Asymmetry))+geom_line()+theme_classic()+ylab('Change vs. No Predator Equilibrium')+scale_color_viridis_d()+facet_grid(variable~.)+xlim(20,30)


###Other sensitivity analyses
sensitivity<-function(parameter=1,
                      from=0,
                      to=1,
                      by=0.01,
                       parms1=parms,
                       ODE=predprey_equations ##choose which model version to run
)
{

  paramsensitivity<-seq(from=from,to=to,by=by)
  js<-numeric(length(paramsensitivity))
  As<-numeric(length(paramsensitivity))
  ps<-numeric(length(paramsensitivity))
  for(i in 1:length(paramsensitivity)){
    

    parms1[parameter]<-paramsensitivity[i]
    rp<-as.data.frame(steady(y=Y0,time=c(0,1e5),func=ODE,parms=parms1,method='runsteady')) #steady function from rootsolve package gets steady state
    js[i]<-rp[1,1]
    As[i]<-rp[2,1]
    ps[i]<-rp[3,1]
  }
  mm<-cbind(js,As,ps,paramsensitivity)
  mmm<-as.data.frame(mm)
  mmm<-melt(mmm,id.vars="paramsensitivity")
  return(mmm)
}

sK<-sensitivity(parameter=8,from=5,to=1000,by=1,ODE=predprey_equations_DD_JD,parms1=parms)
parms

sK
ggplot(sK,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('K')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))
par(mfrow=c(1,1))

s1<-sensitivity(parameter=1,from=0,to=10,by=0.01,ODE=predprey_equations_DD_JD,parms1=parms)

aplot<-ggplot(s1,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('a (Attack rate)')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"),guide=F)

aplot
s2<-sensitivity(parameter=2,from=0,to=10,by=0.01,ODE=predprey_equations_DD_JD, parms1=parms)

cplot<-ggplot(s2,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('c (Conversion rate)')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"),guide=F)

s3<-sensitivity(parameter=3,from=0,to=10,by=0.01,ODE=predprey_equations_DD_JD,parms1=parms)

mplot<-ggplot(s3,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('m (Predator mortality rate')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"),guide=F)

s4<-sensitivity(parameter=4,from=0,to=10,by=0.01,ODE=predprey_equations_DD_JD,parms1=parms)

bplot<-ggplot(s4,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('b (Prey birth rate)')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))

##g
s5<-sensitivity(parameter=5,from=0,to=10,by=0.01,ODE=predprey_equations_DD_JD,parms1=parms)

gplot<-ggplot(s5,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('g (Juvenile growth rate)')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))

#n
s6<-sensitivity(parameter=6,from=0,to=10,by=0.01,ODE=predprey_equations_DD_JD,parms1=parms)

nplot<-ggplot(s6,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('n (Adult death rate)')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))

#K
s7<-sensitivity(parameter=7,from=0,to=10,by=0.01,ODE=predprey_equations_DD_JD,parms1=parms)

oplot<-ggplot(s7,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('o (Juvenile death rate)')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"),guide=F)

#K
s8<-sensitivity(parameter=8,from=1,to=500,by=1,ODE=predprey_equations_DD_JD,parms1=parms)
Kplot<-ggplot(s8,aes(x=paramsensitivity,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density')+xlab('K (Carrying capacity)')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))

plot_grid(aplot,gplot,cplot,bplot,mplot,nplot,oplot,Kplot,nrow=4,ncol=2,labels=c("A","B","C",'D',"E","F","G","H"),rel_widths=c(1,1.3))

###Stability analysis



stability<-function(Pbar=Pbar,
                    Abar=Abar,
                    Jbar=Jbar,
                    parms=parms
  
){
  a=parms[1]; c=parms[2] ; m= parms[3]; b=parms[4]; g=parms[5]; n=parms[6];o=parms[7]; K=parms[8]
  dFdP=c*a-m
  dFdA=0
  dFdJ=c*a*Pbar
  dGdP=0
  dGdA=-n
  dGdJ=g
  dHdP=-a*Jbar
  dHdA=b*(1-Jbar/K-2*Abar/K)
  dHdJ=-b*Abar/K-g-a*Pbar-o
  
  m1<-matrix(data=c(dFdP,dFdA,dFdJ,dGdP,dGdA,dGdJ,dHdP,dHdA,dHdJ),nrow=3,ncol=3)
  LE<-eigen(m1)$values[1]
  if (LE<0){
    print("stable")
  } 
  else {print("unstable")
  }
}




