rm(list=ls())
library(deSolve)
library(rootSolve)
library(reshape2)
library(ggplot2)
library(ggstance)
library(cowplot)
library(tidyverse)
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
parms

#Empirical parameters  --- on a t= 1 day scale
a =.018 ## From field experiment, predation rate per ant per day
c =.1 # Use base value, from generalised trophic conversion
m =0.005 # Predator lives 200 days, 1/200
b = 4.16 #From normal clutch size/adult life span (80/30)
g =0.02  ## 50 days til invulnerable size at Topt 
n =0.008 ### From life span of adults, 1/30
o =0.003## From lab rearing -- death rate per day at Topt
K=1000 ##much more realistic K
parms2 <- c(a,c,m,b,g,n,o,K )

###ODE model
predprey_equations_DD_JD  <- function(t, y, parms) {
  a=parms[1]; c=parms[2] ; m= parms[3]; b=parms[4]; g=parms[5]; n=parms[6];o=parms[7]; K=parms[8]
  A0 <- y[1]; P0 <- y[2]; J0 <- y[3]
  dAdt <- g*J0 - n*A0            #adult ode
  dPdt <- (c*a*J0*P0) - m*P0     #predator ode
  dJdt <- b*A0*(1-(J0+A0)/K) - g*J0 - a*J0*P0 -o*J0  #juvenile ode
  return(list(c(dAdt,dPdt,dJdt)));
}


LV.out <- lsoda(Y0, Time,predprey_equations_DD_JD,parms)
par(las=1,bty="l")
matplot(LV.out[,1],LV.out[,-1], type='l', xlab="time", ylab="density") 


##Temperature response function
taylor<-function(temp,Rm=10,Topt=20,Tsd=1)Rm*exp(-0.5*((temp-Topt)/Tsd)^2)
ts<-seq(20,40,0.1)
gs<-taylor(ts,Topt = 25,Rm=.4,Tsd=5)
plot(gs)
###Stability analysis


eigenv=F
eigenv
stability<-function(Pbar=Pbar,
                    Abar=Abar,
                    Jbar=Jbar,
                    parms=parms,
                    eigenv=FALSE
)
{
  a=parms[1]; c=parms[2] ; m= parms[3]; b=parms[4]; g=parms[5]; n=parms[6];o=parms[7]; K=parms[8]
  dFdP=c*a*Jbar-m
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
  if (eigenv==TRUE){print(LE)} 
  if (Re(LE)<0){
    return("stable")
  } 
  else {return("unstable")
  }
  
  
}

stability(Abar=16.33333,Pbar=16.44444,Jbar=10,parms=parms,eigenv=TRUE)

####Type I, no asymmetry

tempresponse2<-function(ts, ##sequence of temperatures
                        ToptJ,  ##topt for juveniles
                        ToptP,  ##topt for predators
                        TsdP=5,  ##thermal niche breadth
                        TsdA=5,  ##thermal niche breadth
                        RmP=.1, ##max predation rate
                        Rmn=.3, ##max (?) adult death rate
                        Rmm=.1, ##max (?) predator death rate
                        Rmo=.3, ##max (?) juvenile rate
                        RmJ=.4, ##max growth rate
                        parms=parms,
                        ODE=predprey_equations ##choose which model version to run
)
{
  gs<-taylor(ts,Topt = ToptJ,Rm=RmJ,Tsd=TsdA)  ##values of G
  as<-taylor(ts,Topt = ToptP,Rm=RmP,Tsd=TsdP)
  ns<-taylor(ts,Topt = ToptJ,Rm=1-Rmn,Tsd=TsdA)
  ms<-taylor(ts,Topt = ToptP,Rm=1-Rmm,Tsd=TsdP)
  os<-taylor(ts,Topt = ToptJ,Rm=1-Rmo,Tsd=TsdA)
  ns<-1-ns
  ms<-1-ms
  os<-1-os
  
  Js<-numeric(length(ts))
  As<-numeric(length(ts))
  Ps<-numeric(length(ts))
  ss<-numeric(length(ts)) 
  
  for(i in 1:length(as)){
    
    parms1 <- c(a=as[i],c=parms[2],m=ms[i],b=parms[4],g=gs[i],n=ns[i],o=os[i],K=parms[8] )
    rp<-as.data.frame(steady(y=Y0,time=c(0,1e5),func=ODE,parms=parms1,method='runsteady')) #steady function from rootsolve package gets steady state
    As[i]<-rp[1,1]
    Ps[i]<-rp[2,1]
    Js[i]<-rp[3,1]
    ss[i]<-stability(Abar=rp[1,1],Pbar=rp[2,1],Jbar=rp[3,1],parms=parms1)
  }
  out1<-data.frame(As=As,Ps=Ps,Js=Js,stability=ss,Temperature=ts)
  return(out1)
}

###version using empirical values
##Separate Tsds and mortality functions

tempresponse3<-function(ts, ##sequence of temperatures
                        ToptJg=23.42,  ##topt for juvenile
                        ToptJo=17.66,  ##topt for juveniles
                        ToptA=23.42,  ##topt for adults
                        ToptP=23.83,  ##topt for predators
                        TsdP=9.14,  ##thermal niche breadth
                        Tsdg=6.7, 
                        Tsdo=7.9,
                        TsdA=7.9,
                        RmP=0.018, ##max predation rate
                        Rmn=0.008, ##min adult death rate
                        Rmm=0.005, ##min predator death rate
                        Rmo=0.003, ##min  juvenile death rate
                        RmJ=0.02, ##max growth rate
                        parms=parms2,
                        ODE=predprey_equations_DD_JD ##choose which model version to run
)
{
  gs<-taylor(ts,Topt = ToptJg,Rm=RmJ,Tsd=Tsdg)  ##values of G
  as<-taylor(ts,Topt = ToptP,Rm=RmP,Tsd=TsdP)
  ns<-taylor(ts,Topt = ToptJo,Rm=1-Rmn,Tsd=Tsdo)
  ms<-taylor(ts,Topt = ToptP,Rm=1-Rmm,Tsd=TsdP)
  os<-taylor(ts,Topt = ToptJo,Rm=1-Rmo,Tsd=Tsdo)
  
  ns<-1-ns
  ms<-1-ms
  os<-1-os
  
  As<-numeric(length(as))
  Ps<-numeric(length(as))
  Js<-numeric(length(as))
  ss<-numeric(length(as)) 
  
  for(i in 1:length(as)){
    
    parms1 <- c(a=as[i],c=parms[2],m=ms[i],b=parms[4],g=gs[i],n=ns[i],o=os[i],K=parms[8] )
    rp<-as.data.frame(steady(y=Y0,time=c(0,1e5),func=ODE,parms=parms1,method='runsteady')) #steady function from rootsolve package gets steady state
    As[i]<-rp[1,1]
    Ps[i]<-rp[2,1]
    Js[i]<-rp[3,1]
    ss[i]<-stability(Abar=rp[1,1],Pbar=rp[2,1],Jbar=rp[3,1],parms=parms1)
  }
  out1<-data.frame(As=As,Ps=Ps,Js=Js,stability=ss,Temperature=ts)
  return(out1)
}
ts<-seq(15,35,0.1)
##Change pred and prey topt here
same<-tempresponse2(ts,ToptJ=25,ToptP=25,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms,ODE=predprey_equations_DD_JD)

same

str(same)

preym1<-tempresponse2(ts,ToptJ=25,ToptP=23,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
preym1
plot(preym1$Ps)
preyp1<-tempresponse2(ts,ToptJ=25,ToptP=27,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
###Base model plot
preym12<-melt(preym1,id.vars=c("Temperature","stability"))
str(preym12)

ggplot(preym12,aes(x=Temperature,y=value,color=variable,lty=stability))+geom_line()+theme_classic()+ylab('Density ')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+xlim(18,32)+geom_vline(xintercept=25,lty=2)+geom_vline(xintercept=23,lty=2)

preym1plot<-ggplot(preym12,aes(x=Temperature,y=value,color=variable,lty=stability))+geom_line()+theme_classic()+ylab('Density ')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+xlim(18,32)+geom_vline(xintercept=25,lty=2)+geom_vline(xintercept=23,lty=2)+
  annotate('text',x=24,y=49,label=expression(paste(T[opt],' ',Predators)),size=2.5)+
  annotate('text',x=25.75,y=49,label=expression(paste(T[opt],' ',Prey)),size=2.5)+scale_linetype(guide=F)+theme(legend.position = 'none')
preym1plot


preyp12<-melt(preyp1,id.vars=c("Temperature",'stability'))
str(preyp12)
max(preyp12$value[preyp12$variable=='Ps'])
preyp1plot<-ggplot(preyp12,aes(x=Temperature,y=value,color=variable,lty='stability'))+geom_line()+theme_classic()+ylab('Density ')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+xlim(18,32)+ylim(-.1,50)+geom_vline(xintercept=25,lty=2)+geom_vline(xintercept=27,lty=2)+
  annotate('text',x=28,y=49,label=expression(paste(T[opt],' ',Predators)),size=2.5)+
  annotate('text',x=25.75,y=49,label=expression(paste(T[opt],' ',Prey)),size=2.5)+scale_linetype(guide=F)+theme(legend.position = 'none')
preyp1plot


same2<-melt(same,id.vars=c("Temperature","stability"))
str(same2)
sameplot<-ggplot(same2,aes(x=Temperature,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density ')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_vline(xintercept=25,lty=2)+
  annotate('text',x=26.5,y=45,label=expression(paste(T[opt],' ',Predators,' & ' ,Prey)),size=2.5)+scale_linetype(guide=F)+ylim(-.1,50)+xlim(18,32)+theme(legend.position = 'top')

sameplot

###Figure 2
plot_grid(sameplot,preym1plot,preyp1plot,nrow=3,labels=c('A','B','C'),rel_heights = c(1.25,1,1))

##Parameterized versions

ts<-seq(13, 23,0.1)
avirginalis<-tempresponse3(ts=ts)
str(avirginalis)

avirginalis2<-melt(avirginalis,id.vars=c("Temperature",'stability'))
str(avirginalis2)
avirginalis2$stability<-as.factor(avirginalis2$stability)

arrow<-data.frame(x =20,y = 80,xend =22,yend =80,variable="As")
str(arrow)

###Figure 3
ggplot(avirginalis2,aes(x=Temperature,y=value,color=variable))+geom_line()+theme_classic()+ylab('Density ')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_vline(xintercept = 20,lty=3)+geom_vline(xintercept = 22,lty=4)+geom_segment(data=arrow,aes(x=arrow$x,y = arrow$y,xend=arrow$xend,yend=arrow$yend),arrow=arrow(),color="Black")+xlim(13,25)+geom_vline(xintercept=17.66,lty=2)+geom_vline(xintercept=23.8,lty=2)+
  annotate('text',x=24.6,y=400,label=expression(paste(T[opt],' ',Predators)),size=2.5)+
  annotate('text',x=18.2,y=400,label=expression(paste(T[opt],' ',Prey)),size=2.5)

ts<-seq(15,35,0.1)

Y0
Y0[2]<-0
nopred_reponse<-tempresponse2(ts,ToptJ=25,ToptP=25,TsdP=5,TsdA=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms,ODE=predprey_equations_DD_JD)

Y0
Y0[2]<-11



###Show differences by asymmetry level
predp.0<-tempresponse2(ts,ToptJ=25,ToptP=25,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)

predp.25<-tempresponse2(ts,ToptJ=25,ToptP=25.25,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp.5<-tempresponse2(ts,ToptJ=25,ToptP=25.5,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp1<-tempresponse2(ts,ToptJ=25,ToptP=26,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp2<-tempresponse2(ts,ToptJ=25,ToptP=27,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp3<-tempresponse2(ts,ToptJ=25,ToptP=28,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp4<-tempresponse2(ts,ToptJ=25,ToptP=29,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
predp5<-tempresponse2(ts,ToptJ=25,ToptP=30,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)



str(predp.0)
diff0<-predp.0[,1:3]-nopred_reponse[,1:3]
matplot(diff0,type="l",x=ts)
diff.25<-predp.25[,1:3]-nopred_reponse[,1:3]
matplot(diff.25,type="l",x=ts)
diff.5<-predp.5[,1:3]-nopred_reponse[,1:3]
matplot(diff.5,type="l",x=ts)
diff1<-predp1[,1:3]-nopred_reponse[,1:3]
matplot(diff1,type="l",x=ts)
diff2<-predp2[,1:3]-nopred_reponse[,1:3]
matplot(diff2,type="l",x=ts)
diff3<-predp3[,1:3]-nopred_reponse[,1:3]
matplot(diff3,type="l",x=ts)
diff4<-predp4[,1:3]-nopred_reponse[,1:3]
matplot(diff4,type="l",x=ts)
diff5<-predp5[,1:3]-nopred_reponse[,1:3]
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

diff$variable<-factor(diff$variable,levels=c('As','Ps','Js'),labels=c("Adults","Predators","Juveniles"))

diff$variable<-relevel(diff$variable,"Juveniles")
diff$variable<-relevel(diff$variable,"Adults")
diff$variable
library(RColorBrewer)

###Figure 4
ggplot(diff,aes(x=Temperature,y=value,color=Asymmetry))+geom_line()+theme_classic()+ylab('Change vs. No Predator Equilibrium')+ scale_color_brewer(palette = "RdBu",direction = -1)+facet_grid(variable~.)+xlim(20,32.5)


##Figure 5
#initial
A0 = 10
P0 = 11
J0 = 100
Y0 <- c(A0,P0,J0)
out1<-tempresponse2(ts,ToptJ=25,ToptP=27,TsdA=5,TsdP=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)

deltaA<-out1$As[out1$Temperature=='25']-out1$As[out1$Temperature=='23']
deltaJ<-out1$Js[out1$Temperature=='25']-out1$Js[out1$Temperature=='23']
deltaP<-out1$Ps[out1$Temperature=='25']-out1$Ps[out1$Temperature=='23']

deltaA
deltaJ
deltaP
Rma=.1 ##max predation rate
c=0.1
Rmm=.1 ##min predator death rate
b=2
Rmg=.4 ##max growth rate
Rmn=.3 ##min adult death rate
Rmo=.3 ##min juvenile death rate
TsdA=5  ##thermal niche breadth
TsdP=5  ##thermal niche breadth
ToptJ=25  ##topt for juveniles
ToptP=27  ##topt for predators
parms3<-c(Rma,c,Rmm,b,Rmg,Rmn,Rmo,K, TsdA,TsdP,ToptJ,ToptP)
parms<-parms3
parms[10]
sensitivity2<-function(ts=ts,
                      parameter=1,
                      from=0,
                      to=1,
                      by=0.01,
                      parms=parms3,
                      ODE=predprey_equations_DD_JD##choose which model version to run        
){
  ODE=predprey_equations_DD_JD
  paramsensitivity<-seq(from=from,to=to,by=by)
  #paramsensitivity<-seq(from=9,to=10,by=.1)
  dA<-vector(length=length(paramsensitivity))
  dJ<-vector(length=length(paramsensitivity))
  dP<-vector(length=length(paramsensitivity))
 # parameter=10
  for(i in 1:length(paramsensitivity)){
    
    parms[parameter]<-paramsensitivity[i]
    Rma=parms[1]
    c=parms[2]
    Rmm= parms[3]; 
    b=parms[4]; 
    Rmg=parms[5]; 
    Rmn=parms[6];
    Rmo=parms[7]; 
    K=parms[8]
    TsdA=parms[9]
    TsdP=parms[10]
    ToptJ=parms[11]
    ToptP=parms[12]
    parms2<-c(Rma,c,Rmm,b,Rmg,Rmn,Rmo,K)
  out1<-tempresponse2(ts=ts, ##sequence of temperatures
                  ToptJ=ToptJ,  ##topt for juveniles
                  ToptP=ToptP,  ##topt for predators
                  TsdP=TsdP,  ##thermal niche breadth
                  TsdA=TsdA,  ##thermal niche breadth
                  RmP=Rma, ##max predation rate
                  Rmn=Rmn, ##max (?) adult death rate
                  Rmm=Rmm, ##max (?) predator death rate
                  Rmo=Rmo, ##max (?) juvenile rate
                  RmJ=Rmg, ##max growth rate
                  parms=parms2,
                  ODE=predprey_equations_DD_JD)
    deltaA<-out1$As[out1$Temperature=='25']-out1$As[out1$Temperature=='23']
    deltaJ<-out1$Js[out1$Temperature=='25']-out1$Js[out1$Temperature=='23']
    deltaP<-out1$Ps[out1$Temperature=='25']-out1$Ps[out1$Temperature=='23']
    dA[i]<-deltaA
    dJ[i]<-deltaJ
    dP[i]<-deltaP
  }
  out2<-data.frame(dA=dA,dP=dP,dJ=dJ,paramvalue=paramsensitivity)
  out3<-melt(out2,id.vars=c("paramvalue"))
  return(out3)
}



ts=seq(23,25,0.1)


ts
sRma<-sensitivity2(ts=ts,
                    parameter=1,
                    from=0.1,
                    to=0.1,
                    by=0.01,
                    parms=parms3,
                    ODE=predprey_equations_DD_JD
                    
)
plot1<-ggplot(sRma,aes(x=paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab(expression(a[max]))+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'top', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 0.1,lty=3)
plot1

sc<-sensitivity2(ts=ts,
                   parameter=2,
                   from=0.01,
                   to=1,
                   by=0.01,
                   parms=parms3,
                   ODE=predprey_equations_DD_JD
                   
)
plot2<-ggplot(sc,aes(x=paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab('c')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'top', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 0.1,lty=3)
plot2


sRmm<-sensitivity2(ts=ts,
                   parameter=3,
                   from=0.001,
                   to=1,
                   by=0.01,
                   parms=parms3,
                   ODE=predprey_equations_DD_JD
                   
)
plot3<-ggplot(sRmm,aes(x=1-paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab(expression(m[min]))+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'none', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 0.1,lty=3)
plot3

sb<-sensitivity2(ts=ts,
                 parameter=4,
                 from=0.8,
                 to=5,
                 by=0.01,
                 parms=parms3,
                 ODE=predprey_equations_DD_JD
                 
)
plot4<-ggplot(sb,aes(x=paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab('b')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'none', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 2,lty=3)
plot4
sb[sb$value<0.01&sb$value>-0.01,]


sg<-sensitivity2(ts=ts,
                 parameter=5,
                 from=0.1,
                 to=5,
                 by=0.01,
                 parms=parms3,
                 ODE=predprey_equations_DD_JD
                 
)
plot5<-ggplot(sg,aes(x=paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab(expression(g[max]))+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'none', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 0.4,lty=3)
plot5

sRmn<-sensitivity2(ts=ts,
                 parameter=6,
                 from=0.001,
                 to=.95,
                 by=0.01,
                 parms=parms3,
                 ODE=predprey_equations_DD_JD
                 
)
plot6<-ggplot(sRmn,aes(x=1-paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab(expression(n[min]))+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'none', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 0.1,lty=3)
plot6
sRmo<-sensitivity2(ts=ts,
                 parameter=7,
                 from=0.001,
                 to=1,
                 by=0.01,
                 parms=parms3,
                 ODE=predprey_equations_DD_JD
                 
)
plot7<-ggplot(sRmo,aes(x=1-paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab(expression(o[min]))+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'none', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 0.1,lty=3)
plot7
sK<-sensitivity2(ts=ts,
                   parameter=8,
                   from=50,
                   to=1500,
                   by=1,
                   parms=parms3,
                   ODE=predprey_equations_DD_JD
                   
)
plot8<-ggplot(sK,aes(x=paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab('K')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'none', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 100,lty=3)
plot8
sK[sK$value<0.01&sK$value>-0.01,]

stsdA<-sensitivity2(ts=ts,
                    parameter=9,
                    from=1,
                    to=15,
                    by=0.1,
                    parms=parms3,
                    ODE=predprey_equations_DD_JD
                    
)
plot9<-ggplot(stsdA,aes(x=paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab(expression(T[sd]~Prey))+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'none', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 5,lty=3)
plot9
stsdP<-sensitivity2(ts=ts,
                    parameter=10,
                    from=1,
                    to=15,
                    by=0.1,
                    parms=parms3,
                    ODE=predprey_equations_DD_JD
                    
)
plot10<-ggplot(stsdP,aes(x=paramvalue,y=value,color=variable))+geom_line()+theme_classic()+ylab('Change 23-25°C')+xlab(expression(T[sd]~Predators))+scale_color_viridis_d(labels=c("Adults","Predators","Juveniles"),name=NULL)+geom_hline(yintercept =0,lty=2)+theme(legend.position = 'none', axis.title.y = element_text(size = 6))+geom_vline(xintercept = 5,lty=3)
plot10

####Figure 5
plot_grid(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plot10,nrow=5,ncol=2,labels=c("A","B",'C','D','E','F','G','H','I',"J"),rel_heights = c(1.25,1,1,1,1),label_size = 10)





