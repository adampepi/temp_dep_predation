library(deSolve)
library(rootSolve)
library(reshape2)
library(ggplot2)

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





predprey_equations  <- function(t, y, parms) {
  a=parms[1]; c=parms[2] ; m= parms[3]; b=parms[4]; g=parms[5]; n=parms[6]
  A0 <- y[1]; P0 <- y[2]; J0 <- y[3]
  dAdt <- g*J0 - n*A0            #adult ode
  dPdt <- (c*a*J0*P0) - m*P0     #predator ode
  dJdt <- b*A0 - g*J0 - a*J0*P0  #juvenile ode
  return(list(c(dAdt,dPdt,dJdt)));
}
LV.out <- lsoda(Y0, Time,predprey_equations,parms)
par(las=1,bty="l")
matplot(LV.out[,1],LV.out[,-1], type='l', xlab="time", ylab="density") 

predprey_equations_DD  <- function(t, y, parms) {
  a=parms[1]; c=parms[2] ; m= parms[3]; b=parms[4]; g=parms[5]; n=parms[6];o=parms[7]; K=parms[8]
  A0 <- y[1]; P0 <- y[2]; J0 <- y[3]
  dAdt <- g*J0 - n*A0            #adult ode
  dPdt <- (c*a*J0*P0) - m*P0     #predator ode
  dJdt <- b*A0*(1-(J0+A0)/K) - g*J0 - a*J0*P0  #juvenile ode
  return(list(c(dAdt,dPdt,dJdt)));
}
predprey_equations_DD3  <- function(t, y, parms) {
  a=parms[1]; c=parms[2] ; m= parms[3]; b=parms[4]; g=parms[5]; n=parms[6];o=parms[7]; K=parms[8]
  A0 <- y[1]; P0 <- y[2]; J0 <- y[3]
  dAdt <- g*J0*(1-(J0+A0)/K) - n*A0*((J0+A0)/K)            #adult ode
  dPdt <- (c*a*J0*P0) - m*P0     #predator ode
  dJdt <- b*A0*(1-(J0+A0)/K) - g*J0*(1-(J0+A0)/K) - a*J0*P0 - o*J0*((J0+A0)/K) #juvenile ode
  return(list(c(dAdt,dPdt,dJdt)));
}
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
gs<-taylor(ts,Topt = 31,Rm=0.1,Tsd=5)
as<-taylor(ts,Topt = 31,Rm=0.1,Tsd=5)


##Temperature sensitivity analysis
tempresponse<-function(ts, ##sequence of temperatures
                       ToptJ,  ##topt for juveniles
                       ToptP,  ##topt for predators
                       Tsd=6.7,  ##thermal niche breadth
                       RmP=.0133, ##max predation rate
                       RmJ=.02, ##max growth rate
                       parms=parms,
                       ODE=predprey_equations ##choose which model version to run
                       )
  {
gs<-taylor(ts,Topt = ToptJ,Rm=RmJ,Tsd=Tsd)  ##values of G
as<-taylor(ts,Topt = ToptP,Rm=RmP,Tsd=Tsd)

js<-numeric(length(as))
As<-numeric(length(as))
ps<-numeric(length(as))

for(i in 1:length(as)){
  
  parms1 <- c(a=as[i],c=parms[2],m=parms[3],b=parms[4],g=gs[i],n=parms[6],o=parms[7],K=parms[8] )
  rp<-as.data.frame(steady(y=Y0,time=c(0,1e5),func=ODE,parms=parms1,method='runsteady')) #steady function from rootsolve package gets steady state
  js[i]<-rp[1,1]
  As[i]<-rp[2,1]
  ps[i]<-rp[3,1]
}
return(cbind(js,As,ps))
}


os<-taylor(ts,Topt = 30,Rm=1-0.3,Tsd=5)
os<-1-os
plot(os)

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


ts<-seq(15,40,0.1)
rt1<-tempresponse(ts,ToptJ=25,ToptP=30,Tsd=5,RmP=.1,RmJ=.4, parms=parms)
rt1
rt2<-tempresponse(ts,ToptJ=30,ToptP=25,Tsd=5,RmP=.1,RmJ=.4, parms=parms)
rt2

par(mfrow=c(2,1))
matplot(x=ts, y=rt1, type='l', xlab="Temp", ylab="density")
matplot(x=ts, y=rt2, type='l', xlab="Temp", ylab="density")

rt1DD<-tempresponse(ts,ToptJ=25,ToptP=30,Tsd=5,RmP=.1,RmJ=.4, parms=parms,ODE=predprey_equations_DD)
rt1DD
rt2DD<-tempresponse(ts,ToptJ=30,ToptP=25,Tsd=5,RmP=.1,RmJ=.4, parms=parms, ODE=predprey_equations_DD)
rt2DD

rt1DD3<-tempresponse(ts,ToptJ=25,ToptP=30,Tsd=5,RmP=.1,RmJ=.4, parms=parms,ODE=predprey_equations_DD3)
rt1DD3
rt2DD3<-tempresponse(ts,ToptJ=30,ToptP=25,Tsd=5,RmP=.1,RmJ=.4, parms=parms, ODE=predprey_equations_DD3)
rt2DD3


##Change pred and prey topt here
rt1DD<-tempresponse2(ts,ToptJ=30,ToptP=30,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms,ODE=predprey_equations_DD_JD)
rt1DD
rt2DD<-tempresponse2(ts,ToptJ=30,ToptP=30,Tsd=5,RmP=.1,RmJ=.4,Rmo=.3,Rmn=.3,Rmm=.1, parms=parms, ODE=predprey_equations_DD_JD)
rt2DD

###Version with t-dep death rate

par(mfrow=c(2,1))
matplot(x=ts, y=rt1DD, type='l', xlab="Temp", ylab="density",main="Prey Topt=30, Pred Topt=30")
matplot(x=ts, y=rt2DD, type='l', xlab="Temp", ylab="density",main="Prey Topt=30, Pred Topt=30")


##Parameterized versions

ts<-seq(15,30,0.1)
avirginalis<-tempresponse3(ts=ts)
avirginalis1<-as.data.frame(avirginalis)
avirginalis1$Temperature<-ts

avirginalis2<-melt(avirginalis1,id.vars=c("Temperature"))
str(avirginalis2)

ggplot(avirginalis2,aes(x=Temperature,y=value+1,color=variable))+geom_line()+scale_y_log10()+theme_classic()+ylab('Density + 1')+scale_color_viridis_d(name=NULL,labels=c("Adults","Predators","Juveniles"))

###Version with t-dep death rate

par(mfrow=c(1,1))
matplot(x=ts, y=rt1DD+1, type='l', xlab="Temp", ylab="density",log='y')
matplot(x=ts, y=rt1DD+1, type='l', xlab="Temp", ylab="density")

#I reorganized this

matplot(x=ts, y=rt2DD, type='l', xlab="Temp", ylab="density",main="Prey Topt=30, Pred Topt=30")



par(mfrow=c(3,2))

matplot(x=ts, y=rt1, type='l', xlab="Temp", ylab="density", main='No DD, preds higher Topt')
matplot(x=ts, y=rt2, type='l', xlab="Temp", ylab="density",main='No DD, prey higher Topt')
matplot(x=ts, y=rt1DD, type='l', xlab="Temp", ylab="density",main='DD in births, preds higher Topt')
matplot(x=ts, y=rt2DD, type='l', xlab="Temp", ylab="density",main='DD in births, prey higher Topt')
matplot(x=ts, y=rt1DD3, type='l', xlab="Temp", ylab="density",main='DD in all, preds higher Topt')
matplot(x=ts, y=rt2DD3, type='l', xlab="Temp", ylab="density",main='DD in all, prey higher Topt')


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
  return(cbind(js,As,ps,paramsensitivity))
}

sK<-sensitivity(parameter=8,from=5,to=1000,by=1,ODE=predprey_equations_DD,parms1=parms)

sK
par(mfrow=c(1,1))
matplot(x=sK[,4], y=sK[,1:3], type='l', xlab="K", ylab="density")

par(mfrow=c(4,2))

s1<-sensitivity(parameter=1,from=0,to=10,by=0.01,ODE=predprey_equations_DD,parms1=parms)

matplot(x=s1[,4], y=s1[,1:3], type='l', xlab="a", ylab="density")

s2<-sensitivity(parameter=2,from=0,to=10,by=0.01,ODE=predprey_equations_DD,parms1=parms)

matplot(x=s2[,4], y=s2[,1:3], type='l', xlab="c", ylab="density")

s3<-sensitivity(parameter=3,from=0,to=10,by=0.01,ODE=predprey_equations_DD,parms1=parms)

matplot(x=s3[,4], y=s3[,1:3], type='l', xlab="m", ylab="density")


s4<-sensitivity(parameter=4,from=0,to=10,by=0.01,ODE=predprey_equations_DD,parms1=parms)

matplot(x=s4[,4], y=s4[,1:3], type='l', xlab="b", ylab="density")

##g
s5<-sensitivity(parameter=5,from=0,to=10,by=0.01,ODE=predprey_equations_DD,parms1=parms)

matplot(x=s5[,4], y=s5[,1:3], type='l', xlab="g", ylab="density")

#n
s6<-sensitivity(parameter=6,from=0,to=10,by=0.01,ODE=predprey_equations_DD,parms1=parms)

matplot(x=s6[,4], y=s6[,1:3], type='l', xlab="n", ylab="density")

#K
s8<-sensitivity(parameter=8,from=1,to=10000,by=1,ODE=predprey_equations_DD,parms1=parms)

matplot(x=s8[,4], y=s8[,1:3], type='l', xlab="K", ylab="density")

