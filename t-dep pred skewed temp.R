# t-dep pred with skewed temp functions

library(deSolve)
library(rootSolve)

#time vector
Tmax = 1000 # time horizon  
TimeStep = 1 # integration time step
Time <- seq(0, Tmax, by = TimeStep)
#initial
A0 = 10 # initial adults
P0 = 11 # initial predators
J0 = 100 # initial juveniles
Y0 <- c(A0,P0,J0)
#parameters
a =.1 # attack rate (pred)
c =.1 # conversion
m =.1 # predator mortality (predators live longer)
b = 2 # birth rate, must be greater than 1
g =.4 # growth rate (juve)
n =.3 # adult prey mortality
o =.3 
K =100
parms <- c(a,c,m,b,g,n,o,K )

# predator-prey equations ####
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

LV.out <- lsoda(Y0, Time,predprey_equations_DD,parms)
par(las=1,bty="l")
matplot(LV.out[,1],LV.out[,-1], type='l', xlab="time", ylab="density") 


##Temperature response function ####
taylor<-function(temp,Rm=10,Topt=20,Tsd=1)Rm*exp(-0.5*((temp-Topt)/Tsd)^2)
ts<-seq(20,40,0.1)
gs<-taylor(ts,Topt = 31,Rm=0.1,Tsd=5)
as<-taylor(ts,Topt = 31,Rm=0.1,Tsd=5)

## Temperature response function skewed - Gaussian * Gompertz ####
skew <- function(temp,Rm=10,Topt=20,rho,sigma)Rm*exp(-exp((rho*(temp-Topt))-6)-(sigma*((temp-Topt)^2)))
ts<-seq(20,40,0.1)
gs<-skew(ts,Topt = 31,Rm=0.1,rho=0.9,sigma=0.01)
#plot(gs)
as<-skew(ts,Topt = 31,Rm=0.1,rho=0.9,sigma=0.01)
#plot(as)

# Rm, R-max, is height of the curve's peak

##Temperature sensitivity analysis ####
tempresponse<-function(ts, ##sequence of temperatures
                       ToptJ,  ##topt for juveniles
                       ToptP,  ##topt for predators
                       Tsd=5,  ##thermal niche breadth
                       RmP=.1, ##max predation rate
                       RmJ=.4, ##max growth rate
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


par(mfrow=c(3,2))

matplot(x=ts, y=rt1, type='l', xlab="Temp", ylab="density", main='No DD, preds higher Topt')
matplot(x=ts, y=rt2, type='l', xlab="Temp", ylab="density",main='No DD, prey higher Topt')
matplot(x=ts, y=rt1DD, type='l', xlab="Temp", ylab="density",main='DD in births, preds higher Topt')
matplot(x=ts, y=rt2DD, type='l', xlab="Temp", ylab="density",main='DD in births, prey higher Topt')
matplot(x=ts, y=rt1DD3, type='l', xlab="Temp", ylab="density",main='DD in all, preds higher Topt')
matplot(x=ts, y=rt2DD3, type='l', xlab="Temp", ylab="density",main='DD in all, prey higher Topt')