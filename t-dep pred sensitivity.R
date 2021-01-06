library(deSolve)
library(rootSolve)

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
parms <- c(a,c,m,b,g,n )


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

gs<-taylor(ts,Topt = 31,Rm=0.1,Tsd=5)
as<-taylor(ts,Topt = 31,Rm=0.1,Tsd=5)



##Temperature response function
taylor<-function(temp,Rm=10,Topt=20,Tsd=1)Rm*exp(-0.5*((temp-Topt)/Tsd)^2)
ts<-seq(20,40,0.1)
##Temperature sensitivity analysis
tempresponse<-function(ts,ToptJ,ToptP,Tsd=5,RmP=.1,RmJ=.4,parms=parms){
gs<-taylor(ts,Topt = ToptJ,Rm=RmJ,Tsd=Tsd)
as<-taylor(ts,Topt = ToptP,Rm=RmP,Tsd=Tsd)

js<-numeric(length(as))
As<-numeric(length(as))
ps<-numeric(length(as))

for(i in 1:length(as)){
  
  parms1 <- c(a=as[i],c=parms[2],m=parms[3],b=parms[4],g=gs[i],n=parms[6] )
  rp<-as.data.frame(steady(y=Y0,time=c(0,1e5),func=predprey_equations,parms=parms1,method='runsteady')) #steady function from rootsolve package gets steady state
  js[i]<-rp[1,1]
  As[i]<-rp[2,1]
  ps[i]<-rp[3,1]
}
return(cbind(js,As,ps))
}

rt1<-tempresponse(ts,ToptJ=25,ToptP=30,Tsd=5,RmP=.1,RmJ=.4, parms=parms)
rt1
rt2<-tempresponse(ts,ToptJ=30,ToptP=25,Tsd=5,RmP=.1,RmJ=.4, parms=parms)
rt2

par(mfrow=c(2,1))
matplot(x=ts, y=rt1, type='l', xlab="Temp", ylab="density")
matplot(x=ts, y=rt2, type='l', xlab="Temp", ylab="density")
