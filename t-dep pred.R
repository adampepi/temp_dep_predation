#t-dep stage-dep predation 
library(deSolve)

par(mfrow=c(1,1))

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

#condition 2: must be pos
g-n*(n/(b*m)+m+b)
#condition 3: must be pos
n*g+g*b-m*g*b/n-m*g+m*b-n*m

#red = pred
#green = juv
#black = adult


#equilibrium across range of temps using *gaussian function* for a(t) and g(t)
# Rm*exp(-0.5*((T-Topt)/Tsd)^2
# no density dependence
par(mfrow=c(3,2))
# J*=m/ca(t)
# A*=g(t)m/nca(t)
# P*=[g(t)(b/n-1)]/a(t)
T = seq(10,30,1)
Topt = 20
Tsd = 5
### identical ####
at = a*exp(-0.5*((T-Topt)/Tsd)^2) #attack rate as func of temp
gt = g*exp(-0.5*((T-Topt)/Tsd)^2) #growth rate as func of temp

J_hat=m/(c*at)
A_hat=gt*m/(n*c*at)
P_hat=(gt*((b/n)-1))/at
plot(T,J_hat, col="green",type="l", main="No density dependence")
lines(T,A_hat, col="black")
lines(T,P_hat, col="red")

# add density dependence on reproduction
# J*=m/ca(t)
# A*=g(t)m/nca(t)
# P*=bg(t)/na(t)-bg(t)m/Knca(t)2+bg(t)2m/Kn2ca(t)2-g(t)/a(t)
K=200
J_hat=m/(c*at)
A_hat=gt*m/(n*c*at)
P_hat=b*gt/(n*at)-b*gt*m/(K*n*c*at^2)+b*gt^2*m/K*n^2*c*at^2-gt/at

plot(T,J_hat, col="green",type="l",main = "DD on repro")
lines(T,A_hat, col="black")
lines(T,P_hat, col="red")

### growth warmer ####
Topt2 = 25
at = a*exp(-0.5*((T-Topt)/Tsd)^2) #attack rate as func of temp
gt = g*exp(-0.5*((T-Topt2)/Tsd)^2) #growth rate as func of temp
J_hat=m/(c*at)
A_hat=gt*m/(n*c*at)
P_hat=b*gt/(n*at)-b*gt*m/(K*n*c*at^2)+b*gt^2*m/K*n^2*c*at^2-gt/at
plot(T,J_hat, col="green",type="l",main="Growth Topt warmer")
lines(T,A_hat, col="black")
lines(T,P_hat, col="red")

### growth wider ####
Tsd2 = 2*Tsd
at = a*exp(-0.5*((T-Topt)/Tsd)^2) #attack rate as func of temp
gt = g*exp(-0.5*((T-Topt)/Tsd2)^2) #growth rate as func of temp
J_hat=m/(c*at)
A_hat=gt*m/(n*c*at)
P_hat=b*gt/(n*at)-b*gt*m/(K*n*c*at^2)+b*gt^2*m/K*n^2*c*at^2-gt/at
plot(T,J_hat, col="green",type="l",main="Growth Tsd wider")
lines(T,A_hat, col="black")
lines(T,P_hat, col="red")

### growth wider & warmer####
Tsd2 = 2*Tsd
at = a*exp(-0.5*((T-Topt)/Tsd)^2) #attack rate as func of temp
gt = g*exp(-0.5*((T-Topt2)/Tsd2)^2) #growth rate as func of temp
J_hat=m/(c*at)
A_hat=gt*m/(n*c*at)
P_hat=b*gt/(n*at)-b*gt*m/(K*n*c*at^2)+b*gt^2*m/K*n^2*c*at^2-gt/at
plot(T,J_hat, col="green",type="l",main="Growth wider+warmer")
lines(T,A_hat, col="black")
lines(T,P_hat, col="red")

### attack wider ####
Tsd2 = 0.5*Tsd
at = a*exp(-0.5*((T-Topt)/Tsd)^2) #attack rate as func of temp
gt = g*exp(-0.5*((T-Topt)/Tsd2)^2) #growth rate as func of temp
J_hat=m/(c*at)
A_hat=gt*m/(n*c*at)
P_hat=b*gt/(n*at)-b*gt*m/(K*n*c*at^2)+b*gt^2*m/K*n^2*c*at^2-gt/at
plot(T,J_hat, col="green",type="l",main="Attack Tsd wider")
lines(T,A_hat, col="black")
lines(T,P_hat, col="red")

#maybe use rootSolve package function "multiroot" and function "steady"?
