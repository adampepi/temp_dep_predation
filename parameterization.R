##Code to get parameters
library(boot)
library(ggplot2)
library(reshape2)
library(lme4)
library(lsmeans)
library(grid)
library(nlme)

setwd("~/Documents/Research/dissertation/temp dep pred prey modelling/temp_dep_predation")
fullpred<-read.csv('tentpredfull.csv')
tsurv<-read.csv('rearingsurv.csv')
tgrowth1<-read.csv('rearing2018.csv')
tgrowth2<-read.csv('rearing20182.csv')
###Warming tent predation experiment

str(fullpred)
fullpred
glmer1<-glmer(cbind(fullpred$present, fullpred$missing)~treatment*formica.summed+trial+(1|site)+(1|id),
              family=binomial, data=fullpred)

lsmeans(glmer1,specs=1~treatment*formica.summed,type='response')

(1-0.73)*3/15.1/3
##survival*3 caterpillars / number of ants/number of days
#############
###Growth and survival###

str(tgrowth1)
tgrowth2
str(tgrowth2)
#tgrowth1$temp<-as.factor(tgrowth1$temp) # Makes later models not work- but allows
#tgrowth2$temp<-as.factor(tgrowth2$temp)

#Combine dataset
tgrowth1_long <- melt(tgrowth1, id=c("temp", "id",'site',"avg.summer.t",'x.coord','y.coord','state'))
#tgrowth1_long$tempcat<-as.factor(tgrowth1_long$temp)
str(tgrowth1_long)
head(tgrowth1_long)

ggplot(tgrowth1_long,aes(x=variable, y=log(value+1), fill=temp))+ geom_boxplot()

#Combine first and second weighings into separate files
rate1<-tgrowth1[,1:6]
rate2<-tgrowth2[,1:6]
rate1
rate1$rgr1<-tgrowth1$Jul.10/tgrowth1$Jun.22
rate1$rgr2<-tgrowth1$Jul.27/tgrowth1$Jul.10
rate2$rgr1<-tgrowth2$Aug.02/tgrowth2$Jul.19
rate2$rgr2<-tgrowth2$X08.14.2017..dont.use./tgrowth2$Aug.02
rate1$state<-tgrowth1$state
rate2$state<-tgrowth2$state
#Make temp into separate categorical factor
rate1$tempcat<-as.factor(rate1$temp)
rate2$tempcat<-as.factor(rate2$temp)

ggplot(rate1 ,aes(x=tempcat, y=rgr1))+ geom_boxplot()
ggplot(rate1 ,aes(x=tempcat, y=rgr2))+ geom_boxplot()

ggplot(rate2 ,aes(x=tempcat, y=rgr1))+ geom_boxplot()
ggplot(rate2 ,aes(x=tempcat, y=rgr2))+ geom_boxplot()



rate1$rgr1<-log(rate1$rgr1)^1/18
rate2$rgr1<-log(rate2$rgr1)^1/13

str(rate1)
str(rate2)
rate.all<-rbind(rate1,rate2)#Combine both into a single file

head(rate.all)
rate.all[195,] 
rate.all<-rate.all[-195,] #remove ridiculous outlier
str(rate.all)



str(rate.all)
rate.1<-rate.all[,-8] #Omits too many if rgr2 is left in
rate.all1<-na.omit(rate.1) #Remove NAs- gnls can't handle them
str(rate.1)


str(rate.all1)
nls2<-nlsList(rgr1~Rm*exp(-0.5*((temp-Topt)/Tsd)^2)|site,data=rate.all1, 
              start=list(Topt=Topt_start,Tsd=Tsd_start,Rm=Rm_start))
summary(nls2)


str(tsurv)
tsurv



options(max.print = 200)
tsurv
tsurv$t1<-tsurv$Temp-10

logit(tsurv$survival)

ggplot(tsurv, aes(x=Temp, y=survival))+geom_count()+geom_smooth()

glm1<-glm (survival~t1+I(t1^2), family=binomial, data=tsurv)  
summary(glm1)

Topt_start<-20
Tsd_start<-5
Rm_start<-2

nls1<-nls(survival~Rm*exp(-0.5*((Temp-Topt)/Tsd)^2),data=tsurv, 
          start=list(Topt=Topt_start,Tsd=Tsd_start,Rm=Rm_start))
summary(nls1)
#min death rate
(1-0.75299)/81


newdata1<-expand.grid(Temp=seq(10,30,1))

str(newdata1)
newdata1$survival<-0

newdata1$survival<-predict(nls1,newdata=newdata1)
str(newdata1)



ggplot(newdata1, aes(x=Temp,y=survival))+geom_line()+geom_count(data=tsurv)


ggplot(data=tsurv,aes(x=Temp,y=X50mg))+ geom_smooth(method='glm',method.args = list(family = "Gamma"))+geom_count()

glm2<-glm(X50mg~Temp,family="Gamma",data=tsurv)
summary(glm2)
inv(coef(glm2))

n<-data.frame(Temp=23.42)

predict(glm2,n,type="response")
1/48.58318 
