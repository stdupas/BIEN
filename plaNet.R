setwd("/home/dupas/PlaNet/")
setwd("C:/Users/steph/OneDrive/Documents/GitHub/PlaNet")

## PlaNet

## Proof of concept
# ecolink
#     T	H	Pa	Pl	Pe
# T			x	x
# H			x	x
# Pa			x	x	x
# Pl			x	x
# Pe			x

# indlink
#     iT	iH	iPa	iPl	iPe
# iT	 x
# iH			x
# iPa			     x
# iPl			        x
# iPe			            x
setwd("/home/dupas/PlaNet/")
ecoVar=c("Tx","Pr","Pa","Pl","Pe") 
indicVar=c("iT","iPr","iPa","iPl","iPe")

library(raster)

ecoLink <- t(matrix(as.logical(c(0,0,1,1,0,
                  0,0,1,1,0,
                  0,0,1,1,1,
                  0,0,1,1,0,
                  0,0,1,0,0)),nrow=5, ncol =5,dimnames=list(ecoVar,ecoVar)))
dimnames(ecoLink)=list(ecoVar,ecoVar)

ecoLinkTime = t(matrix((c(0,0,1,1,0,
                          0,0,2,1,0,
                          0,0,1,1,1,
                          0,0,1,1,0,
                          0,0,0,0,0)),nrow=5,ncol=5,dimnames=list(ecoVar,ecoVar)))
dimnames(ecoLinkTime)=list(ecoVar,ecoVar)

indicLink <- as.logical(diag(5))
dimnames(indicLink)=list(ecoVar,ecoVar)

# data is a 3 dim array: dim[1] is indicator and ecosystem variables, dim[2] is population, dim[3] is time

setClass("Data",
         contains="array",
         validity=function(object){
           if (length(dim(object))!=3) stop("data should be a 3 dim array, dim[1] is indicator and ecosystem variables, dim[2] is population, dim[3] is time")
         }
         )
#save(dataf,file = "yield.data.RData")
setwd()

load("yield.data.RData")
dataf[,,1]
load(file = "yield.data.RData")
idata <- new("Data",dataf[1:8,c(8,4,2,2,2),])
dimnames(idata) <- list(dimnames(idata)[[1]],c("Tx","Pr","Pa","Pl","Pe"),dimnames(idata)[[3]])
edata <- new("Data",dataf[1:8,c(8,4,2,2,2),])
dimnames(edata) <- list(dimnames(edata)[[1]],c("Tx","Pr","Pa","Pl","Pe"),dimnames(edata)[[3]])
edata[,3,]<-NA
edata[1,4,]<-0;edata[2,4,]<-0.01 # planting offucrs in february
edata[1,4,]<-0;edata[2,4,]<-0.01 # planting offucrs in february
edata[1:2,3,]<-0 # there is no parasite before planting
edata[1:2,5,]<-0 # there is no pesticide before planting
edata[,,1:3]

p=list(PaPa=c(rmax=10,K=20),TxPa=c(min=15,max=30),PrPa=c(min=3),PlPa=c(r=1),PrPl=c(rperPr=.5),TxPl=c(min=10,max=30),PaPl=c(r=-0.03),PlPl=c(r=2.5,K=2,sd=.1),PePa=c(r=0.1),PaPe=c(thr=1.5))
names(p)


#
# Simulating ecosystem history
#

p_PaPa_rmax=10;p_PaPa_K=20
p_TxPa_min=15;p_TxPa_max=30
p_PrPa_min=3
p_PlPa_r=1
p_PrPl_rperPr=.5
p_TxPl_min=10;p_TxPl_max=30
p_PaPl_r=-0.03
p_PlPl_r=2.5;p_PlPl_K=2;p_PlPl_sd=.1
p_PePa_r=0.1
p_PaPe=1.5
p_Pl_sd_r=0.05
p_T_sd=0.3
p_Pe_pDet=.8


# simulating ecosystem 

# using lapply

tmp=aperm(array(unlist(lapply(1:dim(edata)[3], function(k){lapply(3:dim(edata)[1], function (i) {
  Pl=((edata[i-1,1,k]>p_TxPl_min)*0.1)+edata[i-1,4,k]*p_PlPl_r*((((1+edata[i-1,2,k])*p_PrPl_rperPr)*(edata[i-1,1,k]>p_TxPl_min))*
                                                                 (1+(T-p_TxPl_min)/(p_TxPl_max-p_TxPl_min))*(edata[i-1,1,k]<p_TxPl_max))+edata[i-1,3,k]*(p_PaPl_r)
  if(Pl>p_PlPl_K){Pl=2}
  if(Pl<0) {Pl=0}
  Pa=(edata[i-1,3,k]==0)+((!edata[i,5,k])+edata[i,5,k]*p_PePa_r)*edata[i-1,3,k]*p_PaPa_rmax*((edata[i-1,1,k]>p_TxPa_min)*(edata[i-1,1,k]<p_TxPa_max))*(edata[i-1,1,k]-p_TxPa_min)/(p_TxPa_max-p_TxPa_min)*(edata[i-2,2,k]>p_PrPa_min)*(edata[i-1,4,k]*p_PlPa_r)
  c(Pa=rpois(1,Pa*(Pa<p_PaPa_K)+p_PaPa_K*(Pa>=p_PaPa_K)),
    Pl=Pl,Pe=(edata[i-1,3,k]>p_PaPe))
})})),dim=dim(edata)[c(2,1,3)],dimnames=lapply(c(2,1,3),function(i) dimnames(edata)[[i]])),c(2,1,3))

ed=edata
for(k in 1:dim(edata)[3]) for (i in 3:dim(edata)[1]){
  Pl=((edata[i-1,1,k]>p_TxPl_min)*0.1)+edata[i-1,4,k]*p_PlPl_r*((((1+edata[i-1,2,k])*p_PrPl_rperPr)*(edata[i-1,1,k]>p_TxPl_min))*
                                                                  (1+(T-p_TxPl_min)/(p_TxPl_max-p_TxPl_min))*(edata[i-1,1,k]<p_TxPl_max))+edata[i-1,3,k]*(p_PaPl_r)
  if(Pl>p_PlPl_K){Pl=2}
  if(Pl<0) {Pl=0}
  Pa=(edata[i-1,3,k]==0)+((!edata[i,5,k])+edata[i,5,k]*p_PePa_r)*edata[i-1,3,k]*p_PaPa_rmax*((edata[i-1,1,k]>p_TxPa_min)*(edata[i-1,1,k]<p_TxPa_max))*(edata[i-1,1,k]-p_TxPa_min)/(p_TxPa_max-p_TxPa_min)*(edata[i-2,2,k]>p_PrPa_min)*(edata[i-1,4,k]*p_PlPa_r)
  c(Pa=rpois(1,Pa*(Pa<p_PaPa_K)+p_PaPa_K*(Pa>=p_PaPa_K)),
    Pl=Pl,Pe=(edata[i-1,3,k]>p_PaPe))
  ed[i,3:5,k]=c(Pa,Pl,(edata[i-1,3,k]>p_PaPe))
}

Pl=((edata[i-1,1,k]>p_TxPl_min)*0.1)+edata[i-1,4,k]*p_PlPl_r*((((1+edata[i-1,2,k])*p_PrPl_rperPr)*(edata[i-1,1,k]>p_TxPl_min))*
                                                                (1+(T-p_TxPl_min)/(p_TxPl_max-p_TxPl_min))*(edata[i-1,1,k]<p_TxPl_max))+edata[i-1,3,k]*(p_PaPl_r)
if(Pl>p_PlPl_K){Pl=2}
if(Pl<0) {Pl=0}
Pa=(edata[i-1,3,k]==0)+((!edata[i,5,k])+edata[i,5,k]*p_PePa_r)*edata[i-1,3,k]*p_PaPa_rmax*((edata[i-1,1,k]>p_TxPa_min)*(edata[i-1,1,k]<p_TxPa_max))*(edata[i-1,1,k]-p_TxPa_min)/(p_TxPa_max-p_TxPa_min)*(edata[i-2,2,k]>p_PrPa_min)*(edata[i-1,4,k]*p_PlPa_r)
c(Pa=rpois(1,Pa*(Pa<p_PaPa_K)+p_PaPa_K*(Pa>=p_PaPa_K)),
  Pl=Pl,Pe=(edata[i-1,3,k]>p_PaPe))

#
# simulating ecosystem good one
#

for (k in 1:dim(edata)[3]){
  for (i in 3:dim(edata)[1]){
    #Pe
    edata[i,5,k]=edata[i-1,3,k]>p_PaPe
    #Pl
    a=((edata[i-1,1,k]>p_TxPl_min)*0.1)+edata[i-1,4,k]*p_PlPl_r*((((1+edata[i-1,2,k])*p_PrPl_rperPr)*(edata[i-1,1,k]>p_TxPl_min))*
                                                                   (1+(T-p_TxPl_min)/(p_TxPl_max-p_TxPl_min))*(edata[i-1,1,k]<p_TxPl_max))+edata[i-1,3,k]*(p_PaPl_r)
    if(a>p_PlPl_K){a=2}
    if(a<0) {a=0}
    edata[i,4,k]=a
    #Pa
    a=(edata[i-1,3,k]==0)+((!edata[i,5,k])+edata[i,5,k]*p_PePa_r)*edata[i-1,3,k]*p_PaPa_rmax*((edata[i-1,1,k]>p_TxPa_min)*(edata[i-1,1,k]<p_TxPa_max))*(edata[i-1,1,k]-p_TxPa_min)/(p_TxPa_max-p_TxPa_min)*(edata[i-2,2,k]>p_PrPa_min)*(edata[i-1,4,k]*p_PlPa_r)
    edata[i,3,k]=rpois(1,a*(a<p_PaPa_K)+p_PaPa_K*(a>=p_PaPa_K))
  }
}

#
# Calculating probability of ecosystem model
#

dims=dim(edata)
for(k in 1:dims[3]) for (i in 3:dims[1]){

  }

idata=aperm(array(unlist(lapply(1:dim(edata)[3], function(k){lapply(1:dim(edata)[1], function (i) {
  Pl=pgamma(1,shape=edata[i,4,k],rate=p_Pl_var_r)
  if (Pl<0) Pl=0
  c(T=rnorm(1,edata[i,1,k],p_T_sd),Pr=rpois(1,edata[i,2,k]),Pa=rpois(1,edata[i,3,k]),Pl=Pl,Pe=rbinom(1,edata[i,5,k],p_Pe_pDet))
})})),dim=dim(edata)[c(2,1,3)],dimnames=lapply(c(2,1,3),function(i) dimnames(edata)[[i]])),c(2,1,3))



load(file = "edata.RData")
n=length(edata[,1,])

simulate idata from edata

#Pe
  idata[,5,]=(edata[,5,])*rbinom(n,1,p_Pe_pDet)
#Pl
  idata[,4,]=rgamma(n,edata[,4,],rate=p_Pl_var_r)
#Pa
  idata[,3,]=rpois(n,edata[,3,])
#T
  idata[,1,]=rnorm(n,edata[,1,],p_T_sd)
#Pr
  idata[,2,]=rpois(n,edata[,2,])

#for (k in 1:dim(edata)[3]){
#	for (i in 3:dim(edata)[1]){
#Pe
#  if (edata[i,5,k]) idata[i,5,k]=rbinom(1,1,p_Pe_pDet)
#Pl
#  idata[i,4,k]=rgamma(length(edata[,4,]),edata[,4,],rate=p_Pl_var_r)
#rnorm(1,edata[i,4,k],edata[i,4,k]*p_Pl_sd_r)
#  if (idata[i,4,k]<0) idata[i,4,k]=0
#Pa
#  idata[i,3,k]=rpois(1,edata[i,3,k])
#T
#  idata[i,1,k]=rnorm(1,edata[i,1,k],p_T_sd)
#Pr
#  idata[i,2,k]=rpois(1,edata[i,2,k])
#  }
#}

#
# Simulation of indicator data from edata good one
#

idata=aperm(array(unlist(lapply(1:dim(edata)[3], function(k){lapply(1:dim(edata)[1], function (i) {
Pl=rnorm(1,edata[i,4,k],edata[i,4,k]*p_Pl_sd_r)
if (Pl<0) Pl=0
c(T=rnorm(1,edata[i,1,k],p_T_sd),Pr=rpois(1,edata[i,2,k]),Pa=rpois(1,edata[i,3,k]),Pl=Pl,Pe=rbinom(1,edata[i,5,k],p_Pe_pDet))
})})),dim=dim(edata)[c(2,1,3)],dimnames=lapply(c(2,1,3),function(i) dimnames(edata)[[i]])),c(2,1,3))

#
# Simulation of edata from idata
#
# p(i/e)=p(e/i)*p(i)/p(e)
# p(e/i)=p(i/e)*p(e)/p(i)
# p(H)=dgamma3(x,1,5.749)
# P(T)=dnorm(1,10,4)*.55+dnorm(1,21,3)*.45
# P(Pa)=dgamma(x+1,.2,3)
# P(Pl)=rnorm(1,edata[i,4,k],edata[i,4,k]*p_Pl_sd_r)
# plot(dnorm(1:40,10,4)*.6+dnorm(1:40,4)*.4)
# ?dnorm
library(FAdist)
.1*1.5^8
par(mfrow=c(1,2))
hist(idata[,4,])
plot((0:10)/5,dexp((1:11),1.3,1.5))
plot(-5:35,dnorm(-5:35,10,4)*.55+dnorm(-5:35,21,3)*.45)
hist(rpois(length(idata[,2,]),.1))
     plot(1:20,dgamma3(1:20,1,2))
     pgamma3(1,1,5.749)
     sum(edata[,2,]<1)/length(edata[,2,])
     pgamma3(2,1,5.749)
     sum((edata[,2,]<2)&(edata[,2,]>1))/length(idata[,2,])
     ?rgamma3
dimnames(idata)[[2]]
idata=aperm(array(unlist(lapply(1:dim(edata)[3], function(k){lapply(1:dim(edata)[1], function (i) {
  Pl=rnorm(1,edata[i,4,k],edata[i,4,k]*p_Pl_sd_r)
  if (Pl<0) Pl=0
  c(T=rnorm(1,edata[i,1,k],p_T_sd),Pr=rpois(1,edata[i,2,k]),Pa=rpois(1,edata[i,3,k]),Pl=Pl,Pe=rbinom(1,edata[i,5,k],p_Pe_pDet))
})})),dim=dim(edata)[c(2,1,3)],dimnames=lapply(c(2,1,3),function(i) dimnames(edata)[[i]])),c(2,1,3))


#
# Probability of indicator data
#

hist(edata[,4,])
idata
log(exp(1))
pT <- sum(dnorm(idata[,1,],edata[,1,],p_T_sd,log=TRUE))
pH <- sum(dpois(idata[,2,],edata[,2,],log=TRUE)))
pPa <- sum(dpois(idata[,3,],edata[,3,],log=TRUE))
sdPl=edata[,4,]*p_Pl_sd_r
sdPl[sdPl<=0.1]=.1
pPl <- sum(dnorm(x=idata[,4,],mean=edata[,4,],sd=sdPl,log=TRUE))
pPe <- dbinom(idata[,5,],1,prob=edata[,5,]*p_Pe_pDet,log=TRUE)

pIndic <- sum(c(pT = sum(log(dnorm(idata[,1,],edata[,1,],p_T_sd))),
                 pH = sum(log(dpois(idata[,2,],edata[,2,]))),
                 pPa = sum(log(dpois(idata[,3,],edata[,3,]))),
                 sdPl={sdPl=edata[,4,]*p_Pl_sd_r
                 sdPl[sdPl<=0.1]=.1
                 sum(log(dnorm(x=idata[,4,],mean=edata[,4,],sd=sdPl)))},
                 pPe=sum(dbinom(idata[,5,],1,prob=edata[,5,]*p_Pe_pDet,log=TRUE))
))
edata[1:8,,1]
idata[1:8,,1]

p_PaPa_rmax=10
p_PaPa_K=20
p_TxPa_min=15
p_TxPa_max=30
p_PrPa_min=3
p_PlPa_r=1
p_PrPl_rperPr=.5
p_TxPl_min=10
p_TxPl_max=30
p_PaPl_r=-0.03
p_PlPl_r=2.5
p_PlPl_K=2
p_PlPl_sd=.1
p_PePa_r=0.1
p_PaPe=1.5
p_Pl_sd_r=0.05
p_T_sd=0.3
p_Pe_pDet=.8
p_Pe_pFalseDet=.005


#ecosysHistory
p0 = c(p_PaPa_rmax=exp(runif(1,log(5),log(15))),p_PaPa_K=exp(runif(1,log(15),log(25))),
       p_TxPa_min=runif(1,10,20),p_TxPa_max=runif(1,25,35),
       p_PrPa_min=runif(1,1.5,4),p_PlPa_r=runif(1,.7,1.5),
       p_PrPl_rperPr=runif(1,.3,.7),p_TxPl_min=runif(1,5,14),p_TxPl_max=runif(1,26,32),
       p_PaPl_r=exp(runif(1,log(0.015),log(0.045))),
       p_PlPl_r=exp(runif(1,log(1.8),log(3.5))),p_PlPl_K=exp(runif(1,log(1.5),log(3))),p_PlPl_sd=runif(1,.07,.15),
       p_PePa_r=exp(runif(1,log(.06),log(.15))),
       p_PaPe=runif(1,1.3,1.9),
       #ecoindic
       p_Pl_sd_r=exp(runif(1,log(.03),log(.07))),
       p_T_sd=runif(1,.2,.5),
       p_Pe_pDet=runif(1,.6,.99),
       p_Pe_pFalseDet=exp(runif(1,log(.001),log(.02)))
)

#
# Algorithm
#
setwd("/home/dupas/PlaNet/")
ecoVar=c("Tx","Pr","Pa","Pl","Pe") 
indicVar=c("iT","iPr","iPa","iPl","iPe")
load("yield.data.RData")
dataf[,,1]
load(file = "yield.data.RData")

setClass("Data",
         contains="array",
         validity=function(object){
           if (length(dim(object))!=3) stop("data should be a 3 dim array, dim[1] is indicator and ecosystem variables, dim[2] is population, dim[3] is time")
         }
)

idata <- new("Data",dataf[1:8,c(8,4,2,2,2),])
dimnames(idata) <- list(dimnames(idata)[[1]],c("Tx","Pr","Pa","Pl","Pe"),dimnames(idata)[[3]])
edata <- new("Data",dataf[1:8,c(8,4,2,2,2),])
dimnames(edata) <- list(dimnames(edata)[[1]],c("Tx","Pr","Pa","Pl","Pe"),dimnames(edata)[[3]])
edata[,3,]<-NA
edata[1,4,]<-0;edata[2,4,]<-0.01 # planting offucrs in february
edata[1,4,]<-0;edata[2,4,]<-0.01 # planting offucrs in february
edata[1:2,3,]<-0 # there is no parasite before planting
edata[1:2,5,]<-0 # there is no pesticide before planting
edata[,,1:3]

#library(raster)

#edata
#idata



# EXAMPLE OF LEARNING ECOSYSTEM
# simulate ecosystem
# this is a an annual plant pathogen interaction system
# start from the simulation of ecosystem history and indicator data from 
# a true model then assume we have access to indicator data 
# to infer the model using prior.
# The full ecosystem data is simulated from indicator and prior 
# The likelihood of is estimated from ecosystem history
# The posterior is calculated and sampled using metropolis algorithm

setwd("/home/dupas/PlaNet/")
setwd("C:/Users/steph/OneDrive/Documents/GitHub/PlaNet")

# set true parameters

#true parameters
p0=list(PaPa_rmax=10,PaPa_K=20,
     TxPa_min=15,TxPa_max=30,
     PrPa_min=3,
     PlPa_r=1,
     PrPl_rperPr=.5,
     TxPl_min=10,TxPl_max=30,
     PaPl_r=-0.03,
     PlPl_r=2.5,
     PlPl_K=2,
     PlPl_sd=.1,
     PePa_r=0.1,
     PaPe=1.5,
     Pl_var_r=0.05^2,
     T_sd=0.3,
     Pe_pDet=.8,
     Pe_pFalseDet=.005)


# set parameters to true parameters

p_PaPa_rmax=p0["PaPa_rmax"]
p_PaPa_K=p0["PaPa_K"]
p_TxPa_min=p0["TxPa_min"]
p_TxPa_max=p0["TxPa_max"]
p_PrPa_min=p0["PrPa_min"]
p_PlPa_r=p0["PlPa_r"]
p_PrPl_rperPr=p0["PrPl_rperPr"]
p_TxPl_min=p0["TxPl_min"]
p_TxPl_max=p0["TxPl_max"]
p_PaPl_r=p0["PaPl_r"]
p_PlPl_r=p0["PlPl_r"]
p_PlPl_K=p0["PlPl_K"]
p_PlPl_sd=p0["PlPl_sd"]
p_PePa_r=p0["PePa_r"]
p_PaPe=p0["PaPe"]
#ecoindic
p_Pl_var_r=p0["Pl_var_r"]
p_T_sd=p0["T_sd"]
p_Pe_pDet=p0["Pe_pDet"]
p_Pe_pFalseDet=p0["Pe_pFalseDet"]

#
# set time series
#

setwd("/home/dupas/PlaNet/")
ecoVar=c("Tx","Pr","Pa","Pl","Pe") 
indicVar=c("iT","iPr","iPa","iPl","iPe")
load("yield.data.RData")
dataf[,,1]
load(file = "yield.data.RData")

setClass("Data",
         contains="array",
         validity=function(object){
           if (length(dim(object))!=3) stop("data should be a 3 dim array, dim[1] is indicator and ecosystem variables, dim[2] is population, dim[3] is time")
         }
)

idata <- new("Data",dataf[1:8,c(8,4,2,2,2),])
dimnames(idata) <- list(dimnames(idata)[[1]],c("Tx","Pr","Pa","Pl","Pe"),dimnames(idata)[[3]])
edata <- new("Data",dataf[1:8,c(8,4,2,2,2),])
dimnames(edata) <- list(dimnames(edata)[[1]],c("Tx","Pr","Pa","Pl","Pe"),dimnames(edata)[[3]])
edata[,3,]<-NA
edata[1,4,]<-0;edata[2,4,]<-0.01 # planting offucrs in february
edata[1,4,]<-0;edata[2,4,]<-0.01 # planting offucrs in february
edata[1:2,3,]<-0 # there is no parasite before planting
edata[1:2,5,]<-0 # there is no pesticide before planting
edata[,,1:3]

# simulate true edata
for (k in 1:dim(edata)[3]){
  for (i in 3:dim(edata)[1]){
    #Pe
    edata[i,5,k]=edata[i-1,3,k]>p_PaPe
    #Pl
    a=((edata[i-1,1,k]>p_TxPl_min)*0.1)+edata[i-1,4,k]*p_PlPl_r*((((1+edata[i-1,2,k])*p_PrPl_rperPr)*(edata[i-1,1,k]>p_TxPl_min))*
                                                                   (1+(T-p_TxPl_min)/(p_TxPl_max-p_TxPl_min))*(edata[i-1,1,k]<p_TxPl_max))+edata[i-1,3,k]*(p_PaPl_r)
    if(a>p_PlPl_K){a=2}
    if(a<0) {a=0}
    edata[i,4,k]=a
    #Pa
    a=(edata[i-1,3,k]==0)+((!edata[i,5,k])+edata[i,5,k]*p_PePa_r)*edata[i-1,3,k]*p_PaPa_rmax*((edata[i-1,1,k]>p_TxPa_min)*(edata[i-1,1,k]<p_TxPa_max))*(edata[i-1,1,k]-p_TxPa_min)/(p_TxPa_max-p_TxPa_min)*(edata[i-2,2,k]>p_PrPa_min)*(edata[i-1,4,k]*p_PlPa_r)
    edata[i,3,k]=rpois(1,a*(a<p_PaPa_K)+p_PaPa_K*(a>=p_PaPa_K))
  }
}

# save true ecosystem
e0data <-edata
edataTrue <- edata


# simulate true idata
idata=aperm(array(unlist(lapply(1:dim(edata)[3], function(k){lapply(1:dim(edata)[1], function (i) {
  Pl=rgamma(1,edata[i,4,k],edata[i,4,k]*p_Pl_sd_r)
  #if (Pl<0) Pl=0
  c(T=rnorm(1,edata[i,1,k],p_T_sd),Pr=rpois(1,edata[i,2,k]),Pa=rpois(1,edata[i,3,k]),Pl=Pl,Pe=rbinom(1,edata[i,5,k],p_Pe_pDet)+rbinom(1,!edata[i,5,k],p_Pe_pFalseDet))
})})),dim=dim(edata)[c(2,1,3)],dimnames=lapply(c(2,1,3),function(i) dimnames(edata)[[i]])),c(2,1,3))

# save true idata
idataTrue <- idata

Posterior <-  data.frame(runi=0,p_PaPa_rmax=0,p_TxPa_min=0,p_PrPa_min=0,p_PlPa_r=0,p_PrPl_rperPr=0,p_TxPl_min=0,
                         p_PaPl_r=0,p_PlPl_r=0,p_PePa_r=0,p_PaPe=0,p_Pl_sd_r=0,p_T_sd=0,p_Pe_pDet=0,p_Pe_pFalseDet=0,prior=0,posterior=0)

Posterior=Posterior[-1,]
runi=1

# Sampling algorithm

samplePrior <- function(option="prior"){
  if(option=="prior") list(PaPa_rmax=exp(runif(1,log(5),log(15))),
                        PaPa_K=exp(runif(1,log(15),log(25))),
                        TxPa_min=runif(1,10,20),
                        TxPa_max=runif(1,25,35),
                        PrPa_min=runif(1,1.5,4),
                        PlPa_r=runif(1,.7,1.5),
                        PrPl_rperPr=runif(1,.3,.7),
                        TxPl_min=runif(1,5,14),
                        TxPl_max=runif(1,26,32),
                        PaPl_r=-exp(runif(1,log(0.001),log(0.045))),
                        PlPl_r=exp(runif(1,log(1.8),log(3.5))),
                        PlPl_K=exp(runif(1,log(1.5),log(3))),
                        PlPl_sd=runif(1,.07,.15),
                        PePa_r=exp(runif(1,log(.06),log(.15))),
                        PaPe=runif(1,1.3,1.9),
                        #ecoindic
                        Pl_var_r=(exp(runif(1,log(.03),log(.07))))^2,
                        T_sd=runif(1,.2,.5),
                        Pe_pDet=runif(1,.6,.99),
                        Pe_pFalseDet=exp(runif(1,log(.001),log(.02))))
}


getPprime <- function(p,rate=1/20){
list(PaPa_rmax=exp({a={if(as.logical(rbinom(1,1,1/2))) log(p$PaPa_rmax+(15-5)*rate) else log(p$PaPa_rmax-(15-5)*rate)};{ if (a<log(5)) log(5) else if (a>log(15)) log(15) else a}}),
  PaPa_K={a={rbinom(1,1,1/2);if(a==0) a=-1;a}*(log(25)-log(15))*rate;if (p$PaPa_K<15) p$PaPa_K=15;if (p$PaPa_K>15) p$PaPa_K=25;p$PaPa_K},
  TxPa_min={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(20-10)*rate;a=p$TxPa_min+b;if (a<10) a=15;if (a>20) a=20;a},
  TxPa_max={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(35-25)*rate;a=p$TxPa_max+b;if (a<25) a=25;if (a>35) a=35;a},
  PrPa_min={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(4-1.5)*rate;a=p$PrPa_min+b;if (a<1.5) a=1.5;if (a>4) a=4;a},
  PlPa_r={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(1.5-.7)*rate;a=p$PlPa_r+b;if (a<.7) a=.7;if (a>1.5) a=1.5;a},
  PrPl_rperPr={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.7-.3)*rate;a=p$PrPl_rperPr+b;if (a<.3) a=.3;if (a>.7) a=.7;a},
  TxPl_min={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(14-5)*rate;a=p$TxPl_min+b;if (a<5) a=5;if (a>14) a=14;a},
  TxPl_max={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(32-26)*rate;a=p$TxPl_max+b;if (a<26) a=26;if (a>32) a=32;a},
  PaPl_r=-exp({a={if(as.logical(rbinom(1,1,1/2))) log(-p$PaPl_r-(.01-.045)*rate) else log(-p$PaPl_r+(.01-.045)*rate)};{ if (a<log(.01)) log(.01) else if (a>log(.045)) log(.045) else a}}),
  PlPl_r=exp({a={if(as.logical(rbinom(1,1,1/2))) log(p$PlPl_r+(3.5-1.8)*rate) else log(p$PlPl_r-(3.5-1.8)*rate)};{ if (a<log(1.8)) log(1.8) else if (a>log(3.5)) log(3.5) else a}}),
  PlPl_K=exp({a={if(as.logical(rbinom(1,1,1/2))) log(p$PlPl_K+(3-1.5)*rate) else log(p$PlPl_K-(3-1.5)*rate)};{ if (a<log(1.5)) log(1.5) else if (a>log(3)) log(3) else a}}),
  PlPl_sd={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.15-.07)*rate;a=p$PlPl_sd+b;if (a<.07) a=.07;if (a>.15) a=.15;a},
  PePa_r={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.15-.06)*rate;a=p$PePa_r+b;if (a<.06) a=.06;if (a>.15) a=.15;a},
  PaPe={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(1.9-1.3)*rate;a=p$PaPe+b;if (a<1.3) a=1.3;if (a>1.9) a=1.9;a},
  Pl_var_r=exp({a={if(as.logical(rbinom(1,1,1/2))) log(p$Pl_var_r^.5+(.07-.03)*rate) else 2*log(p$Pl_var_r^.5-(.07-.03)*rate)};
  { if (a<2*log(.03)) 2*log(.07) else if (a>2*log(.07)) 2*log(.07) else a}}),
  T_sd={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.5-.2)*rate;a=p$T_sd+b;if (a<.2) a=.2;if (a>.5) a=.5;a},
  Pe_pDet={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.99-.6)*rate;a=p$Pe_pDet+b;if (a<.6) a=.6;if (a>.99) a=.99;a},
  Pe_pFalseDet={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.02-.001)*rate;a=p$Pe_pFalseDet+b;if (a<.001) a=.001;if (a>.02) a=.02;a}
)}


simul_edataFromidata <- function(idata){
  T = rnorm(1,round(idata[,1,]),p["p_T_sd"]))),
  pH = rpois(1,idata[,2,]),
  pPa = rpois(1,idata[,3,]),
  pPl = rgamma(n, shape=idata[,4,], scale = p_Pl_var_r)
  sdPl={sdPl=edata[,4,]*["p_Pl_sd_r"]
  sdPl[sdPl<=0.1]=.1
  sum(log(dnorm(x=idata[,4,],mean=edata[,4,],sd=sdPl)))},
  pPe={
    p1 <- which1 <- which(as.logical(edata[,5,]))
    p1 <- sum(dbinom(idata[,5,][which1],1,["p_Pe_pDet"],log=TRUE))
    p0 <- which0 <- which(!as.logical(edata[,5,]))
    p0 <- sum(dbinom(idata[,5,][which0],1,["p_Pe_pFalseDet"],log=TRUE))
    p0+p1
}
mean(rgamma(1000,2,scale=1/3)) 2/3
var(rgamma(100000,2,3)) #2*(1/3)^2 

# gamma
# mean=shape/rate=shape*scale
# var=shape/(rate^2)=shape*scale^2
# shape=mean^2/var
# scale = var/mean
# rate = 1/scale
# shape = mean*rate = mean*mean/var
# scale=(var/shape)^.5
# shape=shape^.5*mean/var^.5
# shape^.5=mean/var^.5
# scale=var^1.5/mean
# p_pl_sd_r = var.5/mean = 1/shape^.5
# 
# rate = shape/
# rate = p_pl_sd_r^2 = p_pl_var_r
#
# Initialisation of learning


# Sample first prior (or set from true parameters to get to the rightplace at first steps)

p_PaPa_rmax=exp(runif(1,log(5),log(15)))#p_PaPa_rmax=10
p_PaPa_K=exp(runif(1,log(15),log(25)))#p_PaPa_K=20
p_TxPa_min=runif(1,10,20)#p_TxPa_min=15
p_TxPa_max=runif(1,25,35)#p_TxPa_max=30
p_PrPa_min=runif(1,1.5,4)#p_PrPa_min=3
p_PlPa_r=runif(1,.7,1.5)#p_PlPa_r=1
p_PrPl_rperPr=runif(1,.3,.7)#p_PrPl_rperPr=.5
p_TxPl_min=runif(1,5,14)#p_TxPl_min=10
p_TxPl_max=runif(1,26,32)#p_TxPl_max=30
p_PaPl_r=-exp(runif(1,log(0.001),log(0.045)))#p_PaPl_r=-0.03
p_PlPl_r=exp(runif(1,log(1.8),log(3.5)))#p_PlPl_r=2.5
p_PlPl_K=exp(runif(1,log(1.5),log(3)))#p_PlPl_K=2
p_PlPl_sd=runif(1,.07,.15)#p_PlPl_sd=.1
p_PePa_r=exp(runif(1,log(.06),log(.15)))#p_PePa_r=0.1
p_PaPe=runif(1,1.3,1.9)#p_PaPe=1.5
#ecoindic
p_Pl_sd_r=exp(runif(1,log(.03),log(.07)))#p_Pl_sd_r=0.05
p_T_sd=runif(1,.2,.5)#p_T_sd=0.3
p_Pe_pDet=runif(1,.6,.99)#p_Pe_pDet=.8
p_Pe_pFalseDet=exp(runif(1,log(.001),log(.02)))#p_Pe_pFalseDet=.005

#
# Learning loop
#

for (runi in 2:100){# sample prior
print(runi)  
  
  # simulate edata from prior sample
  
  for (k in 1:dim(edata)[3]){
    for (i in 3:dim(edata)[1]){
      #Pe
      edata[i,5,k]=edata[i-1,3,k]>p_PaPe
      #Pl
      a=((edata[i-1,1,k]>p_TxPl_min)*0.1)+edata[i-1,4,k]*p_PlPl_r*((((1+edata[i-1,2,k])*p_PrPl_rperPr)*(edata[i-1,1,k]>p_TxPl_min))*
                                                                     (1+(T-p_TxPl_min)/(p_TxPl_max-p_TxPl_min))*(edata[i-1,1,k]<p_TxPl_max))+edata[i-1,3,k]*(p_PaPl_r)
      if(a>p_PlPl_K){a=2}
      if(a<0) {a=0}
      edata[i,4,k]=a
      #Pa
      a=(edata[i-1,3,k]==0)+((!edata[i,5,k])+edata[i,5,k]*p_PePa_r)*edata[i-1,3,k]*p_PaPa_rmax*((edata[i-1,1,k]>p_TxPa_min)*(edata[i-1,1,k]<p_TxPa_max))*(edata[i-1,1,k]-p_TxPa_min)/(p_TxPa_max-p_TxPa_min)*(edata[i-2,2,k]>p_PrPa_min)*(edata[i-1,4,k]*p_PlPa_r)
      edata[i,3,k]=rpois(1,a*(a<p_PaPa_K)+p_PaPa_K*(a>=p_PaPa_K))
    }
  }
  
  # Calculate idata probability from edata history
  
  pIndic <- sum(c(pT = sum(log(dnorm(round(idata[,1,]),edata[,1,],p_T_sd))),
                  pH = sum(log(dpois(round(idata[,2,]),edata[,2,]))),
                  pPa = {
                    p1 <- which1 <- which(as.logical(edata[,3,]))
                    p1 <- sum(dbinom(idata[,5,][which1],1,p_Pe_pDet,log=TRUE))
                    p1*length(idata[,5,])/length(which1)
                  },
                  sdPl={sdPl=edata[,4,]*p_Pl_sd_r
                  sdPl[sdPl<=0.1]=.1
                  sum(log(dnorm(x=idata[,4,],mean=edata[,4,],sd=sdPl)))},
                  pPe={
                    p1 <- which1 <- which(as.logical(edata[,5,]))
                    p1 <- sum(dbinom(idata[,5,][which1],1,p_Pe_pDet,log=TRUE))
                    p0 <- which0 <- which(!as.logical(edata[,5,]))
                    p0 <- sum(dbinom(idata[,5,][which0],1,p_Pe_pFalseDet,log=TRUE))
                    p0+p1
                  }
  ))
  
  # calculate posterior
  # posterior = likelihood * prior
  
  prior = log(prod(c(p_PaPa_rmax=dunif(log(p_PaPa_rmax),log(5),log(15))/(log(15)-log(5)),
                     p_PaPa_K=1/(log(25)-log(15)),
                     p_TxPa_min=1/(-10+20),
                     p_TxPa_max=1/(-25+35),
                     p_PrPa_min=1/(-1.5+4),
                     p_PlPa_r=1/(-.7+1.5),
                     p_PrPl_rperPr=1/(-.3+.7),
                     p_TxPl_min=1/(5+14),
                     p_TxPl_max=1/(26+32),
                     p_PaPl_r=1/(-log(0.015)+log(0.045)),
                     p_PlPl_r=1/(-log(1.8)+log(3.5)),
                     p_PlPl_K=1/(-log(1.5)+log(3)),
                     p_PlPl_sd=1/(-.07+.15),
                     p_PePa_r=1/(-log(.06)+log(.15)),
                     p_PaPe=1/(-1.3+1.9),
                     #ecoindic
                     p_Pl_sd_r=1/(-log(.03)+log(.07)),
                     p_T_sd=1/(-.2+.5),
                     p_Pe_pDet=1/(-.6+.99),
                     p_Pe_pFalseDet=1/(-.7+.99))))
  
  # since sampled from uniforms props = constant
  # we use onliy likelihood to sample
  
  posterior = pIndic + prior
  Posterior[runi,] <-   c(runi,p_PaPa_rmax,p_TxPa_min,p_PrPa_min,p_PlPa_r,p_PrPl_rperPr,p_TxPl_min,p_PaPl_r,p_PlPl_r,p_PePa_r,p_PaPe,p_Pl_sd_r,p_T_sd,p_Pe_pDet,p_Pe_pFalseDet,prior,posterior)

  # metropolis move
  p_PaPa_rmax=exp({a={if(as.logical(rbinom(1,1,1/2))) log(p_PaPa_rmax+(15-5)/20) else log(p_PaPa_rmax-(15-5)/20)};{ if (a<log(5)) log(5) else if (a>log(15)) log(15) else a}})#exp(runif(1,log(5),log(15)))
  #p_PaPa_K={a={rbinom(1,1,1/2);if(a==0) a=-1;a}*(log(25)-log(15))/20;if (p_PaPa_K<15) p_PaPa_K=15;if (p_PaPa_K>15) p_PaPa_K=25;p_PaPa_K}#exp(runif(1,log(15),log(25)))
  p_TxPa_min={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(20-10)/20;a=p_TxPa_min+b;if (a<10) a=15;if (a>20) a=20;a}#runif(1,10,20)
  p_TxPa_max={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(35-25)/20;a=p_TxPa_max+b;if (a<25) a=25;if (a>35) a=35;a}#runif(1,25,35)
  p_PrPa_min={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(4-1.5)/20;a=p_PrPa_min+b;if (a<1.5) a=1.5;if (a>4) a=4;a}#runif(1,1.5,4)
  p_PlPa_r={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(1.5-.7)/20;a=p_PlPa_r+b;if (a<.7) a=.7;if (a>1.5) a=1.5;a}#runif(1,.7,1.5)
  p_PrPl_rperPr={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.7-.3)/20;a=p_PrPl_rperPr+b;if (a<.3) a=.3;if (a>.7) a=.7;a}#runif(1,.3,.7)
  p_TxPl_min={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(14-5)/20;a=p_TxPl_min+b;if (a<5) a=5;if (a>14) a=14;a}#runif(1,5,14)
  p_TxPl_max={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(32-26)/20;a=p_TxPl_max+b;if (a<26) a=26;if (a>32) a=32;a}#runif(1,26,32)
  p_PaPl_r=-exp({a={if(as.logical(rbinom(1,1,1/2))) log(-p_PaPl_r-(.01-.045)/20) else log(-p_PaPl_r+(.01-.045)/20)};{ if (a<log(.01)) log(.01) else if (a>log(.045)) log(.045) else a}})#-exp(runif(1,log(0.015),log(0.045)))
  p_PlPl_r=exp({a={if(as.logical(rbinom(1,1,1/2))) log(p_PlPl_r+(3.5-1.8)/20) else log(p_PlPl_r-(3.5-1.8)/20)};{ if (a<log(1.8)) log(1.8) else if (a>log(3.5)) log(3.5) else a}})#exp(runif(1,log(1.8),log(3.5)))
  p_PlPl_K=exp({a={if(as.logical(rbinom(1,1,1/2))) log(p_PlPl_K+(3-1.5)/20) else log(p_PlPl_K-(3-1.5)/20)};{ if (a<log(1.5)) log(1.5) else if (a>log(3)) log(3) else a}})#exp(runif(1,log(1.5),log(3)))
  p_PlPl_sd={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.15-.07)/20;a=p_PlPl_sd+b;if (a<.07) a=.07;if (a>.15) a=.15;a}#runif(1,.07,.15)
  p_PePa_r={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.15-.06)/20;a=p_PePa_r+b;if (a<.06) a=.06;if (a>.15) a=.15;a}#exp(runif(1,log(.06),log(.15)))
  p_PaPe={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(1.9-1.3)/20;a=p_PaPe+b;if (a<1.3) a=1.3;if (a>1.9) a=1.9;a}#runif(1,1.3,1.9)
  #ecoindic
  p_Pl_sd_r=exp({a={if(as.logical(rbinom(1,1,1/2))) log(p_Pl_sd_r+(.07-.03)/20) else log(p_Pl_sd_r-(.07-.03)/20)};{ if (a<log(.03)) log(.07) else if (a>log(.07)) log(.07) else a}})#exp(runif(1,log(.03),log(.07)))
    {b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.7-.03)/20;a=p_Pl_sd_r+b;if (a<.03) a=.03;if (a>.7) a=.7;a}
  p_T_sd={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.5-.2)/20;a=p_T_sd+b;if (a<.2) a=.2;if (a>.5) a=.5;a}#runif(1,.2,.5)
  p_Pe_pDet={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.99-.6)/20;a=p_Pe_pDet+b;if (a<.6) a=.6;if (a>.99) a=.99;a}#runif(1,.6,.99)
  p_Pe_pFalseDet={b={b=rbinom(1,1,1/2);if(b==0) b=-1;b}*(.02-.001)/20;a=p_Pe_pFalseDet+b;if (a<.001) a=.001;if (a>.02) a=.02;a}#exp(runif(1,log(.001),log(.02)))
  
  p_PaPa_rmax=pprime["p_PaPa_rmax"]
  p_TxPa_min=pprime["p_TxPa_min"]
  p_PrPa_min=pprime["p_PrPa_min"]
  p_PlPa_r=pprime["p_PlPa_r"]
  p_PrPl_rperPr=pprime["p_PrPl_rperPr"]
  p_TxPl_min=pprime["p_TxPl_min"]
  p_PaPl_r=pprime["p_PaPl_r"]
  p_PlPl_r=pprime["p_PlPl_r"]
  p_PlPl_K=pprime["p_PlPl_K"]
  p_PlPl_sd=pprime["p_PlPl_sd"]
  p_PePa_r=pprime["p_PePa_r"]
  p_PaPe=pprime["p_PaPe"]
  #ecoindic
  p_Pl_sd_r=pprime["p_Pl_sd_r"]
  p_T_sd=pprime["p_T_sd"]
  p_Pe_pDet=pprime["p_Pe_pDet"]
  p_Pe_pFalseDet=pprime["p_Pe_pFalseDet"]
  
  }

Posterior <-  data.frame(p_PaPa_rmax=0,p_TxPa_min=0,p_PrPa_min=0,p_PlPa_r=0,p_PrPl_rperPr=0,p_TxPl_min=0,
                         p_PaPl_r=0,p_PlPl_r=0,p_PePa_r=0,p_PaPe=0,p_Pl_sd_r=0,p_T_sd=0,p_Pe_pDet=0,p_Pe_pFalseDet=0,prior=0,posterior=0)



# 

setClass("parDisFun",
         contains = c("function"),
         slots= c(lengthx="integer",lengthp="integer",xname="character",pnames="character"),
#         prototype = list(function(x,p) rpois(x*p[[1]]),lengthx=1:1,lengthp=1:1,xname="Tx",pnames="TxPa"),
)
# parDisFun are the parametrized distribution functions ]linking ecosystem variables 

setClass("parFunList",
         contains = "list",
         prototype = list(new("parDisFun",function(x,p) rpois(1,p[[1]]),lengthx=1:1,lengthp=1:1,xname="Tx",pnames="TxPa"),
                          new("parDisFun",function(x,p) rnorm(1,p[["PlPl"]][1]+x[["Pl"]]*p[["PlPl"]][2],p[["PlPl"]][3]),lengthx=1:1,lengthp=1:1,xname=c("Pl"),pnames=c("PlPl"))
         ),
         validity = function(object){
           if (any(lapply(object,class)!="parDisFun")) stop("trying to create a funList where not all list components are parDisFun")}
)
setClass("EcoModel",
         contains="Data",
         slots=c(funs="parFunList",ecoLink="matrix",delays="matrix")
)

setClass("prior",
         contains = "function",
         slots = c(priorParam="numeric")
)

pr1=new("prior",runif,priorParam=c(0,10))
sample(pr1)

pr1=new("prior",rnorm,priorParam=c(0,1))
pr2=new("prior",rnorm,priorParam=c(10,5))

setClass("prior",
         contains = "list",
         slot = "FUN",
)


setClass("priorList",
         contains = "list",
         validity = function(object) {
           if (any(lapply(object,class)!="prior")) stop("priorList@funs did not contain list of prior")
         }
)

priorP <- new("priorList",lapply(p,function(p_i) {new("prior",list(PrMin=p_i*0.66,PrMax=p_i*1.5),FUN=runif)}))
priorP <- lapply(p,function(p_i) list(PrMin=p_i*0.66,PrMax=p_i*1.5,FUN=runif))
mn=priorP$PaPl$PrMin
mx=priorP$PaPl$PrMax
priorP$PaPl$PrMin=mx
priorP$PaPl$PrMax=mn

names(priorP)=names(p)
priorP$Papa
x=priorP

setMethod("sample",
          signature = "list",
          definition = function(x){
            l <- lapply(x,function(x) x[[1]])
            for (i in 1:length(l)){for (j in 1:length(l[[i]])) l[[i]][[j]]=x[[i]]$FUN(1,x[[i]][[1]][[j]],x[[i]][[2]][[j]]) }
            l
          })
p1=sample(priorP)

sample(pr2)

x=c(Pa=1.5,Tx=20,Pl=0.5,Pr=2.5,Pe=0)
p=list(PaPa=c(rmax=10,K=20),TxPa=c(min=15,max=30),PrPa=c(min=3),PlPa=c(r=1),PrPl=c(rperPr=.5),TxPl=c(min=10,max=30),PaPl=c(r=-0.03),PlPl=c(r=2.5,K=2,sd=.1),PePa=c(r=0.1),PaPe=c(thr=1.5))
p=new("param",p)


Pafun <- new("parDisFun",function(x,p){
    a=(x["Pa"]==0)+((!x["Pe"])+x["Pe"]*p[["PePa"]]["r"])*x["Pa"]*p[["PaPa"]]["rmax"]*((x["Tx"]>p[["TxPa"]]["min"])*(x["Tx"]<p[["TxPa"]]["max"]))*(x["Tx"]-p[["TxPa"]]["min"])/(p[["TxPa"]]["max"]-p[["TxPa"]]["min"])*(x["Pr"]>p[["PrPa"]]["min"])*(x["Pl"]*p[["PlPa"]])
  rpois(1,a*(a<p[["PaPa"]]["K"])+p[["PaPa"]]["K"]*(a>=p[["PaPa"]]["K"]))
                         },lengthx=as.integer(5),lengthp=as.integer(5),xname=c("Tx","Pr","Pa","Pl","Pe"),pnames=c("TxPa","PrPa","PaPa","PlPa","PePa")) 

#
#Plfun <- new("parDisFun",function(x,p){
#  a = ((x["Tx"]>p[["TxPl"]]["min"])*0.1)+x["Pl"]*p[["PlPl"]]["r"]*(((x["Pr"]>p[["PrPl"]])*(x["Tx"]>p[["TxPl"]]["min"]))*
#      (1+(T-p[["TxPl"]]["min"])/(p[["TxPl"]]["max"]-p[["TxPl"]]["min"]))*(x["Tx"]<p[["TxPl"]]["max"]))+x["Pa"]*(p[["PaPl"]]["r"])
#rnorm(1,a,p[["PlPl"]]["sd"]*(a>0))
Plfun <- new("parDisFun",function(x,p){
  a=((x["Tx"]>p[["TxPl"]]["min"])*0.1)+x["Pl"]*p[["PlPl"]]["r"]*((((1+x["Pr"])*p[["PrPl"]])*(x["Tx"]>p[["TxPl"]]["min"]))*
                                                                 (1+(T-p[["TxPl"]]["min"])/(p[["TxPl"]]["max"]-p[["TxPl"]]["min"]))*(x["Tx"]<p[["TxPl"]]["max"]))+x["Pa"]*(p[["PaPl"]]["r"])
  if(a>p[["PlPl"]]["K"]) {a=2}
  if(a<0) {a=0}
  a
},lengthx=as.integer(4),lengthp=as.integer(4),xname=c("Tx","Pr","Pa","Pl"),pnames=c("TxPl","PrPl","PaPl","PlPl")) 

Pefun <- new("parDisFun",function(x,p){
#  runif(1,0,x["Pa"])>p[["PaPe"]]},lengthx=1:1,lengthp=1:1,xname=c("Pa"),pnames=c("PaPe")
  x["Pa"]>p[["PaPe"]]},lengthx=1:1,lengthp=1:1,xname=c("Pa"),pnames=c("PaPe")
)

iTfun <- new("parDisFun",function(x,p){
  rnorm(1,x["iT"],p["sd"])
},lengthx=1:1,lengthp=1:1,xname=c("Tx"),pnames=c("sd")) 

modlfun =new("parFunList",list(Pe=Pefun,Pa=Pafun,Pl=Plfun))

modl=new("EcoModel",edata,funs=modlfun,ecoLink=ecoLink,delays=ecoLinkTime)
object=modl
modl@funs
names(modl@funs)
dimnames(modl@.Data[,,1])
p
modl

setClass("priorize",
         signature = c("EcoModel","prior")
         definition = function(object){
           
         }
         )

setGeneric(name="simulate",
  def= function(object,p) { return(standardGeneric("simulate"))}
)
object=modl
setClass("simulate",
         signature = c(object="EcoModel",p="list"),
         definition=function(object,p) object*unlist(p))



setMethod("simulate",
         signature=c(object="EcoModel",p="list"),
         definition = function(object,p){
           for (pop in 1:dim(object)[3]){
             for (per in (max(object@delays[,names(object@funs)])+1):dim(object)[1]){
               for (Funi in 1:length(object@funs)){
                 Fun <- object@funs[[Funi]]
                 nameFuni <- names(object@funs)[Funi]
                 #        tmp=
                 #        if (per>tmp) if (any(is.na(object@.Data[(per-tmp):(per-1),Fun@xname,pop]))) stop("data missing for simulation of data coordinates: per=",(per-tmp):(per-1)," dep.var=",Fun@xname," population=",pop) else {
                 tmp2 <- unlist(lapply(which(object@ecoLink[,nameFuni]),FUN=function(i) object@.Data[per-object@delays[i,nameFuni],i,pop]))
                 #          dati <- object@.Data[,,pop]
                 object@.Data[per,nameFuni,pop] <- object@funs[[Funi]](tmp2,p)
                 #print(paste(nameFuni,object@funs[[Funi]](tmp2,p)))
               }
             }
           }
         object
         }
         )


setMethod("probablity",
          signature=c(object="EcoModel",p=list),
          for (pop in 1:dim(object)[3]){
            for (per in (max(object@delays[,names(object@funs)])+1):dim(object)[1]){
              for (Funi in 1:length(object@funs)){
                Fun <- object@funs[[Funi]]
                nameFuni <- names(object@funs)[Funi]
                #        tmp=
                #        if (per>tmp) if (any(is.na(object@.Data[(per-tmp):(per-1),Fun@xname,pop]))) stop("data missing for simulation of data coordinates: per=",(per-tmp):(per-1)," dep.var=",Fun@xname," population=",pop) else {
                tmp2 <- unlist(lapply(which(object@ecoLink[,nameFuni]),FUN=function(i) object@.Data[per-object@delays[i,nameFuni],i,pop]))
                #          dati <- object@.Data[,,pop]
                object@.Data[per,nameFuni,pop] <- object@funs[[Funi]](tmp2,p)
                #print(paste(nameFuni,object@funs[[Funi]](tmp2,p)))
              }
            }
          }
          )

object=modl

sim <- list()
plist=list()

for (i in 2:100){
  plist <- sample(priorP) 
  sim=simulate(modl,plist)
  plsim=list(p=plist,sim=sim)
  save(plsim,file = paste("pAndEcosim",i,".RData",sep=""))
  }


plsim[[1]]
simul <- function(object,p){
  lapply
  Fun <- object@funs[[Funi]]
  nameFuni <- names(object@funs)[Funi]
  #        tmp=
  #        if (per>tmp) if (any(is.na(object@.Data[(per-tmp):(per-1),Fun@xname,pop]))) stop("data missing for simulation of data coordinates: per=",(per-tmp):(per-1)," dep.var=",Fun@xname," population=",pop) else {
  tmp2 <- unlist(lapply(which(object@ecoLink[,nameFuni]),FUN=function(i) object@.Data[per-object@delays[i,nameFuni],i,pop]))
  #          dati <- object@.Data[,,pop]
  object@.Data[per,nameFuni,pop] <- object@funs[[Funi]](tmp2,p)
  
  
simulate
param
class(p)
a= array(1:24,dim=c(2,3,4),dimnames=list(1:2,1:3,1:4))
for ()
fun= object@funs[[1]]
class(fun)
a=NULL;i=0
for (fun in object@funs){i=i+1
  print(fun)
  print(paste("LA CLASSE C'EST", class(fun)))
  print(paste("Le name C'EST", names(fun)))
}
names(object@funs)
Pa=1.5;T=20;Pl=0.5;H=92;Pe=0
p_PaPa=c(rmax=10,K=20);p_TPa=c(min=15,max=30);p_HPa=c(min=90);p_PlPa=c(r=1)
p_HPl=c(min=50);p_TPl=c(min=10,max=30);p_PaPl=c(r=-0.03);p_PlPl=c(r=1.05,sd=.1)
p_PePa=c(r=.1)
# 
a=((!Pe)+Pe*p_PePa)*Pa*p_PaPa[1]*((T>p_TPa[1])*(T<p_TPa[2]))*(T-p_TPa[1])/(p_TPa[2]-p_TPa[1])*(H>p_HPa[1])*(Pl*p_PlPa)
Pa = rpois(1,a*(a<p_PaPa[2])+p_PaPa[2]*(a>=p_PaPa[2]))
# Pa si pas de pesticide pas 
Pl = ((T>p_TPl[1])*0.1)+Pl*p_PlPl[1]*(((H>p_HPl)*(T>p_TPl[1]))*
      (1+(T-p_TPl[1])/(p_TPl[2]-p_TPl[1]))*(T<p_TPl[2]))+Pa*(p_PaPl)
Pl = rnorm(1,Pl*(Pl>0),p_PlPl[2]*Pl*(Pl>0))
Pe = (Pa>=p_PePa)
Pe;Pa;Pl


setClass("model",
         slot=c("parFunList","ecoLink",""))



Model <- function(param,linkH,linkHI,paramx0){
  
}


# Simulate climate data
Njours=900
nbjour_mois <- c(31,28,31,30,31,30,31,31,30,31,30,31)
cumul_jour <- NULL
for (i in 1:12) cumul_jour[i] <- sum(nbjour_mois[1:i])
precmoyparmois <- c(54 ,46 ,50 ,44 ,58 ,56 ,53 ,51 ,56 ,57 ,58, 54)
prec <- rpois(Njours,rep(precmoyparmois/nbjour_mois,nbjour_mois ))
precparmois <- sum(prec[1:cumul_jour[1]])
for (mois in 2:12) precparmois[mois] <- sum(prec[cumul_jour[mois-1]:cumul_jour[mois]])

tmeanmoyparmois <- c(3.3, 4.2, 7.8, 10.8, 14.3, 17.5, 19.4, 19.1, 16.4, 11.6, 7.2, 4.2 )
tmeanmoyjours <- rep(tmeanmoyparmois,nbjour_mois)
tmean <- rnorm(1,tmeanmoyjours,3) 
for (j in 2:Njours) tmean[j] <- (tmean[j-1]+ rnorm(1,0,3) + tmeanmoyjours[{a=j%%365;if(a==0) a=365 ; a}])/2 
tmeanparmois <- mean(tmean[1:cumul_jour[1]])
for (mois in 2:12) tmeanparmois[mois] <- mean(tmean[cumul_jour[mois-1]:cumul_jour[mois]])

plantgrowth <- NULL

for (year in 0:2) {
planted <- FALSE
joursseuil=5
tempseuil=9
precseuil=5
joursseuilP=5
day=1
while (!planted) {planted <- ((sum(prec[(year*365+day-{if (day<=joursseuilP) day-1 else joursseuilP}):(year*365+day)])>precseuil)&(mean(tmean[(year*365):(year*365+day)][(day-{if (day<=joursseuil) day-1 else joursseuil}):day])>tempseuil)&(day>joursseuil)) 
		  day=day+1}
growth=0 while (growth<90) growth = growth +  (sum(prec[(year*365+day-{if (day<=joursseuilP) day-1 else joursseuilP}):(year*365+day)])>precseuil)* tmean(year*365+day)/30
		  }
X2 <- sin((1:Njours)*2*pi/365)*12 + rnorm(Njours,14,5)
tX3<-tX2 <- 1:Njours
plot(tX2,X2)
X3=NULL
X3[1]=0
for(i in (2:Njours)) {
if (i>10) X3[i]=(X3[i-1]+(X2[i-10] + abs(rnorm(1,0,1)))/100) else X3[i]=X3[i-1] + abs(rnorm(1,0,1))/100
if ((i%%365)>100) X3[i]=0
}
par(mfrow=c(1,2))
plot(tX3,X3)
plot(tX2,X2)
a <- sapply(1:5,FUN=function(x) x+1)

F=function(x) 
{
if (x>10) (X3[x-1]+(X2[x-10] + rnorm(1,0,4))/100) 
else (X3[x-1] + rnorm(1,0,4)/100)
}

<- as.data.frame(matrix(runif(800,1,40),nrow=
### Start PlaNet
  
  # continuous-time series


cSTF <- setClass("cSTF",
                                       slots=c(times=c("POSIXct"),coord="matrix"),
                                       validity = function(object){
                                         if (length(object@times)!=nrow(object@coord)) stop("length of @times vector differs from number of rows in @coord matrix")
                                       }
)

setMethod("length",
          signature="cSTF",
          definition=function(x){length(x@times)})

setMethod("dim",
          signature="cSTF",
          definition=function(x){dim(x@coord)})


intCSTS <- setClass("intCSTS",
                       contains="cSTF",
                       slots=c(values="integer"),
                       validity = function(object){
                         if (length(object@values)!=length(object)) stop("length of cSTF vector differs from length of @values list")
                         }
                       )
exiCSTS=new("intCSTS",new("cSTF",times=c(as.POSIXct("2018-10-13 14:28:22",tz="America/Bogota"),as.POSIXct("2018-10-13 14:28:22",tz="America/Bogota")),coord=as.matrix(data.frame(X=c(1.2,1.4),Y=c(3.6,-1)))),values=as.integer(c(4,5)))

charCSTS <- setClass("charCSTS",
                     contains="cSTF",
                     slots=c(values="character"),
                     validity = function(object){
                       if (length(object@values)!=length(object)) stop("length of cSTF vector differs from length of @values list")
                     }
)

excCSTS=new("charCSTS",new("cSTF",times=c(as.POSIXct("2018-10-13 14:28:22 America/Bogota"),as.POSIXct("2018-10-13 14:28:22 America/Bogota")),coord=as.matrix(data.frame(X=c(1.2,1.4),Y=c(3.6,-1)))),values=c("e","r"))

logicCSTS <- setClass("logicCSTS",
                      contains="cSTF",
                      slots=c(values="logical"),
                      validity = function(object){
                        if (length(object@values)!=length(object)) stop("length of cSTF vector differs from length of @values list")
                      }
)

exlCSTS=new("logicCSTS",new("cSTF",times=c(as.POSIXct("2018-10-13 14:28:22 America/Bogota"),as.POSIXct("2018-10-13 14:28:22 America/Bogota")),coord=as.matrix(data.frame(X=c(1.2,1.4),Y=c(3.6,-1)))),values=c(T,F))

numCSTS <- setClass("numCSTS",
                    contains="cSTF",
                    slots=c(values="numeric"),
                    validity = function(object){
                      if (length(object@values)!=length(object)) stop("length of cSTF vector differs from length of @values list")
                    }
)

exnCSTS=new("numCSTS",new("cSTF",times=c(as.POSIXct("2018-10-13 14:28:22 America/Bogota"),as.POSIXct("2018-10-13 14:28:22 America/Bogota")),coord=as.matrix(data.frame(X=c(1.2,1.4),Y=c(3.6,-1)))),values=c(0.88,3.54))

# CSTS : list of Continuous-SpatioTemporal Series

CSTS <- setClass("CSTS",
		      contains="list",
		      validity = function(object){
		       # if (!all(lapply(object,class)%in%c("numCSTS","intCSTS","logicCSTS","charCSTS"))) stop("trying to contrust a TimeSeries object with a list conaining objects other than numCSTS, or intCSTS, or charCSTS or logicCSTS")
		        if (!all(sapply(lapply(object,function(x){colnames(x@coord)}),FUN=identical,colnames(object[[1]]@coord)))) stop("the coordinates colnames differ among time series in CSTS object contruction")
		      }
)

CSTS <- function(X,name=NA){
  #note if X is a matrix, col 1 contains values, col 2 times, and col 3.. coordinates
  #     if X is a list, fisrt element contains values, second elt contains times, and other elements coordinates
  if (any(is.null(names(X)))) stop("requires names for the list as variables names in CSTS constructor")
  if (class(X)== "list") {
    for (i in 1:length(X)) {
      if (!(class(X[[i]])%in%c("numCSTS","intCSTS","logicCSTS","charCSTS"))) {
        X[[i]] = switch(class(X[[i]]),
                        matrix = new("numCSTS",new("cSTF",times=X[i,2],coord=X[i,3:ncol(X)]),values=X[i,1]),
                        list = switch(class(X[[i]][[1]]),
                                      character= new("charCSTS",new("cSTF",times=X[[i]][[2]],coord=(X[[i]][[3]])),values=X[[i]][[1]]),
                                      numeric= new("numCSTS",new("cSTF",times=X[[i]][[2]],coord=(X[[i]][[3]])),values=X[[i]][[1]]),
                                      integer=new("intCSTS",new("cSTF",times=X[[i]][[2]],coord=(X[[i]][[3]])),values=X[[i]][[1]]),
                                      logical= new("logicCSTS",new("cSTF",times=X[[i]][[2]],coord=(X[[i]][[3]])),values=X[[i]][[1]]))
                        )
      }
    }
    new("CSTS",X)
  }
else stop("needs a list as argument")
}


c(as.POSIXct("2018-10-12 20:45:12"),as.POSIXct("2018-10-11 20:45:12"))

x=object=exCSTS <- CSTS(list(landtype=list(values=c("farm","road","wood",'wood'),times=c(as.POSIXct("2018-10-13 14:28:22",tz="America/Bogota"),as.POSIXct("2018-10-12 15:28:24",tz="America/Bogota"),as.POSIXct("2018-10-10 15:28:24",tz="America/Bogota"),as.POSIXct("2018-10-13 14:25:22",tz="America/Bogota")),coord=as.matrix(data.frame(x1=c(1,2,.4,.5),x2=c(2,3,3.4,.5)))),
                    tmean=list(c(10,15,7,8),c(as.POSIXct("2018-10-13 14:28:22",tz="America/Bogota"),as.POSIXct("2018-10-12 15:28:24",tz="America/Bogota"),as.POSIXct("2018-10-10 15:28:24",tz="America/Bogota"),as.POSIXct("2018-10-13 14:25:22",tz="America/Bogota")),as.matrix(data.frame(x1=c(1,2,.4,.5),x2=c(2,3,3.4,.5)))),
                    present=list(c(T,F,F,T),c(as.POSIXct("2018-10-13 14:28:22",tz="America/Bogota"),as.POSIXct("2018-10-12 15:28:24",tz="America/Bogota"),as.POSIXct("2018-10-10 15:28:24",tz="America/Bogota"),as.POSIXct("2018-10-13 14:25:22",tz="America/Bogota")),as.matrix(data.frame(x1=c(1,2,.4,.5),x2=c(2,3,3.4,.5))))))

class(exCSTS)

min(exCSTS[[1]]@times)

setMethod(min,
          signature="CSTS",
          definition=function(x){
            min(as.POSIXct(unlist((lapply(x,FUN=function(xi) as.character(min(xi@times)))))))
          })
min(exCSTS)
setMethod(max,
          signature="CSTS",
          definition=function(x){
            max(as.POSIXct(unlist((lapply(x,FUN=function(xi) as.character(max(xi@times)))))))
          })
max(exCSTS)
setMethod("as.data.frame",
          signature="CSTS",
          definition=function(x){
            df<-as.data.frame(lapply(x,FUN=function(x) x@values))
            vapply(x, function, FUN.VALUE = type, ...)
          })

setMethod(discretize,
          signature="CSTS",
          definition = function(x,unit="day",tZ=Sys.timezone()){
            start = min(x)
            finish= max(x)
            starting=as.POSIXct(switch(unit,
                            second = start,
                            minute= as.Date(start,tz=tZ)+period(hour(start),"hours")+period(minute(start),"minute"),
                            hour= as.Date(start,tz=tZ)+period(hour(start,tz=tZ),"hours"),
                            day=as.Date(start,tz=tZ)+period(1,"hours")-period(1,"hour"),
                            month=as.Date(start,tz=tZ)-day(start,tz=tZ)+period(1,'days'),
                            year=((as.Date(start,tz=tZ)-day(start,tz=tZ)+period(1,"days"))-period(month(start)-1,"month"))
                         ),tz=tZ)
            finished=as.POSIXct(switch(unit,
                            second=finish+period(1,"second"),
                            minute=as.Date(finish)+period(hour(finish),"hours")+period(minute(finish)+1,"minute"),
                            hour=as.Date(finish)+period(hour(finish)+1,"hours"),
                            day=as.Date(finish)+1+period(1,"hours")-period(1,"hours"),
                            month=period(1,"month")+as.Date(finish)-day(finish)+period(1,'days'),
                            year=period(1,"year")+(as.Date(finish)-day(finish)+period(1,"days"))-period(month(finish)-1,"month")
                            ),tz=tZ)
            len <- as.integer(as.period(finished-starting,unit)/period(1,unit))
            tf <- new("timeFrame",starting=starting,finished=as.POSIXct(finished),period=period(1,unit), length=len)
            cutf <- cut(tf)
            listOflist <- list()
            listOflist<-lapply(cutf[-length(cutf)],function(x) {x=list(x);names(x)="times";x})
            dataF <- as.data.frame(matrix(NA,ncol=length(x),nrow=len))
            colnames(dataF) <- names(x)
            for (vari in names(x)) {
              listOflist[[vari]]<-NULL
              for (element in 1:length(x[[vari]])){
                listOflist[[sum(x[[vari]]@times[element]>cutf)]][[vari]]<-append(listOflist[[sum(x[[vari]]@times[element]>cutf)]][[vari]],x[[vari]]@values[element])
                   }
            }
            columnsOfDataFrame <- NULL
              for (vari in names(x)){
                columnsOfDataFrame <- append(columnsOfDataFrame,switch(class(x[[vari]]),
                     intCSTS=c(paste(vari,".n",sep=""),paste(vari,".rep",sep="")),
                     charCSTS=paste(vari,levels(as.factor(x[[vari]]@values)),sep="."),
                     logicCSTS=c(paste(vari,".n",sep=""),paste(vari,".k",sep="")),
                     numCSTS=vari))
                }
            dataF <- data.frame(matrix("",ncol=length(columnsOfDataFrame)))
            colnames(dataF)<- columnsOfDataFrame
            for (time in cutf){
              {
                for(vari in names(x)){
                  for (col in grep(vari,columnsOfDataFrame,value=TRUE)){
                    
                  }
                }
                switch(class(x[[vari]]),
                       intCSTS=,
                       charCSTS={for (state in levels(x[[vari]]@values))
                         
                         },
                       ){
                  
                }
                for (time in cutf){
                  dataF[dataF$times==time,]
                }
              }
              dataF[,vari]
            }
            for (col in colnames(dataF))
            for (elt in 1:length(listOfList)){
              df
            }
          }

              
              referenceTime <- referenceTime - start
            {(referenceTime-start)/period referenceTime <- referenceTime }referenceTime <- referenceTime - start
            timepoints  <- 
            if (referenceTime>=start) starting=referenceTime-ceiling(abs(referenceTime-start))
            if (referenceTime<start) starting=referenceTime+floor(abs(referenceTime-start))
            finished=starting+ceiling(finish-start)
            
          })

x=exCSTS


# discetre-time series #
#######################

library(lubridate)

timeFrame <- setClass("timeFrame",
                          slots=c(period="Period", length="integer"),
                          prototype=prototype(period=period(num = 1,unit = "months"),length=as.integer(6)))


refTimeFrame <- setClass("refTimeFrame",
                      contains="timeFrame",
                      slots = c(starting=c("POSIXct"),finished=c("POSIXct")),
                      prototype = prototype(new("timeFrame",period=period(num = 1,units="seconds"),length=as.integer(6)),starting=as.POSIXct("2011-06-10 08:06:35"),finished=as.POSIXct("2011-06-10 08:06:41")),
                      validity = function(object){
                        if ((object@starting)+object@period*(object@length)!=object@finished) {
                          stop(paste("incompatibility between starting, finished and period values","\n (object@starting)+object@period*(object@length) =",(object@starting)+object@period*(object@length),"\n object@finished =",object@finished))
                          }
                        }
                      )

object=tf6s=new("refTimeFrame",starting=as.POSIXct("2011-06-10 08:06:35"),finished=as.POSIXct("2011-06-10 08:06:41"),period=period(num = 1,units="seconds"),length=as.integer(6))
object=tf6s=new("refTimeFrame",starting=as.POSIXct("2012-06-10 08:06:35"),finished=as.POSIXct("2017-06-10 08:06:35"),period=period(num = 1,units="years"),length=as.integer(5))
object=tf6y=new("refTimeFrame",starting=as.POSIXct("2011-01-01"),finished=as.POSIXct("2017-01-01"),period=period(num = 1,units="years"),length=as.integer(6))

object=new("refTimeFrame",starting=ymd_hms("2011-06-10-08-06-35"),finished=ymd_hms("2011-06-19-08-06-35"),period=period(1,"day"),length=as.integer(9))

setMethod(length,
  signature = "timeFrame",
  definition = function(x){x@length})
length(object)

setMethod(range,
          signature = "refTimeFrame",
          definition = function(x){c(starting=x@starting,finished=x@finished)})
range(object)

setMethod("cut",
          signature="refTimeFrame",
          definition=function(x){
            cuts <- range(x)[1]
            for (i in 1:length(x)){
              cuts = append(cuts,cuts[1]+x@period*i)
            }
          cuts})
cut(object)

# discrete time series DTS
setwd("/home/dupas/PlaNet/")
data <- read.table("TrainingDataSet_Maize.txt")
head(data)

coln <- colnames(data)
coln <- levels(as.factor(unlist(lapply(strsplit(colnames(data),"\\_"),function(x){x[[1]]}))))
?strsplit  
  
dataList <- setClass("dataList",
                     contains="list",
                     validity = function(object){
                       if (!all(lapply(object,class)=="data.frame")) stop("dataList constructor receceived something else than a list of data frame")
                       if (!all(lapply(object,colnames)==colnames(object)[[1]])) stop("colnames of data.frame in dataList should be identicals")
                     })

dataList <- function(x,timeByCol=TRUE,sep="_",listColTag=c("yearHarvest","NUMD","IRR"),timeByRowNColInfo=NULL,connectivity=list(type="temporal",tempVar="yearHarvest",labelVar="NUMD",connectVar="yieldAnomaly",connectRow=9)){
  # timeByRowNcolInfo = list(cols="yearHarvest",colsBy="year",rows="_x",rowsBy="month")
  if (class(x)=="list") return(new("dataList",x))
  if ((class(x)=="data.frame")&(!timeByCol)) return(new("dataList",list(x)))
  if ((class(x)=="data.frame")&(!is.null(timeByRowNColInfo))){
    
  }
  if ((class(x)=="data.frame")&(timeByCol)){
    dataL=list()
    colnlist = strsplit(colnames(x),sep)
    coln <- levels(as.factor(unlist(lapply(strsplit(colnames(x),sep),function(x){x[[1]]}))))
    lastElt = lapply(colnlist,FUN=function(x) x[[length(x)]])
    options(warn=-1)
    times=as.numeric(levels(as.factor(unlist(lastElt[which(!is.na(as.numeric(lastElt)))]))))
    options(warn=0)
    times=times[order(times)]
    for (tag in 1:nrow(x)){
      df0=x[tag,]
      dfn=data.frame(matrix(ncol=length(coln),nrow=length(times),dimnames=list(times,coln)))
      dfn[,coln[which(coln%in%listColTag)]] <- df0[,coln[which(coln%in%listColTag)]]
      for (i in times){
        for (j in coln[!(coln%in%listColTag)]){
          if (paste(j,i,sep="_")%in%colnames(x)) dfn[i,j]=df0[,paste(j,i,sep="_")]
          if (j%in%colnames(x)) dfn[i,j]=df0[i,j]
        }
      }
      dataL[[paste(listColTag,x[tag,listColTag],sep="",collapse="_")]]=dfn
    }
    if (connectivity$typ=="temporal"){
      for (i in 1:length(dataL)){
        previousValue = unlist(lapply(dataL,function(dal){
          if ((dal[1,connectivity$tempVar] == dataL[[i]][1,connectivity$tempVar]-1)&(dal[1,connectivity$labelVar] == dataL[[i]][1,connectivity$labelVar])) dal[connectivity$connectRow,connectivity$connectVar] else NULL}))
        if (is.null(previousValue)) dataL[[i]][,paste("past",connectivity$connectVar,sep="_")]<-NA else dataL[[i]][,paste("past",connectivity$connectVar,sep="_")]<-previousValue
        fivePreviousValue=NULL
        ePreviousValue=append(fivePreviousValue,pv) 
        if (is.null(previousValue)) dataL[[i]][,paste("past",connectivity$connectVar,sep="_")]<-NA else dataL[[i]][,paste("past",connectivity$connectVar,sep="_")]<-previousValue
        }
      }
    }
    return(new("dataList",dataL))
  }
}

dl <- dataList(data,timeByCol=TRUE,sep="_",listColTag=c("yearHarvest","NUMD","IRR"))
save(dl,file = "data.list.RData")
load(file = "data.list.RData")
df <- array(unlist(dl),dim=list(dim(dl[[1]])[1],dim(dl[[1]])[2],length(dl)),dimnames = list(c("JAN","FEB","MAR","APR","MAI","JUN","JUL","AUG","SEP","NOV"),colnames(dl[[1]]),names(dl)))
df[,,1]
?save
save(df,file="yield.data.RData")


x=data.frame(X1=c("a","b","c","d","e","f"),X2=c(1,2,3,4,5,6))
bycol=FALSE
obj=dataList(x=data.frame(X1=c("a","b","c","d","e","f"),X2=c(1,2,3,4,5,6)),timeByCol=FALSE)

setMethod("variable.names",
          signature="dataList",
          definition=function(object){colnames(object[[1]])})

setMethod("nrow",
          signature="dataList",
          definition=function(x) nrow(x[[1]]))

setMethod("dim",          
          signature="dataList",
          definition=function(x) dim(x[[1]]))

variable.names(dl)
nrow(dl)
dim(dl)
length(dl)

dl_TimeFramed <-new("refTimeFrame",starting=as.POSIXct("2018-01-01"),finished=as.POSIXct("2018-11-01"),period=period(36,"minutes")+period(9,"hours")+period(30,"days"),length=as.integer(10))

DTS <- setClass("DTS",
                contains = "refTimeFrame",
                slots=c(variables="character",data="dataList"),
                validity = function(object){
                  if (!all(object@variables%in%variable.names(object@data))) stop("not all variables slot are in data slot colnames")
                  if (length(object)!=nrow(object@data)) stop("timeFrame and number of row in data are different")
                  }
                )

object=dTS <- new("DTS",tf6y,variables=c("X1","X2"),data=dataList(data.frame(X1=c("a","b","c","d","e","f"),X2=c(1,2,3,4,5,6))))
object=dTS <- new("DTS",tf6y,variables=c("X1","X2"),data=data.frame(X1=c("a","b","c","d","e","f"),X2=c(1,2,3,4,5,6)))
object=dTS <- new("DTS",dl_TimeFramed,variables=variable.names(dl),data=dl)

listDTS <- setClass("listDTS",
                    contains="list",
                    )


setMethod("names",
          signature="DTS",
          definition=function(x){
            x@variables
          })

setMethod("as.data.frame",
          signature="DTS",
          definition=function(x){
            cbind(x@data,data.frame(times=cut(x)[-length(cut(x))]))
          })

dfd=as.data.frame(dTS)

setMethod("names",signature="DTS",definition = function(x){x@variables})

names(dTS)

object=eTS <- new("DTS",tf6y,variables=c("E1","E2"),data=data.frame(E1=c(2.45,3.5,4.,6.8,7.0,6.9),E2=as.integer(c(5,5,3,12,67,0))))
names(eTS)
dfe=as.data.frame(eTS)

setGeneric(name = "as.DTS",def = function(x,..){standardGeneric("as.DTS")})

setMethod("as.DTS",
          signature="data.frame",
          definition=function(x,per=NULL,option="bycol"){
            t<-x[,which(lapply(as.list(x),FUN = function(x) class(x)[1])=="POSIXct")]
            df<- x[,-which(lapply(as.list(x),FUN = function(x) class(x)[1])=="POSIXct")]
            starting=min(t)
            t=t[order(t)]
            if(is.null(per)) per = lapply(1:(length(t)-1),function(i){t[i+1]-t[i]})[[which(unlist(lapply(1:(length(t)-1),function(i){t[i+1]-t[i]})) == min(unlist(lapply(1:(length(t)-1),function(i){t[i+1]-t[i]}))))[1]]]
            len = as.integer((max(t)-starting)/as.numeric((per))+1)
            if (per%in%c(difftime("2017-01-01","2016-01-01"),difftime("2018-01-01","2017-01-01"))) per=period(1,"year") else per = as.period(per,unit=unit)
            finished=max(t)+per
            new("DTS",new("timeFrame",starting=starting,finished=finished,period=per,length=len),data=df)
          }
          )



## Ecosystem model structure
# varFunction a character matrix linking each variable of the ecosystem containing XP functions
# 

extendedMatrix <- setClass("extendedMatrix",
                           slots=list)

# XPfun are all the functions with parameters names "x" and "p" where x is a variable , p is a parameter 

xpFun <- setClass("xpFun",
                  contains="function",
                  validity = function(object){if (any(formalArgs(object)!=c("x","p"))) stop("XPFunction was created without 'p' and 'x' as arguments")}
)

proP <- function(x=0,p=0){x*p}
propXP <- new("xpFun",proP(x=0,p=c(p1=0,p2=1)))
propXP(2,3)

pnorM <- new("xpFun",function(x=0,p=c(meaN=2,sigmA=1)){pnorm(x,p[1],p[2])})

pnorM(x=3,p=c(neaN=3,sigmA=4))

# varFunctions is a matrix of functions for the ecosystem graph
varFunctions <- setClass("varFunctions",
                                  contains="list",
                                  validity=function(object){
                                    if (any(lapply(object,class)!="list")) stop("varFunctions should be a list of list")
                                    if (any(lapply(object,length)!=length(object))) stop("varFunctions should be squared")
                                    for (subobject in object) {if (any(lapply(subobject,class)!="xpFun")) stop("varFunctions should be squared list of xpFUN")}
                                  })

vF<-new("varFunctions",list(list(new("xpFun",function(x,p){p*x}),new("xpFun",function(x,p){p[1]*x+p[2]*x^2})),list(new("xpFun",function(x,p){0}),new("xpFun",function(x,p){1}))))
object=vF
vF

setMethod("dim",
          signature = "varFunctions",
          definition = function(x){c(length(x[[1]]),length(x))}
          )
dim(vF)

logimat <- setClass("logimat",
                    contains="matrix",
                    validity=function(object){if (any(lapply(object,class)!="logical")) stop("logimat lentgh should equal its @dim product")}
)

object=G <- new("logimat",matrix(c(T,F,F,T),nr=2))

paramVecMat <- setClass("paramVecMat",
                    contains="list",
                    validity=function(object){
                      if (any(lapply(object,class)!="list")) stop("@p should be a list of list")
                      if (any(lapply(object,length)!=length(object))) stop("@p should be squared")
                      for (subobject in object) {if (any(lapply(subobject,class)!="numeric")) stop("@p should be squared list of numeric vectors")}
                    }
)

setMethod("[","paramVecMat",
          definition = function(x, i, j, ..., drop) {
  x[[i]][[j]]
})

setMethod("dim","paramVecMat",
          definition = function(x) {
            c(length(x),length(x[[1]]))
          })

p=new("paramVecMat",list(list(c(1,2),0),list(0,c(.5,4))))
p[1,2]

dim(p)

r1unif <- function(min=0,max=1){runif(1,mean,sd)}
r1norm <- function(mean=0,sd=1){rnorm(1,mean,sd)}

prior <- setClass("prior",
                  contains="function",
                  slots=c(hyperParam="numeric"),
                  validity = function(object){if (any(names(formals(object))!=names(object@hyperParam))) stop("formals of prior and hyperParam slot names are different")}
)

priorList <- setClass("priorList",
                      contains="list",
                      validity = function(object){if (any(lapply(object,class)!="prior")) stop("priorList should be a list of prior objects")}
                      )

object=new("prior",r1unif,hyperparam=c(min=0,max=1))
object=new('priorList',list(a=object,b=new("prior",r1norm,hyperparam=c(mean=0,sd=1))))
names(object) <- c("b","c")
names(object)

bayesParam <- setClass("bayesParam",
                       slots=c(fun="xpFun",par="numeric",prior="priorList",ecoModel="ecoModel"),
                       validity = funciton(object){
                         if (any(names(object@par)!=names(object@prior))) stop("names of parameter vector and prior list do not coincide")
                         if (any(names(object@par)!=names(object@prior)) stop("names of parameter vector and prior list do not coincide")
                       })

                       

# varFunctions is a matrix of functions with slots 
# @p : parameter vector matrix (list of list) and 
# @Gamma : neighborhood matrix, presented as a 
varFunctionParam <- setClass("varFunctionParam",
                       contains="varFunctions",
                       slot=c(p="list",Gamma="logimat"),
                       validity=function(object){
                         if (any(lapply(object@p,class)!="list")) stop("@p should be a list of list")
                         if (any(lapply(object@p,length)!=length(object@p))) stop("@p should be squared")
                         if (any(dim(object@Gamma)!=dim(object))) stop ("neighborhood matrix and function matrix should be of the same size")
                         if (any(dim(object)!=c(length(object@p),length(object@p)))) stop("dimension of functions lists and dimentions of parameters lists cannot differ")
                         if (any(dim(object)!=c(length(object@p),length(object@p)))) stop("dimension of functions lists and dimentions of parameters lists cannot differ")
                         for (subobject in object@p) {if (any(lapply(subobject,class)!="numeric")) stop("@p should be squared list of numeric vectors")}
                         })

object=   new("varFunctionParam",vF, p=new("paramVecMat",list(list(c(1,2),0),list(0,c(.5,4)))),Gamma=as.logical(c(1,0,0,1))) 


timeFunctions <- setClass("timeFunctions",
                          contains="varFunctions")
timeFunctionParam <- setClass("timeFunctionParam",
                          contains="varFunctionParam")

object=   new("timeFunctions",list(list(new("xpFun",function(x,p){p*x}),new("xpFun",function(x,p){p[1]*x+p[2]*x^2})),list(new("xpFun",function(x,p){0}),new("xpFun",function(x,p){1}))))
object=   new("timeFunctionParameters",list(list(c(1,2),c(.5,4)),list(0,0))) 

#

ecoModel <- setClass("ecoModel",
                     contains=c("ecosysTimeSeries"),
                     slots=c(timeFunctions="list",timeModel="list",varFunctions="list", varModel="list",parameters="list"),
                     validity = function(object){
                       if (!all(levels(as.factor(unlist(object@parents)))%in%names(object))) stop("some parents in ecoModel are not in the ecosystem variables")
                       if (!all(levels(as.factor(unlist(object@parents)))%in%names(object))) stop("some parents in ecoModel are not in the ecosystem variables")
                       if (any(lapply(object@parents,length)!=(lapply(object@timeModel,length)))) stop("some parents in ecoModel are not in the ecosystem variables")
                     }
                       )
                     
                     
                     

3DGamma <- setClass("3DGamma",
                    contains="array",
                    validity=function(object){
                      if (length(dim(object)!=3)) stop("3DGamma should be a 3D array")
                    })
loG <- function(x=0,p=0)
linear <- new("XPFun",function(x=0,p=0){x*p})

varGammaFunctions <- setClass("varGammaFunctions",
                              slots = "list",
                              validity = function(object){
                                if (any(lapply(object,class)!="function")) stop("varGammaFunction is a list of functions, constructor received something else in the list") 
                              }
                              )



ecoGammaFunction <- setClassGammaFunctions

Gamma <- setClass("Gamma",
                  slots=c(GammaFunctions="function")

timeModel <- setClass("timeModel",
                      slots=c(type="character",parameters="numeric"))

data.frame[i,gamma(i)]

variableModel <- setClass("variableModel",
                          slots=c(parents="character",time="list",fun="list")

ecoModelSample <- setClass("ecoModelSample",
                           contains="ecoModel",
                           slots=c(paramValues="numeric"))



setClass("timeModel",
         contains="list",
         validity = function(object){
           lapply(object[["variable"]]
         })

model <- setClass("model",
                  slots=c(variables="character",parents="character",parameters="list", timeModel = "timeModel", residualDistribution ="function"),
                  validity = 
                  )


timeLinks <- setClass("timeLinks",
                      contains="list",
                      validity=function(object){(all(lapply(object,class)=="model"))})

a=c(as.formula(y~ a*x),as.formula(y~1),as.formula(y[1]~x[1]))
names(a)
new("timeLinks",a)

varLinks <- setClass("varLinks",
                      contains="list",
                      validity=function(object){if (!(all(lapply(object,class)=="formula"))) stop("varLinks constructor did not receive alist of formula")})


a=c(as.formula(y~ a*x),as.formula(y~1),as.formula(y[1]~x[1]))
new("varLinks",a)

prior <- setClass("prior",
                  contains = "list",
                  validity = function(object){})


ecoModelSimul <- setClass("ecoModelSimul",
                          contains="ecoModelSample",
                          slots=c(simulations="data.frame"),
                          validity = function(object){
                          if (colnames(simulations)!=)
                          }
                          )

prior <- setClass("pior",
                  )

proportional <- function(param=a,variable=x){param*variable}
exponential <- function(param=a,variable=x){exp(param*variable)}




functions <- setClass("functions",
                      slots=c(parameters="character",)
                      ) 
