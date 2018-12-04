setwd("/home/dupas/PlaNet/")

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
setwd
setwd("C:/Users/steph/OneDrive/Documents/GitHub/PlaNet")

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
names(p)

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

p_Pl_sd_r=0.05
p_T_sd=0.5
p_Pe_pDet=.8

for (k in 1:dim(edata)[3]){
	for (i in 3:dim(edata)[1]){
#Pe
  if (edata[i,5,k]) edata[i,5,k]=rbinom(1,1,p_Pe_pDet)
#Pl
  idata[i,4,k]=rnorm(1,edata[i,4,k],edata[i,4,k]*p_Pl_sd_r)
  if (idata[i,4,k]<0) idata[i,4,k]=0
#Pa
  idata[i,3,k]=rpois(1,edata[i,3,k])
#T
  idata[i,1,k]=rnorm(1,edata[i,1,k],p_T_sd)
#Pr
  idata[i,2,k]=rpois(1,edata[i,2,k])
  }
}
tmp0=c(idata)
tmp=lapply(1:dim(edata)[3], function(k){lapply(3:dim(edata)[1], function (i) {
Pl=rnorm(1,edata[i,4,k],edata[i,4,k]*p_Pl_sd_r)
if (Pl<0) Pl=0
c(T=rnorm(1,edata[i,1,k],p_T_sd),Pr=rpois(1,edata[i,2,k]),Pa=rpois(1,edata[i,3,k]),Pl=Pl,Pe=rbinom(1,edata[i,5,k],p_Pe_pDet))
})})


edata=idata
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
