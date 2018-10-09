## PlaNet

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
  
intValCTS <- setClass("intValCTS",
		      slots=c(values="integer",times="numeric"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @values list")
		      }
)

charValCTS <- setClass("charValCTS",
		      slots=c(values="character",times="numeric"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @values list")
		      }
)

logicValCTS <- setClass("logicValCTS",
		      slots=c(values="logical",times="numeric"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @values list")
		      }
)

numValCTS <- setClass("numValCTS",
		      slots=c(values="numeric",times="numeric"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @values list")
		      }
)


  # CTS : list of continuous-time series

CTS <- setClass("CTS",
		      contains="list",
		      validity = function(object){
		      if (!all(lapply(object,class)%in%c("numValCTS","intValCTS","logicValCTS","charValCTS"))) stop("trying to contrust a TimeSeries object with a list conaining objects other than numericTimeSerie, or discreteTimeSerie, or characterTimeSerie or integerTimeSerie")
		      }
)

CTS <- function(X){
  if (class(X)== "list") {
    for (i in 1:length(X)) {
      if (!(class(X[[i]])%in%c("numValCTS","intValCTS","logicValCTS","charValCTS"))) {
	X[[i]]= switch(class(X[[i]]),
		matrix=new("numericTimeSerie",values=x[,1],times=x[,2]),
		data.frame= switch(class(X[[i]][,1]),
				   character= new("charValCTS",values=X[[i]][,1],times=X[[i]][,2]),
				   numeric= new("numValCTS",values=X[[i]][,1],times=X[[i]][,2]),
				   integer= new("intValCTS",values=X[[i]][,1],times=X[[i]][,2]),
				   logical= new("logicValCTS",values=X[[i]][,1],times=X[[i]][,2]))
		)
      }
    }   				   				   
  new("TimeSeries",X)
  }
else stop("needs a list as argument")
}

library(lubridate)

# disctre-time series

discreteTimeFrame <- setClass("discreteTimeFrame",
                              slots = c(start=c("POSIXct"),end=c("POSIXct"),periodUnit="character",periodLength="integer", tfLength="integer"),
                              #					validity = function(object){
                              #					if (as.integer((object@end-object@start)/object@period)!=(object@length+1)) stop("incompatibility between end, start and period")
                              #					}
)
discreteTimeFrame <- setClass("discreteTimeFrame",
                              slots = c(starts=c("POSIXt"),lasts=c("POSIXct"),period="Period", length="integer"),
                              prototype = prototype("discreteTimeFrame",starts=ymd_hms("2011-06-10-08-06-35"),lasts=ymd_hms("2011-06-10-08-06-40"),period=period(6,"seconds"),length=as.integer(6)),
                              validity = function(object){
                                if ((object@starts)+object@period*(object@length)!=(object@lasts)) stop("incompatibility between starts, lasts and period values")
                              }
)
new("discreteTimeFrame",start=ymd_hms("2011-06-10-08-06-35"),end=ymd_hms("2011-06-10-08-06-40"),periodUnit="seconds",periodLength=as.integer(6),tfLength=as.integer(6))


setMethod(length,
  signature = "discreteTimeFrame",
  definition = function(x){x@length})

dtf<-new("discreteTimeFrame",starts=ymd_hms("2011-06-10-08-06-35"),lasts=ymd_hms("2011-06-10-08-06-41"),period=period(1,"seconds"),length=as.integer(6))


indicatorTimeSeries <- setClass("indicatorTimeSeries",
                                slots=c(indicatorVariables="character",indicatorData="data.frame",times="discreteTimeFrame"),
                                validity = function(object){
                                  if (!all(object@indicatorVariables%in%colnames(object@indicatorData))) stop("not all indicatorVariables slot are in indicatorData slot colnames")
                                  if (length(object@times)!=nrow(object@indicatorData)) stop("discreteTimeFrame and number of row in indicatorData are different")
                                  }
                                )

object=iTS <- new("indicatorTimeSeries",indicatorVariables=c("X1","X2"),indicatorData=data.frame(X1=c("a","b","c","d","e","f"),X2=c(1,2,3,4,5,6)))

ecosysTimeSeries <- setClass("ecosysTimeSeries",
                             contains = "indicatorTimeSeries",
                             slots=c(ecosysVariables="character",ecosysData="data.frame"),
                             validity = function(object){
                               if (!all(object@ecosysVariables%in%colnames(object@ecosysData))) stop("not all indicatorVariables slot are in ecosysData slot colnames")
                               if (length(object@times)!=nrow(object@ecosysData)) stop("discreteTimeFrame and number of row in ecosysData are different")
                               }
                             )

setMethod("names",signature="indicatorTimeSeries",definition = function(x){x@indicatorVariables})
setMethod("names",signature="ecosysTimeSeries",definition = function(x){append(x@indicatorVariables,x@ecosysVariables)})

object=eTS <- new("ecosysTimeSeries",iTS,ecosysVariables=c("E1","E2"),ecosysData=data.frame(E1=c(2.45,3.5,4.,6.8,7.0,6.9),E2=as.integer(c(5,5,3,12,67,0))))

## Ecosystem-Indicator model

timeLinks <- setClass("timeLinks",
                      contains="list",
                      validity=function(object){(all(lapply(object,class)=="formula"))})

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

ecoModel <- setClass("ecoModel",
                     contains=c("ecosysTimeSeries"),
                     slots=c(varLinks="varLinks",timeLinks="timeLinks",parameters="character",prior="prior",posterior="posterior"),
                     )

ecoModelSample <- setClass("ecoModelSample",
                           contains="ecoModel",
                           slots=c(paramValues="numeric"))

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
