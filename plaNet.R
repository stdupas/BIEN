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

integerTimeSerie <- setClass("integerTimeSerie",
		      slots=c(values="integer",times="numeric"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @vaues list")
		      }
)

characterTimeSerie <- setClass("characterTimeSerie",
		      slots=c(values="character",times="numeric"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @vaues list")
		      }
)

logicalTimeSerie <- setClass("logicalTimeSerie",
		      slots=c(values="logical",times="numeric"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @vaues list")
		      }
)

numericTimeSerie <- setClass("numericTimeSerie",
		      slots=c(values="numeric",times="numeric"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @vaues list")
		      }
)

discreteTimeSerie <- setClass("discreteTimeSerie",
		      slots=c(values="integer",times="integer"),
		      validity = function(object){
		      if (length(object@values)!=length(object@times)) stop("length of @times vector differs from length of @vaues list")
		      }
)

TimeSeries <- setClass("TimeSeries",
		      contains="list",
		      validity = function(object){
		      if (!all(lapply(object,class)%in%c("numericTimeSerie","integerTimeSerie","logicalTimeSerie","characterTimeSerie"))) stop("trying to contrust a TimeSeries object with a list conaining objects other than numericTimeSerie, or discreteTimeSerie, or characterTimeSerie or integerTimeSerie")
		      }
)

TimeSeries <- function(X){
  if (class(X)== "list") {
    for (i in 1:length(X)) {
      if (!(class(X[[i]])%in%c("numericTimeSerie","integerTimeSerie","logicalTimeSerie","characterTimeSerie"))) {
	X[[i]]= switch(class(X[[i]]),
		matrix=new("numericTimeSerie",values=x[,1],times=x[,2]),
		data.frame= switch(class(X[[i]][,1]),
				   character= new("characterTimeSerie",values=X[[i]][,1],times=X[[i]][,2]),
				   numeric= new("numericTimeSerie",values=X[[i]][,1],times=X[[i]][,2]),
				   integer= new("integerTimeSerie",values=X[[i]][,1],times=X[[i]][,2]),
				   logical= new("logicalTimeSerie",values=X[[i]][,1],times=X[[i]][,2]))
		)
      }
    }   				   				   
  new("TimeSeries",X)
  }
else stop("needs a list as argument")
}

TS=TimeSeries(list(X1=data.frame(X1,tX1),X2=data.frame(X2,tX2),))

setMethod(
  f = "valuesA",
  signature = "RasterLayer",
  definition = function(object){
    #x=data.frame(variable=na.omit(values(object)))
    select <- !is.na(values(object))
    x=values(object)[select]
    names(x) <- which(select)
    #colnames(x)=names(object)
  x
  }
  )
