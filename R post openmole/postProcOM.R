## 

## libraries
library(tibble)
library(tidyverse)
#library(plyr)
library(dplyr)
library(ggplot2)
library(gganimate)
#library(gifski)
#library(pracma)
theme_set(theme_bw())


## import functions from utils
source("utils_post_proc.R")
# positions_t_Stationnaire_v2
# plot_oneTime_v2
# create_df_v2
# plotTrajectory_v2
# distanceSquared
# nextSimplifiedTrajectory
# simplifiedTrajectory


## Paramètres fixes
# longueur tiges (rayon)
rB = 0.935; rC = 0.935
rH = sqrt(0.50^2+0.05^2); rI = sqrt(0.50^2+0.05^2)
rD = 0.43; rE = 0.38
rF = 0.417; rG = 0.417
# H,I
angleH = atan(0.05/0.5) 
angleI = atan(0.05/0.5)

# position initiale (from angle)
angleIni_B = 0   # from OM
# D,E
angleIni_D = 1
# F,G 
angleIni_F = 2

rayonMax = max(max(rB, max(rH+rD, rH+rE)) , max(rC, max(rI+rF, rI+rG)) )


# 1 Plot Trajectory (extract param from OM res)
# 2 Visualize Measure Results: temps de retour, point singulier, courbure
# 3 obtenir une astroide
# 4 resultat de PSE, OSE, ...


## create and plot trajectory (not from OM)...
Tfinal = 10
deltaT = 0.0005
times = seq(0,Tfinal,by = deltaT )
v1=1.8 ; v2=(-2.66) ; v3=3.14  # from OM
plotTrajectory_v2(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F) 


dir = "One_trajectory_mesures"

###############################################################
###   Plot Trajectory: data frame of points obtained in OM  ###
###############################################################

## import data
res = read.csv(paste0("data/",dir,"/trajectoryElectra.csv"))

## prepare data
df = as.data.frame(res)
df = as_tibble(df)
n = dim(df)[1];n
#df$ind = 1:n
times = df$times


## plot one trajectoire stationnaire
#p <- ggplot(df %>% filter(times <= 0.2)) +
p <- ggplot(df) +  
  expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  geom_point(aes(x=Bx, y=By), size=1) +
  geom_point(aes(x=Cx, y=Cy), size=1) +
  geom_point(aes(x=Dx, y=Dy), size=1, col= "red", shape = 3) +
  geom_point(aes(x=Ex, y=Ey), size=1, col= "yellow") +
  geom_point(aes(x=Fx, y=Fy), size=1, col= "blue") +
  geom_point(aes(x=Gx, y=Gy), size=1, col= "green", shape = 3)
p



## plot trajectoire 1 point
#p <- ggplot(df %>% filter(times <= 0.5)) +
  p <- ggplot(df) +  
  expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  geom_point(aes(x=Dx, y=Dy), size=1) 
p



#######################################
###    TEMPS DE RETOUR SEGMENTS     ###
#######################################

## import data
resTempsRetourSegments = read.csv(paste0("data/",dir,"/retour.csv"))

## prepare data
resTempsRetourSegments = as.data.frame(resTempsRetourSegments)
resTempsRetourSegments = as_tibble(resTempsRetourSegments)

# B
resTempsRetourSegmentsB = resTempsRetourSegments  %>% select(indiceRetourB,tempsRetourB)  %>% drop_na()
nRetourSegmentsB = dim(resTempsRetourSegmentsB)[1];nRetourSegmentsB
#resTempsRetourSegments$ind = 1:nRetourSegmentsB

# D
resTempsRetourSegmentsD = resTempsRetourSegments  %>% select(indiceRetourD,tempsRetourD)  %>% drop_na()
nRetourSegmentsD = dim(resTempsRetourSegmentsD)[1];nRetourSegmentsD
#resTempsRetourSegmentsD$ind = 1:nRetourSegmentsD

# F
resTempsRetourSegmentsF = resTempsRetourSegments %>% select(indiceRetourF,tempsRetourF)  %>% drop_na() 
nRetourSegmentsF = dim(resTempsRetourSegmentsF)[1];nRetourSegmentsF
#resTempsRetourSegmentsF$ind = 1:nRetourSegmentsF


# analysis
# summary B
resTempsRetourSegmentsB %>% group_by(tempsRetourB) %>% summarise(count = n()) 
#resTempsRetourSegmentsB %>% summary()
# summary D
countRetourD = resTempsRetourSegmentsD %>% group_by(tempsRetourD) %>% summarise(count = n()) 
#resTempsRetourSegmentsD %>% summary()
# summary F
resTempsRetourSegmentsF %>% group_by(tempsRetourF) %>% summarise(count = n()) 
#resTempsRetourSegmentsF %>% summary()


# pas la bone notion l'histplot ? ils sont assez regroupés: une valeur par classe
hist(resTempsRetourSegmentsD$tempsRetourD)
hist(diff(resTempsRetourSegmentsD$indiceRetourD))
plot(countRetourD$tempsRetourD,countRetourD$count)

indiceRetourD = resTempsRetourSegmentsD$indiceRetourD

# plot
plot(df$Dx,df$Dy, type = "l", xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax,rayonMax))



nextSimplifyRetour = function(current,res,temp){
  if(isempty(current)){res} else{
    if(size(current)[1] > 1){temp2 = current[1,]} else{temp2 = current[1] }
    l = (size(res)[1])
    if (temp2[1] == temp[1]+1){
      newTemp = temp2
    if(size(current)[1] > 1){newCurrent = current[-1,]} else{newCurrent = c() }       
      if (l>1){
      newRes = rbind( res[1: (l-1), ] , c( res[l,1] , res[l,2] + temp2[2]) )} else{
        newRes = c( res[1] , res[2] + temp2[2]) 
      }
      nextSimplifyRetour(newCurrent,newRes,newTemp)
        }else{
        newTemp = temp2
        if(size(current)[1] > 1){newCurrent = current[-1,]} else{newCurrent = c() }       
        newRes = rbind( res ,  temp2 )
        nextSimplifyRetour(newCurrent,newRes,newTemp)

    }
  }
}


resRetourCount = as.matrix(countRetourD)
current = resRetourCount[2:(dim(resRetourCount)[1]),]

selectRetour = function(resRetourCount){
  res = resRetourCount[1,]
  temp = resRetourCount[1,]
  nextSimplifyRetour(resRetourCount[2:(dim(resRetourCount)[1]),],res,temp)
}

resRetourCount = as.matrix(countRetourD)
resRestourSimpli = selectRetour(resRetourCount)
ll = size(resRestourSimpli)[1] -1
plot(resRestourSimpli[1:ll,1],resRestourSimpli[1:ll,2])

sum(resRestourSimpli[1:ll,2])
sum(countRetourD$count)






# 

# D
plot(res$Dx,res$Dy, type = "l")
selInd = resTempsRetourSegmentsD$indiceRetourD
points(res$Dx[selInd],res$Dy[selInd], col="red")
length(res$Dx)
length(selInd)


# F
plot(res$Fx,res$Fy, type = "l")
selInd = resTempsRetourSegmentsF$indiceRetourF
#selInd = selInd[which(resTempsRetourSegmentsF$tempsRetourF)]
points(res$Fx[selInd],res$Fy[selInd], col="red")
length(res$Fx)
length(selInd)









###############################
###    POINTS SINGULIERS    ###
###############################

## import data
resTimesPointSinguliers = read.csv(paste0("data/",dir,"/pointSinguliers.csv"))

## prepare data
resTimesPointSinguliers = as.data.frame(resTimesPointSinguliers)
resTimesPointSinguliers = as_tibble(resTimesPointSinguliers)


# B
resTimesPointSinguliersB = resTimesPointSinguliers %>% select(timesPointSinguliersB) %>% drop_na()  
nSingulierB = dim(resTimesPointSinguliersD)[1];nSingulierB

plot(df$Bx,df$By, type = "l", xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax,rayonMax))
points(df$Bx[resTimesPointSinguliersB$timesPointSinguliersB],df$By[resTimesPointSinguliersB$timesPointSinguliersB], col="red", cex = 0.5)

# D
resTimesPointSinguliersD = resTimesPointSinguliers %>% select(timesPointSinguliersD) %>% drop_na()  
nSingulierD = dim(resTimesPointSinguliersD)[1];nSingulierD
#resTimesPointSinguliersD$ind = 1:nSingulierD

plot(df$Dx,df$Dy, type = "l", xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax,rayonMax))
points(df$Dx[resTimesPointSinguliersD$timesPointSinguliersD],df$Dy[resTimesPointSinguliersD$timesPointSinguliersD], col="red", cex = 0.5)

# F
resTimesPointSinguliersF = resTimesPointSinguliers %>% select(timesPointSinguliersF) %>% drop_na()  
nSingulierF = dim(resTimesPointSinguliersF)[1];nSingulierF
plot(df$Fx,df$Fy, type = "l", xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax,rayonMax))
points(df$Fx[resTimesPointSinguliersF$timesPointSinguliersF],df$Fy[resTimesPointSinguliersF$timesPointSinguliersF], col="red", cex = 0.5)




#########################
###     COURBURE      ###
#########################

## import data
resCourbures = read.csv(paste0("data/",dir,"/courbures.csv"))

## prepare data
resCourbures = as.data.frame(resCourbures)
resCourbures = as_tibble(resCourbures)

# B
courburesB=resCourbures$courburesB
plot(times,courburesB, type = "l")

# D
courburesD=resCourbures$courburesD
plot(times,courburesD, type = "l")

# F
courburesF=resCourbures$courburesF
plot(times,courburesF, type = "l")


# D
seuilCourbureD = min(courburesD)+(max(courburesD)-min(courburesD))/2
selIndBigCourbureD = which(courburesD > seuilCourbureD)
plot(res$Dx,res$Dy, type = "l")
points(res$Dx[selIndBigCourbureD],res$Dy[selIndBigCourbureD], col="red")

plot(times,courburesD, type = "l")
points(times[selIndBigCourbureD],courburesD[selIndBigCourbureD], col="red")


# F
seuilCourbureF = min(courburesF)+(max(courburesF)-min(courburesF))/2
selIndBigCourbureF = which(courburesF > seuilCourbureF)
plot(res$Fx,res$Fy, type = "l")
points(res$Fx[selIndBigCourbureF],res$Fy[selIndBigCourbureF], col="red")

plot(times,courburesF, type = "l")
points(times[selIndBigCourbureF],courburesF[selIndBigCourbureF], col="red")




##########################################
###     Animation avec gganimate       ###
##########################################

# 1, partial
p <- ggplot(df, aes(x = Dx, y=Dy) ) +
  expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  coord_fixed(ratio=1) +
  geom_point(show.legend = FALSE) 
p
p + transition_time(times) 
  

# 2, complet  
p <- ggplot(df) +
  expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  coord_fixed(ratio=1) +
  geom_point(aes(x=Bx, y=By), col="yellow") +
  geom_point(aes(x=Cx, y=Cy), col="yellow") +
  geom_point(aes(x=Dx, y=Dy), col="yellow") +
  geom_point(aes(x=Ex, y=Ey), col="yellow") +
  geom_point(aes(x=Fx, y=Fy), col="yellow") +
  geom_point(aes(x=Gx, y=Gy), col="yellow")
p
p + transition_time(times) 





















