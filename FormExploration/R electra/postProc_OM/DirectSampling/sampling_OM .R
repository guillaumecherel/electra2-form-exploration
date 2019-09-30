## 

## libraries
library(tibble)
library(tidyverse)
#library(plyr)
library(dplyr)
library(gganimate)
#library(gifski)
library(hypervolume)

## import functions from utils
source("../../utils_electra.R")


## Paramètres fixes
# longueur tiges
rB = 0.935; rC = 0.935
rH = sqrt(0.50^2+0.05^2); rI = sqrt(0.50^2+0.05^2)
rD = 0.43; rE = 0.38
rF = 0.417; rG = 0.417
rayonMax = max(max(rB, max(rH+rD, rH+rE)) , max(rC, max(rI+rF, rI+rG)) )

# angles
angleH = atan(0.05/0.5) 
angleI = atan(0.05/0.5)
angleIni_B = 0
angleIni_D = 0
angleIni_F = 0


# times
Tfinal = 10.0
deltaT = 0.0005
times = seq(0,Tfinal,by = deltaT )

# choix par défaut 
lightB=T ; lightC=T ; lightD=T ; lightE=T ; lightF=T ; lightG=T
v1=1 ; v2=1 ; v3=1



##################################
###  Resultat sampling (histogramme + quelques trajectoires)
##################################

## import dta
resSampling = read.csv("data/resultsDirectSampling3.csv")

## prepare data
df = as.data.frame(resSampling)
df = as_tibble(resSampling)
dim(resSampling)


# histogrammes des colonnes
lapply(1:ncol(df), function(i) hist(df[[i]], xlab= names(df)[i]))


# utils
plotTrajectoryFromDFValues <- function(i,df){
  temp = df[i,]
  plotTrajectorySmart(times, v1=temp$v1,v2=temp$v2,v3=temp$v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B=temp$angleIni_B,angleIni_D=temp$angleIni_D,angleIni_F=temp$angleIni_F, lightB,lightC,lightD,lightE,lightF,lightG)
}


# points singuliers
hist(df$numberPointSinguliersTotal)
max(df$numberPointSinguliersTotal)
selInd = which(df$numberPointSinguliersTotal == max(df$numberPointSinguliersTotal))
plotTrajectoryFromDFValues(selInd,df)

hist(df$meanDiffTempsPointSinguliers)

# densite 
max(df$densiteD)  # changé désormais <1
selInd = which(df$densiteD == max(df$densiteD))[1]
plotTrajectoryFromDFValues(selInd,df)

# Moran
max(df$moranD)
hist(df$moranD)
selInd = which(df$moranD == max(df$moranD))[1]
plotTrajectoryFromDFValues(selInd,df)

min(df$moranD)
selInd = which(df$moranD == min(df$moranD))
plotTrajectoryFromDFValues(selInd,df)


# Courbure
hist(df$meanCourbures)

# Retour
hist(df$nbTempsRetour)
hist((df %>% filter(nbTempsRetour <500))$nbTempsRetour)
selInd = which(df$meanTempsRetour == max(df$meanTempsRetour))
plotTrajectoryFromDFValues(selInd,df)






##################################
###  sampling avec les mesures du PSE
##################################

# points singuliers
resSampling = read.csv("data/resultsDirectSampling_PSE.csv")

dim(resSampling)
## prepare data
df = as.data.frame(resSampling)
df = as_tibble(resSampling)
n = dim(df)[1];n


#df$nbPointsSinguliers
sel0_singulier = which(df$nbPointsSinguliers == 0)
length(sel0_singulier)

# il y a plus de point singulier pour les petites vitesses, mais c'est surement à cause du seuil qui 
# ne depend pas de la vitesse

plot(1:n, abs(df$v1) +  abs(df$v2) +abs(df$v3) ) 
points((1:n)[sel0_singulier], abs(df$v1[sel0_singulier]) +  abs(df$v2[sel0_singulier]) +abs(df$v3[sel0_singulier]), col="red" )   
# faire dépendre le seuil des vitesses / des rayons



# faire un code R pour la densité, voir pourquoi le 0.8 max (visualiser)
# normalement on sait visualiser les points singuliers, les temps de retour (à partie de csv?calculés dans OM), les courbures
















