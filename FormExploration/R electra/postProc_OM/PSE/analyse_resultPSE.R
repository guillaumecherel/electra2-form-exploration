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




########################################
###   l'algo PSE a t'il convergé ?   ###
########################################


dir = "data/resultsPSE_5_3/"
dir = "data/resultsPSE_7/"
files = list.files(dir)
nbFiles = length(files)

# number of points
generation = seq(1000,1000*nbFiles, by = 5000)
generation = seq(200,200*nbFiles, by = 200)

sampleSize = c()
reduceSize = c()

options(scipen=999)
for (i in generation){
  temp <- read.csv(paste0(dir,"population",i,".csv"), header = T)
  sampleSize = c(sampleSize,dim(temp)[1])
  df = as.data.frame(temp)
  #fullSize=dim(resPSE)[1]
  #df = df %>%  filter(nbPointsSinguliers <=10,nbPointsSinguliers >=-10, nbPointsRetour <=10,  -50 <= courbureMoyenne, courbureMoyenne <=50, meanSpeed <=30)
  df = df %>%  filter(nbPointsSinguliersD <=5,nbPointsSinguliersF <=5, courbureMoyenneAbsolueD <=50, courbureMoyenneAbsolueF <=50, indiceMotifD<=2000,indiceMotifF<=2000)
  reduceSize= c(reduceSize,dim(df)[1]) 
}

plot(generation,sampleSize, type = "l")
plot(generation,reduceSize, type = "l")









########################################
###      HSITOGRAMME RESULTAT PSE    ###
########################################

## import dta
resPSE = read.csv("data/resultsPSE_5_3/population801000.csv")
dim(resPSE)

## prepare data
df = as.data.frame(resPSE)
df = as_tibble(resPSE)
df = df %>%  filter(nbPointsSinguliers <=10,nbPointsSinguliers >=-10, nbPointsRetour <=10,  -50 <= courbureMoyenne, courbureMoyenne <=50, meanSpeed <=30)
n = dim(df)[1];n
colnames(df)


# points singuliers
hist(df$nbPointsSinguliers)

# points retour
hist(df$nbPointsRetour)

# densite
max(df$densite)
min(df$densite)
hist(df$densite)

# moran
# max(df$moran)
# max(df$moran[df$moran<1])
# min(df$moran)
# hist(df$moran)
# hist(df$moran[df$moran<1])


# courbure moyenne
max(df$courbureMoyenne)
min(df$courbureMoyenne)
hist(df$courbureMoyenne)


# mean speed
min(df$meanSpeed)
max(df$meanSpeed)
hist(df$meanSpeed)


cor(df$nbPointsSinguliers,df$nbPointsRetour)






