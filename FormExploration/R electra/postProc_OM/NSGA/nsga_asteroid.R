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
Tfinal = 1.0
deltaT = 0.001
times = seq(0,Tfinal,by = deltaT )



##############################
###      nsga2 result      ###
##############################

# times
Tfinal = 1
deltaT = 0.001
times = seq(0,Tfinal,by = deltaT )

# trajectoire cible
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-4 ; v3=1
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD=rB/4,rE,rF,rG,rH=3*rB/4,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)



# résultat
dir = "data/resultsNSGA_astroid"
resNSGA = read.csv(paste0(dir,"/population15000.csv"))
df = as.data.frame(resNSGA)
dim(df)
i=which(df$obj==min(df$obj))[1]


lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=df$v2[i] ; v3=0
angleIni_B = df$angleIni_B[i]; angleIni_D = df$angleIni_D[i]; angleIni_F = df$angleIni_F[i]
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI,angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)



# comparaison (sur le même graphique)
# reference (astroide)
v1=1 ; v2=-4 ; v3=1
res_df_ref = positions_t_Stationnaire(times,v1,v2,v3, rB,rC,rD=rB/4,rE,rF,rG,rH=3*rB/4,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)

# reusltat nsga
v1=1 ; v2=df$v2[i] ; v3=0
angleIni_B = df$angleIni_B[i]; angleIni_D = df$angleIni_D[i]; angleIni_F = df$angleIni_F[i]
res_df_nsga = positions_t_Stationnaire(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI,angleIni_B,angleIni_D,angleIni_F)

# les deux figures sur le même graphique
p <- ggplot(res_df_ref) +
  expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  coord_fixed(ratio=1) +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) 

p = p + geom_point(aes(x=Dx, y=Dy), col= "yellow", size=0.1)
p = p + geom_point(aes(x=res_df_nsga$Dx, y=res_df_nsga$Dy), col= "blue", size=0.1)
p




# si on veut grader dans la trajectoire que des points avec une distamnce minimale entre eux
## import functions from utils local (same directory)
source("utils_nsga.R")
distance = 0.1
resSimpliTraj = simplifiedTrajectory(cbind(res_df_ref$Dx,res_df_ref$Dy), distance)
#mat = cbind(res_df_ref$Dx,res_df_ref$Dy)

plot(res_df_ref$Dx,res_df_ref$Dy, type = "l")
points(res_df_ref$Dx[1],res_df_ref$Dy[1],col="red")

plot(resSimpliTraj[,1],resSimpliTraj[,2], xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax, rayonMax), type = "p", cex = 0.2)
points(resSimpliTraj[1,1],resSimpliTraj[1,2], col="red")


# pour les coins pour lesquels on a l'impression qu'il y a trop de point,
# c'est parcequ'on retire pas mal de points, donc le coin n'est pas exactement la ou on l'imagine
m1=1
m2=14
plot(resSimpliTraj[m1:m2,1],resSimpliTraj[m1:m2,2], xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax, rayonMax), type = "p", cex = 0.5)
#points(resSimpliTraj[1,1],resSimpliTraj[1,2], col="red")
i = which(res_df_ref$Dy == max(res_df_ref$Dy) ) 
points(res_df_ref$Dx[i],res_df_ref$Dy[i],col="blue",cex = 0.3)
points(res_df_ref$Dx[(i-60):(i+60)],res_df_ref$Dy[(i-60):(i+60)],col="red",cex = 0.1)



