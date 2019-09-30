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



####################################
###    Inventaire des figures    ###
####################################

# sans nom 1
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-3/4*v1 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)
plotTrajectorySmart2(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)

# sans nom 1bis
lightB=F ; lightC=F ; lightD=T ; lightE=T ; lightF=F ; lightG=F
v1=1 ; v2=-3/4*v1 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# sans nom 2
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=1 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# sans nom 2 double
lightB=F ; lightC=F ; lightD=T ; lightE=T ; lightF=F ; lightG=F
v1=1 ; v2=1 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# pour être vraiment symetrique, préférer F et G
# sans nom 2 double (symetrique)
lightB=F ; lightC=F ; lightD=F ; lightE=F ; lightF=T ; lightG=T
v1=1 ; v2=0 ; v3=1
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# Adobe Reader
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-3*v1 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# shamrock
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-4*v1 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# sans nom à 4
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=4*v1 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# v1=0
lightB=F ; lightC=F ; lightD=T ; lightE=T ; lightF=T ; lightG=T
v1=0 ; v2=1 ; v3=1
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# v2=0, v3=0
lightB=T ; lightC=F ; lightD=T ; lightE=T ; lightF=T ; lightG=F
v1=1 ; v2=0 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D+pi/2,angleIni_F+pi/4,lightB,lightC,lightD,lightE,lightF,lightG)


# ellipse
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-2 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# sans nom
lightB=F ; lightC=F ; lightD=T ; lightE=T ; lightF=F ; lightG=F
v1=1 ; v2=2 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F+pi/4,lightB,lightC,lightD,lightE,lightF,lightG)


# cercle décentré 
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-1 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# coeur dans cercle
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=1/2 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# sans nom 
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-1/2 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# astroide (il faut forcément changer la géométrie du probème?)
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-4 ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD=rB/4,rE,rF,rG,rH=3*rB/4,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)


# point singuliers
lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=v1*(rH-rD)/rD ; v3=0
v1=1 ; v2=-v1*(rH+rD)/rD ; v3=0
plotTrajectorySmart(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)




###############
###  Test   ###
###############

Tfinal = 1.0
deltaT = 0.001
times = seq(0,Tfinal,by = deltaT )

lightB=F ; lightC=F ; lightD=T ; lightE=F ; lightF=F ; lightG=F
v1=1 ; v2=-pi ; v3=2
v1=1 ; v2=-3 ; v3=2
v1=1 ; v2=0.2 ; v3=2
v1=1 ; v2=3 ; v3=2
plotTrajectory(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)
res = positions_t_Stationnaire(times,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)

plot(res$Dx,res$Dy, cex=0.5, asp=1)
i=1;N=5
points(res$Dx[i:(i+N)],res$Dy[i:(i+N)], col="red")
i=333
points(res$Dx[i:(i+N)],res$Dy[i:(i+N)], col="green")


# par angle
N = 20
rayon = 0.7
i=1
#vec_angle = 2*pi*(0:(N-1))/N
#points(rayon*cos(vec_angle),rayon*sin(vec_angle), col="red")
for (i in 0:(N-1)){
  abline(0,sin(2*pi*i/N)/cos(2*pi*i/N), col="red")
}


tt = 1:1000
r1=0.4
r2=0.9
points(r1*cos(tt),r1*sin(tt),col="green")
points(r2*cos(tt),r2*sin(tt),col="blue")











