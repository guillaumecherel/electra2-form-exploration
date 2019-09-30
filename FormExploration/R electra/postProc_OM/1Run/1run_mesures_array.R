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


# paramètre (depuis OM) 
lightB=T ; lightC=T ; lightD=T ; lightE=T ; lightF=T ; lightG=T
# v1=1.8 ; v2=(-2.66) ; v3=3.14  # from OM
v1=1 ; v2=(-3) ; v3=3.14  # from OM

plotTrajectorySmart(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG) 



## import data
dir = "data/Result_Mesure_Traj"
res = read.csv(paste0(dir,"/trajectoryElectra.csv"))

## prepare data
df = as.data.frame(res)
df = as_tibble(df)
n = dim(df)[1];n
times = df$times


# times
Tfinal = df$times[length(df$times)]  #10.0
deltaT =  diff(df$times)[1]   # 0.0005
times = df$times  #seq(0,Tfinal,by = deltaT )


## import functions from local (same directory)
source("utils_1run_mesures_array.R")


###############################################################
###   Plot Trajectory: data frame of points obtained in OM  ###
###############################################################


## plot one trajectoire stationnaire
p <- ggplot(df) +  
  expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  geom_point(aes(x=Bx, y=By), size=0.1) +
  geom_point(aes(x=Cx, y=Cy), size=.1) +
  geom_point(aes(x=Dx, y=Dy), size=.1, col= "red", shape = 3) +
  geom_point(aes(x=Ex, y=Ey), size=.1, col= "yellow") +
  geom_point(aes(x=Fx, y=Fy), size=.1, col= "blue") +
  geom_point(aes(x=Gx, y=Gy), size=.1, col= "green", shape = 3)
p



plotTrajectoryFromDf(df)



#######################################
###    TEMPS DE RETOUR SEGMENTS     ###
#######################################

## import data
resTempsRetourSegments = read.csv(paste0(dir,"/retour.csv"))

## prepare data
resTempsRetourSegments = as_tibble(resTempsRetourSegments)

# B
resTempsRetourSegmentsB = resTempsRetourSegments  %>% select(indiceRetourB,tempsRetourB)  %>% drop_na()
nRetourSegmentsB = dim(resTempsRetourSegmentsB)[1];nRetourSegmentsB

# D
resTempsRetourSegmentsD = resTempsRetourSegments  %>% select(indiceRetourD,tempsRetourD)  %>% drop_na()
nRetourSegmentsD = dim(resTempsRetourSegmentsD)[1];nRetourSegmentsD

# F
resTempsRetourSegmentsF = resTempsRetourSegments %>% select(indiceRetourF,tempsRetourF)  %>% drop_na() 
nRetourSegmentsF = dim(resTempsRetourSegmentsF)[1];nRetourSegmentsF


# analysis
# summary B
resTempsRetourSegmentsB %>% group_by(tempsRetourB) %>% summarise(count = n()) 
#resTempsRetourSegmentsB %>% summary()

# summary D
countRetourD = resTempsRetourSegmentsD %>% group_by(tempsRetourD) %>% summarise(count = n());countRetourD 
indiceRetourD = resTempsRetourSegmentsD$indiceRetourD
#resTempsRetourSegmentsD %>% summary()

# summary F
resTempsRetourSegmentsF %>% group_by(tempsRetourF) %>% summarise(count = n()) 
#resTempsRetourSegmentsF %>% summary()



# plot
plot(df$Dx,df$Dy, type = "l", xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax,rayonMax), asp=1)


# D
plot(res$Dx,res$Dy, type = "l", asp=1)
selInd = selInd[which(selInd < (1/deltaT))]
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




resRetourCount = as.matrix(countRetourD)
resRetourCountD = resRetourCount
# on regroupe les temps de retour proches
resRestourSimpliD = selectRetour(resRetourCountD)

# verif
#sum(resRestourSimpliD[,2])
#sum(countRetourD$count)


# tracé de temps de retour vs combien de point on ce temps
lg = dim(resRestourSimpliD)[1]
plot(resRestourSimpliD[1:(lg-1),1],resRestourSimpliD[1:(lg-1),2], xlab = "temps de retour", ylab = "nombre de points")

# avant (de regrouper les temps de retour espacés de 1)
#plot(countRetourD$tempsRetourD,countRetourD$count)
# maintenant
#points(resRestourSimpliD[,1],resRestourSimpliD[,2], col="red")
#plot(resRestourSimpliD[,1],resRestourSimpliD[,2])

# length(resRestourSimpliD[,1])
# diff(resRestourSimpliD[,1])



# version liste: on obtient les indices des points concernés (par un même temps de retour)
liste_tempsRetourD = create_liste_indice_tempsRetour(resRetourCountD,resTempsRetourSegmentsD)
resList_retour_D = resRestourSimpli_List = selectRetour_List(liste_tempsRetourD)

# palette de couleur ...
#colorsPalette = colorRampPalette(brewer.pal(9,"YlOrRd")[4:9])(length(resList)) 
colorsPalette = topo.colors(length(resList_retour_D))
#plot(1:length(resList_retour_D),1:length(resList_retour_D), col = colorsPalette)

# trace la trajectoire de D, repère les points de retour, la couleur indique le temps de retour
plot(df$Dx,df$Dy, type = "l", xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax,rayonMax))
for (i in 1:length(resList_retour_D)){
  selIndicesRetour_D = resList_retour_D[[i]][[3]]
  points(df$Dx[selIndicesRetour_D],df$Dy[selIndicesRetour_D], col=colorsPalette[i])
}

selLegende = floor(seq(1,length(resList_retour_D),length.out = 5))
legend("bottomleft",col = colorsPalette[selLegende], legend = c(resList_retour_D[[selLegende[1]]][1],
                                                       resList_retour_D[[selLegende[2]]][1],
                                                       resList_retour_D[[selLegende[3]]][1],
                                                       resList_retour_D[[selLegende[4]]][1],
                                                       resList_retour_D[[selLegende[5]]][1]), pch = c(1,1))




length(resList_retour_D)
i=2
resList_retour_D[[i]][[1]]
resList_retour_D[[i]][[2]]
selInd2 = resList_retour_D[[i]][[3]]
plot(res$Dx,res$Dy, type = "l", asp=1)
points(res$Dx[selInd2],res$Dy[selInd2], col="red")


plot(res$Dx[1:1000],res$Dy[1:1000], type = "l", asp=1)



###############################
###    POINTS SINGULIERS    ###
###############################

## import data
resTimesPointSinguliers = read.csv(paste0(dir,"/pointSinguliers.csv"))

## prepare data
resTimesPointSinguliers = as_tibble(resTimesPointSinguliers)

# B
resTimesPointSinguliersB = resTimesPointSinguliers %>% select(timesPointSinguliersB) %>% drop_na()  
nSingulierB = dim(resTimesPointSinguliersB)[1];nSingulierB
plot(df$Bx,df$By, type = "l", xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax,rayonMax))
points(df$Bx[resTimesPointSinguliersB$timesPointSinguliersB],df$By[resTimesPointSinguliersB$timesPointSinguliersB], col="red", cex = 0.5)

# D
resTimesPointSinguliersD = resTimesPointSinguliers %>% select(timesPointSinguliersD) %>% drop_na()  
nSingulierD = dim(resTimesPointSinguliersD)[1];nSingulierD
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
resCourbures = read.csv(paste0(dir,"/courbures.csv"))

## prepare data
resCourbures = as_tibble(resCourbures)

# B
courburesB=resCourbures$courburesB
plot(times,courburesB, type = "l")

# D
courburesD=resCourbures$courburesD
plot(times,courburesD, type = "l")

#seuilCourbureD = min(courburesD)+(max(courburesD)-min(courburesD))/2
seuilCourbureD_Max = max(courburesD)- 5/100*(max(courburesD)-min(courburesD))
seuilCourbureD_Min = min(courburesD)+ 5/100*(max(courburesD)-min(courburesD))
selIndBigCourbureD_max = which(courburesD > seuilCourbureD_Max)
selIndBigCourbureD_min = which(courburesD < seuilCourbureD_Min)
points(times[selIndBigCourbureD_max],courburesD[selIndBigCourbureD_max], col="red", cex = 0.5)
points(times[selIndBigCourbureD_min],courburesD[selIndBigCourbureD_min], col="blue", cex = 0.5)

plot(res$Dx,res$Dy, type = "l", asp=1)
points(res$Dx[selIndBigCourbureD_max],res$Dy[selIndBigCourbureD_max], col="red", cex = 0.5)
points(res$Dx[selIndBigCourbureD_min],res$Dy[selIndBigCourbureD_min], col="blue", cex = 0.5)

# F
courburesF=resCourbures$courburesF
plot(times,courburesF, type = "l")

seuilCourbureF_Max = max(courburesF)- 5/100*(max(courburesF)-min(courburesF))
seuilCourbureF_Min = min(courburesF)+ 5/100*(max(courburesF)-min(courburesF))
selIndBigCourbureF_max = which(courburesF > seuilCourbureF_Max)
selIndBigCourbureF_min = which(courburesF < seuilCourbureF_Min)
points(times[selIndBigCourbureF_max],courburesF[selIndBigCourbureF_max], col="red", cex = 0.5)
points(times[selIndBigCourbureF_min],courburesF[selIndBigCourbureF_min], col="blue", cex = 0.5)

plot(res$Fx,res$Fy, type = "l", asp=1)
points(res$Fx[selIndBigCourbureF_max],res$Fy[selIndBigCourbureF_max], col="red", cex = 0.5)
points(res$Fx[selIndBigCourbureF_min],res$Fy[selIndBigCourbureF_min], col="blue", cex = 0.5)







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





















