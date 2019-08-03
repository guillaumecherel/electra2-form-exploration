## 

## libraries
library(tibble)
library(tidyverse)
#library(plyr)
library(dplyr)

library(gganimate)
#library(gifski)

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

angleH = atan(0.05/0.5) 
angleI = atan(0.05/0.5)
Tfinal = 10.0
deltaT = 0.0005
times = seq(0,Tfinal,by = deltaT )
rayonMax = max(max(rB, max(rH+rD, rH+rE)) , max(rC, max(rI+rF, rI+rG)) )



## import dta
resSampling = read.csv("data/Sampling/resultsDirectSampling3.csv")

## prepare data
df = as.data.frame(resSampling)
df = as_tibble(resSampling)
dim(resSampling)


# histogrammes des colonnes
lapply(1:ncol(df), function(i) hist(df[[i]], xlab= names(df)[i]))

# utils
plotTrajectoryFromDFValues <- function(i,df){
  temp = df[i,]
  plotTrajectory_v3(times, v1=temp$v1,v2=temp$v2,v3=temp$v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B=temp$angleIni_B,angleIni_D=temp$angleIni_D,angleIni_F=temp$angleIni_F)
}


# points singuliers
hist(df$numberPointSinguliersTotal)
selInd = which(df$numberPointSinguliersTotal == max(df$numberPointSinguliersTotal))
plotTrajectoryFromDFValues(selInd,df)

hist(df$meanDiffTempsPointSinguliers)

# densite 
df$densiteD
max(df$densiteD)
selInd = which(df$densiteD == max(df$densiteD))
plotTrajectoryFromDFValues(selInd,df)

# Moran
hist(df$moranD)
df$moranD
max(df$moranD)
selInd = which(df$moranD == max(df$moranD))
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





# function save trajectory
dirSaveTraj = "ResultSampling/trajectories"
dir.create(dirSaveTraj)

saveTrajectoriesSampling <- function(i,df,times,dirSaveTraj) {
  #i=10 
  temp = df[i,]
  v1 = temp$v1 ; v2 = temp$v2 ; v3 = temp$v3
  angleIni_B = temp$angleIni_B  ; angleIni_D = temp$angleIni_D ; angleIni_F = temp$angleIni_F  
  p = plotTrajectory_v2(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
  
  # round param   
  v1 = round(v1,2)
  v2 = round(v3,2)
  v3 = round(v3,2)
  angleIni_B = round(angleIni_B,2)  
  angleIni_D = round(angleIni_D,2)  
  angleIni_F = round(angleIni_F,2)  
  
  params = paste0("_v1_",v1,"_v2_",v2,"_v3_",v3,"_angleIni_B_",angleIni_B,",_angleIni_D_",angleIni_D,",_angleIni_F_",angleIni_F)
  # save ggplot
  name_plot_file = paste0(dirSaveTraj, "/traj_",i,params,".png")
  ggsave(name_plot_file,p, width = 8, height = 8, dpi = 72)

  }

# i=3
# saveTrajectoriesSampling(i,df,times,dirSaveTraj)

## pour simuler les trajectoires correspond au sampling et les enregeistrer
# map(1:dim(df)[1],function(i){saveTrajectoriesSampling(i,df,times,dirSaveTraj)})


# change la couleur des points, du background
plotTrajectory_v3 <- function(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F) {
  theme_set(theme_bw())
  df = positions_t_Stationnaire_v2(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
  rayonMax = sqrt(df$Bx^2+df$By^2)[1]

  yellow <- rgb(235/255,235/255, 52/255) #235, 235, 52
  jaune2 <- rgb(235/255,229/255, 52/255) #235, 229, 52
  bg <- rgb(4/255,3/255, 51/255) #4, 3, 51
    
  p <- ggplot(df) +
    expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    coord_fixed(ratio=1) +
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_blank(), 
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = bg, colour = bg)) +
    geom_point(aes(x=Bx, y=By), size=.1, col = jaune2, pch=21, fill=yellow, alpha=0.2) +
    geom_point(aes(x=Cx, y=Cy), size=.1, col=jaune2, pch=21, fill=yellow, alpha=0.2) +
    geom_point(aes(x=Dx, y=Dy), size=.1, col= jaune2, pch=21, fill=yellow, alpha=0.2) +
    geom_point(aes(x=Ex, y=Ey), size=.1, col= jaune2, pch=21, fill=yellow, alpha=0.2) +
    geom_point(aes(x=Fx, y=Fy), size=.1, col= jaune2, pch=21, fill=yellow, alpha=0.2) +
    geom_point(aes(x=Gx, y=Gy), size=.1, col= jaune2, pch=21, fill=yellow, alpha=0.2)+
    labs(title = "Trajectoire Stationnaire",
         subtitle = paste0("v1=",v1,", v2=",v2,", v3=",v3,", \n",
                           "angleIni_B=",angleIni_B,", angleIni_D=",angleIni_D,", angleIni_F=",angleIni_F))
  
  p
}

plotTrajectory_v3(times, v1=1,v2=2,v3=0.01, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B=0,angleIni_D=1,angleIni_F=2)
#ggsave("ResultSampling/jaune/plot1.png")
plotTrajectory_v3(times, v1=1,v2=2.14,v3=2, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
#ggsave("ResultSampling/jaune/plot2.png")
plotTrajectory_v3(times, v1=2,v2=2,v3=1, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
#ggsave("ResultSampling/jaune/plot3.png")
plotTrajectory_v3(times, v1=1,v2=-pi,v3=sqrt(2), rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
#ggsave("ResultSampling/jaune/plot4.png")

















