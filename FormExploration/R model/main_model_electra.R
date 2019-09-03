## 

## libraries
library(RColorBrewer)
library(tibble)
library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)
library(gganimate)
library(gifski)
library(pracma)
theme_set(theme_bw())

## import functions from utils
source("utils_model_electra.R")
# fonctions importées:
# positions_t_NONsymetrique
# plot_oneTime
# create_df
# savePlot_index
# save_params_in_file
# savePlot_index_and_paramsInFile
# create_df_persistence
# plot_oneTime_persistence
# scalarProdutRenormalized
# vitesse_rotation_moteur
# positions_t_TRANSITOIRE
# plot_oneTime_TRANSITOIRE
# create_df_TRANSITOIRE
# savePlot_index_TRANSITOIRE
# save_params_in_file_TRANSITOIRE
# savePlot_index_and_paramsInFile_TRANSITOIRE


#########################  Nouveau Stationnaire v2
## Paramètres
# vitesse moteur (tour/ seconde)
v1 = 0.5
v2 = -0.25
v3 = 3

# longueur tiges (rayon)
rB = 4; rC = 4
rH = 2; rI = 2
rD = 1; rE = 1
rF = 1; rG = 1

# position initiale (from angle)
# A,B,C
Ax = 0
Ay = 0
angleIni_B = 0

# H,I
angleH = 0.25 
angleI = 0.25 

# D,E
angleIni_D = 0
# F,G 
angleIni_F = 0


# time
timeStep = 0.005 # sec
#timeStep = 0.02 # sec
Tfinal = 5
times = seq(0,Tfinal, by = timeStep)
t=0.1
t=0


positions_t_Stationnaire_v2(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)

plot_oneTime_v2(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)



##########################################
###   Dépendance condition initiales ? ###
##########################################
# angles par défaut nuls
angleIni_F= 0
#plot_oneTime_v2(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
plotTrajectory_v2(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)

# angles par défaut nuls
angleIni_F= 1
plotTrajectory_v2(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)



timeStep = 0.0005 # sec
Tfinal = 4
times = seq(0,Tfinal, by = timeStep)
plotTrajectory_v2(times, v1=5,v2=4,v3=5, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)



# remplir l'espace (densité) ?
timeStep = 0.005 # sec
Tfinal = 25
times = seq(0,Tfinal, by = timeStep)
plotTrajectory_v2(times, v1,v2,v3=pi, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)





# nsga2 result
v1=5.0
v2=3.0537391026395397
v3=-2.496259000231632
angleIni_B = 0.0
angleIni_D = 0.473719893001516
angleIni_F = 1.0229217740760777


timeStep = 0.005 # sec
Tfinal = 4
times = seq(0,Tfinal, by = timeStep)
df = create_df_v2(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
rayonMax = sqrt(df$Bx^2+df$By^2)[1]

distance = 1
resSimpliTraj = simplifiedTrajectory(cbind(df$Dx,df$Dy), distance)
dim(resSimpliTraj)

n=27
plot(resSimpliTraj[1:n,1],resSimpliTraj[1:n,2], xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax, rayonMax))
points(resSimpliTraj[1,1],resSimpliTraj[1,2], xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax, rayonMax), col="red")

plot(resSimpliTraj[,1],resSimpliTraj[,2], xlim = c(-rayonMax,rayonMax), ylim = c(-rayonMax, rayonMax), type = "l")



#p <- ggplot(df %>% filter(time <0.7)) +  
p <- ggplot(df ) +  
    expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  geom_point(aes(x=Dx, y=Dy), size=1, col= "red", shape = 3) +
  labs(title = "Trajectoire Stationnaire",
       subtitle = paste0("v1=",v1,", v2=",v2,", v3=",v3,", \n",
                         "angleIni_B=",angleIni_B,", angleIni_D=",angleIni_D,", angleIni_F=",angleIni_F))

p





# save a trajectory, true electra param
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

v1=1.8 ; v2=(-2.66) ; v3=3.14  # from OM

# save plots for different times
dirRes = "plotElectra2_1"
dir.create(dirRes)

# time
timeStep = 0.005 # sec
times = seq(0,10, by = timeStep)


#savePlot_index(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
# savePlot_index_and_paramsInFile(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
savePlot_index_v2(times,dirRes,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)

# command line in linux to create a video from images with index (adapt the framerate)
# ffmpeg -framerate 1/0.005 -i plot_%01d.png -crf 15  output1.mp4

















# old modèle

## Paramètres
# vitesse moteur (tour/ seconde)
v1 = 0.5
v2 = -0.25
v3 = 3

# longueur tiges (rayon)
r1B = 4; r1C = 4
r2D = 1; r2E = 1
r3F = 1; r3G = 1

# position initiale (from angle)
# A,B,C
Ax = 0
Ay = 0
angleIni_B = 0

# H,I
alphaH = 2 # >1 (=1 : moteur à la loupiotte)
alphaI = 2

# D,E
angleIni_D = 1
# F,G 
angleIni_F = 2


# intensité lumière

# colors palette
colors = colorRampPalette(brewer.pal(9,"YlOrRd")[4:9])(8) 
# colors
# plot(0,0, xlim = c(1,8), ylim = c(1,8))
# for (i in 1:8){
#   points(i,i, col= colors[i])
# }
oneColor = rgb(255/255,255/255,0/255)


# time
timeStep = 0.005 # sec
timeStep = 0.02 # sec
1/timeStep  # nombre d'image par seconde
times = seq(0,5, by = timeStep)
t=0.1


# Persistence rétinienne
persistenceTime = 0.05
nb_point_persistence = 50



# calcul des positions au temps donné
t=0  # initial conditions
t=0.4
positions_t_NONsymetrique(t,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F)



# trace position au temps donné
t=0  # initial conditions
t=1.45
plot_oneTime(t,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F)





# save plots for different times
dirRes = "plot1"
#dirRes = paste0("plot_","v1=",v1,"_v2=",v2,"_v3=",v3)
dir.create(dirRes)
#times = seq(0,4, by = 0.02)
#savePlot_index(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
# savePlot_index_and_paramsInFile(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)

# command line in linux to create a video from images with index (adapt the framerate)
# ffmpeg -framerate 1/0.02 -i plot_%01d.png -crf 15  output1.mp4


# save current parameters
# dirRes = "params"
# save_params_in_file(dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)








# time
timeStep = 0.005 # sec
times = seq(0,20, by = timeStep)

# param de reférence
df = create_df(times,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
df = create_df(times,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)

# astroide

# influence vitesse: multiplication un même facteur
df = create_df(times,3*v1,3*v2,3*v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
# influence vitesse: seul v2 modifié 
df = create_df(times,v1,(2*v2),v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
df = create_df(times,v1+0.1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
df = create_df(times,v1,(1.9*v2),v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
df = create_df(times,v1,3*v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
df = create_df(times,v1,0,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
df = create_df(times,v1,0,0,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,0,angleIni_F)
# v1=0
df = create_df(times,0,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)

df = create_df(times,v1,0.3,sqrt(2)*v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,0,angleIni_F)


# mesures 
df = df %>% mutate( # distance
                    AD = sqrt( Dx^2 + Dy^2),
                    AF = sqrt( Fx^2 + Fy^2),
                    DF = sqrt( (Fx-Dx)^2 + (Fy-Dy)^2),
                    # distance: pas de suspens
                    HD = sqrt( (Hx-Dx)^2 + (Hy-Dy)^2),
                    HE = sqrt( (Hx-Ex)^2 + (Hy-Ey)^2),
                    FI = sqrt( (Fx-Ix)^2 + (Fy-Iy)^2),
                    GI = sqrt( (Gx-Ix)^2 + (Gy-Iy)^2),
                    # speed: composante tangentielle
                    speedB_Tg = scalarProdutRenormalized(speedBx,speedBy,Bx,By),
                    speedD_Tg = scalarProdutRenormalized(speedDx,speedDy,Dx,Dy),
                    # speed: composante normale
                    speedD_No = scalarProdutRenormalized(speedDx,speedDy,-Dy,Dx),
                    speedB_No = scalarProdutRenormalized(speedBx,speedBy,-By,Bx),
                    # acceleration: composante tangentielle
                    accB_Tg = scalarProdutRenormalized(accBx,accBy,Bx,By),
                    accD_Tg = scalarProdutRenormalized(accDx,accDy,Dx,Dy),
                    # acceleration: composante normale
                    accD_No = scalarProdutRenormalized(accDx,accDy,-Dy,Dx),
                    accB_No = scalarProdutRenormalized(accBx,accBy,-By,Bx)
                    )

# points


# create a dataframe with positions at times prescribed


# astroide
timeStep = 0.005 # sec
times = seq(0,5, by = timeStep)
df = create_df(times,v1,v2=-3/4*v1,v3,r1B,r1C,r1B/8,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)


#p <- ggplot(df %>% filter(times < 1.5)) + coord_fixed(ratio=1) +  
p <- ggplot(df) + coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  #geom_point(aes(x=Bx, y=By), size=1) +
  #geom_point(aes(x=Cx, y=Cy), size=1) +
  geom_point(aes(x=Dx, y=Dy), size=1, col= "red", shape = 3) #+
  #geom_point(aes(x=Ex, y=Ey), size=1, col= "yellow") +
  #geom_point(aes(x=Fx, y=Fy), size=1, col= "blue") +
  #geom_point(aes(x=Gx, y=Gy), size=1, col= "green")
p








#p <- ggplot(df %>% filter(times < 1.5)) + coord_fixed(ratio=1) +  
p <- ggplot(df) + coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  geom_point(aes(x=Bx, y=By), size=1) +
  geom_point(aes(x=Cx, y=Cy), size=1) +
  geom_point(aes(x=Dx, y=Dy), size=1, col= "red", shape = 3) +
  geom_point(aes(x=Ex, y=Ey), size=1, col= "yellow") +
  geom_point(aes(x=Fx, y=Fy), size=1, col= "blue") +
  geom_point(aes(x=Gx, y=Gy), size=1, col= "green")
p


v1*r1B/alphaH - (v1+v2)*r2D


# mesures
#p <- ggplot(df %>% filter(times < 5), aes(x=time)) + 
p <- ggplot(df, aes(x=time)) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) +  
  # distance
  geom_line(aes(y=HD), size=1) +
  geom_line(aes(y=DF), size=1, col= "red") +
  geom_line(aes(y=AD), size=1, col= "blue") #+
  #geom_line(aes(y=AF), size=1, col= "green") +
  # speed
  # geom_line(aes(y=speedB_Tg), size=1, col= "green") +
  # geom_line(aes(y=speedB_No), size=1, col= "orange")+
  #geom_line(aes(y=speedD_Tg), size=1, col= "green") +
  #geom_line(aes(y=speedD_No), size=1, col= "orange")
  # acceleration
  # geom_line(aes(y=accB_Tg), size=1, col= "green") +
  # geom_line(aes(y=accB_No), size=1, col= "orange")+
  # geom_line(aes(y=accD_Tg), size=1, col= "green") +
  # geom_line(aes(y=accD_No), size=1, col= "orange")
p








#########################
###   cas 1: v1 = 0   ###
#########################
times = seq(0,4, by = 0.005)
rayonMax = max(r1B,r1C) + max(r2D,r2E,r3F,r3G)
res4 = create_df(times,v1=0,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)

p <- ggplot(res4 %>% filter(times < 1.5)) + coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  geom_point(aes(x=Bx, y=By), size=1)+
  geom_point(aes(x=Cx, y=Cy), size=1)+
  geom_point(aes(x=Dx, y=Dy), size=1, col= "red") +
  geom_point(aes(x=Fx, y=Fy), size=1, col= "yellow")
p


#############################
###   cas 2: v2, v3 = 0   ###
#############################
times = seq(0,4, by = 0.005)
res5 = create_df(times,v1,v2=0,v3=0,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D+0.2,angleIni_F)
rayonMax = max(r1B,r1C) + max(r2D,r2E,r3F,r3G)
p <- ggplot(res5 %>% filter(times < 1.5)) + coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  geom_point(aes(x=Bx, y=By), size=1)+
  geom_point(aes(x=Cx, y=Cy), size=1)+
  geom_point(aes(x=Dx, y=Dy), size=1, col= "red") +
  geom_point(aes(x=Ex, y=Ey), size=1, col= "green") +
  geom_point(aes(x=Fx, y=Fy), size=1, col= "yellow")
p



########################################
###     Persistence rétinienne       ###
########################################
t=1
plot_oneTime_persistence(t,persistenceTime,nb_point_persistence,10*v1, 10*v2, 10*v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F)




##########################################
###   Dépendance condition initiales ? ###
##########################################





##########################################
###     Animation avec gganimate       ###
##########################################
timeStep = 0.01 # sec
times = seq(0,3, by = timeStep)
df = create_df(times,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)

# 1, partial
p <- ggplot(df, aes(x = Bx, y=By) ) +
  geom_point(show.legend = FALSE) 
p
p + transition_time(time) 
  

# 2, complet  
p <- ggplot(df) +
  geom_point(aes(x=Bx, y=By), col="yellow") +
  geom_point(aes(x=Cx, y=Cy), col="yellow") +
  geom_point(aes(x=Dx, y=Dy), col="yellow") +
  geom_point(aes(x=Ex, y=Ey), col="yellow") +
  geom_point(aes(x=Fx, y=Fy), col="yellow") +
  geom_point(aes(x=Gx, y=Gy), col="yellow") +
  geom_point(aes(x=Hx, y=Hy), col="black")  +
  geom_point(aes(x=Ix, y=Iy), col="black")
p
p + transition_time(time) 


p + transition_states(time,
                  transition_length = 2,
                  state_length = 1)
p



##########################################
###   find different patterns by hand  ###
###         for a single light         ###
##########################################
plot_test_0 <- function(tt,a,b,u,v) {
  x= a*cos(2*pi*tt*u+angleIni_B) + b*cos(2*pi*tt*(u+v)+angleIni_B+angleIni_D) 
  y= a*sin(2*pi*tt*u+angleIni_B) + b*sin(2*pi*tt*(u+v)+angleIni_B+angleIni_D) 
  plot(x,y, type="l")    
}  

plot_test <- function(tt,a,b,u,v) {
  # creates points
  x= a*cos(2*pi*tt*u+angleIni_B) + b*cos(2*pi*tt*(u+v)+angleIni_B+angleIni_D) 
  y= a*sin(2*pi*tt*u+angleIni_B) + b*sin(2*pi*tt*(u+v)+angleIni_B+angleIni_D) 
  res = data_frame(tt,x,y)
  # plot
  rayonMax = a+b
    p <- ggplot(res) + expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    geom_point(aes(x=x, y=y), col="yellow", size = 1) 
  p
}  


#tt = seq(0,200,by=0.001)
#tt = seq(0,500,by=0.001)
tt = seq(0,6,by=0.001)
ttt = seq(0,1,by=0.001)
#tt = seq(0,1,by=0.01)
epsilon = 0.012

a=2
b=1
# u
u=0.5
# v
v=0.5 + epsilon


plot_test(tt,a=2,b=1,u=0.5,v=0.5)
plot_test(tt,a=2,b=1,u=0.5,v=0.5+epsilon)
plot_test(tt,a=2,b=1,u=1,v=0.5)
plot_test(tt,a=2,b=1,u=1.5,v=0.5)
plot_test(tt,a=2,b=1,u=2.5,v=0.5)
plot_test(tt,a=2,b=1,u=1,v=2)
plot_test(tt,a=2,b=1,u=10/3,v=0.5)
plot_test(tt,a=2,b=1,u=pi,v=0.5)
plot_test(tt,a=2,b=1,u=sqrt(2),v=0.5)
plot_test(tt,a=2,b=1,u=4,v=0.5)


# cercles délimitant l'expace si trajectoire dense ?
plot((a+b)*cos(2*pi*ttt),(a+b)*sin(2*pi*ttt), col="red", type = "l")
points((a-b)*cos(2*pi*ttt),(a-b)*sin(2*pi*ttt), col="red", type = "l")

# vitesse u negative
plot_test(tt,a=2,b=1,u=-0.5,v=0.5)
plot_test(tt,a=2,b=1,u=-0.5,v=0.5+epsilon)

# vitesse v negative
plot_test(tt,a=2,b=1,u=0.5,v=-0.5)
plot_test(tt,a=2,b=1,u=0.5,v=-0.5+epsilon)

a*u - b*(u+v)



tt=seq(0,100,by=0.01)
derivX = -2*pi*u*a*sin(2*pi*tt*u+angleIni_B) -2*pi*(u+v)*b*cos(2*pi*tt*(u+v)+angleIni_B+angleIni_D)
plot(tt,derivX, type = "l")








######################################
###       Find loop points         ###
######################################






######################################
###       Non stationary           ###
######################################
# aim: changer la vitesse des moteurs au cours du temps, tenir compte de l'inertie
# https://www.wikimeca.org/index.php/Moteur_%C3%A0_courant_continu

# partie électrique
R = 1 # R est la résistance électrique interne du moteur (Ohm);
Ke = 1 # Ke est la constante de force électromotrice

# partie mécanique
J = 1 # rotor vu comme un volant d'inertie J  
Kc = 1 # Kc est la constante de couple  
f = 1 # f est le coefficient de frottement visqueux.

A =  Kc /(R*f+Kc*Ke)
C = J*R/(R*f+Kc*Ke)


# time
timeStep = 0.005 # sec
timeStep = 0.02 # sec
times = seq(0,5, by = timeStep)

# u: tension imposée = vitesse demandée au moteur

u = function(t){
  1*(t>0)*(t<1)+
  0.5*(t>=1)*(t<2)+
  2*(t>=2)*(t<3)+
  0*(t>=3)  
}

plot(times,u(times), type = "l")

t=1
res = vitesse_rotation_moteur(t,u,A,C)
res = vitesse_rotation_moteur(times,u,A,C)
plot(times,res, type = "l")

# angle
angle_t =  integrate(function(x){vitesse_rotation_moteur(x,u,A,C)},0,t, stop.on.error = F)$value
angle_all = sapply(times, function(a){integrate(function(x){vitesse_rotation_moteur(x,u,A,C)},0,a, stop.on.error = F)$value})
plot(times,angle_all)


# moteur 1
u1=u
A1=A
C1=C
# moteur 2
u2=u
A2=A
C2=C
# moteur 3
u3=u
A3=A
C3=C



t=0
t=1
# positions at time t
positions_t_TRANSITOIRE(t, r1B, r1C, r2D, r2E, r3F, r3G, 
                                   alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                   A1,C1,A2,C2,A3,C3,
                                   u1,u2,u3)


# plot positions at time t
plot_oneTime_TRANSITOIRE(t, r1B, r1C, r2D, r2E, r3F, r3G, 
             alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
             A1,C1,A2,C2,A3,C3,
             u1,u2,u3)



# des u differents
u1 = function(t){
  2
}


u2 = function(t){
  1
}


u3 = function(t){
  3
}


# create data frame of positions % times
df =create_df_TRANSITOIRE(times, r1B, r1C, r2D, r2E, r3F, r3G, 
                      alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                      A1,C1,A2,C2,A3,C3,
                      u1,u2,u3)


p <- ggplot(df) + coord_fixed(ratio=1) + 
  theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
  geom_point(aes(x=Bx, y=By), size=1) +
  geom_point(aes(x=Cx, y=Cy), size=1) +
  geom_point(aes(x=Dx, y=Dy), size=1, col= "red", shape = 3) +
  geom_point(aes(x=Ex, y=Ey), size=1, col= "yellow") +
  geom_point(aes(x=Fx, y=Fy), size=1, col= "blue") +
  geom_point(aes(x=Gx, y=Gy), size=1, col= "green")
p


# res = vitesse_rotation_moteur(times,u1,A1,C1)
# plot(times, res, type = "l")


dirRes = "plot_transitoire_2"
dir.create(dirRes)
times = seq(0,5, by = 0.02)
# start_time <- Sys.time()
# savePlot_index_and_paramsInFile_TRANSITOIRE(times, dirRes, r1B, r1C, r2D, r2E, r3F, r3G, 
#                                             alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
#                                             A1,C1,A2,C2,A3,C3,
#                                             u1,u2,u3)
# end_time <- Sys.time()
# end_time - start_time

# savePlot_index_and_paramsInFile(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)

# command line in linux to create a video from images with index (adapt the framerate)
# ffmpeg -framerate 1/0.02 -i plot_%01d.svg -crf 15  output1.mp4



# test other computation integral
#quadinf(fun, -Inf, Inf, tol=1e-10)
start_time <- Sys.time()
angle_all = sapply(times, function(a){integrate(function(x){vitesse_rotation_moteur(x,u,A,C)},0,a, stop.on.error = F)$value})
end_time <- Sys.time()
end_time - start_time


start_time <- Sys.time()
angle_all = sapply(times, function(a){quadinf(function(x){vitesse_rotation_moteur(x,u,A,C)},0,a)$value})
end_time <- Sys.time()
end_time - start_time






