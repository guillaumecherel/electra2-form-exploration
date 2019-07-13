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



# create a dataframe with positions at times prescribed
df = create_df(times,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)




# save plots for different times
dirRes = "plot1"
#dirRes = paste0("plot_","v1=",v1,"_v2=",v2,"_v3=",v3)
dir.create(dirRes)
#times = seq(0,4, by = 0.02)
#savePlot_index(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
savePlot_index_and_paramsInFile(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)

# command line in linux to create a video from images with index (adapt the framerate)
# ffmpeg -framerate 1/0.02 -i plot_%01d.png -crf 15  output1.mp4


# save current parameters
# dirRes = "params"
# save_params_in_file(dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)








# time
timeStep = 0.005 # sec
times = seq(0,5, by = timeStep)

# param de reférence
df = create_df(times,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
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

df = create_df(times,v1,0.3,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,0,angleIni_F)


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


v1*r1B/alphaH - (v1+2*v2)*r2D


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





##########################################
###   Dépendance condition initiales   ###
##########################################








##  ggplot / animate


p <- ggplot(df, aes(x = Bx, y=By) ) +
  geom_point(show.legend = FALSE) 
p

p + transition_time(time) 
  
  


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

  
  



# test 
tt = seq(0,200,by=0.001)
tt = seq(0,6,by=0.001)
a=2
b=1
u=1 
u=0.5
u=1.5
epsilon = 0.012
v=0.5 
v=0.5 + epsilon
v= 3.4
  
x= a*cos(2*pi*tt*u+angleIni_B) + b*cos(2*pi*tt*(u+v)+angleIni_B+angleIni_D) 
y= a*sin(2*pi*tt*u+angleIni_B) + b*sin(2*pi*tt*(u+v)+angleIni_B+angleIni_D) 
plot(x,y, type="l")  
  

a*u - b*(u+v)


derivX = -2*pi*u*a*sin(2*pi*tt*a*u+angleIni_B) -2*pi*(u+v)*b*cos(2*pi*tt*(u+v)+angleIni_B+angleIni_D)
plot(tt,derivX, type = "l")


plot(tt,cos(0.05*tt)*cos(5*tt), type = "l")



########################################
###     Persistence rétinienne       ###
########################################
t=1
plot_oneTime_persistence(t,persistenceTime,nb_point_persistence,10*v1, 10*v2, 10*v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F)











