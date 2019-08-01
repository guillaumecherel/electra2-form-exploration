## 

# libraries
library(RColorBrewer)
library(tibble)
library(tidyverse)
#library(plyr)
library(dplyr)
library(ggplot2)
library(gganimate)
library(gifski)
library(pracma)
theme_set(theme_bw())





# Positions des points de lumiere en fonction du temps
positions_t_Stationnaire_v2 = function(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F){
# le point H n'est plus sur le segment AC
  # B,C
  Bx = rB * cos(2*pi*v1*t + angleIni_B)    
  By = rB * sin(2*pi*v1*t + angleIni_B)   
  Cx = rC * cos(2*pi*v1*t + angleIni_B +pi)    
  Cy = rC * sin(2*pi*v1*t + angleIni_B +pi)
  
  # H,I
  Hx = rH * cos(2*pi*v1*t + angleIni_B - angleH)    
  Hy = rH * sin(2*pi*v1*t + angleIni_B - angleH)   
  Ix = rI * cos(2*pi*v1*t + angleIni_B + angleI +pi)    
  Iy = rI * sin(2*pi*v1*t + angleIni_B + angleI +pi)

  
  # D,E
  # on ne calcule pas HD en fait, le nom est mal choisit
  HDx = rD  * cos(2*pi*(v2+v1)*t + (angleIni_D + angleIni_B - angleH))
  HDy = rD  * sin(2*pi*(v2+v1)*t + (angleIni_D + angleIni_B - angleH))
  HEx = rE  * cos(2*pi*(v2+v1)*t + (angleIni_D + angleIni_B - angleH) +pi)
  HEy = rE  * sin(2*pi*(v2+v1)*t + (angleIni_D + angleIni_B - angleH) +pi)
  
  Dx = Hx + HDx
  Dy = Hy + HDy
  Ex = Hx + HEx
  Ey = Hy + HEy
  
  # F,G
  IFx = rF  * cos(2*pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+pi) + angleI))
  IFy = rF  * sin(2*pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+pi) + angleI))
  IGx = rG  * cos(2*pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+pi) + angleI) +pi)
  IGy = rG  * sin(2*pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+pi) + angleI) +pi)
  
  Fx = Ix + IFx
  Fy = Iy + IFy
  Gx = Ix + IGx
  Gy = Iy + IGy
  
  
  # speed
  speedBx = -rB * 2*pi*v1* sin(2*pi*v1*t + angleIni_B)
  speedBy = rB * 2*pi*v1* cos(2*pi*v1*t + angleIni_B)
  speedCx = -rC * 2*pi*v1* sin(2*pi*v1*t + angleIni_B+pi)
  speedCy = rC * 2*pi*v1* cos(2*pi*v1*t + angleIni_B+pi)
  
  # D,E
  speedDx = -rH * 2*pi*v1* sin(2*pi*v1*t + angleIni_B -angleH)  - rD *2*pi*(v2+v1)* sin(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH))
  speedDy = rH * 2*pi*v1*  cos(2*pi*v1*t + angleIni_B -angleH)  + rD *2*pi*(v2+v1)* cos(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH))
  speedEx = -rH * 2*pi*v1* sin(2*pi*v1*t + angleIni_B -angleH)  - rE *2*pi*(v2+v1)* sin(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH) +pi)
  speedEy = rH * 2*pi*v1*  cos(2*pi*v1*t + angleIni_B -angleH)  + rE *2*pi*(v2+v1)* cos(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH) +pi)
  
  # F,G
  speedFx = -rI * (2*pi*v1)* sin(2*pi*v1*t + angleIni_B+pi +angleI)   - rF* (2*pi*(v3+v1))* sin(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi +angleI))
  speedFy = rI * (2*pi*v1)* cos(2*pi*v1*t + angleIni_B+pi +angleI)    + rF* (2*pi*(v3+v1))* cos(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi +angleI))
  speedGx = -rI * (2*pi*v1)* sin(2*pi*v1*t + angleIni_B+pi +angleI)   - rG* (2*pi*(v3+v1))* sin(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi +angleI) +pi)
  speedGy = rI * (2*pi*v1)* cos(2*pi*v1*t + angleIni_B+pi +angleI)    + rG* (2*pi*(v3+v1))* cos(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi +angleI) +pi)  
  
  
  #return(c(Bx,By,Cx,Cy,Hx,Hy,Dx,Dy,Ex,Ey,Ix,Iy,Fx,Fy,Gx,Gy,time))
  data.frame(Bx = Bx, By = By, Cx = Cx, Cy = Cy,
                    Hx = Hx, Hy = Hy,
                    Dx = Dx, Dy = Dy, Ex = Ex, Ey = Ey,
                    Ix = Ix, Iy = Iy,
                    Fx = Fx, Fy = Fy, Gx = Gx, Gy = Gy,
                    time = t, 
                    speedBx = speedBx, speedBy = speedBy,
                    speedCx = speedCx, speedCy = speedCy,
                    speedDx = speedDx, speedDy = speedDy,
                    speedEx = speedEx, speedEy = speedEy,
                    speedFx = speedFx, speedFy = speedFy,
                    speedGx = speedGx, speedGy = speedGy) 
}



# plot res at time t
plot_oneTime_v2 = function(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F){
  res = positions_t_Stationnaire_v2(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
  
  rayonMax = max(rB,rC, rH + max(rD,rE), rI + max(rF,rG) )
  p <- ggplot(res) + labs(title= paste0("t = ",t)) +
    geom_point(aes(x=Bx, y=By), col="yellow", size = 3) +
    geom_point(aes(x=Cx, y=Cy), col="yellow", size = 3) +
    geom_point(aes(x=Dx, y=Dy), col="yellow", size = 3) +
    geom_point(aes(x=Ex, y=Ey), col="yellow", size = 3) +
    geom_point(aes(x=Fx, y=Fy), col="yellow", size = 3) +
    geom_point(aes(x=Gx, y=Gy), col="yellow", size = 3) +
    geom_point(aes(x=Hx, y=Hy), col="black")  +
    geom_point(aes(x=Ix, y=Iy), col="black")  +
    geom_point(aes(x=0, y=0), col="black")
  p + expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    coord_fixed(ratio=1) + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
}




# create_df_v2 = function(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F){
#   df = data.frame(Bx = double(), By = double(), Cx = double(), Cy = double(),
#                   Hx = double(), Hy = double(),
#                   Dx = double(), Dy = double(), Ex = double(), Ey = double(),
#                   Ix = double(),
#                   Iy = double(),
#                   Fx = double(), Fy = double(), Gx = double(), Gy = double(),
#                   time = double(),
#                   speedBx = double(), speedBy = double(),
#                   speedCx = double(), speedCy = double(),
#                   speedDx = double(), speedDy = double(),
#                   speedEx = double(), speedEy = double(),
#                   speedFx = double(), speedFy = double(),
#                   speedGx = double(), speedGy = double(),
#                   accBx = double(), accBy = double(),
#                   accCx = double(), accCy = double(),
#                   accDx = double(), accDy = double(),
#                   accEx = double(), accEy = double(),
#                   accFx = double(), accFy = double(),
#                   accGx = double(), accGy = double())
#   
#   for (i in 1:length(times)){
#     t = times[i]
#     res = positions_t_Stationnaire_v2(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
#     df = rbind(df,res)
#   }
#   return(df)
# }



plotTrajectory_v2 <- function(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F) {
  theme_set(theme_bw())
  df = positions_t_Stationnaire_v2(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
  rayonMax = sqrt(df$Bx^2+df$By^2)[1]

  p <- ggplot(df) +
    expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    coord_fixed(ratio=1) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) +
    geom_point(aes(x=Bx, y=By), size=1) +
    geom_point(aes(x=Cx, y=Cy), size=1) +
    geom_point(aes(x=Dx, y=Dy), size=1, col= "red", shape = 3) +
    geom_point(aes(x=Ex, y=Ey), size=0.1, col= "yellow") +
    geom_point(aes(x=Fx, y=Fy), size=1, col= "blue") +
    geom_point(aes(x=Gx, y=Gy), size=0.1, col= "green", shape = 3)+
    labs(title = "Trajectoire Stationnaire",
         subtitle = paste0("v1=",v1,", v2=",v2,", v3=",v3,", \n",
                           "angleIni_B=",angleIni_B,", angleIni_D=",angleIni_D,", angleIni_F=",angleIni_F))

  p
}


distanceSquared = function(v1,v2){
  sqrt((v1[1]-v2[1])^2 +(v1[2]-v2[2])^2)
}


nextSimplifiedTrajectory = function(current,res, distance){
  #current:Vector[(Double,Double)]  => array à deux colonnes
  #res:Vector[(Double,Double)]
  #distance:Double
  if (dim(current)[1]==0) {  #isempty(current)  
    res
  } else {
     temp = current[1,]
     temp2 = apply(current, 1, function(x){distanceSquared(x,temp)})
     temp3 = cumsum(temp2)
     temp4 = which(temp3 > distance)
    if (length(temp4)==0){
      newCurrent = matrix(, nrow = 0, ncol = 2)
      newRes = rbind(res, current[(dim(current)[1]),] )
      nextSimplifiedTrajectory(newCurrent, newRes, distance)
    } else {
      newCurrent =  current[temp4,]
      newRes = rbind(res, current[temp4[1],])
      nextSimplifiedTrajectory(newCurrent, newRes, distance)
    }
  }
}


simplifiedTrajectory = function(mat, distance){
  current = mat
  res =   current[1,]
  nextSimplifiedTrajectory(current,res,distance)
}






















