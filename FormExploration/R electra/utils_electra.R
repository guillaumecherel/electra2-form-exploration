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
positions_t_Stationnaire = function(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F){
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
# pour faire les vidéos
plot_oneTime = function(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F){
  res = positions_t_Stationnaire(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
  
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



# avec titre, grille, axes
plotTrajectory <- function(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,
                                lightB,lightC,lightD,lightE,lightF,lightG) {
  theme_set(theme_bw())
  df = positions_t_Stationnaire(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
  rayonMax = sqrt(df$Bx^2+df$By^2)[1]
  
  p <- ggplot(df) +
    expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    coord_fixed(ratio=1) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) 
  
    if(lightB==T){p = p+geom_point(aes(x=Bx, y=By), col= "yellow", size=0.5)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightC==T){p = p+geom_point(aes(x=Cx, y=Cy), col= "yellow", size=0.5)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightD==T){p = p+geom_point(aes(x=Dx, y=Dy), col= "yellow", size=0.5)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightE==T){p = p+geom_point(aes(x=Ex, y=Ey), col= "yellow", size=0.5)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightF==T){p = p+geom_point(aes(x=Fx, y=Fy), col= "yellow", size=0.5)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightG==T){p = p+geom_point(aes(x=Gx, y=Gy), col= "yellow", size=0.5)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    
    p = p+labs(title = "Trajectoire Stationnaire",
         subtitle = paste0("v1=",round(v1,2),", v2=",round(v2,2),", v3=",round(v3,2) )) #,", \n",
                           #"angleIni_B=",angleIni_B,", angleIni_D=",angleIni_D,", angleIni_F=",angleIni_F))
  
  p
}



# jaune sur fond noir
plotTrajectorySmart <- function(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,
                           lightB,lightC,lightD,lightE,lightF,lightG) {
  theme_set(theme_bw())
  df = positions_t_Stationnaire(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
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
          panel.background = element_rect(fill = bg, colour = bg),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
          ) 
  
  if(lightB==T){p = p+geom_point(aes(x=Bx, y=By), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightC==T){p = p+geom_point(aes(x=Cx, y=Cy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightD==T){p = p+geom_point(aes(x=Dx, y=Dy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightE==T){p = p+geom_point(aes(x=Ex, y=Ey), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightF==T){p = p+geom_point(aes(x=Fx, y=Fy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightG==T){p = p+geom_point(aes(x=Gx, y=Gy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  
  p = p+labs(subtitle = paste0("v1=",round(v1,2),", v2=",round(v2,2),", v3=",round(v3,2))) 
  
  p
}







# jaune sur fond noir, ecriture blanche
plotTrajectorySmart2 <- function(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,
                                lightB,lightC,lightD,lightE,lightF,lightG) {
  theme_set(theme_bw())
  df = positions_t_Stationnaire(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
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
          panel.background = element_rect(fill = bg, colour = bg),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank()
    ) 
  
  if(lightB==T){p = p+geom_path(aes(x=Bx, y=By), colour = yellow, alpha=1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightC==T){p = p+geom_path(aes(x=Cx, y=Cy), colour = yellow, alpha=1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightD==T){p = p+geom_path(aes(x=Dx, y=Dy), colour = yellow, alpha=1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightE==T){p = p+geom_path(aes(x=Ex, y=Ey), colour = yellow, alpha=1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightF==T){p = p+geom_path(aes(x=Fx, y=Fy), colour = yellow, alpha=1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightG==T){p = p+geom_path(aes(x=Gx, y=Gy), colour = yellow, alpha=1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  
  #p = p+labs(subtitle = paste0("v1=",round(v1,2),", v2=",round(v2,2),", v3=",round(v3,2)), color = "white") + theme(plot.subtitle=element_text(color="white"))
  #theme(plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"))
  p
}

























