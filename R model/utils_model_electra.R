## 

# libraries
library(RColorBrewer)
library(tibble)
library(tidyverse)
library(plyr)
library(dplyr)
library(ggplot2)
library(gganimate)
library(gifski)
theme_set(theme_bw())



# Positions des loupiottes en fonction du temps
positions_t_NONsymetrique = function(t,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F){
  # B,C
  Bx = r1B * cos(2*pi*v1*t + angleIni_B)    
  By = r1B * sin(2*pi*v1*t + angleIni_B)   
  Cx = r1C * cos(2*pi*v1*t + angleIni_B+pi)    
  Cy = r1C * sin(2*pi*v1*t + angleIni_B+pi)
  
  # H,I
  Hx = Bx / alphaH 
  Hy = By / alphaH 
  Ix = Cx / alphaI   
  Iy = Cy / alphaI 
  # position de H sur AB: alpha in [1,+infty], =1 => B, =2 milieu de AB, =+infty => A
  
  # D,E
  #r1H = r1B / alphaH 
  HDx = r2D  * cos(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
  HDy = r2D  * sin(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
  HEx = r2E  * cos(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +pi)
  HEy = r2E  * sin(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +pi)
  
  Dx = Hx + HDx
  Dy = Hy + HDy
  Ex = Hx + HEx
  Ey = Hy + HEy
  
  # F,G
  #r1I = r1C / alphaI 
  IFx = r3F  * cos(2*pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+pi)) )
  IFy = r3F  * sin(2*pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+pi)))
  IGx = r3G  * cos(2*pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+pi)) +pi)
  IGy = r3G  * sin(2*pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+pi)) +pi)
  
  Fx = Ix + IFx
  Fy = Iy + IFy
  Gx = Ix + IGx
  Gy = Iy + IGy
  
  
  #################################
  # speed 
  # B,C
  speedBx = -r1B * 2*pi*v1* sin(2*pi*v1*t + angleIni_B)    
  speedBy = r1B * 2*pi*v1* cos(2*pi*v1*t + angleIni_B)   
  speedCx = -r1C * 2*pi*v1* sin(2*pi*v1*t + angleIni_B+pi)    
  speedCy = r1C * 2*pi*v1* cos(2*pi*v1*t + angleIni_B+pi)
  
  # D,E
  speedDx = -r1B/ alphaH * 2*pi*v1* sin(2*pi*v1*t + angleIni_B)    - r2D *2*pi*(v2+v1)* sin(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
  speedDy = r1B/ alphaH * 2*pi*v1*  cos(2*pi*v1*t + angleIni_B)    + r2D *2*pi*(v2+v1)* cos(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
  speedEx = -r1B/ alphaH * 2*pi*v1* sin(2*pi*v1*t + angleIni_B) - r2E *2*pi*(v2+v1)* sin(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +pi)
  speedEy = r1B/ alphaH * 2*pi*v1*  cos(2*pi*v1*t + angleIni_B) + r2E *2*pi*(v2+v1)* cos(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +pi)
  
  # F,G
  speedFx = -r1C/ alphaI * (2*pi*v1)* sin(2*pi*v1*t + angleIni_B+pi)    - r3F* (2*pi*(v3+v1))* sin(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi))
  speedFy = r1C/ alphaI * (2*pi*v1)* cos(2*pi*v1*t + angleIni_B+pi)     + r3F* (2*pi*(v3+v1))* cos(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi))
  speedGx = -r1C/ alphaI * (2*pi*v1)* sin(2*pi*v1*t + angleIni_B+pi) - r3G* (2*pi*(v3+v1))* sin(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi) +pi)
  speedGy = r1C/ alphaI * (2*pi*v1)* cos(2*pi*v1*t + angleIni_B+pi)  + r3G* (2*pi*(v3+v1))* cos(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi) +pi)
  

  
  # acceleration
  # B,C
  accBx = -r1B * (2*pi*v1)^2 * cos(2*pi*v1*t + angleIni_B)    
  accBy = -r1B * (2*pi*v1)^2 * sin(2*pi*v1*t + angleIni_B)   
  accCx = -r1C * (2*pi*v1)^2 * cos(2*pi*v1*t + angleIni_B+pi)    
  accCy = -r1C * (2*pi*v1)^2 * sin(2*pi*v1*t + angleIni_B+pi)
  
  # D,E
  accDx = -r1B/ alphaH * (2*pi*v1)^2* cos(2*pi*v1*t + angleIni_B)    - r2D* (2*pi*(v2+v1))^2* cos(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
  accDy = -r1B/ alphaH * (2*pi*v1)^2* sin(2*pi*v1*t + angleIni_B)    - r2D* (2*pi*(v2+v1))^2* sin(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
  accEx = -r1B/ alphaH * (2*pi*v1)^2* cos(2*pi*v1*t + angleIni_B) - r2E* (2*pi*(v2+v1))^2* cos(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +pi)
  accEy = -r1B/ alphaH * (2*pi*v1)^2* sin(2*pi*v1*t + angleIni_B) - r2E* (2*pi*(v2+v1))^2* sin(2*pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +pi)
  
  # F,G
  accFx = -r1C/ alphaI * (2*pi*v1)^2* cos(2*pi*v1*t + angleIni_B+pi)    - r3F* (2*pi*(v3+v1))^2* cos(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi))
  accFy = -r1C/ alphaI * (2*pi*v1)^2* sin(2*pi*v1*t + angleIni_B+pi)    - r3F* (2*pi*(v3+v1))^2* sin(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi))
  accGx = -r1C/ alphaI * (2*pi*v1)^2* cos(2*pi*v1*t + angleIni_B+pi) - r3G* (2*pi*(v3+v1))^2* cos(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi) +pi)
  accGy = -r1C/ alphaI * (2*pi*v1)^2* sin(2*pi*v1*t + angleIni_B+pi) - r3G* (2*pi*(v3+v1))^2* sin(2*pi*(v3+v1)*t + (angleIni_F+ angleIni_B+pi) +pi)
  
  #return(c(Bx,By,Cx,Cy,Hx,Hy,Dx,Dy,Ex,Ey,Ix,Iy,Fx,Fy,Gx,Gy,time))
  return(data.frame(Bx = Bx, By = By, Cx = Cx, Cy = Cy,
             Hx = Hx, Hy = Hy,
             Dx = Dx, Dy = Dy, Ex = Ex, Ey = Ey,
             Ix = Ix, Iy = Iy,
             Fx = Fx, Fy = Fy, Gx = Gx, Gy = Gy,
             time = t,
             speedBx = speedBx, speedBy = speedBy,
             speedCx = speedCx, speedCy = speedCy,
             speedDx = speedDx, speedDy = speedDy,
             speedEx = speedEx, speedEy = speedEy,
             speedFx = speedFx, speedFy = speedFx,
             speedGx = speedGx, speedGy = speedGy,
             accBx = accBx, accBy = accBy,
             accCx = accCx, accCy = accCy,
             accDx = accDx, accDy = accDy,
             accEx = accEx, accEy = accEy,
             accFx = accFx, accFy = accFx,
             accGx = accGx, accGy = accGy)
         )
}

# t=0
# positions_t_NONsymetrique(t,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F)



# plot res at time t
plot_oneTime = function(t,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F){
  res = positions_t_NONsymetrique(t,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F)
  
  rayonMax = max(r1B,r1C) + max(r2D,r2E,r3F,r3G)
  p <- ggplot(res) + labs(title= paste0("t = ",t)) +
    geom_point(aes(x=Bx, y=By), col="yellow", size = 3) +
    geom_point(aes(x=Cx, y=Cy), col="yellow", size = 3) +
    geom_point(aes(x=Dx, y=Dy), col="yellow", size = 3) +
    geom_point(aes(x=Ex, y=Ey), col="yellow", size = 3) +
    geom_point(aes(x=Fx, y=Fy), col="yellow", size = 3) +
    geom_point(aes(x=Gx, y=Gy), col="yellow", size = 3) #+
    # geom_point(aes(x=Hx, y=Hy), col="black")  +
    # geom_point(aes(x=Ix, y=Iy), col="black")
  #p + xlim(rayonMax, rayonMax) + ylim(rayonMax,rayonMax) 
  p + expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  coord_fixed(ratio=1) + 
  theme(axis.text.y=element_blank(), axis.text.x=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
      #panel.background=element_blank(),
      #panel.border=element_blank(), 
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())
  }

# t=0 # initial conditions
# t=1.45
# plot_oneTime(t,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F)


# create data frame for persistence rétinienne
create_df_persistence = function(t,persistenceTime,nb_point_persistence,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F){
  df = data.frame(Bx = double(), By = double(), Cx = double(), Cy = double(),
                  Hx = double(), Hy = double(),
                  Dx = double(), Dy = double(), Ex = double(), Ey = double(),
                  Ix = double(),
                  Iy = double(),
                  Fx = double(), Fy = double(), Gx = double(), Gy = double(),
                  time = double())
  #times = seq(0,4, by = 0.005)
  
  persistenceTimes = seq(max(0,t-persistenceTime),t,length.out = nb_point_persistence)
  for (current_t in persistenceTimes){
    res = positions_t_NONsymetrique(current_t,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
    # temp_df = data.frame(Bx = res[1], By = res[2], Cx = res[3], Cy = res[4],
    #                      Hx = res[5], Hy = res[6],
    #                      Dx = res[7], Dy = res[8], Ex = res[9], Ey = res[10],
    #                      Ix = res[11], Iy = res[12],
    #                      Fx = res[13], Fy = res[14], Gx = res[15], Gy = res[16],
    #                      time = res[17])
    # df = rbind.fill(df,temp_df)
    df = rbind.fill(df,res)
  }
  return(df)
}


# df = create_df_persistence(t,persistenceTime,nb_point_persistence,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)




# plot res at time t with Persitence rétinienne
# i.e plot the positions of the points during the 50 previous millisecond
plot_oneTime_persistence = function(t,persistenceTime,nb_point_persistence,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F){
  res = create_df_persistence(t,persistenceTime,nb_point_persistence,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
  
  rayonMax = max(r1B,r1C) + max(r2D,r2E,r3F,r3G)
  p <- ggplot(res) + labs(title= paste0("t = ",t, ", persistence = ",persistenceTime, " s")) +
    geom_point(aes(x=Bx, y=By), col="yellow", size = 3) +
    geom_point(aes(x=Cx, y=Cy), col="yellow", size = 3) +
    geom_point(aes(x=Dx, y=Dy), col="yellow", size = 3) +
    geom_point(aes(x=Ex, y=Ey), col="yellow", size = 3) +
    geom_point(aes(x=Fx, y=Fy), col="yellow", size = 3) +
    geom_point(aes(x=Gx, y=Gy), col="yellow", size = 3) #+
  # geom_point(aes(x=Hx, y=Hy), col="black")  +
  # geom_point(aes(x=Ix, y=Iy), col="black")
  #p + xlim(rayonMax, rayonMax) + ylim(rayonMax,rayonMax) 
  p + expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    coord_fixed(ratio=1) + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
          #panel.background=element_blank(),
          #panel.border=element_blank(), 
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
}


# plot_oneTime_persistence(t,persistenceTime,nb_point_persistence,v1, v2, v3, r1B, r1C, r2D, r2E, r3F, r3G, alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F)



# save the result (%time) in a data frame
# function to create data frame
create_df = function(times,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F){
  df = data.frame(Bx = double(), By = double(), Cx = double(), Cy = double(),
                  Hx = double(), Hy = double(),
                  Dx = double(), Dy = double(), Ex = double(), Ey = double(),
                  Ix = double(),
                  Iy = double(),
                  Fx = double(), Fy = double(), Gx = double(), Gy = double(),
                  time = double(),
                  speedBx = double(), speedBy = double(),
                  speedCx = double(), speedCy = double(),
                  speedDx = double(), speedDy = double(),
                  speedEx = double(), speedEy = double(),
                  speedFx = double(), speedFy = double(),
                  speedGx = double(), speedGy = double(),
                  accBx = double(), accBy = double(),
                  accCx = double(), accCy = double(),
                  accDx = double(), accDy = double(),
                  accEx = double(), accEy = double(),
                  accFx = double(), accFy = double(),
                  accGx = double(), accGy = double())
  #times = seq(0,4, by = 0.005)
  for (i in 1:length(times)){
    t = times[i]
    res = positions_t_NONsymetrique(t,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
    # temp_df = data.frame(Bx = res$Bx, By = res$By, Cx = res$Cx, Cy = res$Cy,
    #                      Hx = res$Hx, Hy = res$Hy,
    #                      Dx = res$Dx, Dy = res$Dy, Ex = res$Ex, Ey = res$Ey,
    #                      Ix = res$Ix, Iy = res$Iy,
    #                      Fx = res$Fx, Fy = res$Fy, Gx = res$Gx, Gy = res$Gy,
    #                      speedBx = res$speedBx, speedBy = res$speedBy,
    #                      speedCx = res$speedCx, speedCy = res$speedCy,
    #                      speedDx = res$speedDx, speedDy = res$speedDy,
    #                      speedEx = res$speedEx, speedEy = res$speedEy,
    #                      speedFx = res$speedFx, speedFy = res$speedFx,
    #                      speedGx = res$speedGx, speedGy = res$speedGy,
    #                      accBx = res$accBx, accBy = res$accBy,
    #                      accCx = res$accCx, accCy = res$accCy,
    #                      accDx = res$accDx, accDy = res$accDy,
    #                      accEx = res$accEx, accEy = res$accEy,
    #                      accFx = res$accFx, accFy = res$accFx,
    #                      accGx = res$accGx, accGy = res$accGy,
    #                      time = res$time)
    # df = rbind.fill(df,temp_df)
    df = rbind.fill(df,res)
  }
  return(df)
}

# df = create_df(times,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)



# function to save all plots (%time)
savePlot_index = function(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F){
  df = create_df(times,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
  rayonMax = max(r1B,r1C) + max(r2D,r2E,r3F,r3G)
  for (i in 1:length(times)){
  t=times[i]
  p <- ggplot(df %>% filter(times == t)) + #labs(title= paste0("t = ",t)) +
    geom_point(aes(x=Bx, y=By), col="yellow", size = 3) +
    geom_point(aes(x=Cx, y=Cy), col="yellow", size = 3) +
    geom_point(aes(x=Dx, y=Dy), col="yellow", size = 3) +
    geom_point(aes(x=Ex, y=Ey), col="yellow", size = 3) +
    geom_point(aes(x=Fx, y=Fy), col="yellow", size = 3) +
    geom_point(aes(x=Gx, y=Gy), col="yellow", size = 3) #+
  #geom_point(aes(x=Hx, y=Hy), col="black")  #+
  #geom_point(aes(x=Ix, y=Iy), col="black")
  #p + xlim(rayonMax, rayonMax) + ylim(rayonMax,rayonMax) 
  p <- p + expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    coord_fixed(ratio=1) + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
          #panel.background=element_blank(),
          #panel.border=element_blank(), 
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  # save ggplot
  plot_file_name = paste0(dirRes, "/plot_",i,".svg")
  ggsave(plot_file_name,p, width = 8, height = 8, dpi = 72)
  #dev.off()
  }
}



# save plot
# dirRes = "plotsIndex6"
# dir.create(dirRes)
# times = seq(0,4, by = 0.02)
# savePlot_index(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)



# write parameters used in a txt file: parameters.txt in the directory specified   
save_params_in_file = function(dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F){
  
  # sequence of string with parameter names and values  
  txt <- c( 
    "# vitesse moteur (tour/ seconde)",
    paste0("v1 = ", v1),
    paste0("v2 = ", v2),
    paste0("v3 = ", v3),
    "# rayons",
    paste0("r1B = ", r1B, ", r1C = ", r1C),
    paste0("r2D = ", r2D, ", r2E = ", r2E),
    paste0("r3F = ", r3F, ", r3G = ", r3G),
    "# position initiale (from angle)",
    "A,B,C",
    paste0("angleIni_B = ", angleIni_B),
    "# H,I",
    paste0("alphaH = ", alphaH),
    paste0("alphaI = ", alphaI),
    "# D,E",
    paste0("angleIni_D = ", angleIni_D),
    "# F,G", 
    paste0("angleIni_F = ", angleIni_F))
  
  # save in txt file  
  dir.create(dirRes)
  fileName = "parameters.txt"
  file = paste0(dirRes, "/", fileName)
  writeLines(txt, file)
}

# dirRes = "toto"
# save_params_in_file(dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)


# save plots with index (time) and parameters in a file
savePlot_index_and_paramsInFile = function(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F){
  savePlot_index(times,dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
  save_params_in_file(dirRes,v1,v2,v3,r1B,r1C,r2D,r2E,r3F,r3G,alphaH,alphaI,angleIni_B,angleIni_D,angleIni_F)
}


# le deuxième vecteur est normalisé
scalarProdutRenormalized = function(x1,y1,x2,y2){
  (x1*x2+y1*y2) / sqrt(x2^2+y2^2)
}




# renvoit la vitesse du moteur en fonction de la vitesse demandée
vitesse_rotation_moteur <- function(times,u,A,C) {
  temp = function(t){
    A/C*u(t)*exp(1/C*t)
  }
  
  temp2 = unlist(sapply(times, function(t){integrate(temp,0,t)$value}))
  vitesse = exp(-1/C*times) * temp2
}







# TRANSITOIRE (LES VITESSES CHANGENT) Positions des loupiottes en fonction du temps 
positions_t_TRANSITOIRE = function(t, r1B, r1C, r2D, r2E, r3F, r3G, 
                                   alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                   A1,C1,A2,C2,A3,C3,
                                   u1,u2,u3){
  
  vitesse1 = vitesse_rotation_moteur(t,u1,A1,C1) 
  vitesse2 = vitesse_rotation_moteur(t,u2,A2,C2) 
  vitesse3 = vitesse_rotation_moteur(t,u3,A3,C3) 
  
  # B,C
  Bx = r1B * cos(2*pi*vitesse1*t + angleIni_B)    
  By = r1B * sin(2*pi*vitesse1*t + angleIni_B)   
  Cx = r1C * cos(2*pi*vitesse1*t + angleIni_B+pi)    
  Cy = r1C * sin(2*pi*vitesse1*t + angleIni_B+pi)
  
  # H,I
  Hx = Bx / alphaH 
  Hy = By / alphaH 
  Ix = Cx / alphaI   
  Iy = Cy / alphaI 
  # position de H sur AB: alpha in [1,+infty], =1 => B, =2 milieu de AB, =+infty => A
  
  # D,E
  #r1H = r1B / alphaH 
  HDx = r2D  * cos(2*pi*(vitesse2+vitesse1)*t + (angleIni_D+ angleIni_B))
  HDy = r2D  * sin(2*pi*(vitesse2+vitesse1)*t + (angleIni_D+ angleIni_B))
  HEx = r2E  * cos(2*pi*(vitesse2+vitesse1)*t + (angleIni_D+ angleIni_B) +pi)
  HEy = r2E  * sin(2*pi*(vitesse2+vitesse1)*t + (angleIni_D+ angleIni_B) +pi)
  
  Dx = Hx + HDx
  Dy = Hy + HDy
  Ex = Hx + HEx
  Ey = Hy + HEy
  
  # F,G
  #r1I = r1C / alphaI 
  IFx = r3F  * cos(2*pi*(vitesse3+vitesse1)*t + (angleIni_F+ (angleIni_B+pi)) )
  IFy = r3F  * sin(2*pi*(vitesse3+vitesse1)*t + (angleIni_F+ (angleIni_B+pi)))
  IGx = r3G  * cos(2*pi*(vitesse3+vitesse1)*t + (angleIni_F+ (angleIni_B+pi)) +pi)
  IGy = r3G  * sin(2*pi*(vitesse3+vitesse1)*t + (angleIni_F+ (angleIni_B+pi)) +pi)
  
  Fx = Ix + IFx
  Fy = Iy + IFy
  Gx = Ix + IGx
  Gy = Iy + IGy
  
  
  return(data.frame(Bx = Bx, By = By, Cx = Cx, Cy = Cy,
                    Hx = Hx, Hy = Hy,
                    Dx = Dx, Dy = Dy, Ex = Ex, Ey = Ey,
                    Ix = Ix, Iy = Iy,
                    Fx = Fx, Fy = Fy, Gx = Gx, Gy = Gy,
                    time = t))

}







# plot res at time t TRANSITOIRE
plot_oneTime_TRANSITOIRE = function(t, r1B, r1C, r2D, r2E, r3F, r3G, 
                                    alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                    A1,C1,A2,C2,A3,C3,
                                    u1,u2,u3){
  
  res = positions_t_TRANSITOIRE(t, r1B, r1C, r2D, r2E, r3F, r3G, 
                                     alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                     A1,C1,A2,C2,A3,C3,
                                     u1,u2,u3)
  
  rayonMax = max(r1B,r1C) + max(r2D,r2E,r3F,r3G)
  p <- ggplot(res) + labs(title= paste0("t = ",t)) +
    geom_point(aes(x=Bx, y=By), col="yellow", size = 3) +
    geom_point(aes(x=Cx, y=Cy), col="yellow", size = 3) +
    geom_point(aes(x=Dx, y=Dy), col="yellow", size = 3) +
    geom_point(aes(x=Ex, y=Ey), col="yellow", size = 3) +
    geom_point(aes(x=Fx, y=Fy), col="yellow", size = 3) +
    geom_point(aes(x=Gx, y=Gy), col="yellow", size = 3) #+
  # geom_point(aes(x=Hx, y=Hy), col="black")  +
  # geom_point(aes(x=Ix, y=Iy), col="black")
  #p + xlim(rayonMax, rayonMax) + ylim(rayonMax,rayonMax) 
  p + expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    coord_fixed(ratio=1) + 
    theme(axis.text.y=element_blank(), axis.text.x=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
          #panel.background=element_blank(),
          #panel.border=element_blank(), 
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
}




create_df_TRANSITOIRE = function(times, r1B, r1C, r2D, r2E, r3F, r3G, 
                                 alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                 A1,C1,A2,C2,A3,C3,
                                 u1,u2,u3){
  
  df = data.frame(Bx = double(), By = double(), Cx = double(), Cy = double(),
                  Hx = double(), Hy = double(),
                  Dx = double(), Dy = double(), Ex = double(), Ey = double(),
                  Ix = double(),
                  Iy = double(),
                  Fx = double(), Fy = double(), Gx = double(), Gy = double(),
                  time = double() #,
                  # speedBx = double(), speedBy = double(),
                  # speedCx = double(), speedCy = double(),
                  # speedDx = double(), speedDy = double(),
                  # speedEx = double(), speedEy = double(),
                  # speedFx = double(), speedFy = double(),
                  # speedGx = double(), speedGy = double(),
                  # accBx = double(), accBy = double(),
                  # accCx = double(), accCy = double(),
                  # accDx = double(), accDy = double(),
                  # accEx = double(), accEy = double(),
                  # accFx = double(), accFy = double(),
                  # accGx = double(), accGy = double()
                  )
  #times = seq(0,4, by = 0.005)
  for (i in 1:length(times)){
    t = times[i]
    res = positions_t_TRANSITOIRE(t, r1B, r1C, r2D, r2E, r3F, r3G, 
                                    alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                    A1,C1,A2,C2,A3,C3,
                                    u1,u2,u3)
    df = rbind.fill(df,res)
  }
  return(df)
}



# function to save all plots (%time)
savePlot_index_TRANSITOIRE = function(times, dirRes, r1B, r1C, r2D, r2E, r3F, r3G, 
                                      alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                      A1,C1,A2,C2,A3,C3,
                                      u1,u2,u3){
  
  df = create_df_TRANSITOIRE(times, r1B, r1C, r2D, r2E, r3F, r3G, 
                 alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                 A1,C1,A2,C2,A3,C3,
                 u1,u2,u3)
    
  rayonMax = max(r1B,r1C) + max(r2D,r2E,r3F,r3G)
  for (i in 1:length(times)){
    t=times[i]
    p <- ggplot(df %>% filter(times == t)) + #labs(title= paste0("t = ",t)) +
      geom_point(aes(x=Bx, y=By), col="yellow", size = 3) +
      geom_point(aes(x=Cx, y=Cy), col="yellow", size = 3) +
      geom_point(aes(x=Dx, y=Dy), col="yellow", size = 3) +
      geom_point(aes(x=Ex, y=Ey), col="yellow", size = 3) +
      geom_point(aes(x=Fx, y=Fy), col="yellow", size = 3) +
      geom_point(aes(x=Gx, y=Gy), col="yellow", size = 3) #+
    #geom_point(aes(x=Hx, y=Hy), col="black")  #+
    #geom_point(aes(x=Ix, y=Iy), col="black")
    #p + xlim(rayonMax, rayonMax) + ylim(rayonMax,rayonMax) 
    p <- p + expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
      coord_fixed(ratio=1) + 
      theme(axis.text.y=element_blank(), axis.text.x=element_blank(),
            axis.ticks=element_blank(),
            axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
            #panel.background=element_blank(),
            #panel.border=element_blank(), 
            panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            plot.background=element_blank())
    # save ggplot
    plot_file_name = paste0(dirRes, "/plot_",i,".svg")
    ggsave(plot_file_name,p, width = 8, height = 8, dpi = 72)
    #dev.off()
  }
}




save_params_in_file_TRANSITOIRE = function(dirRes, r1B, r1C, r2D, r2E, r3F, r3G, 
                                           alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                           A1,C1,A2,C2,A3,C3,
                                           u1,u2,u3){
  
  # sequence of string with parameter names and values  
  txt <- c( 
      "# rayons",
    paste0("r1B = ", r1B, ", r1C = ", r1C),
    paste0("r2D = ", r2D, ", r2E = ", r2E),
    paste0("r3F = ", r3F, ", r3G = ", r3G),
    "# position initiale (from angle)",
    "A,B,C",
    paste0("angleIni_B = ", angleIni_B),
    "# H,I",
    paste0("alphaH = ", alphaH),
    paste0("alphaI = ", alphaI),
    "# D,E",
    paste0("angleIni_D = ", angleIni_D),
    "# F,G", 
    paste0("angleIni_F = ", angleIni_F),
    "# paramètres moteurs",
    "# moteur 1",
    paste0("A1 = ", A1),
    paste0("C1 = ", C1),
    "# moteur 2",
    paste0("A2 = ", A2),
    paste0("C2 = ", C2),
    "# moteur 3",
    paste0("A3 = ", A3),
    paste0("C3 = ", C3),
    "# vitesse moteur (tour/ seconde)",
    paste0("fonction vitesse rotation moteur 1: ", paste0(deparse(u1), collapse = " ")),
    paste0("fonction vitesse rotation moteur 2: ", paste0(deparse(u2), collapse = " ")),
    paste0("fonction vitesse rotation moteur 3: ", paste0(deparse(u3), collapse = " "))
    )
  
  # save in txt file  
  dir.create(dirRes)
  fileName = "parameters.txt"
  file = paste0(dirRes, "/", fileName)
  writeLines(txt, file)
}


# dirRes = "toto"
# save_params_in_file_TRANSITOIRE(dirRes, r1B, r1C, r2D, r2E, r3F, r3G, 
#                                 alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
#                                 A1,C1,A2,C2,A3,C3,
#                                 u1,u2,u3)




# save plots with index (time) and parameters in a file  TRANSITOIRE
savePlot_index_and_paramsInFile_TRANSITOIRE = function(times, dirRes, r1B, r1C, r2D, r2E, r3F, r3G, 
                                                       alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                                       A1,C1,A2,C2,A3,C3,
                                                       u1,u2,u3){

  savePlot_index_TRANSITOIRE(times, dirRes, r1B, r1C, r2D, r2E, r3F, r3G, 
                                              alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                              A1,C1,A2,C2,A3,C3,
                                              u1,u2,u3) 
  
  save_params_in_file_TRANSITOIRE(dirRes, r1B, r1C, r2D, r2E, r3F, r3G, 
                                   alphaH, alphaI, angleIni_B, angleIni_D, angleIni_F,
                                   A1,C1,A2,C2,A3,C3,
                                   u1,u2,u3)
  }


