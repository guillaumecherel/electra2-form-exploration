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




savePlot_index = function(times,dirRes,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,
                          lightB,lightC,lightD,lightE,lightF,lightG){
  
  # df = positions_t_Stationnaire_v2(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
  # rayonMax = max( max(rB,rC) ,  max(rH+rD,rH+rE,rI+rF,rI+rG) ) 
  # for (i in 1:length(times)){
  #   t=times[i]
  #   p <- ggplot(df %>% filter(times == t)) + #labs(title= paste0("t = ",t)) +
  #     geom_point(aes(x=Bx, y=By), col="yellow", size = 3) +
  #     geom_point(aes(x=Cx, y=Cy), col="yellow", size = 3) +
  #     geom_point(aes(x=Dx, y=Dy), col="yellow", size = 3) +
  #     geom_point(aes(x=Ex, y=Ey), col="yellow", size = 3) +
  #     geom_point(aes(x=Fx, y=Fy), col="yellow", size = 3) +
  #     geom_point(aes(x=Gx, y=Gy), col="yellow", size = 3) #+
  #   p <- p + expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
  #     coord_fixed(ratio=1) + 
  #     theme(axis.text.y=element_blank(), axis.text.x=element_blank(),
  #           axis.ticks=element_blank(),
  #           axis.title.x=element_blank(), axis.title.y=element_blank(),legend.position="none",
  #           panel.grid.major=element_blank(),
  #           panel.grid.minor=element_blank(),
  #           plot.background=element_blank())
  #   # save ggplot
  #   plot_file_name = paste0(dirRes, "/plot_",i,".svg")
  #   ggsave(plot_file_name,p, width = 8, height = 8, dpi = 72)
  #   #dev.off()
  #}
  
  theme_set(theme_bw())
  rayonMax = sqrt(df$Bx^2+df$By^2)[1]
  
  for ( t in times){
    # juste le i ieme
    df = positions_t_Stationnaire(t, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
    p <- ggplot(df) +
      expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
      coord_fixed(ratio=1) +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) 
    
    if(lightB==T){p = p+geom_point(aes(x=Bx, y=By), col= "yellow", size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightC==T){p = p+geom_point(aes(x=Cx, y=Cy), col= "yellow", size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightD==T){p = p+geom_point(aes(x=Dx, y=Dy), col= "yellow", size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightE==T){p = p+geom_point(aes(x=Ex, y=Ey), col= "yellow", size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightF==T){p = p+geom_point(aes(x=Fx, y=Fy), col= "yellow", size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    if(lightG==T){p = p+geom_point(aes(x=Gx, y=Gy), col= "yellow", size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= "white", size=0.5)} 
    
    p = p+labs(  # title = "Trajectoire Stationnaire",
               subtitle = paste0("v1=",round(v1,2),", v2=",round(v2,2),", v3=",round(v3,2) )) #,", \n",
    #"angleIni_B=",angleIni_B,", angleIni_D=",angleIni_D,", angleIni_F=",angleIni_F))

    # save ggplot
    plot_file_name = paste0(dirRes, "/plot_",i,".svg")
    ggsave(plot_file_name,p, width = 8, height = 8, dpi = 72)
    #dev.off()
    
  }
}


save_params_in_file = function(dirRes,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,
                               lightB,lightC,lightD,lightE,lightF,lightG){
  
  # sequence of string with parameter names and values  
  txt <- c( 
    "# vitesse moteur (tour/ seconde)",
    paste0("v1 = ", v1),
    paste0("v2 = ", v2),
    paste0("v3 = ", v3),
    "# rayons",
    paste0("rB = ", rB, ", rC = ", rC),
    paste0("rD = ", rD, ", rE = ", rE),
    paste0("rF = ", rF, ", rG = ", rG),
    "# position initiale (from angle)",
    paste0("angleIni_B = ", angleIni_B),
    "# H,I",
    paste0("alphaH = ", alphaH),
    paste0("alphaI = ", alphaI),
    "# D,E",
    paste0("angleIni_D = ", angleIni_D),
    "# F,G", 
    paste0("angleIni_F = ", angleIni_F),
    # lumieres on / of
    paste0("lightB = ", lightB),
    paste0("lightC = ", lightC),
    paste0("lightD = ", lightD),
    paste0("lightE = ", lightE),
    paste0("lightF = ", lightF),
    paste0("lightG = ", lightG),
    )
  
  # save in txt file  
  dir.create(dirRes)
  fileName = "parameters.txt"
  file = paste0(dirRes, "/", fileName)
  writeLines(txt, file)
}



# save plots with index (time) and parameters in a file
savePlot_index_and_paramsInFile = function(times,dirRes,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,
                                           lightB,lightC,lightD,lightE,lightF,lightG){
  savePlot_index(times,dirRes,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)
  save_params_in_file(dirRes,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)
}




# save plots for different times
dirRes = "plotElectra"
dir.create(dirRes)

#savePlot_index(times,dirRes,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)
savePlot_index_and_paramsInFile(times,dirRes,v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)

# command line in linux to create a video from images with index (adapt the framerate)
# ffmpeg -framerate 1/0.005 -i plot_%01d.png -crf 15  output1.mp4






