## 

## libraries
library(ggplot2)

# plot trajectoire stationnaire
plotTrajectory <- function(Bx,By,Cx,Cy,Dx,Dy,Ex,Ey,Fx,Fy,Gx,Gy,
                           v1,v2,v3,
                           angleIni_B,angleIni_D,angleIni_F) {
  theme_set(theme_bw())
  df = data.frame(Bx=Bx,By=By,Cx=Cx,Cy=Cy,Dx=Dx,Dy=Dy,Ex=Ex,Ey=Ey,Fx=Fx,Fy=Fy,Gx=Gx,Gy=Gy)
  rayonMax = sqrt(df$Bx^2+df$By^2)[1]
  
  p <- ggplot(df) +  
    expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
    coord_fixed(ratio=1) + 
    theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) + 
    geom_point(aes(x=Bx, y=By), size=1) +
    geom_point(aes(x=Cx, y=Cy), size=1) +
    geom_point(aes(x=Dx, y=Dy), size=1, col= "red", shape = 3) +
    geom_point(aes(x=Ex, y=Ey), size=1, col= "yellow") +
    geom_point(aes(x=Fx, y=Fy), size=1, col= "blue") +
    geom_point(aes(x=Gx, y=Gy), size=1, col= "green", shape = 3)+
         labs(title = "Trajectoire Stationnaire",
              subtitle = paste0("v1=",v1,", v2=",v2,", v3=",v3,", \n",
                                "angleIni_B=",angleIni_B,", angleIni_D=",angleIni_D,", angleIni_F=",angleIni_F))
  
  # save ggplot
  #plot_file_name = paste0(dirRes, "/plot_",i,".svg")
  plot_file_name = "plot.png"
  ggsave(plot_file_name,p, width = 8, height = 8, dpi = 72)
  
}





# test
#   v1=1;v2=2.4444;v3=1.5
#   angleIni_B=0;angleIni_D=1.2;angleIni_F=3.14
#   theme_set(theme_bw())
#   df = data.frame(Bx=c(1,2),By=c(3,4))
#   rayonMax = sqrt(df$Bx^2+df$By^2)[1]
#   p <- ggplot(df) +
#     expand_limits(x=c(-rayonMax,rayonMax), y=c(-rayonMax, rayonMax)) +
#     coord_fixed(ratio=1) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_blank() ) +
#     geom_point(aes(x=Bx, y=By), size=1) +
#     labs(title = "Trajectoire Stationnaire",
#          subtitle = paste0("v1=",v1,", v2=",v2,", v3=",v3,", \n",
#                            "angleIni_B=",angleIni_B,", angleIni_D=",angleIni_D,", angleIni_F=",angleIni_F))
# p









