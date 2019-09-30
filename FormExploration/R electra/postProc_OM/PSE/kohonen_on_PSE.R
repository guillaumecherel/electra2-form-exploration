## 

# lien utile pour tableau graphe ggplot
# https://github.com/thomasp85/patchwork#patchwork

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
 Tfinal = 5.0
# deltaT = 0.001
times = seq(0,Tfinal,by = deltaT )

# choix par défaut 
lightB=T ; lightC=T ; lightD=T ; lightE=T ; lightF=T ; lightG=T
v1=1 ; v2=1 ; v3=1  





## import functions from utils local (same directory)
source("utils_kohonen_on_PSE.R")

## import data
resPSE = read.csv("data/resultsPSE_5_3/population284000.csv")
resPSE = read.csv("data/PSE_6_1/population28600.csv")
dim(resPSE)

## prepare data
df = as_tibble(resPSE)
dim(resPSE)
n=dim(resPSE)[1];n
colnames(df)
#df = df %>%  filter(nbPointsSinguliers <=10, nbPointsRetour <=10,  -50 <= courbureMoyenne, courbureMoyenne <=50, meanSpeed <=30)
df = df %>%  filter(meanSpeed <=30, indiceMotifD<2000,indiceMotifF<2000)
dim(df)

# test pour voir le nombre de les pettes vitesses
df = df %>%  filter(v1 < -0.1 | v1 > 0.1)
df = df %>%  filter(v2 < -0.1 | v2 > 0.1)  
df = df %>%  filter(v3 < -0.1 | v3 > 0.1)
df = df %>%  filter(v2 < -0.1 | v2 > 0.1, v3 < -0.1 | v3 > 0.1)
dim(df)


df = df %>%  filter(v1_PSE < -0.1 | v1_PSE > 0.1)
df = df %>%  filter(v2_PSE < -0.1 | v2_PSE > 0.1)  
df = df %>%  filter(v3_PSE < -0.1 | v3_PSE > 0.1)


min(df$densite)
max(df$densite)

hist(df$indiceMotifD)
hist(df$indiceMotifF)

hist(df$densite)


df2 = df %>% select(nbPointsSinguliers,nbPointsRetour,densite,courbureMoyenne,meanSpeed)
df2 = df %>% select(densite,meanSpeed,indiceMotifD,indiceMotifF)
df2 = as.matrix(df2)
df2 = scale(df2)
dim(df2)


# som (kohonen)
library(kohonen) #SOM
set.seed(100)
tailleGrid = 3
tailleGrid = 6
carte <-som(df2,grid=somgrid(tailleGrid,tailleGrid,"hexagonal"), rlen=5000)
nbPointsGrille = carte$grid$xdim * carte$grid$ydim


plot(carte)

##summary
#print(summary(carte))
##architecture of the grid
#print(carte$grid)
# count plot
plot(carte,type="count",palette.name=degrade.bleu)
# commentaire de http://eric.univ-lyon2.fr/~ricco/tanagra/fichiers/fr_Tanagra_Kohonen_SOM_R.pdf
# Dans l’idéal, la répartition devrait être assez homogène.
# La taille de la carte doit être réduite s’il y a de nombreuses cellules vides. 
# A contrario,  elle  doit  être  augmentée  si  des  zones  de très forte densité apparaissent (Lynn, 2014).

#noeud (entier entre 1 et taille de grille) d’appartenance des observations
print(carte$unit.classif)
length(carte$unit.classif)

#nombre d’observations affectés à chaque noeud
# même info  que count plot
nb <-table(carte$unit.classif)
as.matrix(nb)
#print(nb)


#check if there are empty nodes
print(length(nb))
nbPointsGrille


# exemples de codebook (représentant de classe), en vrai le plus proche des jeux de param qu'on a qui en est le plus proche
i=3 # noeud
jj = findClosestRepresentantInDf(i,carte,df2)
resRepr = df[jj,]; resRepr

plot_traj_representantNode(resRepr)

i=25
df[findClosestRepresentantInDf(i,carte,df2),]
plot_traj_representantNode(df[findClosestRepresentantInDf(i,carte,df2),])



# save trajectory of almost codebook
dirRes = "resultsPlotKohonen"
#dir.create(dirRes)
#savePlot(i,carte,df2,dirRes)
#lapply(1:nbPointsGrille, function(i){savePlot(i,carte,df2,dirRes)})


# create plot of trajectory of almost codebook
createPlot(i,carte,df2,df)
p = createPlot(i,carte,df2,df)
p




liste_plots = lapply(1:nbPointsGrille, function(i){createPlot(i,carte,df2,df) })
#liste_plots[[1]]


#install.packages("devtools")
#devtools::install_github("thomasp85/patchwork")
devtools::install_github("steveharoz/patchwork")
# liste_plots = liste de plots (obtenus avec createPlot)
library(patchwork)
# crée un beau tableau avec les plots dans la liste

pp = reduce(liste_plots, function(a,b){a+b}) + plot_layout(ncol = 6)
pp = reduce(liste_plots, function(a,b){a+b}) + plot_layout(ncol = 3)

#ggsave("tableauPlot_kohonen.png",pp, width = 15, height = 15, dpi = 400) #dpi = 72)
#ggsave("tableauPlot_kohonen_PSE_6_1_2.png",pp, width = 15, height = 15, dpi = 400) #dpi = 72)

bg <- rgb(4/255,3/255, 51/255)
pp = pp + plot_theme(background = bg)
ggsave("tableauPlot_kohonen_PSE_6_1_2_2.png",pp, width = 15, height = 15, dpi = 400) #dpi = 72)

#plot_theme(background = "gray92")



# tracer tous les résultats

createTrajectorySmart <- function(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,
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
  
  if(lightB>0.5){p = p+geom_point(aes(x=Bx, y=By), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightC>0.5){p = p+geom_point(aes(x=Cx, y=Cy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightD>0.5){p = p+geom_point(aes(x=Dx, y=Dy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightE>0.5){p = p+geom_point(aes(x=Ex, y=Ey), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightF>0.5){p = p+geom_point(aes(x=Fx, y=Fy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  if(lightG>0.5){p = p+geom_point(aes(x=Gx, y=Gy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)} 
  
  p = p+labs(subtitle = paste0("v1=",round(v1,2),", v2=",round(v2,2),", v3=",round(v3,2))) # + theme(plot.subtitle=element_text(color="white"))
  #theme(plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"))
  #p
}


i=1
N = dim(df)[1]
N=20
dirRes="plotResPSE/"
dir.create(dirRes)
for (i in 1:N) {
 
  p = createTrajectorySmart(times, df[i,]$v1_PSE,df[i,]$v2_PSE,df[i,]$v3_PSE, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,df[i,]$lightB_double,df[i,]$lightC_double,df[i,]$lightD_double,df[i,]$lightE_double,df[i,]$lightF_double,df[i,]$lightG_double)
  ggsave(paste0(dirRes,"plot_",i,".png"),p, width = 10, height = 10, dpi = 200) #dpi = 72)
}


length(which((abs(df$v1_PSE)<0.1) + (abs(df$v2_PSE)<0.1) + (abs(df$v3_PSE)<0.1) >= 1))








