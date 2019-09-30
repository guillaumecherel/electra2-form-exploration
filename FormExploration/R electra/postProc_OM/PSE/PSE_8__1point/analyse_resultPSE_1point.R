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
source("../../../utils_electra.R")


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


# choix par défaut 
lightB=T ; lightC=T ; lightD=T ; lightE=T ; lightF=T ; lightG=T
v1=1 ; v2=1 ; v3=1




########################################
###   l'algo PSE a t'il convergé ?   ###
########################################


dir = "data/resultsPSE_8_1/"
files = list.files(dir)
nbFiles = length(files)

# number of points
generation = seq(1000,1000*nbFiles, by = 1000)

sampleSize = c()
reduceSize = c()

options(scipen=999)
for (i in generation){
  temp <- read.csv(paste0(dir,"population",i,".csv"), header = T)
  sampleSize = c(sampleSize,dim(temp)[1])
  df = as.data.frame(temp)
  #fullSize=dim(resPSE)[1]
  df = df %>%  filter(nbPointsSinguliersD <=10,courbureMoyenneAbsolueD <=50, indiceMotifD<=2000, nbPointsRetourD<=10)
  reduceSize= c(reduceSize,dim(df)[1]) 
}

plot(generation,sampleSize, type = "l")
plot(generation,reduceSize, type = "l")





##############################
###        Kohonen         ###
##############################

## import functions from utils local (same directory)
source("utils_kohonen_on_PSE_1_point.R")

## import data
resPSE = read.csv("data/resultsPSE_8_1/population113000.csv")
dim(resPSE)

## prepare data
df = as_tibble(resPSE)
dim(resPSE)
n=dim(resPSE)[1];n
colnames(df)
#df = df %>%  filter(nbPointsSinguliers <=10, nbPointsRetour <=10,  -50 <= courbureMoyenne, courbureMoyenne <=50, meanSpeed <=30)
df = df %>%  filter(nbPointsSinguliersD <=10,courbureMoyenneAbsolueD <=50, indiceMotifD<=2000, nbPointsRetourD<=10)
dim(df)

# test pour voir le nombre de les pettes vitesses
#df = df %>%  filter(v2 < -0.2 | v2 > 0.2)  
df = df %>%  filter(v2_PSE < -0.1 | v2_PSE > 0.1)  
dim(df)

# densite
hist(df$densiteD)
# point retour
df %>% group_by(nbPointsRetourD) %>% summarise(count = n())
# point singulier
df %>% group_by(nbPointsSinguliersD) %>% summarise(count = n())
# courbure 
hist(df$courbureMoyenneAbsolueD)
# indice motif
hist(df$indiceMotifD)


df2 = df %>% select(densiteD,nbPointsRetourD,nbPointsSinguliersD,courbureMoyenneAbsolueD,indiceMotifD)
df2 = df %>% select(nbPointsRetourD,nbPointsSinguliersD,courbureMoyenneAbsolueD,indiceMotifD)
df2 = as.matrix(df2)
df2 = scale(df2)
dim(df2)


# som (kohonen)
library(kohonen) #SOM
set.seed(100)
tailleGrid = 3
tailleGrid = 6
carte <-som(df2,grid=somgrid(tailleGrid,tailleGrid,"hexagonal"), rlen=5000)
carte <-som(df2,grid=somgrid(3,4,"hexagonal"), rlen=5000)
nbPointsGrille = carte$grid$xdim * carte$grid$ydim


plot(carte)

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
i=1 # noeud
jj = findClosestRepresentantInDf(i,carte,df2)
resRepr = df[jj,]; resRepr

plot_traj_representantNode(resRepr)

i=13
df[findClosestRepresentantInDf(i,carte,df2),]
plot_traj_representantNode(df[findClosestRepresentantInDf(i,carte,df2),])



# save trajectory of almost codebook
dirRes = "resultsPlotKohonen_onePoint"
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
ggsave("tableauPlot_kohonen_PSE_8_1_onePoint3_sans_densite.png",pp, width = 15, height = 15, dpi = 400) #dpi = 72)

#plot_theme(background = "gray92")






# tracer/ enregistrer tous les résultats

createTrajectorySmart <- function(times, v1,v2,v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F) {
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
  

  lightB = F
  lightC = F
  lightD = T
  lightE = F
  lightF = F
  lightG = F
  
  if(lightB ==T){p = p+geom_point(aes(x=Bx, y=By), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)}
  if(lightC ==T){p = p+geom_point(aes(x=Cx, y=Cy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)}
  if(lightD ==T){p = p+geom_point(aes(x=Dx, y=Dy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)}
  if(lightE ==T){p = p+geom_point(aes(x=Ex, y=Ey), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)}
  if(lightF ==T){p = p+geom_point(aes(x=Fx, y=Fy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)}
  if(lightG ==T){p = p+geom_point(aes(x=Gx, y=Gy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)} else{p = p+geom_point(aes(x=0, y=0), col= bg, size=0.1)}
  

  
    p = p+labs(subtitle = paste0("v1=",round(v1,2),", v2=",round(v2,2),", v3=",round(v3,2))) # + theme(plot.subtitle=element_text(color="white"))
  #theme(plot.subtitle=element_text(size=18, hjust=0.5, face="italic", color="black"))
  #p
}



i=1
N = dim(df)[1]
#N=3
# times
Tfinal = 10.0
deltaT = 0.0005
times = seq(0,Tfinal,by = deltaT )
# dir
dirRes=paste0("plotResPSE_onePoint_T=",Tfinal,"/")
dir.create(dirRes)

for (i in 1:N) {
  
  p = createTrajectorySmart(times, 1,df[i,]$v2_PSE,0, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F)
  ggsave(paste0(dirRes,"plot_",i,".png"),p, width = 10, height = 10, dpi = 200) #dpi = 72)
}


df[25,]
df[27,]




########################
# à mesures fixées (sauf densité), regrouper celles de densité différentes et ne garder que la minimale?

temp = df %>% group_by(nbPointsSinguliersD,courbureMoyenneAbsolueD,indiceMotifD,nbPointsRetourD) %>% summarise(count = n())


plot.histogramm(df$densiteD)
min(df$densiteD)
max(df$densiteD)


