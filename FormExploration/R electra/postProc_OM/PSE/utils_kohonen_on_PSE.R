## 

# lien utile pour tableau graphe ggplot
# https://github.com/thomasp85/patchwork#patchwork


#jeu de couleurs pour les nœuds de la carte
degrade.bleu <-function(n){
  return(rgb(0,0.4,1,alpha=seq(0,1,1/n)))
}

distanceVectors = function(a,b){
  sum(abs(a-b))
} 



# function that gives the param set in the data that is the best (% distance on outputMesure) = closest 
# to the codebook (représentant) of node i (i in 1: grid size som)
findClosestRepresentantInDf = function(i, carte, df2){
  selInd_i = which(carte$unit.classif == i)
  #length(selInd_i)
  codebokk_i = carte$codes[[1]][i,]
  
  #dim(df2[selInd_i,])
  #df2[selInd_i,]
  temp = df2[selInd_i[i],]
  #temp = df2[selInd_i[1],] %>% select(nbPointsSinguliers,nbPointsRetour,densite,courbureMoyenne,meanSpeed)
  #temp - codebokk_i

  distanceVectors = function(a,b){
    sum(abs(a-b))
  } 
  
  resCompa  = apply(df2[selInd_i,], 1, function(a){distanceVectors(a,codebokk_i)})
  # le codebook n'est pas parmis les noeuds
  # on prend donc le point qui en est le plus proche
  selId_Repr = which(resCompa == min(resCompa))  # l'indice dans les points du noeud
  
  selId_Repr_2 = selInd_i[selId_Repr] # l'indice dans df2 / df
  selId_Repr_2
}


# jj = findClosestRepresentantInDf(3, carte, df2)
# resRepr = df[jj,]




# plot la trajectoire du representant
plot_traj_representantNode = function(resRepr){
  lightB = if(resRepr$lightB_double>=0.5){T}else{F}
  lightC = if(resRepr$lightC_double>=0.5){T}else{F}
  lightD = if(resRepr$lightD_double>=0.5){T}else{F}
  lightE = if(resRepr$lightE_double>=0.5){T}else{F}
  lightF = if(resRepr$lightF_double>=0.5){T}else{F}
  lightG = if(resRepr$lightG_double>=0.5){T}else{F}
  
  # Tfinal = 3
  # deltaT = 0.0005
  # times = seq(0,Tfinal,by = deltaT )
  
  if (is.null(resRepr$v1) ){ 
  plotTrajectorySmart2(times, resRepr$v1_PSE,resRepr$v2_PSE,resRepr$v3_PSE, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)  
  }else{
  plotTrajectorySmart2(times, resRepr$v1,resRepr$v2,resRepr$v3, rB,rC,rD,rE,rF,rG,rH,rI, angleH,angleI, angleIni_B,angleIni_D,angleIni_F,lightB,lightC,lightD,lightE,lightF,lightG)
  }
}


#plot_traj_representantNode(resRepr)


# function to save the plot of the (alomost) codebook
# rq: df2 et df c'est presque pareil (dans df2 on selectionne des colonne de df + df2 matrice et pas dataframe pour les besoin de som)
savePlot = function(i, carte, df2,df,dirRes){
  p = plot_traj_representantNode(df[findClosestRepresentantInDf(i,carte,df2),])  
  
  plot_file_name = paste0(dirRes, "/plot_",i,".png")
  ggsave(plot_file_name,p, width = 8, height = 8, dpi = 200) #dpi = 72)
}




# create a plot (no save)
createPlot = function(i,carte,df2,df){
  
  resRepr = df[findClosestRepresentantInDf(i,carte,df2),]
  p=plot_traj_representantNode(resRepr)  
  p
}












