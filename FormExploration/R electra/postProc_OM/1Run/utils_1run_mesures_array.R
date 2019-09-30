##

# fonction pour tracer la trajectoir quand on a un dataframe (par exemple crée par OM)
# et qu'on ne connait pas forcément les vitesses
plotTrajectoryFromDf = function(df){
  yellow <- rgb(235/255,235/255, 52/255) 
  jaune2 <- rgb(235/255,229/255, 52/255) 
  bg <- rgb(4/255,3/255, 51/255) #
  
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
  
  p = p+
  geom_point(aes(x=Bx, y=By), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)  +
  geom_point(aes(x=Cx, y=Cy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)  +
  geom_point(aes(x=Dx, y=Dy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)  +
  geom_point(aes(x=Ex, y=Ey), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)  +
  geom_point(aes(x=Fx, y=Fy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)  +
  geom_point(aes(x=Gx, y=Gy), col = jaune2, pch=21, fill=yellow, alpha=0.2, size=0.1)  
  
  #p = p+labs(subtitle = paste0("v1=",round(v1,2),", v2=",round(v2,2),", v3=",round(v3,2))) 
  p
  
}

#plotTrajectoryFromDf(df)





#######################################
###    TEMPS DE RETOUR SEGMENTS     ###
#######################################

# but: regrouper les temps de retour qui sont espacés de 1 pas de temps

# a faire sur countRetourD = resTempsRetourSegmentsD %>% group_by(tempsRetourD) %>% summarise(count = n()) 
# car les temps de retour (pas les indices) sont décalés de 1 parfois => les regrouper
nextSimplifyRetour = function(current,res,temp){
  if(isempty(current)){res} else{
    if(size(current)[1] > 1){temp2 = current[1,]} else{temp2 = current[1] }
    l = (size(res)[1])
    if (temp2[1] == temp[1]+1){
      newTemp = temp2
      if(size(current)[1] > 1){newCurrent = current[-1,]} else{newCurrent = c() }
      if (l>1){
        newRes = rbind( res[1: (l-1), ] , c( res[l,1] , res[l,2] + temp2[2]) )} else{
          newRes = c( res[1] , res[2] + temp2[2])
        }
      nextSimplifyRetour(newCurrent,newRes,newTemp)
    }else{
      newTemp = temp2
      if(size(current)[1] > 1){newCurrent = current[-1,]} else{newCurrent = c() }
      newRes = rbind( res ,  temp2 )
      nextSimplifyRetour(newCurrent,newRes,newTemp)
      
    }
  }
}


#resRetourCount = as.matrix(countRetourD)
#current = resRetourCount[2:(dim(resRetourCount)[1]),]

selectRetour = function(resRetourCount){
  res = resRetourCount[1,]
  temp = resRetourCount[1,]
  res2 = nextSimplifyRetour(resRetourCount[2:(dim(resRetourCount)[1]),],res,temp)
  res2[1:(dim(res2)[1]-1),]
}

#resRetourCount = as.matrix(countRetourD)
#resRestourSimpli = selectRetour(resRetourCount)




########### version liste + on recupere les indices des points concernés (pas juste leur nombre)


# retrouver les points associés au temps de retour dans la version simplifiée (rassemblée)
# la fonction crée une liste de liste.
# le premier niveau est de taille lg = dim(resRetourCount)[1] (points non regroupés)
# pour chaque élément de cette première liste, à nouveau une liste dont le premier element est la valeur du temps de retour (durée en pas de temps)
# le deuxieme le nombre de points concernés, le dernier est un vecteur d'indice des points auquels ont eu lieu ce temps de retour
# il est donc de taille le deuxieme element de la liste
create_liste_indice_tempsRetour = function(resRetourCount,resTempsRetourSegments){
  ll = list()
  lg = dim(resRetourCount)[1]
  for (i in 1:lg){
    temp2 = list()
    temp2[[1]] = resRetourCount[i,1]
    temp2[[2]] = resRetourCount[i,2]
    tt = resRetourCount[i,1]
    selInd_tt = which(resTempsRetourSegments$tempsRetourD == tt)
    selIndicesRetour_D = resTempsRetourSegments$indiceRetourD[selInd_tt]
    temp2[[3]] = selIndicesRetour_D
    ll[[i]] = temp2
  }
  ll
}



#liste_tempsRetourD = create_liste_indice_tempsRetour(resRetourCountD,resTempsRetourSegmentsD)


# la même fonction qu'avant pour regrouper, adaptée aux listes
# current et res : list, temp: int (le temps de retour en cours)
nextSimplifyRetourList = function(current,res,temp){
  if(isempty(current)){res} else{
    temp2 = current[[1]] # liste
    l = length(res)
    if (temp2[[1]] == temp[[1]][[1]] + 1){
      newTempRetour = temp2[[1]]
      if(length(current) > 1){newCurrent = current[-1]} else{newCurrent = list() }       
      if (l>1){
        newRes = c( res[1: (l-1)] , list( list( res[[l]][[1]], res[[l]][[2]] + temp2[[2]], c(res[[l]][[3]],  temp2[[3]]) ))) 
      } else{
        newRes = list( list(res[[1]][[1]], res[[1]][[2]] + temp2[[2]],  c(res[[1]][[3]], temp2[[3]]) ) )    #c( res[[1]][1] , res[[1]][2] + newTemp) 
      }
      nextSimplifyRetourList(newCurrent,newRes,newTempRetour)
    }else{
      newTempRetour = temp2[[1]]
      if(length(current) > 1){newCurrent = current[-1]} else{newCurrent = list() }       
      newRes = c( res ,  list(temp2) )
      nextSimplifyRetourList(newCurrent,newRes,newTempRetour)
    }
  }
}

# pour tester
# current = newCurrent
# res = newRes
# temp = newTempRetour



# ll = ll[1:5]
selectRetour_List = function(liste_tempsRetour){
  res = liste_tempsRetour[1]
  temp = liste_tempsRetour[1]
  current = liste_tempsRetour[2:length(liste_tempsRetour)]
  nextSimplifyRetourList(current,res,temp)
}



#resList = resRestourSimpli_List = selectRetour_List(liste_tempsRetourD)

















