## 


# utile ? si on veut gerder des points à une distance fixe (environ) dans la trajectoire

distanceSquared = function(v1,v2){
  sqrt((v1[1]-v2[1])^2 + (v1[2]-v2[2])^2)
}

distanceSquared2 = function(v1x,v1y,v2x,v2y){
  sqrt((v1x-v2x)^2 + (v1y-v2y)^2)
}

successiveDistance = function(vx,vy){
  n = length(vx)
  v1x = vx[1:(n-1)]
  v1y = vy[1:(n-1)]
  v2x = vx[2:n]
  v2y = vy[2:n]
  distanceSquared2(v1x,v1y,v2x,v2y)
}


nextSimplifiedTrajectory = function(current,res, distance){
  #current:Vector[(Double,Double)]  => array à deux colonnes
  #res:Vector[(Double,Double)]
  #distance:Double
  if (dim(current)[1]==0) {  #isempty(current)
    res
  } else {
     #temp = current[1,]
     #temp2 = apply(current, 1, function(x){distanceSquared(x,temp)})
     temp2 = successiveDistance(current[,1],current[,2])
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



#distance = 1
#resSimpliTraj = simplifiedTrajectory(cbind(res_df_ref$Dx,res_df_ref$Dy), distance)






