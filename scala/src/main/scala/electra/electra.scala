// TODO
// (pas urgent) faire une fonction transitoire pour le nouveau model (adaptation)


package electra


import electra.NumericalIntegration._
import electra.Model._
import electra.MesurePSE._
import electra.Mesure.{arrayDensitySquare, sumVectorOfArraysOfArrays, _}
import org.apache.commons.math3._
import org.apache.commons.math3.stat.StatUtils

import scala.math._
import org.openmole.spatialdata._
import org.openmole.spatialdata.grid.measures.GridMorphology
import org.openmole.spatialdata.points.measures._ //SpatStat


object Electra extends App {

  ////////////////////////////////
  //    1 POSITION AT TIME T
  ////////////////////////////////


  // 1 position: Stationnary
/*
    val t = 1.0
    val res = dynamicStationnary()(DefaultValuesParameterModel.v1, DefaultValuesParameterModel.v2, DefaultValuesParameterModel.v3)(t)
    println(res)
*/


  // 1 position: Transitoire, avec next
/*
  val fixedParametersModel = new FixedParametersModel()
  //val parametersTransitoire = new ParametersTransitoire()
  val parametersTransitoire = new ParametersTransitoire(DefaultValuesParameterModel.u1_sin,
    DefaultValuesParameterModel.u2_sin,DefaultValuesParameterModel.u3_sin)
  val Nmax = 100
  val deltaT = 0.01
  val stepsTransitoryNext = HiddenParameters.stepsTransitoryNext
  val initialConditions = new InitialConditions()
  val res = simuTransitoire(fixedParametersModel)(parametersTransitoire )(Nmax,deltaT,initialConditions,stepsTransitoryNext)
  println(res)
  //println(res.Ex)
*/


  // 1 position: Transitoire, sans next
/*
  val Nmax = 100
  val deltaT = 0.01
  val stepsTransitory = HiddenParameters.stepsTransitory
  val t = Nmax*deltaT
  val res2 =  dynamicTransitoire(stepsTransitory)()(t,DefaultValuesParameterModel.u1_sin, DefaultValuesParameterModel.u2_sin, DefaultValuesParameterModel.u3_sin)
  println(res2)
  //println(res2.Ex)
*/


  ////////////////////////////////
  //    TRAJECTOIRE
  ////////////////////////////////


  // Trajectoire stationnaire
/*
    val T = 0.02
    val deltaT = 0.01
    val res = dynamicTrajectoryStationnary()()(T,deltaT)
    //println(res)
    val res2 = convertResultStationnary(res)
    println(res2)
    //println(res2.Bx)
*/




  // Trajectoire transitoire (next)
/*
  val fixedParametersModel = new FixedParametersModel()
  val parametersTransitoire = new ParametersTransitoire()
  val Nmax = 2
  val deltaT = 0.01
  val initialConditions = new InitialConditions()
  val stepsTransitoryNext = HiddenParameters.stepsTransitoryNext

  val res = simuTrajectoireTransitoire(fixedParametersModel)(parametersTransitoire)(Nmax,deltaT,initialConditions,stepsTransitoryNext)
  print(res)
*/





  ////////////////////////////////
  //    MESURES
  ////////////////////////////////

  // créer trajectoire stationnaire

/*
  val T = 10.0
  val deltaT = 0.001
  val res = dynamicTrajectoryStationnary(rB=1.0,rH=0.5,rD=0.25)(v1=2.0,v2=(-6.0))(T,deltaT)
  //val res = dynamicTrajectoryStationnary()(v1 = 5.0, v2 = 4.0, v3 = Pi)(T, deltaT)
  //val res = dynamicTrajectoryStationnary(angleIni_B = 0.0,angleIni_D = 0.0,angleIni_F = 0.0)(v1 = -5.0, v2 = -5.0, v3 = -5.0)(T, deltaT)
  val res2 = convertResultStationnary(res)
*/


  ////////////////////////////////
  //    MESURES. points singuliers
  ////////////////////////////////

/*
  val T = 10.0
  val deltaT = 0.001
  val res = dynamicTrajectoryStationnary(rB=1.0,rH=0.5,rD=0.25)(v1=2.0,v2=(-6.0))(T,deltaT)
  val res2 = convertResultStationnary(res)
*/

/*
  val seuilPointSingulier = 4
  val numberPointSinguliersD = countSingularPoints(res2.speedDx,res2.speedDy,seuilPointSingulier)
  //println(numberPointSinguliersD)

  val timesPointSinguliersD = timesOfSingularPoints(res2.speedDx,res2.speedDy,seuilPointSingulier)
  //println(timesPointSinguliersD.mkString(" "))
  //println(timesPointSinguliersD.length)

  val anglesPointSinguliersD = angleOfSingularPoints(res2.speedDx,res2.speedDy,seuilPointSingulier)
  //println(anglesPointSinguliersD.mkString(" "))
  //println(anglesPointSinguliersD.length)
*/


  ////////////////////////////////
  //    MESURES. densité
  ////////////////////////////////

/*
    val T = 10.0
    val deltaT = 0.001
    val res = dynamicTrajectoryStationnary()(v1 = 5.0, v2 = 4.0, v3 = Pi)(T, deltaT)
    val res2 = convertResultStationnary(res)
*/


/*
    val N = 50  // pas de la subdiviion du carré
    val xmax = maxSquareForDensity()
    val xmin = -xmax
    val ymax = xmax
    val ymin = xmin
    val rMax = xmax
    val rMin = minSquareForDensity()


    val densiteB =  pointsDensity(res2.Bx,res2.By,xmin,xmax,ymin,ymax,N,rMin,rMax)
    println(densiteB)

    val totalDensity = densityAllTrajectories(res2,xmin,xmax,ymin,ymax,N,rMin,rMax)
    println(totalDensity)
*/



  ////////////////////////////////
  //    MESURES. Courbure
  ////////////////////////////////

/*
      val T = 10.0
      val deltaT = 0.001
      val res = dynamicTrajectoryStationnary()(v1 = 5.0, v2 = 4.0, v3 = Pi)(T, deltaT)
      val res2 = convertResultStationnary(res)
*/


/*
  // point B
  val courburesB = courbure(res2.speedBx,res2.speedBy,res2.accBx,res2.accBy)
  //println(courburesB)

  // point D
  val courburesD = courbure(res2.speedDx,res2.speedDy,res2.accDx,res2.accDy)
  //println(courburesD)

  // à tester sur une courbe qui à une ligne droite
  val seuilCourbureLigneDroite = 3
  //val indicesLigneDroite = detectStraitLine(courburesB,seuilCourbureLigneDroite)
  //print(indicesLigneDroite)

*/


  ////////////////////////////////
  //    MESURES. loop points: temps de premier retour
  ////////////////////////////////

/*
  val T = 10.0
  val deltaT = 0.001
  val res = dynamicTrajectoryStationnary()(v1 = 5.0, v2 = 4.0, v3 = Pi)(T, deltaT)
  val res2 = convertResultStationnary(res)
*/

/*
  val seuilRetour = max(math.floor(1/ (deltaT * DefaultValuesParameterModel.v1)).toInt - 10 ,0)

  val (indicesRetourSegmentsB,tempsRetourSegmentsB,anglesRetourSegmentsB) = indiceEtTempsPremierRetourEtAngle(res2.Bx,res2.By)
  val (indicesRetourSegmentsD,tempsRetourSegmentsD,anglesRetourSegmentsD) = indiceEtTempsPremierRetourEtAngle(res2.Dx,res2.Dy)
  val (indicesRetourSegmentsF,tempsRetourSegmentsF,anglesRetourSegmentsF) = indiceEtTempsPremierRetourEtAngle(res2.Fx,res2.Fy)
  //println(indicesRetourSegmentsD.mkString(""))
  //println(tempsRetourSegmentsD.mkString(""))
  //println(anglesRetourSegmentsD.mkString(""))

  //println(mean(tempsRetourAllTrajectories(res2,1000).map(_.toDouble).toArray).getOrElse(0))

*/




  ////////////////////////////////
  //    MESURES. Moran
  ////////////////////////////////

/*
  val T = 10.0
  val deltaT = 0.001
  val res = dynamicTrajectoryStationnary()(v1 = 5.0, v2 = 4.0, v3 = Pi)(T, deltaT)
  val res2 = convertResultStationnary(res)
*/


/*
  val N = 50  // pas de la subdiviion du carré
  val xmax = maxSquareForDensity()
  val xmin = -xmax
  val ymax = xmax
  val ymin = xmin


  // D
  val moran_directD = moranDirectTraj(res2.Dx,res2.Dy,xmin,xmax,ymin,ymax,N)
  val moranD = moranTraj(res2.Dx,res2.Dy,xmin,xmax,ymin,ymax,N)
  println(moran_directD)
  println(moranD)

  // B
  val moran_directB = moranDirectTraj(res2.Bx,res2.By,xmin,xmax,ymin,ymax,N)
  val moranB = moranTraj(res2.Bx,res2.By,xmin,xmax,ymin,ymax,N)
  println(moran_directB)
  println(moranB)

  // F
  val moran_directF = moranDirectTraj(res2.Fx,res2.Fy,xmin,xmax,ymin,ymax,N)
  val moranF = moranTraj(res2.Fx,res2.Fy,xmin,xmax,ymin,ymax,N)
  println(moran_directF)
  println(moranF)
*/



  ////////////////////////////////
  //   Comparaison trajectoires sans aspect temporel (trace)
  ////////////////////////////////

/*
  val rayon = 2.3
  val nbPointsCircle = 100
  val (x_circle,y_circle) = createCircle(rayon,nbPointsCircle)
  //println(x_circle)
  //println(y_circle)

  val distance = 1.0
  val res3 = simplifiedTrajectory(x_circle.toVector,y_circle.toVector,distance)
  println(res3.length)
  //println(res3.map(x=>x._1))
  //println(res3.map(x=>x._2))
  //val res4 = distanceTwoTrajectories(x_circle,y_circle,res2.Dx.toArray,res2.Dy.toArray,distance)
  //println(res4)

*/

  ////////////////////////////////
  //    Persistence stationnaire
  ////////////////////////////////

/*
  val t = 1.0
  val resPrsistenceStationnaire = figurePersistenceStationnaireAtTimeT(time=t)()()
  println(resPrsistenceStationnaire)
*/








  ////////////////////////////////
  //    PAS TOUTES LES LUMIRES ALLUMÉES
  ////////////////////////////////


  val T = 10.0
  val deltaT = 0.001
  val res = dynamicTrajectoryStationnary(rB=1.0,rH=0.5,rD=0.25)(v1=2.0,v2=(-6.0))(T,deltaT)
  val res2 = convertResultStationnary(res)
  val res3 = trajectoryLightChoice(res2)(lightB = false)
  //val res3 = trajectoryLightChoice(res2)()
  //val res3 = trajectoryLightChoice(res2)(lightB = false,lightC = false,lightD = false,lightE = false,lightF = false,lightG = false)
  //println(res3.Cx)


  // Mesures PSE (pas toutes les lumières)


  // points singuliers
/*
  val seuilPointSingulier = 4
  val numberPointSinguliersB = countSingularPoints(res2.speedBx,res2.speedBy,seuilPointSingulier)
  println(numberPointSinguliersB)

  println(countSingularPointsAllTrajectories(res3,seuilPointSingulier))
*/




  // point retour (boucle) [le seuil est mis pour éviter de compter les points des trajectoires périodiques]
/*
  val seuilRetour = max(math.floor(1/ (deltaT * DefaultValuesParameterModel.v1)).toInt - 10 ,0)
  //println(numberRetour(res3.Bx,res3.By,seuilRetour))
  println(numberRetourAllTrajectories(res3,seuilRetour))
*/




  // densite (dans carre)
/*
  val N = 50  // Nombre de points pour la subdiviion du segment pour le carré
  val xmax = maxSquareForDensity(rB=1.0,rH=0.5,rD=0.25)  // être raccord avec la dynamique de la trajectoire
  val xmin = -xmax
  val ymax = xmax
  val ymin = xmin

  val densiteB =  pointsDensitySquare(res3.Bx,res3.By,xmin,xmax,ymin,ymax,N)
  println(densiteB)

  println(densityAllTrajectoriesSquare(res3,xmin,xmax,ymin,ymax,N))
*/




  // Moran
/*
  val N = 50  // Nombre de points pour la subdiviion du segment pour le carré
  val xmax = maxSquareForDensity(rB=1.0,rH=0.5,rD=0.25)  // être raccord avec la dynamique de la trajectoire
  val xmin = -xmax
  val ymax = xmax
  val ymin = xmin

  val moranD = moranTraj(res3.Dx,res3.Dy,xmin,xmax,ymin,ymax,N)
  println(moranD)

  println( moranAllTrajectories(res3,xmin,xmax,ymin,ymax,N) )
*/



  // courbure moyenne

  println(meanCourbure(res3.speedBx,res3.speedBy,res3.accBx,res2.accBy))
  println(meanCourbureAllTrajectories(res3))


}











case class DynamicalCurrentStateTransitory(time : Double, Bx: Double, By:Double, Cx:Double=0.0, Cy:Double=0.0, Dx:Double=0.0, Dy:Double=0.0,
                                           Ex:Double=0.0, Ey:Double=0.0, Fx:Double=0.0, Fy:Double=0.0, Gx:Double=0.0, Gy:Double=0.0, Hx:Double=0.0, Hy:Double=0.0,
                                           Ix:Double=0.0, Iy:Double=0.0,
                                           HDx:Double=0.0, HDy:Double=0.0, HEx:Double=0.0, HEy:Double=0.0,
                                           IFx:Double=0.0, IFy:Double=0.0, IGx:Double=0.0, IGy:Double=0.0,
                                           angle1:Double=0.0, angle2:Double=0.0, angle3:Double=0.0,
                                           vitesseMoteur1:Double=0.0, vitesseMoteur2:Double=0.0, vitesseMoteur3:Double=0.0,
                                           currentIntegralMoteur1 :Double=0.0, currentIntegralMoteur2 :Double=0.0, currentIntegralMoteur3 :Double=0.0
                               )



case class DynamicalCurrentStateStationary(time : Double=0.0, Bx: Double=0.0, By:Double=0.0, Cx:Double=0.0, Cy:Double=0.0, Dx:Double=0.0, Dy:Double=0.0,
                                           Ex:Double=0.0, Ey:Double=0.0, Fx:Double=0.0, Fy:Double=0.0, Gx:Double=0.0, Gy:Double=0.0, Hx:Double=0.0, Hy:Double=0.0,
                                           Ix:Double=0.0, Iy:Double=0.0,
                                           HDx:Double=0.0, HDy:Double=0.0, HEx:Double=0.0, HEy:Double=0.0,
                                           IFx:Double=0.0, IFy:Double=0.0, IGx:Double=0.0, IGy:Double=0.0,
                                           speedBx: Double=0.0, speedBy:Double=0.0, speedCx:Double=0.0, speedCy:Double=0.0, speedDx:Double=0.0, speedDy:Double=0.0,
                                           speedEx:Double=0.0, speedEy:Double=0.0, speedFx:Double=0.0, speedFy:Double=0.0, speedGx:Double=0.0, speedGy:Double=0.0,
                                           accBx: Double=0.0, accBy:Double=0.0, accCx:Double=0.0, accCy:Double=0.0, accDx:Double=0.0, accDy:Double=0.0,
                                           accEx:Double=0.0, accEy:Double=0.0, accFx:Double=0.0, accFy:Double=0.0, accGx:Double=0.0, accGy:Double=0.0
                                         )




case class TrajectoryStationnary (time:Vector[Double], Bx:Vector[Double], By:Vector[Double], Cx:Vector[Double], Cy:Vector[Double],
                                          Dx:Vector[Double], Dy:Vector[Double],Ex:Vector[Double], Ey:Vector[Double], Fx:Vector[Double],
                                          Fy:Vector[Double], Gx:Vector[Double], Gy:Vector[Double],
                                  speedBx:Vector[Double], speedBy:Vector[Double], speedCx:Vector[Double], speedCy:Vector[Double], speedDx:Vector[Double],
                                  speedDy:Vector[Double], speedEx:Vector[Double], speedEy:Vector[Double], speedFx:Vector[Double], speedFy:Vector[Double],
                                  speedGx:Vector[Double], speedGy:Vector[Double],
                                  accBx:Vector[Double], accBy:Vector[Double], accCx:Vector[Double], accCy:Vector[Double], accDx:Vector[Double],
                                  accDy:Vector[Double], accEx:Vector[Double], accEy:Vector[Double], accFx:Vector[Double], accFy:Vector[Double],
                                  accGx:Vector[Double], accGy:Vector[Double]
                                   )




case class FixedParametersModel(rB :Double=FixedParameterModel.rB, rC:Double=FixedParameterModel.rC, rD:Double=FixedParameterModel.rD,
                                rE:Double=FixedParameterModel.rE, rF:Double=FixedParameterModel.rF,  rG:Double=FixedParameterModel.rG,
                                rH:Double=FixedParameterModel.rH, rI:Double=FixedParameterModel.rI,
                                A1:Double=FixedParameterModel.A1 ,C1:Double=FixedParameterModel.C1,
                                A2:Double=FixedParameterModel.A2,C2:Double=FixedParameterModel.C2,
                                A3:Double=FixedParameterModel.A3,C3:Double=FixedParameterModel.C3)

case class ParametersTransitoire(u1:Double=>Double = DefaultValuesParameterModel.u1,
                                 u2:Double=>Double = DefaultValuesParameterModel.u2,
                                 u3:Double=>Double = DefaultValuesParameterModel.u3)

case class InitialConditions(angleIni_B:Double = DefaultValuesParameterModel.angleIni_B,
                           angleIni_D:Double = DefaultValuesParameterModel.angleIni_D,
                           angleIni_F:Double = DefaultValuesParameterModel.angleIni_F)



object Model {

  ////////////////////////////////
  //    STATIONNAIRE
  ////////////////////////////////


  def dynamicStationnary(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD,
                         rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                         rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI,
                         angleH:Double = FixedParameterModel.angleH, angleI:Double = FixedParameterModel.angleI,
                         angleIni_B:Double = DefaultValuesParameterModel.angleIni_B, angleIni_D:Double = DefaultValuesParameterModel.angleIni_D,
                         angleIni_F:Double = DefaultValuesParameterModel.angleIni_F)(v1:Double=DefaultValuesParameterModel.v1, v2:Double=DefaultValuesParameterModel.v2, v3:Double=DefaultValuesParameterModel.v3)(t:Double) = {

    // Position
    // B,C
    val Bx = rB * cos(2*Pi*v1*t + angleIni_B)
    val By = rB * sin(2*Pi*v1*t + angleIni_B)
    val Cx = rC * cos(2*Pi*v1*t + angleIni_B+Pi)
    val Cy = rC * sin(2*Pi*v1*t + angleIni_B+Pi)

    // H,I
    val Hx = rH * cos(2*Pi*v1*t + angleIni_B - angleH)
    val Hy = rH * sin(2*Pi*v1*t + angleIni_B - angleH)
    val Ix = rI * cos(2*Pi*v1*t + angleIni_B + angleI +Pi)
    val Iy = rI * sin(2*Pi*v1*t + angleIni_B + angleI +Pi)

    // D,E
    val HDx = rD  * cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH))
    val HDy = rD  * sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH))
    val HEx = rE  * cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH) +Pi)
    val HEy = rE  * sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH) +Pi)

    val Dx = Hx + HDx
    val Dy = Hy + HDy
    val Ex = Hx + HEx
    val Ey = Hy + HEy

    // F,G
    val IFx = rF  * cos(2*Pi*(v3+v1)*t + (angleIni_F+ (angleIni_B +angleI +Pi)) )
    val IFy = rF  * sin(2*Pi*(v3+v1)*t + (angleIni_F+ (angleIni_B +angleI +Pi)))
    val IGx = rG  * cos(2*Pi*(v3+v1)*t + (angleIni_F+ (angleIni_B +angleI +Pi)) +Pi)
    val IGy = rG  * sin(2*Pi*(v3+v1)*t + (angleIni_F+ (angleIni_B +angleI +Pi)) +Pi)

    val Fx = Ix + IFx
    val Fy = Iy + IFy
    val Gx = Ix + IGx
    val Gy = Iy + IGy

    // speed
    // B,C
    val speedBx = -rB * 2*Pi*v1* sin(2*Pi*v1*t + angleIni_B)
    val speedBy = rB * 2*Pi*v1* cos(2*Pi*v1*t + angleIni_B)
    val speedCx = -rC * 2*Pi*v1* sin(2*Pi*v1*t + angleIni_B+Pi)
    val speedCy = rC * 2*Pi*v1* cos(2*Pi*v1*t + angleIni_B+Pi)

    // D,E
    val speedDx = -rH * 2*Pi*v1* sin(2*Pi*v1*t + angleIni_B -angleH)  - rD *2*Pi*(v2+v1)* sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH))
    val speedDy = rH * 2*Pi*v1*  cos(2*Pi*v1*t + angleIni_B -angleH)  + rD *2*Pi*(v2+v1)* cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH))
    val speedEx = -rH * 2*Pi*v1* sin(2*Pi*v1*t + angleIni_B -angleH)  - rE *2*Pi*(v2+v1)* sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH) +Pi)
    val speedEy = rH * 2*Pi*v1*  cos(2*Pi*v1*t + angleIni_B -angleH)  + rE *2*Pi*(v2+v1)* cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH) +Pi)

    // F,G
    val speedFx = -rI * (2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B+Pi +angleI)   - rF* (2*Pi*(v3+v1))* sin(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi +angleI))
    val speedFy = rI * (2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B+Pi +angleI)    + rF* (2*Pi*(v3+v1))* cos(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi +angleI))
    val speedGx = -rI * (2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B+Pi +angleI)   - rG* (2*Pi*(v3+v1))* sin(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi +angleI) +Pi)
    val speedGy = rI * (2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B+Pi +angleI)    + rG* (2*Pi*(v3+v1))* cos(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi +angleI) +Pi)


    ////////////////////////////////
    // acceleration
    // B,C
    val accBx = -rB * (2*Pi*v1)*(2*Pi*v1) * cos(2*Pi*v1*t + angleIni_B)
    val accBy = -rB * (2*Pi*v1)*(2*Pi*v1) * sin(2*Pi*v1*t + angleIni_B)
    val accCx = -rC * (2*Pi*v1)*(2*Pi*v1) * cos(2*Pi*v1*t + angleIni_B+Pi)
    val accCy = -rC * (2*Pi*v1)*(2*Pi*v1) * sin(2*Pi*v1*t + angleIni_B+Pi)

    // D,E
    val accDx = -rH * (2*Pi*v1)*(2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B -angleH)  - rD* (2*Pi*(v2+v1))*(2*Pi*(v2+v1))* cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH))
    val accDy = -rH * (2*Pi*v1)*(2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B -angleH)  - rD* (2*Pi*(v2+v1))*(2*Pi*(v2+v1))* sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH))
    val accEx = -rH * (2*Pi*v1)*(2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B -angleH)  - rE* (2*Pi*(v2+v1))*(2*Pi*(v2+v1))* cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH) +Pi)
    val accEy = -rH * (2*Pi*v1)*(2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B -angleH)  - rE* (2*Pi*(v2+v1))*(2*Pi*(v2+v1))* sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B -angleH) +Pi)

    // F,G
    val accFx = -rI * (2*Pi*v1)*(2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B+Pi +angleI)  - rF* (2*Pi*(v3+v1))*(2*Pi*(v3+v1))* cos(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi +angleI))
    val accFy = -rI * (2*Pi*v1)*(2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B+Pi +angleI)  - rF* (2*Pi*(v3+v1))*(2*Pi*(v3+v1))* sin(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi +angleI))
    val accGx = -rI * (2*Pi*v1)*(2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B+Pi +angleI)  - rG* (2*Pi*(v3+v1))*(2*Pi*(v3+v1))* cos(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi +angleI) +Pi)
    val accGy = -rI * (2*Pi*v1)*(2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B+Pi +angleI)  - rG* (2*Pi*(v3+v1))*(2*Pi*(v3+v1))* sin(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi +angleI) +Pi)

    //(Bx,By,Cx,Cy,Dx,Dy,Ex,Ey,Fx,Fy,Gx,Gy,Hx,Hy,Ix,Iy)

    new DynamicalCurrentStateStationary(time = t, Bx=Bx, By=By, Cx=Cx, Cy=Cy, Dx=Dx, Dy=Dy,
      Ex=Ex, Ey=Ey, Fx=Fx, Fy=Fy, Gx=Gx, Gy=Gy, Hx=Hx, Hy=Hy,
      Ix=Ix, Iy=Iy,
      HDx=HDx, HDy=HDy, HEx=HEx, HEy=HEy,
      IFx=IFx, IFy=IFy, IGx=IGx, IGy=IGy,
      speedBx=speedBx, speedBy=speedBy, speedCx=speedCx, speedCy=speedCy, speedDx=speedDx, speedDy=speedDy,
      speedEx=speedEx, speedEy=speedEy, speedFx=speedFx, speedFy=speedFy, speedGx=speedGx, speedGy=speedGy,
      accBx=accBx, accBy=accBy, accCx=accCx, accCy=accCy, accDx=accDx, accDy=accDy,
      accEx=accEx, accEy=accEy, accFx=accFx, accFy=accFy, accGx=accGx, accGy=accGy
    )
  }


  def dynamicTrajectoryStationnary(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD,
                                   rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                                   rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI,
                                   angleH:Double = FixedParameterModel.angleH, angleI:Double = FixedParameterModel.angleI,
                                   angleIni_B:Double = DefaultValuesParameterModel.angleIni_B, angleIni_D:Double = DefaultValuesParameterModel.angleIni_D,
                                   angleIni_F:Double = DefaultValuesParameterModel.angleIni_F)(v1:Double=DefaultValuesParameterModel.v1, v2:Double=DefaultValuesParameterModel.v2, v3:Double=DefaultValuesParameterModel.v3)(T:Double, deltaT:Double):Vector[DynamicalCurrentStateStationary] = {
    val N = floor(T / deltaT).toInt
    val vecTimes = (0 until N).map(x => x * deltaT) :+ T
    val res = vecTimes.map(x => dynamicStationnary(rB,rC,rD,rE,rF,rG,rH,rI,angleH,angleI,angleIni_B,angleIni_D,angleIni_F)(v1,v2,v3)(x)).toVector
    // type de res: Vector[DynamicalCurrentStateStationary]
    res
  }






  def trajectoryLightChoice(res:TrajectoryStationnary)(lightB :Boolean = DefaultValuesParameterModel.lightB, lightC :Boolean = DefaultValuesParameterModel.lightC,
                                                       lightD :Boolean = DefaultValuesParameterModel.lightD, lightE :Boolean = DefaultValuesParameterModel.lightE,
                                                       lightF :Boolean = DefaultValuesParameterModel.lightF, lightG :Boolean = DefaultValuesParameterModel.lightG) = {

    new TrajectoryStationnary(
      time = res.time,
      Bx = if (lightB) res.Bx else Vector.empty,
      By = if (lightB) res.By else Vector.empty,
      Cx = if (lightC) res.Cx else Vector.empty,
      Cy = if (lightC) res.Cy else Vector.empty,
      Dx = if (lightD) res.Dx else Vector.empty,
      Dy = if (lightD) res.Dy else Vector.empty,
      Ex = if (lightE) res.Ex else Vector.empty,
      Ey = if (lightE) res.Ey else Vector.empty,
      Fx = if (lightF) res.Fx else Vector.empty,
      Fy = if (lightF) res.Fy else Vector.empty,
      Gx = if (lightG) res.Gx else Vector.empty,
      Gy = if (lightG) res.Gy else Vector.empty,
      speedBx = if (lightB) res.speedBx else Vector.empty,
      speedBy = if (lightB) res.speedBy else Vector.empty,
      speedCx = if (lightC) res.speedCx else Vector.empty,
      speedCy = if (lightC) res.speedCy else Vector.empty,
      speedDx = if (lightD) res.speedDx else Vector.empty,
      speedDy = if (lightD) res.speedDy else Vector.empty,
      speedEx = if (lightE) res.speedEx else Vector.empty,
      speedEy = if (lightE) res.speedEy else Vector.empty,
      speedFx = if (lightF) res.speedFx else Vector.empty,
      speedFy = if (lightF) res.speedFy else Vector.empty,
      speedGx = if (lightG) res.speedGx else Vector.empty,
      speedGy = if (lightG) res.speedGy else Vector.empty,
      accBx = if (lightB) res.accBx else Vector.empty,
      accBy = if (lightB) res.accBy else Vector.empty,
      accCx = if (lightC) res.accCx else Vector.empty,
      accCy = if (lightC) res.accCy else Vector.empty,
      accDx = if (lightD) res.accDx else Vector.empty,
      accDy = if (lightD) res.accDy else Vector.empty,
      accEx = if (lightE) res.accEx else Vector.empty,
      accEy = if (lightE) res.accEy else Vector.empty,
      accFx = if (lightF) res.accFx else Vector.empty,
      accFy = if (lightF) res.accFy else Vector.empty,
      accGx = if (lightG) res.accGx else Vector.empty,
      accGy = if (lightG) res.accGy else Vector.empty)
  }





  // convert the result of dynamicTrajectoryStationnary from  Vector[DynamicalCurrentStateStationary]   to TrajectoryStationnary
  def convertResultStationnary(res:Vector[DynamicalCurrentStateStationary])={

    val time = res.map(_.time)

    // Positions
    val Bx = res.map(_.Bx)
    val By = res.map(_.By)
    val Cx = res.map(_.Cx)
    val Cy = res.map(_.Cy)

    val Dx = res.map(_.Dx)
    val Dy = res.map(_.Dy)
    val Ex = res.map(_.Ex)
    val Ey = res.map(_.Ey)

    val Fx = res.map(_.Fx)
    val Fy = res.map(_.Fy)
    val Gx = res.map(_.Gx)
    val Gy = res.map(_.Gy)

    // speed
    val speedBx = res.map(_.speedBx)
    val speedBy = res.map(_.speedBy)
    val speedCx = res.map(_.speedCx)
    val speedCy = res.map(_.speedCy)

    val speedDx = res.map(_.speedDx)
    val speedDy = res.map(_.speedDy)
    val speedEx = res.map(_.speedEx)
    val speedEy = res.map(_.speedEy)

    val speedFx = res.map(_.speedFx)
    val speedFy = res.map(_.speedFy)
    val speedGx = res.map(_.speedGx)
    val speedGy = res.map(_.speedGy)


    // acc
    val accBx = res.map(_.accBx)
    val accBy = res.map(_.accBy)
    val accCx = res.map(_.accCx)
    val accCy = res.map(_.accCy)

    val accDx = res.map(_.accDx)
    val accDy = res.map(_.accDy)
    val accEx = res.map(_.accEx)
    val accEy = res.map(_.accEy)

    val accFx = res.map(_.accFx)
    val accFy = res.map(_.accFy)
    val accGx = res.map(_.accGx)
    val accGy = res.map(_.accGy)


    new TrajectoryStationnary(time = time, Bx = Bx,By = By,Cx = Cx,Cy = Cy,Dx = Dx,Dy = Dy,
      Ex = Ex,Ey = Ey,Fx = Fx,Fy = Fy,Gx = Gx,Gy = Gy,
      speedBx = speedBx,speedBy = speedBy,speedCx = speedCx,speedCy = speedCy,speedDx = speedDx,speedDy = speedDy,
      speedEx = speedEx,speedEy = speedEy,speedFx = speedFx,speedFy = speedFy,speedGx = speedGx,speedGy = speedGy,
      accBx = accBx,accBy = accBy,accCx = accCx,accCy = accCy,accDx = accDx,accDy = accDy,
      accEx = accEx,accEy = accEy,accFx = accFx,accFy = accFy,accGx = accGx,accGy = accGy)

  }




  ////////////////////////////////
  //    TRANSITOIRE
  ////////////////////////////////


  def vitesse_rotation_moteur(stepsTransitory:Int)(u: Double =>Double ,A:Double,C:Double)(t:Double)= {
    def temp(t:Double):Double ={
      A / C * u(t) * exp(t/C)
    }

    val vitesse = exp(-t/C ) * NumericalIntegration.integrate(temp,0.0,t,(ceil(abs(t))+1).toInt*stepsTransitory,simpson)
    vitesse
  }


  def calcul_angle(stepsTransitory:Int)(u: Double =>Double ,A:Double,C:Double)(t:Double) ={
    def vit(t:Double):Double = {Model.vitesse_rotation_moteur(stepsTransitory)(u,A,C)(t)}
    //val angle = NumericalIntegration.integrate((x:Double) => x*x,0.0,t,ceil(t).toInt*HiddenParameters.steps,simpson)
    val angle = NumericalIntegration.integrate(vit,0.0,t,(ceil(abs(t))+1).toInt*stepsTransitory,simpson)
    //vitesse(2.0)
    angle
  }



  def dynamicTransitoire(stepsTransitory:Int=HiddenParameters.stepsTransitory)(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD ,
                         rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                         rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI,
                         angleIni_B:Double = DefaultValuesParameterModel.angleIni_B, angleIni_D:Double = DefaultValuesParameterModel.angleIni_D,
                         angleIni_F:Double = DefaultValuesParameterModel.angleIni_F,
                         A1:Double = FixedParameterModel.A1,C1:Double=FixedParameterModel.C1,
                         A2:Double = FixedParameterModel.A2,C2:Double=FixedParameterModel.C2,
                         A3:Double = FixedParameterModel.A3,C3:Double=FixedParameterModel.C3)
                        (t:Double, u1:Double=>Double, u2:Double=>Double, u3:Double=>Double) = {

    // angle
    val angle1_t = calcul_angle(stepsTransitory)(u1,A1,C1)(t)
    val angle2_t = calcul_angle(stepsTransitory)(u2,A2,C2)(t)
    val angle3_t = calcul_angle(stepsTransitory)(u3,A3,C3)(t)

    // Position
    // B,C
    val Bx = rB * cos(2 * Pi * angle1_t + angleIni_B)
    val By = rB * sin(2 * Pi * angle1_t + angleIni_B)
    val Cx = rC * cos(2 * Pi * angle1_t + angleIni_B + Pi)
    val Cy = rC * sin(2 * Pi * angle1_t + angleIni_B + Pi)

    // H,I
    val Hx = rH * cos(2 * Pi * angle1_t + angleIni_B)
    val Hy = rH * sin(2 * Pi * angle1_t + angleIni_B)
    val Ix = rI * cos(2 * Pi * angle1_t + angleIni_B + Pi)
    val Iy = rI * sin(2 * Pi * angle1_t + angleIni_B + Pi)

    // D,E
    val HDx = rD * cos(2 * Pi * (angle2_t + angle1_t)  + (angleIni_D + angleIni_B))
    val HDy = rD * sin(2 * Pi * (angle2_t + angle1_t)  + (angleIni_D + angleIni_B))
    val HEx = rE * cos(2 * Pi * (angle2_t + angle1_t)  + (angleIni_D + angleIni_B) + Pi)
    val HEy = rE * sin(2 * Pi * (angle2_t + angle1_t)  + (angleIni_D + angleIni_B) + Pi)

    val Dx = Hx + HDx
    val Dy = Hy + HDy
    val Ex = Hx + HEx
    val Ey = Hy + HEy

    // F,G
    val IFx = rF * cos(2 * Pi * (angle3_t + angle1_t)  + (angleIni_F + (angleIni_B + Pi)))
    val IFy = rF * sin(2 * Pi * (angle3_t + angle1_t)  + (angleIni_F + (angleIni_B + Pi)))
    val IGx = rG * cos(2 * Pi * (angle3_t + angle1_t)  + (angleIni_F + (angleIni_B + Pi)) + Pi)
    val IGy = rG * sin(2 * Pi * (angle3_t + angle1_t)  + (angleIni_F + (angleIni_B + Pi)) + Pi)

    val Fx = Ix + IFx
    val Fy = Iy + IFy
    val Gx = Ix + IGx
    val Gy = Iy + IGy

    //(Bx,By,Cx,Cy,Dx,Dy,Ex,Ey,Fx,Fy,Gx,Gy,Hx,Hy,Ix,Iy)

    new DynamicalCurrentStateTransitory(time = t,Bx = Bx, By= By, Cx = Cx, Cy= Cy,
      Hx = Hx, Hy = Hy, Ix = Ix, Iy = Iy,
      HDx = HDx, HDy = HDy, HEx = HEx, HEy = HEy,
      Dx = Dx, Dy = Dy, Ex = Ex, Ey = Ey,
      IFx = IFx, IFy = IFy, IGx = IGx, IGy = IGy,
      Fx = Fx, Fy = Fy, Gx = Gx, Gy = Gy
    )
  }


  ////////////////////////////////
  //  TRANSITOIRE, NEXT
  ////////////////////////////////

  // fonctions utiles pour nextStepDynamic

  // calcul vitesse moteur
  def integrandVitesseMoteur(A:Double,C:Double,u:Double=>Double)(t:Double):Double ={
    A / C * u(t) * exp(t/C)
  }

  def slaveIntegralNewVitesseMoteur(stepsTransitoryNext:Int)(currentIntegral:Double,t:Double,deltaT:Double,f:Double=>Double) = {currentIntegral+ NumericalIntegration.integrate(f,t,t+deltaT,stepsTransitoryNext,simpson)}

  def computeNewMotorSpeed(C:Double, intergralTerm:Double,x:Double) = {exp(-x/C)* intergralTerm}


  // calcul angle
  def integrandAngle(stepsTransitoryNext:Int)(C:Double,integrandVitesseMoteur:Double=>Double,t:Double)(tt: Double) = { exp(-tt/C)* NumericalIntegration.integrate(integrandVitesseMoteur,t,tt,stepsTransitoryNext,simpson)}

  def computeNewAngle(stepsTransitoryNext:Int)(currentIntegralMotor:Double,t:Double,deltaT:Double, C:Double, integrandAngle:Double=>Double) = {currentIntegralMotor * C * (exp(-t/C)-exp(-(t+deltaT)/C)) + NumericalIntegration.integrate(integrandAngle,t,t+deltaT,stepsTransitoryNext,simpson)  }





  def nextStepDynamic(fixedParametersModel: FixedParametersModel)(parametersTransitoire: ParametersTransitoire)
                     (dynmicalCurrentState: DynamicalCurrentStateTransitory, deltaT:Double, stepsTransitoryNext:Int)={


    val t = dynmicalCurrentState.time


    // vitesse moteur
    // moteur 1
    val newTempIntegralVitesseMotor1 = slaveIntegralNewVitesseMoteur(stepsTransitoryNext)(dynmicalCurrentState.currentIntegralMoteur1,t,deltaT,integrandVitesseMoteur(fixedParametersModel.A1,fixedParametersModel.C1,parametersTransitoire.u1))
    val newVitesseMotor1 = computeNewMotorSpeed(fixedParametersModel.C1,newTempIntegralVitesseMotor1,t+deltaT)
    // moteur 2
    val newTempIntegralVitesseMotor2 = slaveIntegralNewVitesseMoteur(stepsTransitoryNext)(dynmicalCurrentState.currentIntegralMoteur2,t,deltaT,integrandVitesseMoteur(fixedParametersModel.A2,fixedParametersModel.C2,parametersTransitoire.u2))
    val newVitesseMotor2 = computeNewMotorSpeed(fixedParametersModel.C2,newTempIntegralVitesseMotor2,t+deltaT)
    // moteur 3
    val newTempIntegralVitesseMotor3 = slaveIntegralNewVitesseMoteur(stepsTransitoryNext)(dynmicalCurrentState.currentIntegralMoteur3,t,deltaT,integrandVitesseMoteur(fixedParametersModel.A3,fixedParametersModel.C3,parametersTransitoire.u3))
    val newVitesseMotor3 = computeNewMotorSpeed(fixedParametersModel.C3,newTempIntegralVitesseMotor3,t+deltaT)


    // angles
    // angle 1
    val tempAngle1 = computeNewAngle(stepsTransitoryNext)(dynmicalCurrentState.currentIntegralMoteur1,t,deltaT,fixedParametersModel.C1,
      integrandAngle(stepsTransitoryNext)(fixedParametersModel.C1,integrandVitesseMoteur(fixedParametersModel.A1,fixedParametersModel.C1,parametersTransitoire.u1),t))
    val newAngle1 = dynmicalCurrentState.angle1 + tempAngle1
    // angle 2
    val tempAngle2 = computeNewAngle(stepsTransitoryNext)(dynmicalCurrentState.currentIntegralMoteur2,t,deltaT,fixedParametersModel.C2,
      integrandAngle(stepsTransitoryNext)(fixedParametersModel.C2,integrandVitesseMoteur(fixedParametersModel.A2,fixedParametersModel.C2,parametersTransitoire.u2),t))
    val newAngle2 = dynmicalCurrentState.angle2 + tempAngle2
    // angle 3
    val tempAngle3 = computeNewAngle(stepsTransitoryNext)(dynmicalCurrentState.currentIntegralMoteur3,t,deltaT,fixedParametersModel.C3,
      integrandAngle(stepsTransitoryNext)(fixedParametersModel.C3,integrandVitesseMoteur(fixedParametersModel.A3,fixedParametersModel.C3,parametersTransitoire.u3),t))
    val newAngle3 = dynmicalCurrentState.angle3 + tempAngle3



    // Positions
    // B
    val newBx = dynmicalCurrentState.Bx * cos(2*Pi*tempAngle1) - dynmicalCurrentState.By * sin(2*Pi*tempAngle1)
    val newBy = dynmicalCurrentState.By * cos(2*Pi*tempAngle1) + dynmicalCurrentState.Bx * sin(2*Pi*tempAngle1)

    // C
    val newCx = dynmicalCurrentState.Cx * cos(2*Pi*tempAngle1) - dynmicalCurrentState.Cy * sin(2*Pi*tempAngle1)
    val newCy = dynmicalCurrentState.Cy * cos(2*Pi*tempAngle1) + dynmicalCurrentState.Cx * sin(2*Pi*tempAngle1)

    // H,I
    val newHx = dynmicalCurrentState.Hx * cos(2*Pi*tempAngle1) - dynmicalCurrentState.Hy * sin(2*Pi*tempAngle1)
    val newHy = dynmicalCurrentState.Hy * cos(2*Pi*tempAngle1) + dynmicalCurrentState.Hx * sin(2*Pi*tempAngle1)
    val newIx = dynmicalCurrentState.Ix * cos(2*Pi*tempAngle1) - dynmicalCurrentState.Iy * sin(2*Pi*tempAngle1)
    val newIy = dynmicalCurrentState.Iy * cos(2*Pi*tempAngle1) + dynmicalCurrentState.Ix * sin(2*Pi*tempAngle1)

    // D,E
    val newHDx = dynmicalCurrentState.HDx * cos(2*Pi*(tempAngle1+tempAngle2)) - dynmicalCurrentState.HDy * sin(2*Pi*(tempAngle1+tempAngle2))
    val newHDy = dynmicalCurrentState.HDy * cos(2*Pi*(tempAngle1+tempAngle2)) + dynmicalCurrentState.HDx * sin(2*Pi*(tempAngle1+tempAngle2))
    val newHEx = dynmicalCurrentState.HEx * cos(2*Pi*(tempAngle1+tempAngle2)) - dynmicalCurrentState.HEy * sin(2*Pi*(tempAngle1+tempAngle2))
    val newHEy = dynmicalCurrentState.HEy * cos(2*Pi*(tempAngle1+tempAngle2)) + dynmicalCurrentState.HEx * sin(2*Pi*(tempAngle1+tempAngle2))

    val newDx = newHx + newHDx
    val newDy = newHy + newHDy
    val newEx = newHx + newHEx
    val newEy = newHy + newHEy

    // F,G
    val newIFx = dynmicalCurrentState.IFx * cos(2*Pi*(tempAngle1+tempAngle3)) - dynmicalCurrentState.IFy * sin(2*Pi*(tempAngle1+tempAngle3))
    val newIFy = dynmicalCurrentState.IFy * cos(2*Pi*(tempAngle1+tempAngle3)) + dynmicalCurrentState.IFx * sin(2*Pi*(tempAngle1+tempAngle3))
    val newIGx = dynmicalCurrentState.IGx * cos(2*Pi*(tempAngle1+tempAngle3)) - dynmicalCurrentState.IGy * sin(2*Pi*(tempAngle1+tempAngle3))
    val newIGy = dynmicalCurrentState.IGy * cos(2*Pi*(tempAngle1+tempAngle3)) + dynmicalCurrentState.IGx * sin(2*Pi*(tempAngle1+tempAngle3))

    val newFx = newIx + newIFx
    val newFy = newIy + newIFy
    val newGx = newIx + newIGx
    val newGy = newIy + newIGy


    new DynamicalCurrentStateTransitory(time = t+deltaT,Bx = newBx, By=newBy, Cx = newCx, Cy=newCy,
      Hx = newHx, Hy = newHy, Ix = newIx, Iy = newIy,
      HDx = newHDx, HDy = newHDy, HEx = newHEx, HEy = newHEy,
      Dx = newDx, Dy = newDy, Ex = newEx, Ey = newEy,
      IFx = newIFx, IFy = newIFy, IGx = newIGx, IGy = newIGy,
      Fx = newFx, Fy = newFy, Gx = newGx, Gy = newGy,
      angle1=newAngle1, angle2=newAngle2, angle3=newAngle3,
      vitesseMoteur1= newVitesseMotor1, vitesseMoteur2=newVitesseMotor2, vitesseMoteur3=newVitesseMotor3,
      currentIntegralMoteur1= newTempIntegralVitesseMotor1, currentIntegralMoteur2=newTempIntegralVitesseMotor2,currentIntegralMoteur3=newTempIntegralVitesseMotor3
    )


  }



  def slaveSimuTransitoire(fixedParametersModel: FixedParametersModel)(parametersTransitoire: ParametersTransitoire)
                          (iter:Int,Nmax:Int,deltaT:Double,stepsTransitoryNext:Int, res:DynamicalCurrentStateTransitory):DynamicalCurrentStateTransitory = {

    if (iter == Nmax) res
    else slaveSimuTransitoire(fixedParametersModel)(parametersTransitoire)(iter + 1, Nmax, deltaT, stepsTransitoryNext, nextStepDynamic(fixedParametersModel)(parametersTransitoire)(res,deltaT,stepsTransitoryNext))
  }


  def slaveSimuTrajectoireTransitoire(fixedParametersModel: FixedParametersModel)(parametersTransitoire: ParametersTransitoire)
                          (iter:Int,Nmax:Int,deltaT:Double,stepsTransitoryNext:Int, res:Vector[DynamicalCurrentStateTransitory]):Vector[DynamicalCurrentStateTransitory] = {

    if (iter == Nmax) res
    else slaveSimuTrajectoireTransitoire(fixedParametersModel)(parametersTransitoire)(iter + 1, Nmax, deltaT, stepsTransitoryNext, res :+ nextStepDynamic(fixedParametersModel)(parametersTransitoire)(res.last,deltaT,stepsTransitoryNext))
  }


  def createInitialState(fixedParametersModel: FixedParametersModel)(initialConditions: InitialConditions):DynamicalCurrentStateTransitory={
    // Position
    // B,C
    val Bx = fixedParametersModel.rB * cos(initialConditions.angleIni_B)
    val By = fixedParametersModel.rB * sin(initialConditions.angleIni_B)
    val Cx = fixedParametersModel.rC * cos(initialConditions.angleIni_B+Pi)
    val Cy = fixedParametersModel.rC * sin(initialConditions.angleIni_B+Pi)

    // H,I
    val Hx = fixedParametersModel.rH * cos(initialConditions.angleIni_B)
    val Hy = fixedParametersModel.rH * sin(initialConditions.angleIni_B)
    val Ix = fixedParametersModel.rI * cos(initialConditions.angleIni_B+Pi)
    val Iy = fixedParametersModel.rI * sin(initialConditions.angleIni_B+Pi)

    // D,E
    val HDx = fixedParametersModel.rD  * cos((initialConditions.angleIni_D+ initialConditions.angleIni_B))
    val HDy = fixedParametersModel.rD  * sin((initialConditions.angleIni_D+ initialConditions.angleIni_B))
    val HEx = fixedParametersModel.rD  * cos((initialConditions.angleIni_D+ initialConditions.angleIni_B) +Pi)
    val HEy = fixedParametersModel.rD  * sin((initialConditions.angleIni_D+ initialConditions.angleIni_B) +Pi)

    val Dx = Hx + HDx
    val Dy = Hy + HDy
    val Ex = Hx + HEx
    val Ey = Hy + HEy

    // F,G
    val IFx = fixedParametersModel.rF  * cos((initialConditions.angleIni_F+ (initialConditions.angleIni_B+Pi)) )
    val IFy = fixedParametersModel.rF  * sin((initialConditions.angleIni_F+ (initialConditions.angleIni_B+Pi)))
    val IGx = fixedParametersModel.rG  * cos((initialConditions.angleIni_F+ (initialConditions.angleIni_B+Pi)) +Pi)
    val IGy = fixedParametersModel.rG  * sin((initialConditions.angleIni_F+ (initialConditions.angleIni_B+Pi)) +Pi)

    val Fx = Ix + IFx
    val Fy = Iy + IFy
    val Gx = Ix + IGx
    val Gy = Iy + IGy

    new DynamicalCurrentStateTransitory(time =0, Bx=Bx, By=By, Cx=Cx, Cy=Cy, Dx=Dx, Dy=Dy,
                                    Ex=Ex, Ey=Ey, Fx=Fx, Fy=Fy, Gx=Gx, Gy=Gy, Hx=Hx,Hy=Hy,
                                    Ix=Ix, Iy=Iy,
                                    HDx=HDx, HDy=HDy, HEx=HEx, HEy=HEy,
                                    IFx=IFx,IFy=IFy,IGx=IGx,IGy=IGy,
                                    angle1=initialConditions.angleIni_B, angle2= initialConditions.angleIni_D,angle3=initialConditions.angleIni_F,
                                    vitesseMoteur1=0.0, vitesseMoteur2=0.0, vitesseMoteur3=0.0,
                                    currentIntegralMoteur1=0.0, currentIntegralMoteur2=0.0, currentIntegralMoteur3=0.0)
  }



  // just the result at time t
  def simuTransitoire(fixedParametersModel: FixedParametersModel)(parametersTransitoire: ParametersTransitoire)
                     (Nmax:Int,deltaT:Double, initialConditions:InitialConditions,stepsTransitoryNext:Int):DynamicalCurrentStateTransitory = {

    val initialState = createInitialState(fixedParametersModel)(initialConditions)
    slaveSimuTransitoire(fixedParametersModel)(parametersTransitoire)(0, Nmax, deltaT, stepsTransitoryNext, initialState)

  }


  // all the trajectory: from 0 sec  to Nmax*deltaT sec  by deltaT
  def simuTrajectoireTransitoire(fixedParametersModel: FixedParametersModel)(parametersTransitoire: ParametersTransitoire)
                     (Nmax:Int,deltaT:Double, initialConditions:InitialConditions,stepsTransitoryNext:Int):Vector[DynamicalCurrentStateTransitory] = {

    val initialState = createInitialState(fixedParametersModel)(initialConditions)
    slaveSimuTrajectoireTransitoire(fixedParametersModel)(parametersTransitoire)(0, Nmax, deltaT, stepsTransitoryNext, Vector(initialState))

  }




  ////////////////////////////////
  //    Persistance rétinienne at time t
  ////////////////////////////////

  def figurePersistenceStationnaireAtTimeT(time:Double,dureePersistence:Double=DefaultValuesParameterModel.dureePersistence,nbPoints:Int=DefaultValuesParameterModel.nbPointsPersistence)(
    rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD ,
    rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
    rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI,
    angleH:Double = FixedParameterModel.angleH, angleI:Double = FixedParameterModel.angleI,
    angleIni_B:Double = DefaultValuesParameterModel.angleIni_B, angleIni_D:Double = DefaultValuesParameterModel.angleIni_D,
    angleIni_F:Double = DefaultValuesParameterModel.angleIni_F)(v1:Double=DefaultValuesParameterModel.v1, v2:Double=DefaultValuesParameterModel.v2, v3:Double=DefaultValuesParameterModel.v3)={

    val vecTimes = (0 until nbPoints).map(x => time + x * dureePersistence / nbPoints )  :+ (time+dureePersistence)
    val res = vecTimes.map(x => dynamicStationnary(rB,rC,rD,rE,rF,rG,rH,rI,angleIni_B,angleIni_D,angleIni_F)(v1,v2,v3)(x)).toVector
    res
    // type de res: Vector[DynamicalCurrentStateStationary]

  }







  ////////////////////////////////
  //    ODE solver
  ////////////////////////////////

  /*
  def integrate(f: (Double, Vector[Double]) => Vector[Double])(t0: Double, dt: Double, counter: Int, ysol: List[Vector[Double]]): List[Vector[Double]] = {
    def multiply(v: Vector[Double], s: Double) = v.map(_ * s)
    def divide(v: Vector[Double], s: Double) = v.map(_ / s)
    def add(vs: Vector[Double]*) = {
      def add0(v1: Vector[Double], v2: Vector[Double]) = (v1 zip v2).map { case(a, b) => a + b }
      vs.reduceLeft(add0)
    }

    val yn = ysol.head

    if (counter > 0) {
      val dy1 = multiply(f(t0, yn), dt)
      val dy2 = multiply(f(t0 + dt / 2, add(yn, divide(dy1, 2))), dt)
      val dy3 = multiply(f(t0 + dt / 2, add(yn, divide(dy2, 2))), dt)
      val dy4 = multiply(f(t0 + dt, add(yn, dy3)), dt)
      val y = add(yn, divide(add(dy1, multiply(dy2, 2), add(multiply(dy3, 2), dy4)), 6))::ysol
      val t = t0 + dt
      integrate(f)(t, dt, counter - 1, y)
    } else ysol.reverse
  }
*/




}




object NumericalIntegration {
  // adapted from https://rosettacode.org/wiki/Numerical_integration#Scala
  def leftRect(f:Double=>Double, a:Double, b:Double)=f(a)
  def midRect(f:Double=>Double, a:Double, b:Double)=f((a+b)/2)
  def rightRect(f:Double=>Double, a:Double, b:Double)=f(b)
  def trapezoid(f:Double=>Double, a:Double, b:Double)=(f(a)+f(b))/2
  def simpson(f:Double=>Double, a:Double, b:Double)=(f(a)+4*f((a+b)/2)+f(b))/6;

  /*
  def fn1(x:Double)=x*x*x
  def fn2(x:Double)=1/x
  def fn3(x:Double)=x
  */

  type Method = (Double=>Double, Double, Double) => Double
  def integrate(f:Double=>Double, a:Double, b:Double, steps:Int, m:Method)={
    val delta:Double=(b-a)/steps
    def ff(a:Double,b:Double,steps:Int)(h:Int)= {a+h*(b-a)/steps}
    val subdiv= List.tabulate(steps)(ff(a,b,steps))
    delta*subdiv.foldLeft(0.0)((s,x) => s+m(f, x, x+delta))
  }

}



object HiddenParameters{
  val stepsTransitory = 1000 // for numerical integration
  val stepsTransitoryNext = 30 // for numerical integration
}



object DefaultValuesParameterModel {

  // initial positions (angle)
  val  angleIni_B = 0.0
  val  angleIni_D = 1.0
  val  angleIni_F = 2.0


  // Persistence rétinienne stationnaire
  val dureePersistence = 1/25.toDouble
  val nbPointsPersistence = 50

  // vitesse (stationnaire)
  val v1 = 2
  val v2 = 1
  val v3 = 3

  // vitesse (transitoire) constante
  def u1(x: Double)={2.0}
  def u2(x: Double)={1.0}
  def u3(x: Double)={3.0}

  // vitesse (transitoire) sinuosidales
  def u1_sin(x: Double)={2.0*cos(x)}
  def u2_sin(x: Double)={1.0*cos(x)}
  def u3_sin(x: Double)={3.0*sin(x)}

  // Pour savoir si la lumière est alumée
  val lightB = true
  val lightC = true
  val lightD = true
  val lightE = true
  val lightF = true
  val lightG = true

}


object FixedParameterModel {

  // rayon
  val rB = 0.935
  val rC = 0.935
  val rD = 0.43
  val rE = 0.38
  val rF = 0.417
  val rG = 0.417
  val rH = sqrt(pow(0.50,2)+pow(0.05,2))
  val rI = sqrt(pow(0.50,2)+pow(0.05,2))
  val angleH = atan(0.05/0.5)
  val angleI = atan(0.05/0.5)

  // moteur 1
  // partie électrique
  val R1 = 1 // R est la résistance électrique interne du moteur (Ohm);
  val Ke1 = 1 // Ke est la constante de force électromotrice

  // partie mécanique
  val J1 = 1 // rotor vu comme un volant d'inertie J
  val Kc1 = 1 // Kc est la constante de couple
  val f1 = 1 // f est le coefficient de frottement visqueux.

  // moteur 2
  val R2 = 1
  val Ke2 = 1
  val J2 = 1
  val Kc2 = 1
  val f2 = 1

  // moteur 3
  val R3 = 1
  val Ke3 = 1
  val J3 = 1
  val Kc3 = 1
  val f3 = 1

  def computeA(Kc: Double, Ke: Double, R: Double, f: Double) = {
    Kc / (R * f + Kc * Ke)
  }

  def computeC(J: Double, Kc: Double, Ke: Double, R: Double, f: Double) = {
    J * R / (R * f + Kc * Ke)
  }

  val A1 = computeA(Kc1, Ke1, R1, f1)
  val C1 = computeC(J1, Kc1, Ke1, R1, f1)

  val A2 = computeA(Kc2, Ke2, R2, f2)
  val C2 = computeC(J2, Kc2, Ke2, R2, f2)

  val A3 = computeA(Kc3, Ke3, R3, f3)
  val C3 = computeC(J3, Kc3, Ke3, R3, f3)

}





  object Mesure{

    ////////////////////////////////
    //    STAT
    ////////////////////////////////

    def diff(v: Array[Int])= {
      if (v.isEmpty || v.tail.isEmpty) {
        Array[Int]()
      } else {
        (v zip v.tail).map { case (a, b) => b - a }
      }
    }

        def mean(v: Array[Double]): Option[Double]={
          if (v.isEmpty) None
          else Some(v.sum / v.length)
        }

        def variance(v:Array[Double]): Option[Double]={
          mean(v).flatMap{
            m => (mean(v.map{x=> pow(x-m,2)}))
          }
        }

        def momentOrdre3(v:Array[Double]): Option[Double] = {
          mean(v).flatMap{m =>
            variance(v).map(sqrt).flatMap{sd =>
              if(sd == 0) Some(0)
              else mean(v.map{x : Double => pow((x-m)/ sd,3)})
            }
          }
        }




    // general utils
        def calculAngle(x:Double,y:Double)= {
          val temp = math.atan(y / x)
          if (x >= 0 && y >= 0) temp
          else if (x < 0 && y >= 0) temp + Pi
          else if (x >= 0 && y < 0) 2 * Pi + temp
          else temp + Pi
        }





        ////////////////////////////////
        //    POINT SINGULIER
        ////////////////////////////////


        def distance(a: Double, b: Double)=  {sqrt( a*a + b*b)}

        def slaveFindSingularPoints(vec1:Vector[Double],vec2:Vector[Double], seuil:Double) ={
          vec1.zip(vec2).zipWithIndex.collect{ case((a,b),c) if distance(a,b)<seuil => (a,b,distance(a,b),c)}
          // le resultat c'est un vecteur de  (Doulble,Double,Double,Int)
        }



        // selon le seuil choisit sur la vitesse, on peut avoir plusieurs points singulier proche, alors qu'il n'y en a qu'un, on sélectionne celui
        // qui a la plus petite vitesse
        def nextSelectSingularPoints(vec:Vector[(Double,Double,Double,Int)],res:Vector[(Double,Double,Double,Int)], temp:(Double,Double,Double,Int)) : Vector[(Double,Double,Double,Int)]= {

          if (vec.isEmpty) { res
          } else {

            val temp2 = vec(0)

            if (temp2._4 == temp._4 + 1) {
              if (temp2._3 < res.last._3) {
                val newRes = res.dropRight(1) :+ temp2 // replace the last element with the one that has lower distance
                val newTemp = temp2
                val newVec = vec.tail
                nextSelectSingularPoints(newVec, newRes, newTemp)
              } else {
                // if temp2._3 >= res.last._3
                val newRes = res
                val newTemp = temp2
                val newVec = vec.tail
                nextSelectSingularPoints(newVec, newRes, newTemp)
              }

            } else {
              // if temp2._4 != temp._4 +1 ie start new consecutive sequence of index
              val newRes = res :+ temp2
              val newTemp = temp2
              val newVec = vec.tail
              nextSelectSingularPoints(newVec, newRes, newTemp)
            }

          }
        }


        def SelectSingularPoints(vec:Vector[(Double,Double,Double,Int)])={
          val res = Vector(vec(0))
          val temp = vec(0)
          nextSelectSingularPoints(vec.tail,res,temp)
        }


        // fonctions à executer sur une trajectoire:
        def countSingularPoints(vec1:Vector[Double],vec2:Vector[Double], seuil:Double) = {
          // return an integer
          val temp = slaveFindSingularPoints(vec1,vec2,seuil)
          if (temp.isEmpty){0} else {
            SelectSingularPoints(temp).length
          }
        }


        def timesOfSingularPoints(vec1:Vector[Double],vec2:Vector[Double], seuil:Double): Array[Int] = {
          // return a vector of integer (indices)
          val temp = slaveFindSingularPoints(vec1,vec2,seuil)
          if (temp.isEmpty){ Array[Int]()} else {
            SelectSingularPoints(temp).map(_._4).toArray
          }
        }


        def angleOfSingularPoints(vec1:Vector[Double],vec2:Vector[Double], seuil:Double): Array[Double] = {
          // return a vector of integer (indices)
          val temp = slaveFindSingularPoints(vec1, vec2, seuil)
          if (temp.isEmpty) {
            Array[Double]()
          } else {
            SelectSingularPoints(temp).map { case (a, b, c, d) => calculAngle(a, b) }.toArray
          }
        }




    // all trajectories
        def countSingularPointsAllTrajectories(res:TrajectoryStationnary, seuil:Double) = {
          val nbSingularPointsB = countSingularPoints(res.speedBx, res.speedBy,seuil)
          val nbSingularPointsC = countSingularPoints(res.speedCx, res.speedCy,seuil)
          val nbSingularPointsD = countSingularPoints(res.speedDx, res.speedDy,seuil)
          val nbSingularPointsE = countSingularPoints(res.speedEx, res.speedEy,seuil)
          val nbSingularPointsF = countSingularPoints(res.speedFx, res.speedFy,seuil)
          val nbSingularPointsG = countSingularPoints(res.speedGx, res.speedGy,seuil)

          nbSingularPointsB + nbSingularPointsC  + nbSingularPointsD + nbSingularPointsE + nbSingularPointsF + nbSingularPointsG
        }



        def angleSingularPointsAllTrajectories(res:TrajectoryStationnary, seuil:Double) = {

          val angleSingularPointsB = angleOfSingularPoints(res.speedBx, res.speedBy,seuil)
          val angleSingularPointsC = angleOfSingularPoints(res.speedCx, res.speedCy,seuil)
          val angleSingularPointsD = angleOfSingularPoints(res.speedDx, res.speedDy,seuil)
          val angleSingularPointsE = angleOfSingularPoints(res.speedEx, res.speedEy,seuil)
          val angleSingularPointsF = angleOfSingularPoints(res.speedFx, res.speedFy,seuil)
          val angleSingularPointsG = angleOfSingularPoints(res.speedGx, res.speedGy,seuil)

          angleSingularPointsB ++ angleSingularPointsC ++ angleSingularPointsD ++ angleSingularPointsE ++ angleSingularPointsF ++ angleSingularPointsG

        }


        def diffTimesSingularPointsAllTrajectories(res:TrajectoryStationnary, seuil:Double) = {
          val diffTimesSingularPointsB = diff(timesOfSingularPoints(res.speedBx, res.speedBy,seuil))
          val diffTimesSingularPointsC = diff(timesOfSingularPoints(res.speedCx, res.speedCy,seuil))
          val diffTimesSingularPointsD = diff(timesOfSingularPoints(res.speedDx, res.speedDy,seuil))
          val diffTimesSingularPointsE = diff(timesOfSingularPoints(res.speedEx, res.speedEy,seuil))
          val diffTimesSingularPointsF = diff(timesOfSingularPoints(res.speedFx, res.speedFy,seuil))
          val diffTimesSingularPointsG = diff(timesOfSingularPoints(res.speedGx, res.speedGy,seuil))

          diffTimesSingularPointsB ++ diffTimesSingularPointsC ++ diffTimesSingularPointsD ++ diffTimesSingularPointsE ++ diffTimesSingularPointsF ++ diffTimesSingularPointsG
        }







        ////////////////////////////////
        //    DENSITE
        ////////////////////////////////

        // utils
        def sum2Arrays(a:Array[Double],b:Array[Double])={
          a.zip(b).map(x=>x._1+x._2)
        }

        def sum2ArraysOfArrays(a:Array[Array[Double]],b:Array[Array[Double]])={
          a.zip(b).map(x => sum2Arrays(x._1,x._2))
        }

        def nextSumVectorOfArraysOfArrays(v:Vector[Array[Array[Double]]],acc:Array[Array[Double]]): Array[Array[Double]]= {
          if (v.isEmpty){acc} else{
            val temp = v(0)
            val newV = v.tail
            val newAcc = sum2ArraysOfArrays(acc,temp)
            nextSumVectorOfArraysOfArrays(newV,newAcc)
          }

        }

        def sumVectorOfArraysOfArrays(v:Vector[Array[Array[Double]]]) = {
          val acc = v(0)
          nextSumVectorOfArraysOfArrays(v.tail,acc)
        }


        def maxSquareForDensity(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD ,
                                    rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                                    rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI)= {
              max(max(rB, max(rH+rD, rH+rE)) , max(rC, max(rI+rF, rI+rG)) )
            }


        def minSquareForDensity(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD ,
                                    rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                                    rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI)= {
              min(min(rH-rD, rH-rE) , min(rI-rF, rI-rG) )
            }



        // trouver les indices (i,j) du point dans la grille
        // subdivsion des axes en N segments, donc grille de N^2 rectangles: indices entre 0 et N-1
        def findCaseForPoint(point:(Double,Double), xmin:Double,xmax:Double,ymin:Double,ymax:Double, N:Int) = {
          val deltaX = (xmax-xmin)/N
          val indiceX = if(point._1 == xmax){N-1} else {floor( (point._1 - xmin)/deltaX ).toInt }

          val deltaY = (ymax-ymin)/N
          val indiceY = if(point._2 == ymax){N -1} else {floor( (point._2 - ymin)/deltaY ).toInt }

          (indiceX,indiceY)

        }


        def arrayDensitySquare(vec1:Vector[Double],vec2:Vector[Double],xmin:Double,xmax:Double,ymin:Double,ymax:Double,N:Int) = {

          val temp = vec1.zip(vec2)

          var array = Array.ofDim[Double](N,N)  // initialisé à 0

          for (point <- temp) {

            val (i, j) = findCaseForPoint((point._1, point._2), xmin, xmax, ymin, ymax, N)
            array(i)(j) += 1.0
          }

          //println( array.map(x => x.sum).sum ) // on retrouve le nombre de points, ouf!
          //println(convertResultStationnary(res).time.length)
          array
          // to print the array array.map(_.mkString(" ")).mkString("\n")
        }


        def nbPointGrilleDansZoneAtteignable(xmin:Double,xmax:Double,ymin:Double,ymax:Double,N:Int,rMin:Double,rMax:Double)={
          val deltaX = (xmax-xmin)/N
          val deltaY = (ymax-ymin)/N
          val posX = (0 to N).map(x=> xmin + x*deltaX)
          val posY = (0 to N).map(x=> ymin + x*deltaY)

          var count = 0
          for ( x <-  posX){
            for ( y <- posY){
              if  ( distance(x,y) >= rMin  &  distance(x,y) <= rMax ) { count += 1}
            }
          }
          count
        }


        // fonction à executer sur une trajectoire:
        def pointsDensity(vec1:Vector[Double], vec2:Vector[Double], xmin:Double, xmax:Double, ymin:Double, ymax:Double, N:Int, rMin:Double, rMax:Double) = {

          // nombre de cases non vides dans l'array
          val array = arrayDensitySquare(vec1,vec2,xmin,xmax,ymin,ymax,N)
          val nbCasesNonEmpty = array.map(x => x.filter(_ != 0).length).sum

          // la ratio d'occupation:
          //(nbCasesNonEmpty.toFloat / (N * N)).toDouble   // dans un carre
          (nbCasesNonEmpty.toFloat / nbPointGrilleDansZoneAtteignable(xmin,xmax,ymin,ymax,N,rMin,rMax) ).toDouble  // couronne (à l'intérieur de 2 cercles)
        }


        def pointsDensitySquare(vec1:Vector[Double], vec2:Vector[Double], xmin:Double, xmax:Double, ymin:Double, ymax:Double, N:Int) = {

          // nombre de cases non vides dans l'array
          val array = arrayDensitySquare(vec1,vec2,xmin,xmax,ymin,ymax,N)
          val nbCasesNonEmpty = array.map(x => x.filter(_ != 0).length).sum

          // la ratio d'occupation:
          (nbCasesNonEmpty.toFloat / (N * N)).toDouble   // dans un carre

        }




        // all trajectories
        def densityAllTrajectories(res:TrajectoryStationnary, xmin:Double, xmax:Double, ymin:Double, ymax:Double, N:Int, rMin:Double, rMax:Double) = {

          val arrayB = arrayDensitySquare(res.Bx,res.By,xmin,xmax,ymin,ymax,N)
          val arrayC = arrayDensitySquare(res.Cx,res.Cy,xmin,xmax,ymin,ymax,N)
          val arrayD = arrayDensitySquare(res.Dx,res.Dy,xmin,xmax,ymin,ymax,N)
          val arrayE = arrayDensitySquare(res.Ex,res.Ey,xmin,xmax,ymin,ymax,N)
          val arrayF = arrayDensitySquare(res.Fx,res.Fy,xmin,xmax,ymin,ymax,N)
          val arrayG = arrayDensitySquare(res.Gx,res.Gy,xmin,xmax,ymin,ymax,N)

          val temp = Vector(arrayB,arrayC,arrayD,arrayE,arrayF,arrayG)
          val array = sumVectorOfArraysOfArrays(temp)


          // nombre de cases non vides dans l'array (sommes des array)
          val nbCasesNonEmpty = array.map(x => x.filter(_ != 0).length).sum


          // la ratio d'occupation:
          //(nbCasesNonEmpty.toFloat / (N * N)).toDouble
          ( nbCasesNonEmpty.toFloat  / nbPointGrilleDansZoneAtteignable(xmin,xmax,ymin,ymax,N,rMin,rMax) ).toDouble

        }


        def densityAllTrajectoriesSquare(res:TrajectoryStationnary, xmin:Double, xmax:Double, ymin:Double, ymax:Double, N:Int) = {

          val arrayB = arrayDensitySquare(res.Bx,res.By,xmin,xmax,ymin,ymax,N)
          val arrayC = arrayDensitySquare(res.Cx,res.Cy,xmin,xmax,ymin,ymax,N)
          val arrayD = arrayDensitySquare(res.Dx,res.Dy,xmin,xmax,ymin,ymax,N)
          val arrayE = arrayDensitySquare(res.Ex,res.Ey,xmin,xmax,ymin,ymax,N)
          val arrayF = arrayDensitySquare(res.Fx,res.Fy,xmin,xmax,ymin,ymax,N)
          val arrayG = arrayDensitySquare(res.Gx,res.Gy,xmin,xmax,ymin,ymax,N)

          val temp = Vector(arrayB,arrayC,arrayD,arrayE,arrayF,arrayG)
          val array = sumVectorOfArraysOfArrays(temp)


          // nombre de cases non vides dans l'array (sommes des array)
          val nbCasesNonEmpty = array.map(x => x.filter(_ != 0).length).sum


          // la ratio d'occupation:
          (nbCasesNonEmpty.toFloat / (N * N)).toDouble

        }






        ////////////////////////////////
        //   RETOUR D'UN POINT (boucle)
        ////////////////////////////////

    /*

        def slavePremierRetour(point:(Double,Double,Int), vect:Vector[(Double,Double,Int)], seuil:Double) = {
          val temp = vect.filter(_._3 >= point._3)
          // on ne garde que les éléments du vecteur dont le point a un
          //  indice supérieur qu le point de référence (qui sera aussi dans le vecteur)
          val res2 = comparePointsToAllElemnts2(point,temp) // calcule les distances, on a toujours l'indice  vector[(Double,Int)] et le temps de retour

          val res3 = res2.filter(_._1 < seuil)


          val indices = res3.map(_._2)
          val a = indices.tail
          val b = indices.dropRight(1)   // tous sauf le dernier
          val res4 = a.zip(b).collect{ case(x,y) if (x-y)>1  => x } // permet de ne pas avoir les points qui collent notre point de référence
          //res4

          if (res4.length >=1 ){ (point,res4(0)-point._3) } else { (point,0)}
          // on retourne le point et le temps de retour associé (s'il y en a un , sinon le point et 0)

        }


        def premierRetour(v1:Vector[Double], v2:Vector[Double], seuil:Double) = {
          // l'indice du point qui est concerné par le retour, et le nombre de pas de temps pour y revenir
          val temp = v1.zip(v2)
          val res3 = temp.zipWithIndex.map( x =>  (x._1._1, x._1._2, x._2) ) // pour avoir le format (x:Double,y:Double,indice:Int)
          val res4 = res3.map(x => slavePremierRetour(x,res3, seuil))
          res4
          // le format du retour: un vecteur (nombre de composante = nb de pas de temps = nb de positions du point): pour chaque position
          // renvoit  ((x,y,i),t)  x,y position du point pour lequel il y a un retour, i l'indice auquel il est parcouru, et t le nombre de pas de temps pour le retour
        }


        def indiceEtTempsPremierRetour(v1:Vector[Double], v2:Vector[Double], seuil:Double) : (Vector[Int],Vector[Int]) = {
          val resPremierRetour = premierRetour(v1,v2,seuil)
          //val res4 = resPremierRetour .filter(_ != () )  // retire les vecteur vide
          //val res4 = resPremierRetour.filterNot(_ == ())
          val res4 = resPremierRetour.filterNot(_._2 ==0 )
          val tempsRetour = res4.map{ case( ((a,b,c),d) )  => d}  //  le temps qu'on met pour revenir au point concerné
          val indicesRetour = res4.map{ case( ((a,b,c),d) )  => c}   // l'indice du point concerné
          (indicesRetour,tempsRetour)
        }

*/
        // version segments Guillaume

        // utils
        def distanceLoopPoint( x: (Double,Double), y: (Double,Double) )=  {sqrt( (x._1 - y._1)*(x._1 - y._1) + (x._2 - y._2)*(x._2 - y._2)  )}

        // on garde l'indice du point de référence
        def distanceLoopPoint2( x: (Double,Double,Int), y: (Double,Double,Int) )=  {(sqrt( (x._1 - y._1)*(x._1 - y._1) + (x._2 - y._2)*(x._2 - y._2) ) , x._3, y._3-x._3)}

        def comparePointsToAllElemnts(point:(Double,Double), vect: Vector[(Double,Double)])={
          vect.map(x => distanceLoopPoint(x,point))
        }

        def comparePointsToAllElemnts2(point:(Double,Double,Int), vect: Vector[(Double,Double,Int)])={
          vect.map(x => distanceLoopPoint2(x,point))
        }



        def intersection(A:(Double,Double),B:(Double,Double),C:(Double,Double),D:(Double,Double)) : Boolean = {
          val ABx = B._1 - A._1
          val ABy = B._2 - A._2
          val CDx = D._1 - C._1
          val CDy = D._2 - C._2
          val det = -ABx * CDy + ABy * CDx

          if (det == 0) {
            false
          } else {
            val k = (-CDy * (C._1 - A._1) + CDx * (C._2 - A._2)) / det
            val t = (-ABy * (C._1 - A._1) + ABx * (C._2 - A._2)) / det
            if ((0 <= k) & (k <= 1) & (0 <= t) & (t <= 1)  ) {
              true
            } else {
              false
            }
          }

        }


        def nextRetour(trajCurrent : Vector[((Double,Double),Int)], res:Vector[(Double,Double,Int,Int)]) : Vector[(Double,Double,Int,Int)] = {
          if (trajCurrent.isEmpty || trajCurrent.tail.isEmpty || trajCurrent.tail.tail.isEmpty) {
            res
          } else {
            val segments = trajCurrent.dropRight(1).zip(trajCurrent.tail) // le deuxieme vecteur a forcément une intersection avec le premier
            val temp = segments(0)
            val segments2 = segments.tail.tail
            val resIntersect = segments2.map { x => (intersection(temp._1._1, temp._2._1, x._1._1, x._2._1), x._1._1._1,x._1._1._2 , x._1._2 )}.filter(_._1 == true)
            if (resIntersect.isEmpty) {
              val newTrajCurrent = trajCurrent.tail
              val newRes = res
              nextRetour(newTrajCurrent, newRes)
            } else {
              val newTrajCurrent = trajCurrent.tail
              val newRes = res :+ ( resIntersect(0)._2 , resIntersect(0)._3, trajCurrent(0)._2, resIntersect(0)._4 - trajCurrent(0)._2)
              nextRetour(newTrajCurrent, newRes)
            }
          }
        }





        def retour(traj: Vector[(Double,Double)]) = {
          val trajCurrent = traj.zipWithIndex
          val res = Vector[(Double,Double,Int,Int)]()
          nextRetour(trajCurrent,res)
        }



        // fonctions à executer sur une trajectoire:
        // Le seuil c'est pour ne pas avoir la période, pour les trajectoires périodiques, mais les boucle locales
        def numberRetour(v1:Vector[Double], v2:Vector[Double], seuil:Int)={
          val resPremierRetour = retour( v1.zip(v2) )
          val nbRetour= resPremierRetour.map(x=>x._4).length
          nbRetour
        }


        def indiceRetour(v1:Vector[Double], v2:Vector[Double], seuil:Int)={

          val resPremierRetour = retour( v1.zip(v2) ).filter(_._4 < seuil)
          val indiceRetour= resPremierRetour.map(x=>x._3)
          indiceRetour
        }


        def tempsRetour(v1:Vector[Double], v2:Vector[Double], seuil:Int)={

          val resPremierRetour = retour( v1.zip(v2) ).filter(_._4 < seuil)
          val tempsRetour= resPremierRetour.map(x=>x._4)
          tempsRetour
        }


        def angleRetour(v1:Vector[Double], v2:Vector[Double], seuil:Int)={

          val resPremierRetour = retour( v1.zip(v2) ).filter(_._4 < seuil)
          val angleRetour = resPremierRetour .map(x=>x._1). zip(resPremierRetour .map(x=>x._2)).map{ x => calculAngle(x._1,x._2)}
          angleRetour
        }




        def indiceEtTempsPremierRetourEtAngle(v1:Vector[Double], v2:Vector[Double]) : (Array[Int],Array[Int],Array[Double]) = {
          val resPremieretourSegments = retour( v1.zip(v2) )
          val indicesRetour  = resPremieretourSegments.map(x=>x._3) // l'indice du point concerné => sa position
          val tempsRetour= resPremieretourSegments.map(x=>x._4)   // le temps de retour du point concerné
          val angleRetour =  resPremieretourSegments.map(x=>x._1). zip(resPremieretourSegments.map(x=>x._2)).map{ x => calculAngle(x._1,x._2)}
          (indicesRetour.toArray,tempsRetour.toArray,angleRetour.toArray)
        }



        // all trajectories
        def resPremieretour(v1:Vector[Double], v2:Vector[Double]) = {retour( v1.zip(v2) ) }


        def numberRetourAllTrajectories(res:TrajectoryStationnary, seuil:Int)={

          val resRetourB = resPremieretour(res.Bx,res.By).filter(_._4 < seuil)
          val resRetourC = resPremieretour(res.Cx,res.Cy).filter(_._4 < seuil)
          val resRetourD = resPremieretour(res.Dx,res.Dy).filter(_._4 < seuil)
          val resRetourE = resPremieretour(res.Ex,res.Ey).filter(_._4 < seuil)
          val resRetourF = resPremieretour(res.Fx,res.Fy).filter(_._4 < seuil)
          val resRetourG = resPremieretour(res.Gx,res.Gy).filter(_._4 < seuil)

          val nbRetourB= resRetourB.map(x=>x._4).length
          val nbRetourC= resRetourC.map(x=>x._4).length
          val nbRetourD= resRetourD.map(x=>x._4).length
          val nbRetourE= resRetourE.map(x=>x._4).length
          val nbRetourF= resRetourF.map(x=>x._4).length
          val nbRetourG= resRetourG.map(x=>x._4).length

          (nbRetourB + nbRetourC + nbRetourD + nbRetourE + nbRetourF + nbRetourG )
        }


        def tempsRetourAllTrajectories(res:TrajectoryStationnary, seuil:Int)={

              val resRetourB = resPremieretour(res.Bx,res.By).filter(_._4 < seuil)
              val resRetourD = resPremieretour(res.Dx,res.Dy).filter(_._4 < seuil)
              val resRetourF = resPremieretour(res.Fx,res.Fy).filter(_._4 < seuil)

              val tempsRetourB= resRetourB.map(x=>x._4)
              val tempsRetourD= resRetourD.map(x=>x._4)
              val tempsRetourF= resRetourF.map(x=>x._4)

              (tempsRetourB ++ tempsRetourD ++ tempsRetourF)
            }



        def angleRetourAllTrajectories(res:TrajectoryStationnary, seuil:Int)={

          val resRetourB = resPremieretour(res.Bx,res.By).filter(_._4 < seuil)
          val resRetourD = resPremieretour(res.Dx,res.Dy).filter(_._4 < seuil)
          val resRetourF = resPremieretour(res.Fx,res.Fy).filter(_._4 < seuil)

          val angleRetourB= resRetourB.map(x=>x._1). zip(resRetourB.map(x=>x._2)).map{ x => calculAngle(x._1,x._2)}
          val angleRetourD= resRetourD.map(x=>x._1). zip(resRetourD.map(x=>x._2)).map{ x => calculAngle(x._1,x._2)}
          val angleRetourF= resRetourF.map(x=>x._1). zip(resRetourF.map(x=>x._2)).map{ x => calculAngle(x._1,x._2)}

          angleRetourB ++ angleRetourD ++ angleRetourF
        }









        ////////////////////////////////
        //    COURBURE de la l'arc paramétré
        ////////////////////////////////


        // formule: x' y'' - y' x'' / ( x'^2 + y' ^2 )^(3/2)
        def slaveCourbure(sx:Double, sy:Double, acx:Double, acy:Double)={
          ((sx*acy) - (sy*acx)) / pow( sx*sx + sy*sy , 3/2)
        }



        def courbure(speedX:Vector[Double], speedY:Vector[Double], accX:Vector[Double], accY:Vector[Double])={

          val temp = speedX.zip(speedY).zip(accX).zip(accY).map( x  => (x._1._1._1,x._1._1._2,x._1._2,x._2))  // pour avoir autre format
          temp.map(x=>  slaveCourbure(x._1,x._2,x._3,x._4))
        }


        // all trajectories
        def courbureAllTrajectories(res:TrajectoryStationnary)={
          val courburesB = courbure(res.speedBx,res.speedBy,res.accBx,res.accBy)
          val courburesC = courbure(res.speedCx,res.speedCy,res.accCx,res.accCy)
          val courburesD = courbure(res.speedDx,res.speedDy,res.accDx,res.accDy)
          val courburesE = courbure(res.speedEx,res.speedEy,res.accEx,res.accEy)
          val courburesF = courbure(res.speedFx,res.speedFy,res.accFx,res.accFy)
          val courburesG = courbure(res.speedGx,res.speedGy,res.accGx,res.accGy)

          courburesB ++ courburesC ++ courburesD ++ courburesE ++ courburesF ++ courburesG

        }


    ////////// pas encore utilisé mais pour voir si il y a des lignes droites à partir de la courbure

    // a partir du zip des courbures et index selectionnés par seuil (ici on n'utilise que les index)
    // on repère les point consécutifs qui sont dans ce vecteur
    // on garde un vector de (Int,Int) le premier Int et le premier indice de la suite (valeurs d'indice consécutives)
    // le second est la longeur de cette suite
    def nextSlaveDetectStraitLine(v:Vector[(Double,Int)],acc:Vector[(Int,Int)],temp:(Double,Int)):Vector[(Int,Int)]={
      if (v.isEmpty ){ acc
      } else {
        val temp2 = v(0)
        if (temp2._2 == temp._2+1){
          val newV = v.tail
          val newTemp = temp2
          val newAcc = acc.dropRight(1) :+ (acc.last._1,acc.last._2 +1)
          nextSlaveDetectStraitLine(newV,newAcc,newTemp)
        } else {
          val newV = v.tail
          val newTemp = temp2
          val newAcc = acc :+ (temp2._2,1)
          nextSlaveDetectStraitLine(newV,newAcc,newTemp)
        }
      }
    }




    def detectStraitLine(courbures:Vector[Double], seuil:Double)={
      val temp = courbures.zipWithIndex.filter(_._1 < seuil)
      // temp
      if (temp.isEmpty){
        Vector[(Int,Int)]()
      }else {
        val newV = temp.tail
        val newTemp = temp(0)
        val newAcc = Vector((temp(0)._2, 1))
        nextSlaveDetectStraitLine(newV, newAcc, newTemp)
      }
    }









    ////////////////////////////////
        //    MORAN
        ////////////////////////////////

        // pour avoir un array de array, pour un point

        def convertFromMoran(v1:Vector[Double], v2:Vector[Double])={
          v1.zip(v2).map(x => Array(x._1,x._2) ).toArray
        }

        //GridMorphology.moranDirect()  // grid
        //Spatstat.moran()


        def moranDirectTraj(vec1:Vector[Double],vec2:Vector[Double],xmin:Double,xmax:Double,ymin:Double,ymax:Double,N:Int)={
          val temp = arrayDensitySquare(vec1,vec2,xmin,xmax,ymin,ymax,N) //: Array[Array[Double]]
          GridMorphology.moranDirect(temp)
        }


        def moranTraj(vec1:Vector[Double],vec2:Vector[Double],xmin:Double,xmax:Double,ymin:Double,ymax:Double,N:Int)={
          if (vec1.isEmpty){1.0} else {
            val temp = arrayDensitySquare(vec1, vec2, xmin, xmax, ymin, ymax, N)
            GridMorphology.moran(temp)
          }
        }


        def moranAllTrajectories(res:TrajectoryStationnary, xmin:Double, xmax:Double, ymin:Double, ymax:Double, N:Int)= {

          if (res.Bx.isEmpty && res.Cx.isEmpty && res.Dx.isEmpty && res.Ex.isEmpty && res.Fx.isEmpty && res.Gx.isEmpty) {
            1.0
          } else {
            val arrayB = arrayDensitySquare(res.Bx, res.By, xmin, xmax, ymin, ymax, N)
            val arrayC = arrayDensitySquare(res.Cx, res.Cy, xmin, xmax, ymin, ymax, N)
            val arrayD = arrayDensitySquare(res.Dx, res.Dy, xmin, xmax, ymin, ymax, N)
            val arrayE = arrayDensitySquare(res.Ex, res.Ey, xmin, xmax, ymin, ymax, N)
            val arrayF = arrayDensitySquare(res.Fx, res.Fy, xmin, xmax, ymin, ymax, N)
            val arrayG = arrayDensitySquare(res.Gx, res.Gy, xmin, xmax, ymin, ymax, N)

            val temp = Vector(arrayB, arrayC, arrayD, arrayE, arrayF, arrayG)
            val array = sumVectorOfArraysOfArrays(temp)
            GridMorphology.moran(array)
          }
        }



        //println(res3.map(_.mkString(" ")).mkString("\n"))


        ////////////////////////////////
        //    utils pour comparer une trajectoire simuée et une désirée
        ////////////////////////////////

        // on souhaite s'affranchir du temps, on s'intéresse seulement à la trace, l'image à la fin de la simu
        // on de garde sur un trajectoire que des points distant d'un distance fixée

        def cumSum[A](xs: Seq[A])(implicit num: Numeric[A]): Seq[A] = {
          xs.tail.scanLeft(xs.head)(num.plus)
        }



        def nextSimplifiedTrajectory (current:Vector[(Double,Double)],res:Vector[(Double,Double)], distance:Double):Vector[(Double,Double)] = {
          if (current.isEmpty) {
            res
          } else {
            val temp = current(0)
            val temp2 = comparePointsToAllElemnts(temp, current)
            val temp3 = current.zip(cumSum(temp2))
            val temp4 = temp3.filter(_._2 > distance)
            if (temp4.isEmpty){
              val newCurrent = Vector[(Double,Double)]()
              val newRes = res :+ temp3.last._1
              nextSimplifiedTrajectory(newCurrent, newRes, distance)
            } else {
              val newCurrent = temp4.map(x=>x._1)
              val newRes = res :+ temp4(0)._1
              nextSimplifiedTrajectory(newCurrent, newRes, distance)
            }
          }
        }


        def simplifiedTrajectory(v1:Vector[Double], v2:Vector[Double], distance:Double)={
          val current = v1.zip(v2)
          val res = Vector(current(0))
          nextSimplifiedTrajectory(current,res,distance)
        }

        def createCircle(rayon:Double, nbPoints:Int)={
          val temp = (0 to nbPoints).map(x=> (x*2*Pi/nbPoints)).toArray
          (temp.map(x=> rayon*cos(x)), temp.map(x=> rayon*sin(x)) )
        }


        // une distance simulée, l'autre construite, pas nécessairement de la même longueur
        def slaveDistanceTwoTrajectories(v1:Vector[(Double,Double)],v2:Vector[(Double,Double)])={
          //val v1 = u_x.zip(u_y)
          //val v2 = v_x.zip(v_y)
          val minSize = min(v1.length,v2.length)

          val v1Bis = v1.take(minSize)
          val v2Bis = v2.take(minSize)
          v1Bis.zip(v2Bis).map(x=>distanceLoopPoint(x._1,x._2)).sum

        }



        def distanceTwoTrajectories(u_x:Array[Double],u_y:Array[Double],v_x:Array[Double],v_y:Array[Double],distance:Double)= {
          val u =simplifiedTrajectory(u_x.toVector,u_y.toVector,distance)
          val v =simplifiedTrajectory(v_x.toVector,v_y.toVector,distance)
          slaveDistanceTwoTrajectories(u,v)
        }



      }




// on ne garde que les mesures nécessaires pour PSE
object MesurePSE{


  def meanCourbure(speedX:Vector[Double],speedY:Vector[Double],accX:Vector[Double],accY:Vector[Double]): Double ={
    mean(courbure(speedX,speedY,accX,accY).toArray).getOrElse(0.0)
  }

  def meanCourbureAllTrajectories(res:TrajectoryStationnary): Double ={
    mean(courbureAllTrajectories(res).toArray).getOrElse(0.0)
  }


}













