package electra


//import electra.Model.integrate
import java.awt.geom.Point2D

import electra.NumericalIntegration._
import electra.Model._
import electra.Mesure._

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

  // trajectoire stationnaire


  val T = 2.0
  val deltaT = 0.0005
  val res = dynamicTrajectoryStationnary()(v1=2.0,v2=(-6.0))(T,deltaT)
  val res2 = convertResultStationnary(res)



  ////////////////////////////////
  //    MESURES. points singuliers
  ////////////////////////////////

 /*
  val seuilPointSingulier = 4
  val numberPointSinguliersD = countSingularPoints(res2.speedDx,res2.speedDy,seuilPointSingulier)
  println(numberPointSinguliersD)

  val timesPointSinguliersD = timesOfSingularPoints(res2.speedDx,res2.speedDy,seuilPointSingulier)
  println(timesPointSinguliersD)
*/


  ////////////////////////////////
  //    MESURES. densité
  ////////////////////////////////

/*
    val N = 30  // pas de la subdiviion du carré
    val xmax = maxSquareForDensity()
    val xmin = -xmax
    val ymax = xmax
    val ymin = xmin


    val densiteB =  pointsDensitySquare(res2.Bx,res2.By,xmin,xmax,ymin,ymax,N)
    //println(densiteB)

    val totalDensity = allTrajectoriesDensitySquare(res2,xmin,xmax,ymin,ymax,N)
    //println(totalDensity)
*/




  ////////////////////////////////
  //    MESURES. Courbure
  ////////////////////////////////

/*
  // point B
  val courburesB = courbure(res2.speedBx,res2.speedBy,res2.accBx,res2.accBy)
  //println(courburesB)
  //println(mean(courburesB))

  // point D
  val courburesD = courbure(res2.speedDx,res2.speedDy,res2.accDx,res2.accDy)
  //println(courburesD)
  // pourquoi c'est presque constant ?


  //println(res2.speedDx)

  //val seuilCourbure = 13.0
  //val indicesLigneDroite = detectStraitLine(courburesB,seuilCourbure)
  //print(indicesLigneDroite)
*/




  ////////////////////////////////
  //    MESURES. loop points: temps de premier retour
  ////////////////////////////////

/*
    // ici on ne garde que les temps de premier retour (et pas les points qui y sont associés)
    // voir la fonction premierRetour pour avoir toutes ces infos
    val seuilLoop = 0.01
    val tempsRetourB = tempsPremierRetour(res2.Bx,res2.By,seuilLoop)
    //println(tempsRetourB)
    val tempsRetourD = tempsPremierRetour(res2.Dx,res2.Dy,seuilLoop)
    println(tempsRetourD)
*/





  ////////////////////////////////
  //    MESURES. Moran
  ////////////////////////////////

/*
    val piMoranB = convertFromMoran(res2.Bx,res2.Bx)
    //println(piMoranB .map(_.mkString(" ")).mkString("\n"))

    val xForMoran = List.fill(piMoranB.length)(1.0).toArray
    //println(xForMoran.mkString(""))

    val resMoranB = Spatstat.moran(piMoranB,xForMoran)
    println(resMoranB)
    // NaN ?

    // test avec des autres valeurs
    val piTestMoran = Array( Array(1.0,2.0),Array(3.0,4.0))
    val xMoran = List.fill(2)(1.0).toArray
    val resMoran = Spatstat.moran(piTestMoran,xMoran)
    //println(resMoran)
*/

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

  def dynamicStationnary(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD ,
                         rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                         rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI,
                         angleIni_B:Double = DefaultValuesParameterModel.angleIni_B, angleIni_D:Double = DefaultValuesParameterModel.angleIni_D,
                         angleIni_F:Double = DefaultValuesParameterModel.angleIni_F)(v1:Double=DefaultValuesParameterModel.v1, v2:Double=DefaultValuesParameterModel.v2, v3:Double=DefaultValuesParameterModel.v3)(t:Double) = {

    // Position
    // B,C
    val Bx = rB * cos(2*Pi*v1*t + angleIni_B)
    val By = rB * sin(2*Pi*v1*t + angleIni_B)
    val Cx = rC * cos(2*Pi*v1*t + angleIni_B+Pi)
    val Cy = rC * sin(2*Pi*v1*t + angleIni_B+Pi)

    // H,I
    val Hx = rH * cos(2*Pi*v1*t + angleIni_B)
    val Hy = rH * sin(2*Pi*v1*t + angleIni_B)
    val Ix = rI * cos(2*Pi*v1*t + angleIni_B+Pi)
    val Iy = rI * sin(2*Pi*v1*t + angleIni_B+Pi)

    // D,E
    val HDx = rD  * cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
    val HDy = rD  * sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
    val HEx = rE  * cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +Pi)
    val HEy = rE  * sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +Pi)

    val Dx = Hx + HDx
    val Dy = Hy + HDy
    val Ex = Hx + HEx
    val Ey = Hy + HEy

    // F,G
    val IFx = rF  * cos(2*Pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+Pi)) )
    val IFy = rF  * sin(2*Pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+Pi)))
    val IGx = rG  * cos(2*Pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+Pi)) +Pi)
    val IGy = rG  * sin(2*Pi*(v3+v1)*t + (angleIni_F+ (angleIni_B+Pi)) +Pi)

    val Fx = Ix + IFx
    val Fy = Iy + IFy
    val Gx = Ix + IGx
    val Gy = Iy + IGy

    ////////////////////////////////
    // speed
    // B,C
    val speedBx = -rB * 2*Pi*v1* sin(2*Pi*v1*t + angleIni_B)
    val speedBy = rB * 2*Pi*v1* cos(2*Pi*v1*t + angleIni_B)
    val speedCx = -rC * 2*Pi*v1* sin(2*Pi*v1*t + angleIni_B+Pi)
    val speedCy = rC * 2*Pi*v1* cos(2*Pi*v1*t + angleIni_B+Pi)

    // D,E
    val speedDx = -rH * 2*Pi*v1* sin(2*Pi*v1*t + angleIni_B)  - rD *2*Pi*(v2+v1)* sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
    val speedDy = rH * 2*Pi*v1*  cos(2*Pi*v1*t + angleIni_B)  + rD *2*Pi*(v2+v1)* cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
    val speedEx = -rH * 2*Pi*v1* sin(2*Pi*v1*t + angleIni_B)  - rE *2*Pi*(v2+v1)* sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +Pi)
    val speedEy = rH * 2*Pi*v1*  cos(2*Pi*v1*t + angleIni_B)  + rE *2*Pi*(v2+v1)* cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +Pi)

    // F,G
    val speedFx = -rI * (2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B+Pi)   - rF* (2*Pi*(v3+v1))* sin(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi))
    val speedFy = rI * (2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B+Pi)    + rF* (2*Pi*(v3+v1))* cos(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi))
    val speedGx = -rI * (2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B+Pi)   - rG* (2*Pi*(v3+v1))* sin(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi) +Pi)
    val speedGy = rI * (2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B+Pi)    + rG* (2*Pi*(v3+v1))* cos(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi) +Pi)


    ////////////////////////////////
    // acceleration
    // B,C
    val accBx = -rB * (2*Pi*v1)*(2*Pi*v1) * cos(2*Pi*v1*t + angleIni_B)
    val accBy = -rB * (2*Pi*v1)*(2*Pi*v1) * sin(2*Pi*v1*t + angleIni_B)
    val accCx = -rC * (2*Pi*v1)*(2*Pi*v1) * cos(2*Pi*v1*t + angleIni_B+Pi)
    val accCy = -rC * (2*Pi*v1)*(2*Pi*v1) * sin(2*Pi*v1*t + angleIni_B+Pi)

    // D,E
    val accDx = -rH * (2*Pi*v1)*(2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B)  - rD* (2*Pi*(v2+v1))*(2*Pi*(v2+v1))* cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
    val accDy = -rH * (2*Pi*v1)*(2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B)  - rD* (2*Pi*(v2+v1))*(2*Pi*(v2+v1))* sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B))
    val accEx = -rH * (2*Pi*v1)*(2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B)  - rE* (2*Pi*(v2+v1))*(2*Pi*(v2+v1))* cos(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +Pi)
    val accEy = -rH * (2*Pi*v1)*(2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B)  - rE* (2*Pi*(v2+v1))*(2*Pi*(v2+v1))* sin(2*Pi*(v2+v1)*t + (angleIni_D+ angleIni_B) +Pi)

    // F,G
    val accFx = -rI * (2*Pi*v1)*(2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B+Pi)  - rF* (2*Pi*(v3+v1))*(2*Pi*(v3+v1))* cos(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi))
    val accFy = -rI * (2*Pi*v1)*(2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B+Pi)  - rF* (2*Pi*(v3+v1))*(2*Pi*(v3+v1))* sin(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi))
    val accGx = -rI * (2*Pi*v1)*(2*Pi*v1)* cos(2*Pi*v1*t + angleIni_B+Pi)  - rG* (2*Pi*(v3+v1))*(2*Pi*(v3+v1))* cos(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi) +Pi)
    val accGy = -rI * (2*Pi*v1)*(2*Pi*v1)* sin(2*Pi*v1*t + angleIni_B+Pi)  - rG* (2*Pi*(v3+v1))*(2*Pi*(v3+v1))* sin(2*Pi*(v3+v1)*t + (angleIni_F+ angleIni_B+Pi) +Pi)

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



  // compute a trajctory in stationary mode
  def dynamicTrajectoryStationnary(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD ,
                                   rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                                   rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI,
                                   angleIni_B:Double = DefaultValuesParameterModel.angleIni_B, angleIni_D:Double = DefaultValuesParameterModel.angleIni_D,
                                   angleIni_F:Double = DefaultValuesParameterModel.angleIni_F)(v1:Double=DefaultValuesParameterModel.v1, v2:Double=DefaultValuesParameterModel.v2, v3:Double=DefaultValuesParameterModel.v3)(T:Double, deltaT:Double):Vector[DynamicalCurrentStateStationary] = {
    val N = floor(T / deltaT).toInt
    val vecTimes = (0 until N).map(x => x * deltaT) :+ T
    val res = vecTimes.map(x => dynamicStationnary(rB,rC,rD,rE,rF,rG,rH,rI,angleIni_B,angleIni_D,angleIni_F)(v1,v2,v3)(x)).toVector
    // type de res: Vector[DynamicalCurrentStateStationary]
    res
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

  /*
  def print(f:Double=>Double, a:Double, b:Double, steps:Int)={
    println("rectangular left   : %f".format(integrate(f, a, b, steps, leftRect)))
    println("rectangular middle : %f".format(integrate(f, a, b, steps, midRect)))
    println("rectangular right  : %f".format(integrate(f, a, b, steps, rightRect)))
    println("trapezoid          : %f".format(integrate(f, a, b, steps, trapezoid)))
    println("simpson            : %f".format(integrate(f, a, b, steps, simpson)))
  }

  print(fn1, 0, 1, 100)
  println("------")
  print(fn2, 1, 100, 1000)
  println("------")
  print(fn3, 0, 5000, 5000000)
  println("------")
  print(fn3, 0, 6000, 6000000)
  */


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


  // Persistence rétinienne
  val persistenceTime = 0.05
  val nb_point_persistence = 50

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

}

object FixedParameterModel {

  // rayon
  val rB = 4.0
  val rC = 4.0
  val rD = 1.0
  val rE = 1.0
  val rF = 1.0
  val rG = 1.0
  val rH = 2
  val rI = 2

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
    // ctrl b pour avoir le code


    ////////////////////////////////
    //    POINT SINGULIER
    ////////////////////////////////


    def distance(a: Double, b: Double)=  {sqrt( a*a + b*b)}

    def slaveFindSingularPoints(vec1:Vector[Double],vec2:Vector[Double], seuil:Double)={
      vec1.zip(vec2).zipWithIndex.collect{ case((a,b),c) if distance(a,b)<seuil => (a,b,distance(a,b),c)}
      // le resultat c'est un vecteur de triplets (Doulble,Double,Double,Int)
    }


    def countSingularPoints(vec1:Vector[Double],vec2:Vector[Double], seuil:Double) = {
      // return an integer
      val temp = slaveFindSingularPoints(vec1,vec2,seuil)
      SelectSingularPoints(temp).length

    }


    def timesOfSingularPoints(vec1:Vector[Double],vec2:Vector[Double], seuil:Double) = {
      // return a vector of integer (indices)
      val temp = slaveFindSingularPoints(vec1,vec2,seuil)
      SelectSingularPoints(temp).map(_._4)
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






    ////////////////////////////////
    //    DENSITE
    ////////////////////////////////


    def maxSquareForDensity(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD ,
                            rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                            rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI)= {
      max(max(rB, max(rH+rD, rH+rE)) , max(rC, max(rI+rF, rI+rG)) )
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

      var array = Array.ofDim[Int](N,N)  // initialisé à 0

      for (point <- temp) {

        val (i, j) = findCaseForPoint((point._1, point._2), xmin, xmax, ymin, ymax, N)
        array(i)(j) += 1
      }

      //println( array.map(x => x.sum).sum ) // on retrouve le nombre de points, ouf!
      //println(convertResultStationnary(res).time.length)
      array
      // to print the array array.map(_.mkString(" ")).mkString("\n")
    }


    def pointsDensitySquare(vec1:Vector[Double],vec2:Vector[Double],xmin:Double,xmax:Double,ymin:Double,ymax:Double,N:Int) = {

      // nombre de cases non vides dans l'array
      val array = arrayDensitySquare(vec1,vec2,xmin,xmax,ymin,ymax,N)
      val nbCasesNonEmpty = array.map(x => x.filter(_ != 0).length).sum

      // la ratio d'occupation:
      (nbCasesNonEmpty.toFloat / (N * N)).toDouble
    }


    def sum2Arrays(a:Array[Int],b:Array[Int])={
      a.zip(b).map(x=>x._1+x._2)
    }

    def sum2ArraysOfArrays(a:Array[Array[Int]],b:Array[Array[Int]])={
      a.zip(b).map(x => sum2Arrays(x._1,x._2))
    }

    def nextSumVectorOfArraysOfArrays(v:Vector[Array[Array[Int]]],acc:Array[Array[Int]]): Array[Array[Int]]= {
      if (v.isEmpty){acc} else{
        val temp = v(0)
        val newV = v.tail
        val newAcc = sum2ArraysOfArrays(acc,temp)
        nextSumVectorOfArraysOfArrays(newV,newAcc)
      }

    }

    def sumVectorOfArraysOfArrays(v:Vector[Array[Array[Int]]]) = {
      val acc = v(0)
      nextSumVectorOfArraysOfArrays(v.tail,acc)
    }

/*
    // to test the sums
    val a = Array(1.0,2.0)
    val b = Array(3.0,4.0)
    //println(sum2Arrays(a,b).mkString(""))

    val c = Array(a,a)
    val d = Array(b,a)
    //println(sum2ArraysOfArrays(c,d).map(_.mkString(" ")).mkString("\n") )

    //println(sumVectorOfArraysOfArrays(Vector(c,d)).map(_.mkString(" ")).mkString("\n"))
*/

    def allTrajectoriesDensitySquare(res:TrajectoryStationnary,xmin:Double,xmax:Double,ymin:Double,ymax:Double,N:Int) = {

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
    //   BOUCLE ou presque (TEMPS DE PREMIER RETOUR)
    ////////////////////////////////



    def distanceLoopPoint( x: (Double,Double), y: (Double,Double) )=  {sqrt( (x._1 - y._1)*(x._1 - y._1) + (x._2 - y._2)*(x._2 - y._2)  )}

    // on garde l'indice du point de référence
    def distanceLoopPoint2( x: (Double,Double,Int), y: (Double,Double,Int) )=  {(sqrt( (x._1 - y._1)*(x._1 - y._1) + (x._2 - y._2)*(x._2 - y._2) ) , x._3, y._3-x._3)}

    def comparePointsToAllElemnts(point:(Double,Double), vect: Vector[(Double,Double)])={
      vect.map(x => distanceLoopPoint(x,point))
    }

    def comparePointsToAllElemnts2(point:(Double,Double,Int), vect: Vector[(Double,Double,Int)])={
      vect.map(x => distanceLoopPoint2(x,point))
    }




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


    def tempsPremierRetour(v1:Vector[Double], v2:Vector[Double], seuil:Double) : Vector[Int] = {
      val resPremierRetour = premierRetour(v1,v2,seuil)
      //val res4 = resPremierRetour .filter(_ != () )  // retire les vecteur vide
      //val res4 = resPremierRetour.filterNot(_ == ())
      val res4 = resPremierRetour.filterNot(_._2 ==0 )
      res4.map{ case( ((a,b,c),d) )  => d}

    }






    ////////////////////////////////
    //    COURBURE de la l'arc paramétré
    ////////////////////////////////

    def mean(v: Vector[Double])={
      (v.sum / v.length)
    }


    // formule: x' y'' - y' x'' / ( x'^2 + y' ^2 )^(3/2)
    def slaveCourbure(sx:Double, sy:Double, acx:Double, acy:Double)={
      ((sx*acy) - (sy*acx)) / pow( sx*sx + sy*sy , 3/2)
    }

    def courbure(speedX:Vector[Double], speedY:Vector[Double], accX:Vector[Double], accY:Vector[Double])={

      val temp = speedX.zip(speedY).zip(accX).zip(accY).map( x  => (x._1._1._1,x._1._1._2,x._1._2,x._2))  // pour avoir autre format
      temp.map(x=>  slaveCourbure(x._1,x._2,x._3,x._4))
    }




    // a partir du zip des courbures et index selectionnés par seuil (ici on n'utilise que les index)
    // on repère les point consécutifs qui sont dans ce vecteur
    // on garde un vector de (Int,Int) le premier Int et le premier indice de la suite (valeurs d'indice consécutives)
    // le second est la longeur de cette suite
    def nextSlaveDetectStraitLine(v:Vector[(Double,Int)],acc:Vector[(Int,Int)],temp:(Double,Int)):Vector[(Int,Int)]={
      if (v.isEmpty){ acc
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

      val newV = temp.tail
      val newTemp = temp(0)
      val newAcc = Vector((temp(0)._2,1))
      nextSlaveDetectStraitLine(newV,newAcc,newTemp)

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



  }






///////////////////  old Loop points

/*
def slaveFindLoopPoints(vec1:Vector[Double],vec2:Vector[Double], seuil:Double)={

  val temp = vec1.zip(vec2)

  // vector de vector avec toutes les distances
  val res = temp.map(y => comparePointsToAllElemnts(y,temp))

  // on ne garde que les distances < seuil
  //val res2 = res.map(_.filter( _< seuil))
  val res2 = res.map(_.zipWithIndex.collect{ case((a,b)) if a<seuil => (a,b)})

  res2
  // pour chaque vecteur, au moins une distance nulle (avec sois même), on la retire (ou condition sur la différence entre les temps?)

  //  temp.collect{ case((a,b),c) if distance(a,b)<seuil => (a,b,c)}
  // le resultat c'est un vecteur de triplets (Doulble,Double,Int)

  // pour voir ce que ca donne
  //val res3 = res2(0)
  //res3
}
*/


/*
  def timeLoopPoints(res:Vector[DynamicalCurrentStateStationary], seuil:Double) = {

    val res2 = convertResultStationnary(res)
    // Loop for B

    val res3 = slaveFindLoopPoints(res2.Bx,res2.By, 1.0) //(0)  // first point of B
    //val indices = res3.map(_._2)




    // Loop for F
    //val res3 = slaveFindLoopPoints(res2.Fx,res2.Fy, seuil) //(0)  // 0: first point (is a loop?)
    // res3  //  de type Vector[(Doublle,Doublle)]

    val indices = res3.map(_._2)


    val a = indices.tail
    // tous sauf le dernier
    val b = indices.dropRight(1)
    //println(b)
    val res4 = a.zip(b).collect{ case(x,y) if (x-y)>1  => x }
    //res4

    // pour chaque positions de la trajctoire d'un point, on compte le nombre de loop qu'il y a
    // si la figure est périodique, il y en a au moins 1





    def trieLoopPoint(res:Vector[(Double,Int)])={
      val indices = res.map(_._2)
      val a = indices.tail
      // tous sauf le dernier
      val b = indices.dropRight(1)
      Vector(a(0)) ++ a.zip(b).collect{ case(x,y) if (x-y)>1  => x }
      // on ne garde que les indices qui ne sont pas voisins (continuité) ie 1 indice par classe d'indice voisins (le premier)
      // au moins un élément car on est à distance 0 de nous même
    }

    //trieLoopPoint(res3(0))
    res3.map(x => trieLoopPoint(x))
*/






/*
val indices = res3.map(_._2)
print(indices)

// tous sauf le premier
val a = indices.tail
//println(a)
// tous sauf le dernier
val b = indices.dropRight(1)
//println(b)
val seuiTimeLoop : Int = 5
//val res4 = a.zip(b).map( x => x._1 -x._2).filter(_ > 1 )
val res4 = a.zip(b).collect{ case(x,y) if (x-y)>1  => x } // .filter(_ > seuiTimeLoop )
// collect{ case((a,b),c) if distance(a,b)<seuil => (a,b,c)}
println(res4)

}
*/
