package electra

import com.sun.org.apache.xerces.internal.impl.xpath.XPath.Step
//import electra.Model.integrate
import electra.NumericalIntegration._
import electra.Model._
import scala.math._

object Electra extends App {

  /*
  val a=0
  val b=1
  val steps=10000
  def f(t:Double) = { t*t }
  val res = integrate(f,a,b,steps,simpson)
  */


  /*
  val res = dynamicStationnary()(1.0,DefaultValuesParameterModel.v1, DefaultValuesParameterModel.v2, DefaultValuesParameterModel.v3)
  println(res)


  val t= 1.0
  val v = vitesse_rotation_moteur(DefaultValuesParameterModel.u1,1.0,1.0)(t)
  println(v)


  val angle1_t = calcul_angle(DefaultValuesParameterModel.u1,1.0,1.0)(t)
  println(angle1_t)


  val res2 =  dynamicTransitoire()(t,DefaultValuesParameterModel.u1, DefaultValuesParameterModel.u2, DefaultValuesParameterModel.u3)
  println(res2)
  */

  // avec next
  val fixedParametersModel = new FixedParametersModel()
  val parametersTransitoire = new ParametersTransitoire()
  val Nmax = 100
  val deltaT = 0.01
  val initialConditions = new InitialConditions()
  val initialState = createInitialState(fixedParametersModel)(initialConditions)
  val res = simu(fixedParametersModel)(parametersTransitoire)(Nmax,deltaT,initialState)

  println(res)
  println( sqrt((res.Dx)*(res.Dx) + (res.Dy)*(res.Dy)))


  // sans next
  val t = Nmax*deltaT
  val res2 =  dynamicTransitoire()(t,DefaultValuesParameterModel.u1, DefaultValuesParameterModel.u2, DefaultValuesParameterModel.u3)
  println(res2)

  //
  /*
  println(res.Dx)
  println(res2.Dx)
  println(res.Dy)
  println(res2.Dy)
  println(sqrt((res.Dx)*(res.Dx) + (res.Dy)*(res.Dy)))
  println(sqrt((res2.Dx)*(res2.Dx) + (res2.Dy)*(res2.Dy)))
  */

/*
  println(res.Hx)
  println(res2.Hx)
  println(res.Hy)
  println(res2.Hy)
  println(res.HDx)
  println(res2.HDx)
  println(res.HDy)
  println(res2.HDy)
*/

}



case class DynmicalCurrentState(time : Double, Bx: Double, By:Double, Cx:Double=0.0, Cy:Double=0.0, Dx:Double=0.0, Dy:Double=0.0,
                                Ex:Double=0.0, Ey:Double=0.0, Fx:Double=0.0, Fy:Double=0.0, Gx:Double=0.0, Gy:Double=0.0, Hx:Double=0.0,Hy:Double=0.0,
                                Ix:Double=0.0, Iy:Double=0.0,
                                HDx:Double=0.0, HDy:Double=0.0, HEx:Double=0.0, HEy:Double=0.0,
                                IFx:Double=0.0, IFy:Double=0.0, IGx:Double=0.0, IGy:Double=0.0,
                                angle1:Double=0.0, angle2:Double=0.0,angle3:Double=0.0,
                                vitesseMoteur1:Double=0.0, vitesseMoteur2:Double=0.0, vitesseMoteur3:Double=0.0
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
                         angleIni_F:Double = DefaultValuesParameterModel.angleIni_F)(t:Double,v1:Double, v2:Double, v3:Double) = {

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

    (Bx,By,Cx,Cy,Dx,Dy,Ex,Ey,Fx,Fy,Gx,Gy,Hx,Hy,Ix,Iy)
  }


  ////////////////////////////////
  //    TRANSITOIRE
  ////////////////////////////////


  def vitesse_rotation_moteur(u: Double =>Double ,A:Double,C:Double)(t:Double)= {
    def temp(t:Double):Double ={
      A / C * u(t) * exp(t/C)
    }

    val vitesse = exp(-t/C ) * NumericalIntegration.integrate(temp,0.0,t,(ceil(abs(t))+1).toInt*HiddenParameters.steps,simpson)
    vitesse
  }


  def calcul_angle(u: Double =>Double ,A:Double,C:Double)(t:Double) ={
    def vit(t:Double):Double = {Model.vitesse_rotation_moteur(u,A,C)(t)}
    //val angle = NumericalIntegration.integrate((x:Double) => x*x,0.0,t,ceil(t).toInt*HiddenParameters.steps,simpson)
    val angle = NumericalIntegration.integrate(vit,0.0,t,(ceil(abs(t))+1).toInt*HiddenParameters.steps,simpson)
    //vitesse(2.0)
    angle
  }



  def dynamicTransitoire(rB: Double = FixedParameterModel.rB, rC: Double = FixedParameterModel.rC, rD:Double = FixedParameterModel.rD ,
                         rE:Double = FixedParameterModel.rE, rF:Double = FixedParameterModel.rF, rG:Double = FixedParameterModel.rG,
                         rH:Double = FixedParameterModel.rH, rI:Double = FixedParameterModel.rI,
                         angleIni_B:Double = DefaultValuesParameterModel.angleIni_B, angleIni_D:Double = DefaultValuesParameterModel.angleIni_D,
                         angleIni_F:Double = DefaultValuesParameterModel.angleIni_F,
                         A1:Double = FixedParameterModel.A1,C1:Double=FixedParameterModel.C1,
                         A2:Double = FixedParameterModel.A2,C2:Double=FixedParameterModel.C2,
                         A3:Double = FixedParameterModel.A3,C3:Double=FixedParameterModel.C3)
                        (t:Double, u1:Double=>Double, u2:Double=>Double, u3:Double=>Double) = {

    // angle
    val angle1_t = calcul_angle(u1,A1,C1)(t)
    val angle2_t = calcul_angle(u2,A2,C2)(t)
    val angle3_t = calcul_angle(u3,A3,C3)(t)

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

    new DynmicalCurrentState(time = t,Bx = Bx, By= By, Cx = Cx, Cy= Cy,
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


  def nextStepDynamic(fixedParametersModel: FixedParametersModel)(parametersTransitoire: ParametersTransitoire)
                     (dynmicalCurrentState: DynmicalCurrentState, deltaT:Double)={


    val t = dynmicalCurrentState.time

    // vitesse moteur
    // 1
    def tempVitesseMoteur1 (t:Double):Double ={
      fixedParametersModel.A1 / fixedParametersModel.C1 * parametersTransitoire.u1(t) * exp(t/fixedParametersModel.C1)
    }

    //val newVitesseMoteur1 = dynmicalCurrentState.vitesseMoteur1 + NumericalIntegration.integrate(tempVitesseMoteur1,max((t-deltaT),0),t,HiddenParameters.steps,simpson)

    // 2
    def tempVitesseMoteur2 (t:Double):Double ={
      fixedParametersModel.A2 / fixedParametersModel.C2 * parametersTransitoire.u2(t) * exp(t/fixedParametersModel.C2)
    }

    //val newVitesseMoteur2 = dynmicalCurrentState.vitesseMoteur2 + NumericalIntegration.integrate(tempVitesseMoteur2,max((t-deltaT),0),t,HiddenParameters.steps,simpson)

    // 3
    def tempVitesseMoteur3 (t:Double):Double ={
      fixedParametersModel.A3 / fixedParametersModel.C3 * parametersTransitoire.u3(t) * exp(t/fixedParametersModel.C3)
    }

    //val newVitesseMoteur3 = dynmicalCurrentState.vitesseMoteur3 + NumericalIntegration.integrate(tempVitesseMoteur3,max((t-deltaT),0),t,HiddenParameters.steps,simpson)


    // angle
    // angle 1
    def tempAngle1_1(tt: Double) = { exp(-tt/FixedParameterModel.C1)* NumericalIntegration.integrate(tempVitesseMoteur1,t,tt,HiddenParameters.steps,simpson)}
    //val tempAngle1_2 = fixedParametersModel.C1 * (exp(-t/fixedParametersModel.C1)-exp(-(t+deltaT)/fixedParametersModel.C1)) * newVitesseMoteur1 + NumericalIntegration.integrate(tempAngle1_1,t,t+deltaT,HiddenParameters.steps,simpson)
    val tempAngle1_2 = fixedParametersModel.C1 * (exp(-t/fixedParametersModel.C1)-exp(-(t+deltaT)/fixedParametersModel.C1)) * NumericalIntegration.integrate(tempVitesseMoteur1,0,t,HiddenParameters.steps,simpson)  + NumericalIntegration.integrate(tempAngle1_1,t,t+deltaT,HiddenParameters.steps,simpson)

    val newAngle1 = dynmicalCurrentState.angle1 + tempAngle1_2


    // angle 2
    def tempAngle2_1(tt: Double) = { exp(-tt/FixedParameterModel.C2)* NumericalIntegration.integrate(tempVitesseMoteur2,t,tt,HiddenParameters.steps,simpson)}
    val tempAngle2_2 = fixedParametersModel.C2 * (exp(-t/fixedParametersModel.C2)-exp(-(t+deltaT)/fixedParametersModel.C2)) * NumericalIntegration.integrate(tempVitesseMoteur2,0,t,HiddenParameters.steps,simpson)  + NumericalIntegration.integrate(tempAngle2_1,t,t+deltaT,HiddenParameters.steps,simpson)
    val newAngle2 = dynmicalCurrentState.angle2 + tempAngle2_2


    // angle 3
    def tempAngle3_1(tt: Double) = { exp(-tt/FixedParameterModel.C3)* NumericalIntegration.integrate(tempVitesseMoteur3,t,tt,HiddenParameters.steps,simpson)}
    val tempAngle3_2 = fixedParametersModel.C3 * (exp(-t/fixedParametersModel.C3)-exp(-(t+deltaT)/fixedParametersModel.C3)) * NumericalIntegration.integrate(tempVitesseMoteur3,0,t,HiddenParameters.steps,simpson)  + NumericalIntegration.integrate(tempAngle3_1,t,t+deltaT,HiddenParameters.steps,simpson)
    val newAngle3 = dynmicalCurrentState.angle3 + tempAngle3_2

    // Positions
    // B
    val newBx = dynmicalCurrentState.Bx * cos(2*Pi*tempAngle1_2) - dynmicalCurrentState.By * sin(2*Pi*tempAngle1_2)
    val newBy = dynmicalCurrentState.By * cos(2*Pi*tempAngle1_2) + dynmicalCurrentState.Bx * sin(2*Pi*tempAngle1_2)

    // C
    val newCx = dynmicalCurrentState.Cx * cos(2*Pi*tempAngle1_2) - dynmicalCurrentState.Cy * sin(2*Pi*tempAngle1_2)
    val newCy = dynmicalCurrentState.Cy * cos(2*Pi*tempAngle1_2) + dynmicalCurrentState.Cx * sin(2*Pi*tempAngle1_2)

    // H,I
    val newHx = dynmicalCurrentState.Hx * cos(2*Pi*tempAngle1_2) - dynmicalCurrentState.Hy * sin(2*Pi*tempAngle1_2)
    val newHy = dynmicalCurrentState.Hy * cos(2*Pi*tempAngle1_2) + dynmicalCurrentState.Hx * sin(2*Pi*tempAngle1_2)
    val newIx = dynmicalCurrentState.Ix * cos(2*Pi*tempAngle1_2) - dynmicalCurrentState.Iy * sin(2*Pi*tempAngle1_2)
    val newIy = dynmicalCurrentState.Iy * cos(2*Pi*tempAngle1_2) + dynmicalCurrentState.Ix * sin(2*Pi*tempAngle1_2)

    // D,E
    val newHDx = dynmicalCurrentState.HDx * cos(2*Pi*(tempAngle1_2+tempAngle2_2)) - dynmicalCurrentState.HDy * sin(2*Pi*(tempAngle1_2+tempAngle2_2))
    val newHDy = dynmicalCurrentState.HDy * cos(2*Pi*(tempAngle1_2+tempAngle2_2)) + dynmicalCurrentState.HDx * sin(2*Pi*(tempAngle1_2+tempAngle2_2))
    val newHEx = dynmicalCurrentState.HEx * cos(2*Pi*(tempAngle1_2+tempAngle2_2)) - dynmicalCurrentState.HEy * sin(2*Pi*(tempAngle1_2+tempAngle2_2))
    val newHEy = dynmicalCurrentState.HEy * cos(2*Pi*(tempAngle1_2+tempAngle2_2)) + dynmicalCurrentState.HEx * sin(2*Pi*(tempAngle1_2+tempAngle2_2))

    val newDx = newHx + newHDx
    val newDy = newHy + newHDy
    val newEx = newHx + newHEx
    val newEy = newHy + newHEy

    // F,G
    val newIFx = dynmicalCurrentState.IFx * cos(2*Pi*(tempAngle1_2+tempAngle3_2)) - dynmicalCurrentState.IFy * sin(2*Pi*(tempAngle1_2+tempAngle3_2))
    val newIFy = dynmicalCurrentState.IFy * cos(2*Pi*(tempAngle1_2+tempAngle3_2)) + dynmicalCurrentState.IFx * sin(2*Pi*(tempAngle1_2+tempAngle3_2))
    val newIGx = dynmicalCurrentState.IGx * cos(2*Pi*(tempAngle1_2+tempAngle3_2)) - dynmicalCurrentState.IGy * sin(2*Pi*(tempAngle1_2+tempAngle3_2))
    val newIGy = dynmicalCurrentState.IGy * cos(2*Pi*(tempAngle1_2+tempAngle3_2)) + dynmicalCurrentState.IGx * sin(2*Pi*(tempAngle1_2+tempAngle3_2))

    val newFx = newIx + newIFx
    val newFy = newIy + newIFy
    val newGx = newIx + newIGx
    val newGy = newIy + newIGy


    new DynmicalCurrentState(time = t+deltaT,Bx = newBx, By=newBy, Cx = newCx, Cy=newCy,
      Hx = newHx, Hy = newHy, Ix = newIx, Iy = newIy,
      HDx = newHDx, HDy = newHDy, HEx = newHEx, HEy = newHEy,
      Dx = newDx, Dy = newDy, Ex = newEx, Ey = newEy,
      IFx = newIFx, IFy = newIFy, IGx = newIGx, IGy = newIGy,
      Fx = newFx, Fy = newFy, Gx = newGx, Gy = newGy
      )


  }



  def slaveSimu(fixedParametersModel: FixedParametersModel)(parametersTransitoire: ParametersTransitoire)
  (iter:Int,Nmax:Int,deltaT:Double, res:DynmicalCurrentState):DynmicalCurrentState = {

    if (iter == Nmax) res
    else slaveSimu(fixedParametersModel)(parametersTransitoire)(iter + 1, Nmax, deltaT, nextStepDynamic(fixedParametersModel)(parametersTransitoire)(res,deltaT))
  }


  def createInitialState(fixedParametersModel: FixedParametersModel)(initialConditions: InitialConditions):DynmicalCurrentState={
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

    new DynmicalCurrentState(time =0, Bx=Bx, By=By, Cx=Cx, Cy=Cy, Dx=Dx, Dy=Dy,
                                    Ex=Ex, Ey=Ey, Fx=Fx, Fy=Fy, Gx=Gx, Gy=Gy, Hx=Hx,Hy=Hy,
                                    Ix=Ix, Iy=Iy,
                                    HDx=HDx, HDy=HDy, HEx=HEx, HEy=HEy,
                                    IFx=IFx,IFy=IFy,IGx=IGx,IGy=IGy,
                                    angle1=initialConditions.angleIni_B, angle2= initialConditions.angleIni_D,angle3=initialConditions.angleIni_F,
                                    vitesseMoteur1=0.0, vitesseMoteur2=0.0, vitesseMoteur3=0.0)
  }



  def simu(fixedParametersModel: FixedParametersModel)(parametersTransitoire: ParametersTransitoire)
          (Nmax:Int,deltaT:Double, initialState:DynmicalCurrentState):DynmicalCurrentState = {

    slaveSimu(fixedParametersModel)(parametersTransitoire)(0, Nmax, deltaT, initialState)

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
 val steps = 1000 // for numerical integration
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

  // vitesse (transitoire)
  def u1(double: Double)={2.0}
  def u2(double: Double)={1.0}
  def u3(double: Double)={3.0}

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


  }