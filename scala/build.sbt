import sbt._

name := "scala"

version := "0.1"

scalaVersion := "2.12.8"

//libraryDependencies +=


enablePlugins(SbtOsgi)

libraryDependencies += "org.apache.commons" % "commons-math3" % "3.6.1"

OsgiKeys.exportPackage := Seq("electra.*")

OsgiKeys.importPackage := Seq("*;resolution:=optional")

//OsgiKeys.privatePackage := Seq("*")
OsgiKeys.privatePackage := Seq("!scala.*,!java.*,!monocle.*,!META-INF.*.RSA,!META-INF.*.SF,!META-INF.*.DSA,META-INF.services.*,META-INF.*,*")

OsgiKeys.requireCapability := """osgi.ee;filter:="(&(osgi.ee=JavaSE)(version=1.8))""""

resolvers += "osgeo" at "http://download.osgeo.org/webdav/geotools"
libraryDependencies += "org.openmole.library" %% "spatialdata" % "0.2"



/*
lazy val bundle = Project("electra", file("electra")) enablePlugins(SbtOsgi) settings(
  scalaVersion := "2.12.8",
  OsgiKeys.exportPackage := Seq("electra.*;-split-package:=merge-first"),
  OsgiKeys.importPackage := Seq("*;resolution:=optional"),
  OsgiKeys.privatePackage := Seq("*"),
  OsgiKeys.requireCapability := """osgi.ee;filter:="(&(osgi.ee=JavaSE)(version=1.8))"""
)

enablePlugins (SbtOsgi)

*/