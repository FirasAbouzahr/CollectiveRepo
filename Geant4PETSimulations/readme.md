# PET Simulations Using Geant4 Monte Carlo Simulations

This directory includes all data files needed to run our groups canonical PET simulations in Geant4. This is a collective work from many group members in Dr. Lang's high-energy physics group. My main contributions are via novel geometries, materials, simulated sources, and some other logistical source codes as well. Here's how to compile and run this code: 

* Start by making a directory, say directory PETSimulation, and putting this entire directory inside of it.
```
$:PETSimulation  pwd 
/path/to/PETSimulation

```
Geant4PETSimulations should be located inside of it now. 
```
$:PETSimulation  ls 
Geant4PETSimulations

```
Now create a build directory in the same main directory PETSimulation and move into it. 

```
$:PETSimulation mkdir build && cd build 

```
Now, we need to source Geant4 installation. See my other tutorial on this. 

```
$:PETSimulation source /path/to/geant4_install/bin/geant4.sh

```
Now, we compile and make it
