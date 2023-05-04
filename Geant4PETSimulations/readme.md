# PET Simulations Using Geant4

This directory includes all C++ codes needed to run our group's canonical Geant4-based PET simulations. This is a collective work from many group members in Dr. Lang's high-energy physics group. My main contributions are via novel geometries, materials, simulated sources, and some other contributions to more logistical source codes.

Inside of Geant4PETSimulations/source/geometry you'll want to choose a single geometry script and move into the source directory with the other .cc scripts. You will also only want to choose one of the PETmain.cc files to keep. PETmainMPI.cc is only for running in parallel with MPI, it also requires some additonal changes to the other files that I will cover later. I personally use this only when running large geometries a full body PET scanner which can contain 5000+ individual crystals. 

Note the following requirements to actually run this simulation:
* Geant4 10.7.4 or older (in versions greater than 11, some pieces of this code of have been depreciated)
* cmake
* xquartz or openGL (check Geant4 documentation for other alternatives too)
* for mac: XCode and installation of XCode Command Line Tools 

Now, here's how to compile and run this code: 

* Start by saving this directory while keeping the above comments in mind (like choosing only one geometry .cc script to move to Geant4PETSimulations/source/geometry) and move into into another diretory, say we call it PETSimulation. 

```
$:PETSimulation pwd 
/path/to/PETSimulation
```
* Geant4PETSimulations should be located inside of it now. 
```
$:PETSimulation ls 
Geant4PETSimulations
```
* Now create a build directory in the same main directory PETSimulation and move there. 

```
$:PETSimulation mkdir build && cd build 
```
* We need to source Geant4 now.

```
$:build source /path/to/geant4_install/bin/geant4.sh
```
* Compile with cmake and then run make:

```
$:PETSimulation cmake GEANT4_DIR=/path/to/geant4_install/lib/Geant4-10.7.4 /path/to/Geant4PETSimulations #and now you should see an output like:
CMake Deprecation Warning at CMakeLists.txt:3 (cmake_minimum_required):
  Compatibility with CMake < 2.8.12 will be removed from a future version of
  CMake.

  Update the VERSION argument <min> value or use a ...<max> suffix to tell
  CMake that the project does not need compatibility with older versions.


-- The C compiler identification is AppleClang 14.0.0.14000029
-- The CXX compiler identification is AppleClang 14.0.0.14000029
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- Detecting CXX compiler ABI info
-- Detecting CXX compiler ABI info - done
-- Check for working CXX compiler: /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/c++ - skipped
-- Detecting CXX compile features
-- Detecting CXX compile features - done
-- Found EXPAT: /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX13.1.sdk/usr/lib/libexpat.tbd (found suitable version "2.4.8", minimum required is "2.4.8") 
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD
-- Performing Test CMAKE_HAVE_LIBC_PTHREAD - Success
-- Found Threads: TRUE  
-- Found X11: /usr/X11R6/include   
-- Looking for XOpenDisplay in /usr/X11R6/lib/libX11.dylib;/usr/X11R6/lib/libXext.dylib
-- Looking for XOpenDisplay in /usr/X11R6/lib/libX11.dylib;/usr/X11R6/lib/libXext.dylib - found
-- Looking for gethostbyname
-- Looking for gethostbyname - found
-- Looking for connect
-- Looking for connect - found
-- Looking for remove
-- Looking for remove - found
-- Looking for shmat
-- Looking for shmat - found
-- Looking for IceConnectionNumber in ICE
-- Looking for IceConnectionNumber in ICE - found
-- Found XQuartzGL: /usr/X11R6/include  
-- Configuring done
-- Generating done
-- Build files have been written to: /Users/feef/Desktop/PETSimulation/build

$:build make -j4 
```
* I'll skip the long output that follows from running make. Finally, after all of this, run the command ./PETSim while still in the same build directory and your simulation will start!
