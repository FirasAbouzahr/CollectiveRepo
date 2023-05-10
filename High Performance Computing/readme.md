# High Performance Computing Resources

This directory serves as a hub for resources I have developed for using TACC's (Texas Advanced Computing Center) super computers 
to run computational expensive tasks. For instance, detector simulations in Geant4 can require an immense amount of computational power, 
especially when simulating complex geometries with many events.

# simulationSerialJob.py
- this script is used as a means to generate 2 files needed to Geant4 simulations in parallel on TACC's Stampede2 (or any other for that matter)
- supercomputer. The first file is the launcher file which is a script that contains a line for each executable that you wish to run simultaneously.
