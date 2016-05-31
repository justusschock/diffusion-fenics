This repository is for solving Partial Differential Equations (PDEs) in a numerical way.

# Dependencies #
This repository is based on [FEniCS-Project (1.6.0 or above)](http://fenicsproject.org/)
Essentially needed parts of the FEniCs-Project are dolfin, UFL, FFC and UFC.

For the build-process [CMake (2.8.15or above)](https://cmake.org/) and [GNU make (4.2 or above)](http://www.gnu.org/software/make/) are recommed.

#Set Up#
To setup the project make sure, all dependencies are installed and working well.
Clone this repository into the folder, you wish to work at and open a terminal in directory's src-folder.

In this terminal with 
```
#!terminal
ffc -l dolfin *.ufl
```
you can update the c++-headerfiles by compiling the .ufl-files.

After doing this, change the terminal's working directory to the directory containing this repository and run 
```
#!terminal
cmake ./
```
to generate the makefiles for your system.

Compiling this application can be done by running 
```
#!terminal
make
```
in the home-directory of this project.

After compiling, the application can be started by running
```
#!terminal
./diffusion-fenics
```
in the build-directory.
To start this application with specified arguments please take a lokk at the [FEniCs Manual](https://launchpadlibrarian.net/84116499/fenics-manual-2011-10-31.pdf)