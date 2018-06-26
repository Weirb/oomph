# Fourier decomposed Helmholtz equation with oomph-lib

This repository contains the code for my dissertation on solving the Fourier decomposed Helmholtz equation with oomph-lib.

## Building

With oomph-lib installed on your system in the directory $OOMPH_HOME, clone the repository into the directory $OOMPH_HOME/user_drivers/billy.
Then, return to $OOMPH_HOME and run
   ./autogen.sh 
and wait for the command to complete.
After this, enter first the billy/src directory and run
   make
then enter billy/driver and run
   make

The program can then be run by typing ./main

## To do list

0. Fix the Newton solver iteration error
0. Remove the PML/fix the issues
0. Add an exact solution for the problem
0. Cleanup/comment the code

Alt-list
0. oscillating_sphere.cc
0. Copy the pml_helmholtz_equation code and build in directory, don't create lib
0. Copy across the Fourier decomposed stuff
