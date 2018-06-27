# Fourier decomposed Helmholtz equation with oomph-lib

This repository contains the code for my dissertation on solving the Fourier decomposed Helmholtz equation with oomph-lib.

## Building

Clone the directory into $OOMPH_HOME/user_drivers and run $OOMPH_HOME/autogen.sh.
After this has finished running, enter the cloned directory and run `make`.
Finally, the program can then be run by typing `./main`.

## To do list

0. Figure out how to use multigrid to solve my mesh
0. Add an exact solution for the problem
0. Cleanup/comment the code
