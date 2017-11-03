# advanced-cmp-midsem

This branch implements the averaged density Hartree-Fock method to find SCF corrections to a system with multiple gaussian wells at fixed separation.

# Building and Running

The potential, its parameters and other constants are placed in [configuration.hpp](configuration.hpp) file.  Any modification of fixed paramters can be performed within the [configuration.hpp](configuration.hpp) file.

To build the program compile using "**g++ -std=c++14 -o HF hartree-fock.cpp configuration.hpp -llapack -llapacke**".

The program takes **2 arguments** with **main()**. First argument specifies 100 x well-separation, and the second argument is the number of wells. Therefore to execute the program with _N_  separated by _x_, run **./HF x*100 N**. For example, for finding the correction for 5 wells with well-separation of 0.35 unit, run _./HF 35 5_.

# System Requirements:
 1. C++ (g++ -std=c++14)
 2. Eigen Library
 3. Lapack and Lapacke

# Installation Instruction for Eigen
Please download Eigen tarball from [here](http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2) . Follow the instruction manual from http://eigen.tuxfamily.org/dox/GettingStarted.html#title2 regarding installation of Eigen.

(**TLDR for Linux and Mac OS X users:** download and extract the tar ball, then create a symlink in /usr/local/include to the "Eigen" folder.)
