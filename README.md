# Hartree-fock-scf

This project is a solution to the mid-semester examination of course-<a href="http://www.niser.ac.in/sps/course/p451-advanced-solid-state-physics" target="_blank">P451</a>, Advanced Solid State Physics, offered during Fall, 2017 in School of Physical Sciences, NISER, Bhubaneswar. This project, with continued improvement, seeks to develop a generalized routine for various implementations of Hartree-Fock like mean-field techniques.

This branch implements the averaged density Hartree-Fock method to find SCF corrections to a system with given potential and number of electrons.

# Building and Running

The potential, its parameters and other constants are defined in configuration.hpp file.  To change any fixed paramters, such as the number of electrons or the density of grid points etc., the corresponding variables can be modified within the configuration.hpp file.

## Defining the potential

The potential is defined in the [**double V (double x)**](https://github.com/bineet-dash/hartree-fock/blob/1e6b9e808b0ea9da7f681aaa5f67ecbac277f2cf/configuration.hpp#L18) function in the [configuration.hpp](configuration.hpp) file. Modify the definition as required, but _do not_ add other arguments to the function. According to the nature of the potential, modify the upper limit, lower limit and the density of grid points in the [configuration.hpp](configuration.hpp) file.

## Defining the parameters in the potential
Extra parameters, if needed in the definition of potential,
can be declared in the [configuration.hpp](configuration.hpp) file. If these variables need to be used in the main program (e.g. if needed as an user input), then declare these variables as extern in [common_globals.hpp](common_globals.hpp).

Suppose you need an extra parameter *A*  in the definition of potential, which has to be obtained from user input. Then
declare **double A**  in the beginning of [configuration.hpp](configuration.hpp) and add **extern double A** at the end of [common_globals.hpp](common_globals.hpp). Now you can use the variable *A* in the main() as a global variable without declaring it anymore.

To build the program compile using "g++ -std=c++14 -o HF hartree-fock.cpp configuration.hpp -llapack -llapacke".

# System Requirements:
 1. C++ (g++ -std=c++14)
 2. Eigen Library
 3. Lapack and Lapacke


# Installation Instruction for Eigen
Please download Eigen tarball from [here](http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2) . Follow the instruction manual from http://eigen.tuxfamily.org/dox/GettingStarted.html#title2 regarding installation of Eigen.

(**TLDR for Linux and Mac users:** download and extract the tar ball, then create a symlink in /usr/local/include to the "Eigen" folder.)
