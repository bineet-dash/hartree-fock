# advanced-cmp-midsem

This project is a solution to the mid-semester examination of course-[P451](http://www.niser.ac.in/sps/course/p451-advanced-solid-state-physics ), Advanced Solid State Physics, offered during Fall, 2017 in School of Physical Sciences, NISER, Bhubaneswar.

This project, with continued improvement, seeks to develop a generalized routine for various implementations of Hartree-Fock like mean-field techniques. I've already implemented Closed-Shell Roothan method, Density averaged Hartree-Fock, and Independent Density Hartree-Fock methods.

# System Requirements:
 1. C++ (g++ -std=c++14)
 2. Eigen Library
 3. Lapack and Lapacke


# Installation Instruction for Eigen
Please download Eigen tarball from [here](http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2) . Follow the instruction manual from http://eigen.tuxfamily.org/dox/GettingStarted.html#title2 regarding installation of Eigen.

(**TLDR for Linux and Mac OS X users:** download and extract the tar ball, then create a symlink in /usr/local/include to the "Eigen" folder.)

The codes for Roothan have some issues as of now. The working codes for Hartree-Fock methods are placed in hartree-fock directory.

To build the program compile using "g++ -std=c++14 -o HF hartree-fock.cpp hermite_polynomial.cpp -llapack -llapacke".
