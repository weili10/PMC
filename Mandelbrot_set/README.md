The goal of this project is using MPI and OpenMP to parallelize a program for computing the number of points in the Mandelbrot Set, and make it run as fast as possible on a 4x32 core machine.

Firstly, I tried to optimize the given serial algorithm by using the tips in Wiki1, and by the symmetry property of the Mandelbrot Set. 

Then, I use MPI and OpenMP to parallelize the program on the 4x32 core machine. 

Finally, I implemented cache method.
