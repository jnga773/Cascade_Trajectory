# Cascade_Trajectory

Fortran 90 code for calculating the auto-correlations of a three-level atom. 

cross_correlation.f90 is the actual code to be worked on.

plot_correlation_traj.py is a python code using numpy and matplotlib.pyplot to plot the results.

## Benchmark test case

For the NeSI project we work in the BMtest directory, in cc.f90. This code prints out an
output file that can be compared with a reference file to ensure the results have not changed.

Check the README in BMtest.

## Building the code

CMake and a Fortran compiler (e.g. `gfortran` or `ifort`) is required to build the code.

### Mahuika

From this directory run:

```
module load CMake intel
mkdir build && cd build  # the directory can be called anything you like

# build with default options
cmake ..

# compile
make VERBOSE=1
```

To run the shorter test case:

```
ctest -R short -V
```

There is a longer test case that is used for benchmarking.
