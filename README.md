# Cascade_Trajectory

Fortran 90 code for calculating the auto-correlations of a three-level atom. 

cross_correlation.f90 is the actual code to be worked on.

plot_correlation_traj.py is a python code using numpy and matplotlib.pyplot to plot the results.

## Benchmark test case

For the NeSI project we work in the BMtest directory, in cc.f90. This code prints out an
output file that can be compared with a reference file to ensure the results have not changed.

There are two versions of the benchmark case; the original version is cc.f90 and the matrix
version is cc_matrix.f90. Both give the same results.

Check the README in BMtest too.

## Building the code

CMake and a Fortran compiler (e.g. `gfortran` or `ifort`) is required to build the code.

### Mahuika

From this directory run:

```
module load CMake intel
mkdir build && cd build  # the directory can be called anything you like

# build with default options
cmake ..

# compile using Intel compiler
FC=ifort make VERBOSE=1
```

To run the shorter test cases:

```
ctest -R short -V
```

There are two shorter test cases, one for the original version and one for the matrix
version.

To list all available test cases run: `ctest -N`. Then `ctest -R <regular-expression> -V`
can be used to run tests that match the given regular expression.

There is a longer test case that is used for benchmarking, e.g.:

```
ctest -R long -V
```
