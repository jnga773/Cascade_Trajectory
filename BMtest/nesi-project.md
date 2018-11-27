# NeSI Project

Keep track of work done in the NeSI project.

## CMake

Adding the CMake build system makes it easier to build the code without having
to remember different command line options.

On Mahuika load the CMake and compiler modules:

```
module load CMake intel
```

Then make a build directory (can be named anything you like) and build the code
in there, e.g.:

```
mkdir build
cd build

# configure
cmake .. -DDOUBLE_PRECISION=ON

# compile
make
```

There is a short and long test case; to run the short run:

```
ctest -R short -V
```

The executables are in the build directory and also in the BMtest subdirectory
underneath the build directory.

## Benchmarks

Unless mentioned otherwise, the benchmarks are run using the command
(starting from the build directory):

```
./BMtest/cctest 25000  # test original version
./BMtest/cctest_matrix 25000  # test matrix version
```

These are equivalent to running:

```
ctest -R "long$"
ctest -R long-matrix
```

except these also check results are correct.

## Initial timings

Initial timings and timings for a version with some of the calculations in the
integration section moved outside the timestep loop, using the Intel compiler
on Mahuika.

| Case                          | Time (s)    |
|-------------------------------|-------------|
| Initial, ifort                | 254.7 ± 0.0 |
| Precalc, ifort                | 229.8 ± 0.1 |


## Standardising precision

* Remove `(KIND=8)` from `REAL` and `COMPLEX` variable declarations
* Use command line options when building the code manually:
  - default if you don't specify any extra option is single precision
  - if you specify `-r8` with ifort or `-fdefault-real-8` with gfortran you
    will have double precision
* With the CMake build:
  - the default is to build using double precision
  - to build with single precision use: `cmake -DDOUBLE_PRECISION=OFF`

Note, timings will probably be different now.

## Benchmarks after precision change

Timings on Mahuika with Intel compiler (note this is the original version apart
from precision change). Single precision was configured with `cmake -DDOUBLE_PRECISION=OFF`.

| Precision    | Timings (s)   |
|--------------|---------------|
| Double       | 243.8 ± 1.2   |
| Single       | 221.3 ± 1.3   |

Approximately 9.2 % reduction in run time from using single precision.

## Matrix version

There is also "matrix" version, which constructs a matrix at the beginning
and uses that for the integration, instead of recalculating a number of
values. Need to add benchmarks.

| Version      | Timings (s)   |
|--------------|---------------|
| Original     | 243.8 ± 1.2   |
| Matrix       | 240.6 ± 1.0   |

Initial results show the matrix version has similar performance. Need to look
closer at the code. In both cases there are optimisations that can be made,
such as taking out the `i * dt`. Will work on these next.





## Todo

* Finish optimising the current code doing a single run
* add scripts for running array jobs for many shorter runs at once
* scripts for parameter sweeps?

