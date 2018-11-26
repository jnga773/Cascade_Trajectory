# NeSI Project

Keep track of work done in the NeSI project.

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
* Use command line options instead
  - default if you don't specify any extra option is single precision
  - if you specify `-r8` with ifort or `-fdefault-real-8` with gfortran you
    will have double precision

Note, timings will probably be different now.

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

