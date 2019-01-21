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

#### Intel Compiler

From this directory run:

```
module load CMake intel
mkdir build-intel && cd build-intel  # the directory can be called anything you like

# build with default options
FC=ifort cmake ..

# compile using Intel compiler
make VERBOSE=1
```

#### Cray compiler

From this directory run:

```
module load CMake PrgEnv-cray craype-broadwell
mkdir build-cray && cd build-cray  # the directory can be called anything you like

# build with default options
FC=ftn cmake ..

# compile using Intel compiler
make VERBOSE=1
```

#### GNU compiler

From this directory run:

```
module load CMake intel
mkdir build-gnu && cd build-gnu  # the directory can be called anything you like

# build with default options
FC=gfortran cmake ..

# compile using Intel compiler
make VERBOSE=1
```

#### Running tests

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

## Array job

There are scripts to make it easier to run many copies of the same simulation
as an array job on Mahuika.

Here is an example of how to run an array job:

```
# configure and build the code like normal
mkdir array-job-test
cd array-job-test
module load PrgEnv-cray CMake
FC=ftn cmake ..
make -j

# edit the params.nml file if required (e.g. to update run time)
# each simulation will use the same input file
nano params.nml

# look at the options in the array job script
./run-array.py --help

# submit 20 jobs (with default memory and time limits, 500 MB and 1 hour)
./run-array.py --name=arraytest 20

# check the queue
squeue -u $USER
```

The output folder was printed to standard output when running the Python script
and starts with "arraytest.XXXXX" where XXXXX is a timestamp and "arraytest"
was the argument passed to `--name`. Each job in the array has its own directory
where its output is stored.

This runs the `cctest_matrix` executable by default but this can be changed, for
example to the `cctest` executable:

```
cmake . -DARRAY_JOB_EXE=BMtest/cctest
```

### Parameter sweeps

Sweeping over values for a single parameter can also be done using the array
job script. Similar to the above, this is an example of a parameter sweep for
the `omega` parameter, with values from 1.0 to 10.0 by 1.0, and submitting 10
jobs per value:

```
# configure and build the code like normal
mkdir parameter-sweep-test
cd parameter-sweep-test
module load PrgEnv-cray CMake
FC=ftn cmake ..
make -j

# edit the params.nml file if required (e.g. to update run time)
# each simulation will use the same input file with the parameter
# sweep variable modified by the python script
nano params.nml

# look at the options in the array job script
./run-array.py --help

# submit 10 jobs per value in the sweep
./run-array.py --name=sweeptest --sweep=omega --sweep=min=1.0 --sweep-max=10.0 --sweep-step=1 10

# check the queue
squeue -u $USER
```
