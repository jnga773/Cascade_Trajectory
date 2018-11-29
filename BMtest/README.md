# Benchmark case

On Mahuika we can run the benchmark using the Slurm scripts in this directory.

Running:

```
sbatch run-benchmark-ifort.sl
```

will generate an output directory ("benchmark-XXXXXX-ifort") and output file
("slurm-XXXXXX.out"). The build and results are in the directory. Timings are
in the output file. Run:

```
python analyse-benchmark.py slurm-XXXXXX.out
```

to get the timings (replacing XXXXXX with the slurm job ID).

Tests can be run using the CTest in the CMake build folder - see the main
README.md file for details.
