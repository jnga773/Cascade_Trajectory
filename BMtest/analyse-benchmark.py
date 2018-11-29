#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import math

import numpy as np


if len(sys.argv) != 2:
    print("Usage: ./analyse-benchmark.py <SLURM_OUTPUT_FILE>")
    sys.exit(1)

fn = sys.argv[1]
print(fn)

# load lines from slurm output file
with open(fn) as fh:
    lines = fh.readlines()

# take lines that start with "real"
lines = [line for line in lines if line.startswith("real")]
print(len(lines))

# extract times and convert to seconds
times = []
for line in lines:
    tstr = line.split()[1]
    array = tstr.rstrip("s").split("m")
    runtime = float(array[0]) * 60.0 + float(array[1])
    times.append(runtime)
times = np.asarray(times)
print(times)

# compute mean and standard error
mean = np.mean(times)
stderr = np.std(times) / math.sqrt(float(len(times)))
print("%.1f ± %.1f s" % (mean, stderr))
