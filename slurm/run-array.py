#!/usr/bin/env python

"""
Script to run Cascade Trajectory in array job.

"""
from __future__ import print_function
import os
import sys
import argparse
import subprocess
import datetime


SLURM_SCRIPT = "@CMAKE_BINARY_DIR@/_run-array.sl"


def main():
    # parse arguments
    parser = argparse.ArgumentParser(description="Run Cascade Trajectory array job")
    parser.add_argument('array_size', type=int, help="How many jobs to submit in the array")
    parser.add_argument('--account', default=None, help="NeSI project code (only required if not using your default project)")
    parser.add_argument('--name', default="CascadeTrajectory", help='Job name for the queue and naming output directory (default="CascadeTrajectory")')
    parser.add_argument('--time', default="01:00:00", help='Time per job, for the queue (default is "01:00:00", i.e. 1 hour)')
    parser.add_argument('--mem', default="500M", help='Memory per job, for the queue (default is "500", i.e. 500 MB)')
    args = parser.parse_args()

    print("Submitting job name: {}".format(args.name))
    print("Array size: {}".format(args.array_size))

    # need to make a directory for running the simulation in
    # based on current date and time to make it unique
    timestamp = datetime.datetime.now().strftime("%y%m%dT%H%M%S")
    output_dir = args.name + "." + timestamp

    if os.path.exists(output_dir):
        print("Error: output directory exists: '{}'".format(output_dir))
        sys.exit(1)

    # create the main output directory
    print("Output directory is: '{}'".format(output_dir))
    os.makedirs(output_dir)
    os.chdir(output_dir)

    # make an output directory for each job
    for i in range(1, args.array_size + 1):
        os.mkdir("%05d" % i)

    # build the submit command
    cmd = ["sbatch", "--array=1-{}".format(args.array_size)]
    if args.account is not None:
        cmd.append("--account={}".format(args.account))
    cmd.append("--time={}".format(args.time))
    cmd.append("--mem={}".format(args.mem))
    cmd.append("--job-name={}".format(args.name))
    cmd.append(SLURM_SCRIPT)
    print(cmd)

    # submit the job
    output = subprocess.check_output(cmd)
    print(output.strip())


if __name__ == "__main__":
    main()
