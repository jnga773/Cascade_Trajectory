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
import shutil


SLURM_SCRIPT = "@CMAKE_BINARY_DIR@/_run-array.sl"


def main():
    # parse arguments
    parser = argparse.ArgumentParser(description="Run Cascade Trajectory array job")
    parser.add_argument('array_size', type=int, help="How many jobs to submit (for each sweep value if doing a parameter sweep)")
    parser.add_argument('--account', default=None, help="NeSI project code (only required if not using your default project)")
    parser.add_argument('--name', default="CascadeTrajectory", help='Job name for the queue and naming output directory (default="CascadeTrajectory")')
    parser.add_argument('--time', default="01:00:00", help='Time per job, for the queue (default is "01:00:00", i.e. 1 hour)')
    parser.add_argument('--mem', default="500M", help='Memory per job, for the queue (default is "500", i.e. 500 MB)')
    parser.add_argument('--input', default='params.nml', help='Input namelist file (default="params.nml")')
    parser.add_argument('--sweep', default=None, help='Name of variable in nml file to sweep over (default=disabled)')
    parser.add_argument('--sweep-min', type=float, default=1.0, help="First value in the parameter sweep (default=1.0)")
    parser.add_argument('--sweep-max', type=float, default=40.0, help="Last value in the parameter sweep (default=40.0)")
    parser.add_argument('--sweep-step', type=float, default=1.0, help="Interval between values in the parameter sweep (default=1.0)")
    parser.add_argument('--dryrun', action="store_true", default=False, help="Do everything except submit the job (for debugging)")
    args = parser.parse_args()

    print("Submitting job name: {}".format(args.name))
    print("Array size: {}".format(args.array_size))

    # full path to input file
    input_file = os.path.abspath(args.input)
    print("Input file: {}".format(input_file))

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

    # not parameter sweep
    if args.sweep is None:
        # number of jobs
        num_jobs = args.array_size

        # make an output directory for each job and copy namelist
        for i in range(1, args.array_size + 1):
            dirname = "%05d" % i
            os.mkdir(dirname)
            shutil.copy(input_file, dirname)

    # parameter sweep
    else:
        sweep_size = int((args.sweep_max - args.sweep_min) / args.sweep_step + 1)
        num_jobs = args.array_size * sweep_size
        print("Sweep for {} from {} to {} by {} ({} values)".format(args.sweep, args.sweep_min,
                                                                    args.sweep_max, args.sweep_step,
                                                                    sweep_size))

        # read in namelist
        with open(input_file) as fh:
            nml_lines = fh.readlines()

        # check the parameter exists and find the line number
        match_count = 0
        for i, line in enumerate(nml_lines):
            if line.strip().startswith("{} =".format(args.sweep)):
                match_count += 1
                sweep_line = i
        if match_count != 1:
            print("There should be exactly one line in the nml file for the sweep parameter", match_count)
            sys.exit(1)

        # make output directory for each job and write modified namelist file
        array_index = 1
        for i in range(sweep_size):
            sweep_value = args.sweep_min + i * args.sweep_step
            for j in range(args.array_size):
                # create directory for this task
                dirname = "%05d" % array_index
                array_index += 1
                os.mkdir(dirname)

                # modify nml file parameter value
                nml_lines[sweep_line] = "    %s = %s\n" % (args.sweep, str(sweep_value))

                # write the file
                with open(os.path.join(dirname, "params.nml"), "w") as fh:
                    fh.write("".join(nml_lines))

    # build the submit command
    cmd = ["sbatch", "--array=1-{}".format(num_jobs)]
    if args.account is not None:
        cmd.append("--account={}".format(args.account))
    cmd.append("--time={}".format(args.time))
    cmd.append("--mem={}".format(args.mem))
    cmd.append("--job-name={}".format(args.name))
    cmd.append(SLURM_SCRIPT)
    print("Array job size is", num_jobs)
    print(cmd)

    # submit the job
    if not args.dryrun:
        output = subprocess.check_output(cmd)
        print(output.strip())


if __name__ == "__main__":
    main()
