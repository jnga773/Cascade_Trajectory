#!/bin/bash -e
# Don't edit this file, look at run_array.py instead
#SBATCH --partition=large
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --output=%5a/slurm-%A-%5a.out

# program to run (set from CMake)
EXECUTABLE=@ARRAY_JOB_EXECUTABLE@

# check we are in an array job
if [[ -z "${SLURM_ARRAY_TASK_ID}" ]]; then
    echo "Must be submitted as an array job"
    exit 1
fi

# change to directory for this element
mydir=$(printf "%05d" ${SLURM_ARRAY_TASK_ID})
cd "${mydir}"

# run
srun ${EXECUTABLE}
