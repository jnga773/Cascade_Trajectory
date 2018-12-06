#!/bin/bash -e
#SBATCH --job-name=xcor
#SBATCH --account=nesi99999
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=500M
#SBATCH --partition=large,prepost

# set the Fortran compiler
export FC=ifort

# max_time parameter
MAX_TIME=25000

# Note: Assumes you submit the job from the BMtest directory

SCRATCH_DIR="/nesi/nobackup/${SLURM_JOB_ACCOUNT}/${USER}/scratch-${SLURM_JOB_ID}"
mkdir -p ${SCRATCH_DIR}

module load intel/2018b
module list

echo "Working in ${SCRATCH_DIR}"
cd ${SCRATCH_DIR}

# build code first
srun cmake ${SLURM_SUBMIT_DIR}/..
srun make VERBOSE=1

# check the short test passes
srun ctest --output-on-error -R short

# run three times to get average
echo "Running matrix 1..."
time srun ./BMtest/cctest_matrix ${MAX_TIME}
echo "Running matrix 2..."
time srun ./BMtest/cctest_matrix ${MAX_TIME}
echo "Running matrix 3..."
time srun ./BMtest/cctest_matrix ${MAX_TIME}

cd "${SLURM_SUBMIT_DIR}"
mv ${SCRATCH_DIR} ./benchmark-matrix-${SLURM_JOB_ID}-${FC}
