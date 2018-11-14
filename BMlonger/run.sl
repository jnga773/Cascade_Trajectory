#!/bin/bash -e
#SBATCH --job-name=xcor
#SBATCH --account=nesi99999
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --mem=100M
#SBATCH --partition=large,prepost

SCRATCH_DIR="/nesi/nobackup/${SLURM_JOB_ACCOUNT}/${USER}/scratch_${SLURM_JOB_ID}"
mkdir -p ${SCRATCH_DIR}

module load intel/2018b
module list

make -f makefile.ifort clean
make -f makefile.ifort
cp cc ${SCRATCH_DIR}/
cd ${SCRATCH_DIR}
echo "Running in ${SCRATCH_DIR}"

echo "Running 1..."
time srun ./cc
echo "Running 2..."
time srun ./cc
echo "Running 3..."
time srun ./cc

cd "${SLURM_SUBMIT_DIR}"
mv ${SCRATCH_DIR} ./run-${SLURM_JOB_ID}
