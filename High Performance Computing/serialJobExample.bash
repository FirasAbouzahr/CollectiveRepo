#!/bin/bash
#SBATCH -J serialJob
#SBATCH -o PETSimulation
#SBATCH -N launcher.o%j
#SBATCH -n 1
#SBATCH -p 20
#SBATCH -t development
#SBATCH -A PET

module load launcher

export LAUNCHER_WORKDIR=/home1/08038/firas/sparse/build
export LAUNCHER_JOB_FILE=launcher

${LAUNCHER_DIR}/paramrun
