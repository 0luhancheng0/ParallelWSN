#!/bin/bash
#SBATCH --ntasks=21
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lche0021@student.monash.edu
#SBATCH --time=00:20:00
make wsn
OMP_NUM_THREADS=1
mpirun -n 21 -x OMP_NUM_THREADS=$(OMP_NUM_THREADS) ./wsn
