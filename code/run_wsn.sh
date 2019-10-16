#!/bin/bash
#SBATCH --ntasks=21
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lche0021@student.monash.edu

hostname
make build
make run OMP_NUM_THREADS=4
