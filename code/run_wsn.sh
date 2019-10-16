#!/bin/bash
#SBATCH --ntasks=21
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lche0021@student.monash.edu
#SBATCH --time=00:20:00
hostname
module load mpip
make build
make run
