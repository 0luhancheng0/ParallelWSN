#!/bin/bash
#SBATCH --ntasks=21
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lche0021@student.monash.edu
#SBATCH --time=00:20:00
hostname
module load mpip
make run MESSAGE_LEN=32768 X=4 Y=5 N_ITERATION=2000 N_BIT_RAND=3 INTERVAL=10
