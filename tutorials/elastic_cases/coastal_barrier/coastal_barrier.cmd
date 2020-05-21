#!/bin/bash
#parallel job using 56 proccesors and runs for 15 hours
#SBATCH -N 2 # node count
#SBATCH --ntasks-per-node=28 #number of threads desired
#SBATCH -t 15:30:00

# Load openmpi environment
module load openmpi/gcc

srun elasticHBIF -parallel 2>&1 1>log

