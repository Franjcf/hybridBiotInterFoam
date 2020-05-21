#!/bin/bash
#parallel job using 4 proccesors and runs for 1 hour (max)
#SBATCH -N 1 # node count
#SBATCH --ntasks-per-node=1 #number of threads desired
#SBATCH -t 0:10:00
# Load openmpi environment
module load openmpi/gcc

srun elasticHBIF 2>&1 1>log

