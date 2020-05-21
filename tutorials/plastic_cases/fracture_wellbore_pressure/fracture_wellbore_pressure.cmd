#!/bin/bash
#parallel job using 28 proccesors and runs for 2.5 hour (max)
#SBATCH -N 1 # node count
#SBATCH --ntasks-per-node=28 #number of threads desired
#SBATCH -t 2:35:00

# Load openmpi environment
module load openmpi/gcc

srun plasticHBIF -parallel 2>&1 1>log

