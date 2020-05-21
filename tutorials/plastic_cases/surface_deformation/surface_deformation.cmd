#!/bin/bash
#parallel job using 20 proccesors and runs for 4.5 hours
#SBATCH -N 1 # node count
#SBATCH --ntasks-per-node=20 #number of threads desired
#SBATCH -t 4:30:00

# Load openmpi environment
module load openmpi/gcc

srun plasticHBIF -parallel 2>&1 1>log

