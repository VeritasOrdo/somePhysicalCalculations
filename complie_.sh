#!/bin/bash

#SBATCH --partition=vip_39 --nodes=1 --ntasks-per-node=1 --cpus-per-task=128 --mem=1G --job-name=compile --output=compile.out --error=compile.err

#envirnoment
module purge
module load gcc/9.3.0
module load mpi/openmpi/4.1.1-gcc9.3.0_new 

# openmp
export OMP_NUM_THREADS=128

# run the compiled code
./DrawRadiationPlots/DrawRadiationPlots.out

# End of file
