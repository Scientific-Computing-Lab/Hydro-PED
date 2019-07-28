#!/bin/bash
#SBATCH -N 1
#SBATCH -o out_Vectorized
#SBATCH --exclusive

#export OMP_SCHEDULE=dynamic
#export OMP_NUM_THREADS=118
export  MIC_USE_2MB_BUFFERS=64k
./poro.exe
