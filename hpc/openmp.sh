#!/bin/bash
#SBATCH --nodes=1
#SBATCH --output=./slurm/openmp/%j.out
#SBATCH --partition=batch
#SBATCH --qos=small
#SBATCH --time=1:00

# USAGE: `sbatch --ntasks=<threads> --export=A=<filename>,B=<filename>,F=<filename>`

export OMP_NUM_THREADS=$SLURM_NTASKS
./bin/openmp $A $B $F
#./bin/openmp matrices/A_n1e6_d4.mtx matrices/B_n1e6_d4.mtx matrices/F_n1e6_d4.mtx