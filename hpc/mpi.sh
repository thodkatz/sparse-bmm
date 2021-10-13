#!/bin/bash
#SBATCH --output=./slurm/mpi/%j.out
#SBATCH --partition=htc
#SBATCH --qos=small
#SBATCH --time=1:00

# USAGE: `sbatch --ntasks=<threads> --export=A=<filename>,B=<filename>,F=<filename>`

srun -n $SLURM_NTASKS ./bin/openmpi $A $B $F
