#!/bin/bash
#SBATCH --output=./slurm/hybrid/%j.out
#SBATCH --partition=htc
#SBATCH --qos=small
#SBATCH --time=1:00

# USAGE: `sbatch --ntasks=<threads> --export=A=<filename>,B=<filename>,F=<filename>`

export OMP_NUM_THREADS=$SLURM_NTASKS
srun -n $SLURM_NTASKS  ./bin/hybrid $A $B $F
