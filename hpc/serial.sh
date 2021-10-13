#!/bin/bash
#SBATCH --output=./slurm/serial/%j.out
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --qos=small
#SBATCH --time=10:00

./bin/serial $A $B $F