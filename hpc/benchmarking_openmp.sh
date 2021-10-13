#!/bin/bash

# go to root project
cd ../

module load gcc openmpi

matrices=("n1e6_d4.mtx" "n2e6_d4.mtx" "n3e6_d4.mtx" "n5e6_d4.mtx")
path="matrices/"

#echo "tasks,type,matrix,time" > log.csv

# openmp
echo "OPENMP"
make openmp

for j in ${matrices[@]}; do
    A=$path"A_"$j
    B=$path"B_"$j
    F=$path"F_"$j
    for task in 1 2 4 6 8 10 12 14 16 18 20; do
        echo "task: "$task
        echo $A
        sbatch --ntasks=$task --export=A=$A,B=$B,F=$F openmp.sh
    done
done