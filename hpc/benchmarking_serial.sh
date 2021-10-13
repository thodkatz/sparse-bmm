#!/bin/bash

# go to root project
cd ../

module load gcc

matrices=("n1e6_d4.mtx" "n2e6_d4.mtx" "n3e6_d4.mtx" "n5e6_d4.mtx")
path="matrices/"

make serial

# serial
echo "SERIAL"
for i in ${matrices[@]}; do
    A=$path"A_"$i
    B=$path"B_"$i
    F=$path"F_"$i

    sbatch --export=A=$A,B=$B,F=$F serial.sh
done