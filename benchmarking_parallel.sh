#!/bin/bash

#matrices=("n5e5_d4.mtx" "n1e6_d4.mtx" "n2e6_d4.mtx" "n4e6_d4.mtx" "n5e6_d4.mtx")
matrices=("n1e6_d2.mtx")
path="matrices/"

echo "tasks,type,matrix,time" > log.csv

# openmp, openmpi
types=("openmp" "openmpi")
echo "OPENMP,OPENMPI"
for i in ${types[@]}; do
    if [[ $i == "openmp" ]]; then
        make openmp
    else 
        make openmpi
    fi

    for j in ${matrices[@]}; do
        A=$path"A_"$j
        B=$path"B_"$j
        F=$path"F_"$j
        for task in 1 2 4 6 8 10 12 14 16 18 20; do
            echo "type: "$i", task:"$task
            if [[ $i == "openmp" ]]; then
                export OMP_NUM_THREADS=$task
                echo $A
                ./bin/openmp $A $B $F
            else
                mpirun --oversubscribe -n $(($task+1)) ./bin/openmpi $A $B $F
            fi
        done
    done
done

# hybrid
echo "HYBRID"
make hybrid
for j in ${matrices[@]}; do
    A=$path"A_"$j
    B=$path"B_"$j
    F=$path"F_"$j
    for task in 1 2 4 6; do
        echo $task
        export OMP_NUM_THREADS=$task-1
        mpirun --oversubscribe -n $(($task+1)) ./bin/openmpi $A $B $F
    done
done
