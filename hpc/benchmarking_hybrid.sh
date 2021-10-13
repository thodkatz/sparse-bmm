module load gcc openmpi

# go to root project
cd ../

matrices=("n1e6_d4.mtx" "n2e6_d4.mtx" "n3e6_d4.mtx" "n5e6_d4.mtx")
path="matrices/"

echo "HYBRID"
make hybrid
for j in ${matrices[@]}; do
    A=$path"A_"$j
    B=$path"B_"$j
    F=$path"F_"$j
    for task in 1 2 3 4 5 6; do
    #for task in 3; do
        echo $task
        #export OMP_NUM_THREADS=$task
        #mpirun --oversubscribe -n $(($task+1)) ./bin/openmpi $A $B $F
        sbatch --ntasks=$(($task+1)) --cpus-per-task=$task --export=A=$A,B=$B,F=$F,OMP_NUM_THREADS=$task hpc/hybrid.sh
        #sbatch --nodes=$(($task+1)) --ntasks-per-node=1 --cpus-per-task=$task --export=A=$A,B=$B,F=$F --exclude=cn[4-5] hybrid.sh
    done
done
