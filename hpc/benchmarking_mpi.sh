module load gcc openmpi

# go to root project
cd ../

matrices=("n1e6_d4.mtx" "n2e6_d4.mtx" "n3e6_d4.mtx" "n5e6_d4.mtx")
path="matrices/"

echo "OPENMPI"
make openmpi

for j in ${matrices[@]}; do
    A=$path"A_"$j
    B=$path"B_"$j
    F=$path"F_"$j
    for task in 1 2 4 6 8 10 12 14 16 18 19; do
    #for task in 20; do
        echo "task: "$task
        echo $A
        sbatch --ntasks=$(($task+1)) --export=A=$A,B=$B,F=$F mpi.sh
    done
done
