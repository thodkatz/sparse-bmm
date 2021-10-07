#!/usr/bin/bash

PROJECT=/home/tkatz/repos/sparse-bmm
MATRICES=$PROJECT/matrices

#target=("serial" "openmp")
target=("hybrid")
export OMP_NUM_THREADS=8

for i in ${target[@]}; do
    figlet Build | lolcat
    cd ../ && make $i
    
    A="A.mtx"
    B="B.mtx"
    F="F.mtx"
    
    #A="com-Youtube/com-Youtube.mtx"
    #A="pli/pli.mtx"
    #A="jgl009/jgl009.mtx"
    #A="mycielskian3/mycielskian3.mtx"
    #B=$A
    #F=$A
    
    FILEA=$MATRICES/$A
    FILEB=$MATRICES/$B
    FILEF=$MATRICES/$F
    
    echo -e "\nFull file path matrix A: $FILEA"
    echo -e "\nFull file path matrix B: $FILEB"
    echo -e "\nFull file path matrix F: $FILEF"
    
    figlet Run | lolcat
    if [[ "$i" == "hybrid" ]]; then
        mpirun --oversubscribe -n 2 ./bin/$i $FILEA $FILEB $FILEF
    else
        ./bin/$i $FILEA $FILEB $FILEF
    fi
    
    figlet Validation | lolcat
    #cd test/ && python spgemm.py $FILEA $FILEB $FILEF && cd ../
    cd test/ && python spgemm.py $FILEA $FILEB $FILEF
    
done
