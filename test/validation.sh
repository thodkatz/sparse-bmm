#!/usr/bin/bash

PROJECT=/home/tkatz/repos/sparse-bmm
MATRICES=$PROJECT/matrices

figlet Build | lolcat
cd ../ && make serial

matrices=("jgl009/jgl009.mtx" "mycielskian3/mycielskian3.mtx" "mycielskian10/mycielskian10.mtx")
#matrices=("roadNet-PA/roadNet-PA.mtx")
#matrices=("mycielskian10/mycielskian10.mtx")
#matrices=("Stanford/Stanford.mtx")
#matrices=("G47/G47.mtx")
#matrices=("wing_nodal/wing_nodal.mtx")
#matrices=("biplane-9/biplane-9.mtx")
#matrices=("pli/pli.mtx")
#matrices=("dblp-2010/dblp-2010.mtx")

for i in "${matrices[@]}"
do
    FILE=$MATRICES/$i
    echo -e "\nFull file path: $FILE"
    pwd

    # to be compatible with the matrix market reading of the C++ program we assume that symmetric is general
    sed -i '/symmetric/ s//general/g' $FILE

    figlet Run | lolcat
    ./bin/serial $FILE

    figlet Validation | lolcat
    cd test/ && python spgemm.py $FILE && cd ../
done
