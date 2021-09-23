#!/usr/bin/bash

# to be compatible with the matrix market reading of the C++ program we assume that symmetric is general
PROJECT=~/repos/sparse-bmm
FILE=$PROJECT/$1

sed -i '/symmetric/ s//general/g' $FILE

figlet Build | lolcat
cd ../ && make serial

figlet Run | lolcat
./bin/serial $FILE

figlet Validation | lolcat
cd test/ && python spgemm.py $FILE
