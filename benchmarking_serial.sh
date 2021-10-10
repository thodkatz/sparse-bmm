#!/bin/bash

# serial
echo "SERIAL"
for i in ${matrices[@]}; do
    A=$path"A_"$i
    B=$path"B_"$i
    F=$path"F_"$i

    make serial
    ./bin/serial $A $B $F
done
