#!/bin/bash
nx=("7" "8" "9")
nl=("10" "20" "40")
for i in ${nx[*]}; do for j in ${nl[*]}; do
    echo "nx:"$i "nl:"$j
    ./multilayer $i $j
done; done;