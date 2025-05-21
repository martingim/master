#!/bin/bash

for i in "7" "8" "9"; do for j in "10" "20" "40"; do
    echo "nx:"$i "nl:"$j
    ./multilayer -L $i -nl $j
done; done;