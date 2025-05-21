#!/bin/bash

r="1"
for n_threads in "1" "2" "3" "4" "5" "6" "7" "8"
do
    for L in "10"
    do
        #echo "./piston -L ${L} -r ${r} |& tee -a LEVEL${L}_r${r}_log.txt"
        ./piston -L "${L}" -r "${r}" --nthreads "${n_threads}" |& tee  "LEVEL${L}_r${r}_n${n_threads}_log_5s.txt"
        #./piston -L "${L}" -r "${r}" |& tee  "run11/LEVEL${L}_r${r}_log.txt"
    done
done