#!/bin/bash

SEQ=list

N_LIST=$(cat $SEQ | tr '\n' ' ')

rm -f out/* err/*

for N in $N_LIST ; do

	echo "Running $N"

	#../../cpic $F > out/$N 2>err/$N
	#srun  --cpu-bind=v,cores -c 32 -n 1 ../../cpic $F > out/$N 2>err/$N
	sbatch jobs/$N

done
