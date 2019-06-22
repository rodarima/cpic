#!/bin/bash

SEQ=list

N_LIST=$(cat $SEQ | tr '\n' ' ')

rm -f out/* err/*

for N in $N_LIST ; do
	echo "Running $N"

	sbatch jobs/$N
done
