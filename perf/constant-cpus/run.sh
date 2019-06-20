#!/bin/bash

SEQ=list

rm -f out/* err/*

PROC_LIST=$(cat $SEQ | tr '\n' ' ')

for NP in $PROC_LIST ; do

	NCPUS=$((48/$NP))

	echo Running with $NP process, with $NCPUS CPUs per process
	srun  --slurmd-debug=4 --cpu-bind=v,cores -c $NCPUS -n $NP ../../cpic base.conf > out/$NP 2>err/$NP
	#srun -c $NCPUS -n 1 echo $NP
done
