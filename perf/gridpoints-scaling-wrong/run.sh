#!/bin/bash

rm -f out/* err/*

for F in jobs/* ; do
	echo "Running $F"

	sbatch $F

	read nop
done
