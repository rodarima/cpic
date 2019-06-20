#!/bin/bash

rm out/* err/*

for F in conf/* ; do

	echo "Running $F"

	N=$(basename $F)

	../../cpic $F > out/$N 2>err/$N

done
