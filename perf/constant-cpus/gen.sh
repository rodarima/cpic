#!/bin/bash

SEQ=list

rm -f jobs/*

while read NCPUS ; do

	NP=$((32/$NCPUS))

	sed "s/@NPROCS@/$NP/g;s/@NCPUS@/$NCPUS/g" job.sh > "jobs/$NCPUS"

done < $SEQ
