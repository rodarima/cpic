#!/bin/bash

SEQ=list
OUT=csv/time.csv

rm $OUT

while read NP ; do

	LINE=$(cat out/$NP | grep '^stats ' | sed 's/[^ ]*=//g;s/stats //g' | tail -1)

	echo $NP $LINE >> $OUT

done < $SEQ

python regression.py
