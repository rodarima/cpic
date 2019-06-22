#!/bin/bash

SEQ=list
OUT=time.csv

rm $OUT

while read N ; do

	LINE=$(cat out/$N | grep '^stats ' | sed 's/[^ ]*=//g;s/stats //g' | tail -1)

	echo $(($N)) $LINE >> $OUT

done < $SEQ
