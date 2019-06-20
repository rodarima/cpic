#!/bin/bash

if [ $# != 1 ]; then
	echo "Usage $0 CONFIG-FILE"
	exit
fi

FILE=$1
SEQ=list

rm conf/*



while read NP ; do

	sed "s/\<points = \[.*\]/points = [$NP, $NP]/g" $FILE > "conf/$NP"

done < $SEQ
