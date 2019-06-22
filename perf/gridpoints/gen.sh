#!/bin/bash

SEQ=list

rm -f jobs/* conf/*

while read NP ; do

	sed "s/\<points = \[.*\]/points = [$NP, $NP]/g" base.conf > "conf/$NP"
	sed "s/@PARAM@/$NP/g" job.sh > "jobs/$NP"

done < $SEQ
