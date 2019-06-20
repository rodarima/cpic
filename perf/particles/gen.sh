#!/bin/bash

SEQ=list

while read NP ; do

	sed "s/\<particles = .*/particles = $NP/g" base.conf > "conf/$NP"
	sed "s/@PARAM@/$NP/g" job.sh > "jobs/$NP"

done < $SEQ
