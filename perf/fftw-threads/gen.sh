#!/bin/bash

SEQ=list

rm -f jobs/* conf/*

while read N ; do

	sed "s/@NT@/$N/g" base.conf > "conf/$N"
	sed "s/@NT@/$N/g" job.sh > "jobs/$N"

done < $SEQ
