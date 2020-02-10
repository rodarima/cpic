#!/bin/bash

printf "%d %d %d\n" \
	$(getconf LEVEL1_DCACHE_LINESIZE) \
	$(getconf LEVEL2_CACHE_LINESIZE) \
	$(getconf LEVEL3_CACHE_LINESIZE)
