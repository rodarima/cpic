#!/bin/bash

echo "THE WORKER.SH BEGINS" >> log

#valgrind --log-file=log-vg ./worker >> log
#valgrind ./worker
./worker
