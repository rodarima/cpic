#CC=OMPI_CC=mcc mpicc
CC=I_MPI_CC=mcc mpiicc

CFLAGS=-g --ompss-2

#CFLAGS+=$(shell mpicc --showme:compile)
#LDFLAGS+=$(shell mpicc --showme:link)

#CFLAGS+=-fsanitize=address -fno-omit-frame-pointer

LDLIBS+=-lm
LDLIBS+=-ltampi-c
#LDLIBS+=-lompitrace

all: test

test: perf.o

run: test
	mpirun -n 2 --map-by NUMA:PE=2 ./test

clean:
	rm -f test *.o
