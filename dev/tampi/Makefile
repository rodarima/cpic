#CC=OMPI_CC=mcc mpicc
CC=mcc

CFLAGS=--ompss-2 -ltampi-c

CFLAGS+=$(shell mpicc --showme:compile)
LDFLAGS+=$(shell mpicc --showme:link)

all: test

clean:
	rm -f test
