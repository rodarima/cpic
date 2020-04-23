CC=mpicc
CFLAGS=-Wall
LDLIBS=-lfftw3_mpi -lfftw3 -lm

example: example.o perf.o

clean:
	rm -f *.o example

.PHONY: clean
