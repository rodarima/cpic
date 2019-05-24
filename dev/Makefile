CC=gcc

CFLAGS:=-g -pg
LDLIBS:=-lm -lstdc++ -lGL -lmgl2

CFLAGS+=`pkg-config --cflags glfw3`
LDLIBS+=`pkg-config --libs glfw3`

CFLAGS+=`pkg-config --cflags fftw3`
LDLIBS+=`pkg-config --libs-only-L fftw3` -lfftw3_mpi -lfftw3

CFLAGS+=`mpicc --showme:compile`
LDLIBS+=`mpicc --showme:link`


#BIN=mgl glvid fftmpi
BIN=fftmpi

all: $(BIN)

fftmpi: fftmpi.c

clean:
	rm -rf *.o $(BIN)
