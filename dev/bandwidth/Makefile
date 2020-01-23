# Common
CFLAGS+=-O3
#CFLAGS+=-O0
#CFLAGS+=-O2
CFLAGS+=-g -Wall
#CFLAGS+=-masm=intel

CFLAGS+=-DUSE_VECTOR_256

CFLAGS+=-I../../src/
#CFLAGS+=-fno-omit-frame-pointer
CFLAGS+=-march=native -mllvm -force-vector-width=8

CFLAGS+=-Rpass-analysis=loop-vectorize

# Default compiler
CC=clang
#CC=mcc --cc=clang
#CFLAGS+=--ompss-2

LDLIBS+=-lm

all: main

main: ../../src/perf.o

clean:
	rm -f *.o main
