# Common
#CFLAGS+=-O3
CFLAGS+=-O0
#CFLAGS+=-O2
CFLAGS+=-g -Wall
#CFLAGS+=-masm=intel

CFLAGS+=-DUSE_VECTOR_256

CFLAGS+=-I../../src/
#CFLAGS+=-fno-omit-frame-pointer
CFLAGS+=-march=core-avx2
#CFLAGS+=-fsanitize=address

# Default compiler
CC=clang

LDLIBS+=-lm

all: main

main: ../../src/perf.o

clean:
	rm -f *.o main
