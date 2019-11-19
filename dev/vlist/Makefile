# Common
#CFLAGS=-O3
CFLAGS=-g -Wall

CC=c99

# GCC
#CC=gcc
#CFLAGS+=-fopenmp
#CFLAGS+=-fopenmp-simd
#CFLAGS+=-fopt-info-vec-all
#CFLAGS+=-march=native
#CFLAGS+=-mtune=native

# Intel
#CC=icc
#CFLAGS+=-axCORE-AVX512
#CFLAGS+=-xCOMMON-AVX512 -qopt-report=5

#CC=mcc
#CFLAGS+=--ompss-2 -k --v
#CFLAGS+=-fsanitize=address -fno-omit-frame-pointer

LDLIBS+=-lm

all: test

test: test.c

clean:
	rm -f test *.o

