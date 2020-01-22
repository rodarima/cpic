# Common
CFLAGS+=-O3
#CFLAGS+=-O0
#CFLAGS+=-O2
CFLAGS+=-g -Wall
#CFLAGS+=-masm=intel


CFLAGS+=-I../../src/
#CFLAGS+=-fno-omit-frame-pointer

# Default compiler
#CC=clang
CC=mcc --cc=clang
CFLAGS+=--ompss-2

LDLIBS+=-lm

all: main

main: ../../src/perf.o
