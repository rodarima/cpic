#CC=gcc
#CC=mcc -k
#CC=/usr/bin/mcc
CC=clang
CFLAGS+=-I../../src
CFLAGS+=-march=core-avx2
#CFLAGS+=-O3

all: main main.s

main.s: main.c
	$(CC) $(CFLAGS) -S -o $@ $^

main: main.s
	$(CC) -g $(CFLAGS) -o $@ $^

clean:
	rm -f *.o *.s main
