CC=mcc
CFLAGS=-g --ompss-2

CFLAGS+=-fsanitize=address -fno-omit-frame-pointer

LDLIBS+=-lm

all: test

test: test.c

clean:
	rm -f test *.o
