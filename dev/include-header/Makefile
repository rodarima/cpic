#CC=gcc
CC=mcc -k
#CC=/usr/bin/mcc
#CC=clang

all: main main2

main: main.c obj.o

main2.c: main.c
	gcc main.c -E main2.c

main2: main2.c obj.o

clean:
	rm *.o main main2
