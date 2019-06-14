#MCC=mcc
#MCC_CFLAGS=--ompss-2
#CC=gcc
#LDLIBS=-lnanos6 -ldl /usr/lib/nanos6-main-wrapper.o
#
#BIN=main
#
#all: $(BIN)
#
#%.mcc.c: %.c
#	mcc $(MCC_CFLAGS) -y $^ -o $@
#
#%.mcc.o: %.mcc.c
#	$(COMPILE.c) $(OUTPUT_OPTION) $<
#
#main: main.mcc.o
#	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@
#
#clean:
#	rm -f $(BIN) *.o

CC=mcc
CFLAGS=--ompss-2

BIN=main

all: $(BIN)

main: main.c

clean:
	rm -f $(BIN) *.o
