CC=clang
CFLAGS=-g -Wall -std=c11 -O0 -fPIE -march=core-avx2 -fno-omit-frame-pointer -fstack-protector-all
LDFLAGS=-pthread -Wl,-rpath -Wl,/usr/lib/openmpi -Wl,--enable-new-dtags -L/usr/lib/openmpi
LDLIBS+=-lmpi -lm -lconfig -lfftw3_threads -lfftw3_mpi -lfftw3

OMPSS2=-fompss-2
#OMPSS2=--ompss-2
CFLAGS+=$(OMPSS2)
LDFLAGS+=$(OMPSS2)

#CC=mcc

all: fft

fft: fft.o perf.o

fft.o: fft.c

clean:
	rm -f *.o fft

.PHONY: clean
