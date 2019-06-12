CC=gcc
MCC_CFLAGS=--ompss-2
CFLAGS=-fPIC
LDLIBS=-lnanos6 -ldl /usr/lib/nanos6-main-wrapper.o

BIN=libbug.so main

#CFLAGS+=-D_MCC -D_MERCURIUM -D_OMPSS_2=1 -include nanos6.h
#OBJ+=/usr/lib/nanos6-main-wrapper.o
#obj+=/usr/lib/nanos6-main-wrapper.o
#LDFLAGS+=-L/usr/lib -Wl,-z,lazy -Xlinker -rpath -Xlinker /usr/lib
#LDLIBS+=-lnanos6 -ldl /usr/lib/nanos6-main-wrapper.o
#GEN+=$(src)

all: $(BIN)

libbug.so: libbug.mcc.o
	$(CC) -shared $(CFLAGS) -lnanos6 $^ -o $@

%.mcc.c: %.c
	mcc $(MCC_CFLAGS) -y $^ -o $@

%.mcc.o: %.mcc.c
	$(COMPILE.c) $(OUTPUT_OPTION) $<

main: main.o libbug.so
	$(CC) $(CFLAGS) $(LDLIBS) $^ -o $@

clean:
	rm -f $(BIN) *.o
