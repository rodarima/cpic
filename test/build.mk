module:=test

bin:=

bin+=$(module)/cyclotron
#bin+=$(module)/fft

all-tests:=$(addsuffix .test, $(bin))

src:=$(wildcard $(module)/*.c)
obj:=$(subst .c,.o,$(src))

test_ldlibs:=-lm -lconfig -lfftw3 -lgsl -lgslcblas -lGL -lmgl2
test_cflags:=-g -pthread -Wall

test_cflags+=`pkg-config --cflags glfw3`
test_ldlibs+=`pkg-config --libs glfw3`

$(module)/cyclotron: $(module)/cyclotron.o cpic.a
	$(CC) $(CFLAGS) $(test_cflags) $(LDLIBS) $(test_ldlibs) $^ -o $@

$(module)/fft: $(module)/fft.o cpic.a
	$(CC) $(CFLAGS) $(test_cflags) $(LDLIBS) $(test_ldlibs) $^ -o $@

.PHONY: test

test: $(bin)
	@for bin in $(bin); do $$bin || exit 1; done

# Add to main rules
SRC += $(src)
BIN += $(bin)
