module:=test

src:=$(wildcard $(module)/*.c)
tests:=$(subst .c,.test,$(src))
obj:=$(subst .c,.o,$(src))

test_ldlibs:=-lm -lconfig -lfftw3 -lgsl -lgslcblas -lGL -lmgl2
test_cflags:=-g -pthread -Wall

test_cflags+=`pkg-config --cflags glfw3`
test_ldlibs+=`pkg-config --libs glfw3`

%.test: %.o cpic.a
	$(CC) $(CFLAGS) $(test_cflags) $(LDLIBS) $(test_ldlibs) $^ -o $@

.PHONY: test

test: $(tests)
	@for f in $(tests); do $$f || exit 1; done

# Add to main rules
SRC += $(src)
BIN += $(tests)
