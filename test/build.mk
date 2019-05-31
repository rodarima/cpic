module:=test

test_src:=$(wildcard $(module)/*.c)
tests:=$(subst .c,.test,$(test_src))
test_obj:=$(subst .c,.o,$(test_src))

test_ldlibs:=-lm -lconfig
test_cflags:=-g -pthread -Wall

test_cflags+=$(src_cflags)
test_ldlibs+=$(src_ldlibs)

#test_ldlibs+=-lfftw3 -lgsl -lgslcblas -lGL -lmgl2
#test_cflags+=`pkg-config --cflags glfw3`
#test_ldlibs+=`pkg-config --libs glfw3`
#
#test_cflags+=`mpicc --showme:compile`
#test_ldlibs+=`mpicc --showme:link`

%.test: %.o cpic.a
	$(CC) $(CFLAGS) $(test_cflags) $(LDLIBS) $(test_ldlibs) $^ -o $@

.PHONY: test

test: $(tests)
	@for f in $(tests); do $$f || exit 1; done

# Add to main rules
SRC += $(test_src)
BIN += $(tests)
