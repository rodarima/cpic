module:=src

bin:=cpic cpic.a

src:=$(wildcard $(module)/*.c)
obj:=$(subst .c,.o,$(src))

src_lib:=$(filter-out $(module)/cpic.c,$(src))
obj_lib:=$(subst .c,.o,$(src))

src_ldlibs:=-lm -lconfig -lfftw3 -lgsl -lgslcblas -lGL -lmgl2
src_cflags:=-g -pthread -Wall

src_cflags+=`pkg-config --cflags glfw3`
src_ldlibs+=`pkg-config --libs glfw3`

cpic: $(obj)
	$(CC) $(CFLAGS) $(src_cflags) $(LDFLAGS) $(LDLIBS) $(src_ldlibs) $^ -o $@

cpic.a: $(obj_lib)
	ar rcs $@ $^

# Add to main rules
SRC += $(src)
BIN += $(bin)
