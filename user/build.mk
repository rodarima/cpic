module:=user

user_src:=$(wildcard $(module)/*.c)
user_bin:=$(subst .c,.bin,$(user_src))
user_obj:=$(subst .c,.o,$(user_src))

user_ldlibs:=-lm -lconfig -lfftw3 -lgsl -lgslcblas -lGL -lmgl2
user_cflags:=-g -pthread -Wall

user_cflags+=`pkg-config --cflags glfw3`
user_ldlibs+=`pkg-config --libs glfw3`

user_cflags+=`mpicc --showme:compile`
user_ldlibs+=`mpicc --showme:link`

%.bin: %.o cpic.a
	$(CC) $(CFLAGS) $(user_cflags) $(LDLIBS) $(user_ldlibs) $^ -o $@

.PHONY: user

# Add to main rules
SRC += $(user_src)
BIN += $(user_bin)
