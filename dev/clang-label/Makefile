COMPILERS=gcc clang
BIN:=
BIN+=$(addprefix main., $(COMPILERS))

all: $(BIN)

%.gcc: %.c
	gcc $^ -o $@

%.clang: %.c
	clang $^ -o $@

clean:
	rm -f $(BIN)
