COMPILERS=gcc icc clang mcc imcc
BIN:=
BIN+=$(addprefix main., $(COMPILERS))
BIN+=$(addprefix task., $(COMPILERS))

all: $(BIN)

%.gcc: %.c
	gcc -Ofast -ffast-math $^ -o $@

%.icc: %.c
	icc -no-ftz $^ -o $@

%.clang: %.c
	clang -Ofast -ffast-math $^ -o $@

%.mcc: %.c
	mcc --ompss-2 -ffast-math $^ -o $@

%.imcc: %.c
	imcc --v --ompss-2 -Xcompiler -Ofast -fp-model=fast $^ -o $@

run: $(BIN)
	@for EXE in $(BIN) ; do \
		echo -en "$$EXE\\t";\
		./$$EXE ; \
	done

clean:
	rm -f $(BIN)
