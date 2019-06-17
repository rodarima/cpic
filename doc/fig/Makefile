all: comm-particles.pdf deadlock-particles.pdf

%.pdf: %.dot
	dot -Tpdf -o $@ $^
