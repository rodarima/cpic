BIN=tap.pdf slides.pdf
#BIN=thesis.pdf slides.pdf

CHAPTERS=$(wildcard ch/*.tex)

all: $(BIN)

thesis.pdf: thesis.tex $(CHAPTERS)
	latexmk -latexoption=-shell-escape $< -pdf

%.pdf: %.tex
	latexmk -latexoption=-shell-escape $< -pdf
