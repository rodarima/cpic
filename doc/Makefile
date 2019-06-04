BIN=thesis.pdf

CHAPTERS=$(wildcard ch/*.tex)

all: $(BIN)

thesis.pdf: thesis.tex $(CHAPTERS)
	latexmk -latexoption=-shell-escape $< -pdf
