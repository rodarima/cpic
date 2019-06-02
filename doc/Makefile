BIN=thesis.pdf

CHAPTERS=$(wildcard chapter/*.tex)

all: $(BIN)

thesis.pdf: thesis.tex $(CHAPTERS)
	latexmk -latexoption=-shell-escape $< -pdf
