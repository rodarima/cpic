CSS=document.css
OUT=html
FLAGS+=--css-ref=$(CSS)
FLAGS+=--no-split

all: html

html: manual.texi
	makeinfo $(FLAGS) --html $^ -o $(OUT)

.PHONY: html
