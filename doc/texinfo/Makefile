CSS=document.css
FLAGS+=--css-ref=$(CSS)
FLAGS+=--no-split

all: index.html

index.html: manual.texi
	makeinfo $(FLAGS) --html $^ -o $@
	sed -i '/<title>/a<link rel="icon" type="image/png" href="/logo16.png"/>' $@

manual.pdf: manual.texi
	makeinfo $(FLAGS) --pdf $^ -o $@

upload:	index.html
	rsync -va ./ cpic-web:/home/project-web/cpic-doc/htdocs/

.PHONY: html
