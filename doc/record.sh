#!/bin/bash

DIR="$(dirname $0)"
PDF="$DIR/thesis.pdf"
LOG="$DIR/working.log"


WORDS=$(pdftotext $PDF - | wc -w)
PERCENT=$(printf "$WORDS/200\n" | calc -p)


printf "$(date +%s)\t$WORDS\t$PERCENT\n" >> "$LOG"

cd "$DIR"

python regression.py
