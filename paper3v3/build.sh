#!/usr/bin/env bash
# Build Paper 3 v2 in place (self-contained; this is what Overleaf does).
set -e; cd "$(dirname "$0")"
latexmk -pdf -interaction=nonstopmode -halt-on-error BBCN_paper3v3.tex
latexmk -c
echo "-> BBCN_paper3v3.pdf"
