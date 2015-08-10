#!/bin/bash
# A script to compile the PhD Thesis - Krishna Kumar 
# Distributed under GPLv2.0 License

compile="compile";
clean="clean";

if test -z "$2"
then
if [ $1 = $clean ]; then
	echo "Cleaning please wait ..."
	rm -f *~
	rm -rf *.aux
	rm -rf *.bbl
	rm -rf *.blg
	rm -rf *.d
	rm -rf *.fls
	rm -rf *.ilg
	rm -rf *.ind
	rm -rf *.toc*
	rm -rf *.lot*
	rm -rf *.lof*
	rm -rf *.log
	rm -rf *.idx
	rm -rf *.out*
	rm -rf *.nlo
	rm -rf *.nls
	rm -rf $filename.pdf
	rm -rf $filename.ps
	rm -rf $filename.dvi
	rm -rf *#*

	rm -rf Abstract/*.aux
	rm -rf Acknowledgement/*.aux
	rm -rf Declaration/*.aux
	rm -rf Dedication/*.aux
	rm -rf Appendix1/*.aux
	rm -rf Appendix2/*.aux

	rm -rf Chapter1/*.aux
	rm -rf Chapter2/*.aux
	rm -rf Chapter3/*.aux
	rm -rf Chapter4/*.aux
	rm -rf Chapter5/*.aux
	rm -rf Chapter6/*.aux

	rm -rf *.backup
	rm -rf Abstract/*.backup
	rm -rf Acknowledgement/*.backup
	rm -rf Declaration/*.backup
	rm -rf Dedication/*.backup
	rm -rf Appendix1/*.backup
	rm -rf Appendix2/*.backup

	rm -rf Chapter1/*.backup
	rm -rf Chapter2/*.backup
	rm -rf Chapter3/*.backup
	rm -rf Chapter4/*.backup
	rm -rf Chapter5/*.backup
	rm -rf Chapter6/*.backup

	echo "Cleaning complete!"
	exit
else
	echo "Shell script for compiling the PhD Thesis"
	echo "Usage: sh ./compile-thesis.sh [OPTIONS] [filename]"
	echo "[option]  compile: Compiles the PhD Thesis"
	echo "[option]  clean: removes temporary files no filename required"
	exit
fi
fi

filename=$2;

if [ $1 = $clean ]; then
	echo "Cleaning please wait ..."
	rm -f *~
	rm -rf *.aux
	rm -rf *.bbl
	rm -rf *.blg
	rm -rf *.d
	rm -rf *.fls
	rm -rf *.ilg
	rm -rf *.ind
	rm -rf *.toc*
	rm -rf *.lot*
	rm -rf *.lof*
	rm -rf *.log
	rm -rf *.idx
	rm -rf *.out*
	rm -rf *.nlo
	rm -rf *.nls
	rm -rf $filename.pdf
	rm -rf $filename.ps
	rm -rf $filename.dvi
	rm -rf *#* 
	echo "Cleaning complete!"
	exit
elif [ $1 = $compile ]; then
	echo "Compiling your PhD Thesis...please wait...!"
	pdflatex -interaction=nonstopmode $filename.tex
	bibtex $filename.aux 	
	makeindex $filename.aux
	makeindex $filename.idx
	makeindex $filename.nlo -s nomencl.ist -o $filename.nls
	pdflatex -interaction=nonstopmode $filename.tex
	makeindex $filename.nlo -s nomencl.ist -o $filename.nls
	pdflatex -interaction=nonstopmode $filename.tex
	echo "Success!"
	exit
fi


if test -z "$3"
then
	exit
fi
