#!/usr/bin/env bash 

idxfiles=`find . -name "*.index*"`
[ -z "$idxfiles" ] && echo "you do not have indexed reads for nanopolish! Don't you want to extract current intensity values?" || rm $idxfiles
 
csvfiles=`find . -name  "*.csv"`
[ -z "$csvfiles" ] && echo "you have not made any predictions! Don't you want to give it a go?"  || rm $csvfiles

pdffiles=`find . -name "*.pdf"`
[ -z "$pdffiles" ] && echo  "you have not plotted prediction results! Don't you want to visualize your results?" || rm $pdffiles

alntsv=`find . -name "*eventalign.tsv.gz"`
[ -z "$alntsv" ] && echo "you have not performed nanopolish eventalign! Should give it a go unless you have no interest in current intensity information!" || rm $alntsv

evndir=`find . -name "*eventalign.*_events.collapsed" -type d`
[ -z "$evndir" ] && echo  "you do not have collapsed current intensity table! Aren't you interested?" || rm -r $evndir 

rm *.+.evn.tbl
