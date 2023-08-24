#!/bin/bash
# input format is only FASTA
start=$1
echo " Start prodigal "
prodigal -i $1 -a output_1.faa -f gbk -o output_1.gbk
echo " Finish prodigal"
python3 true_gbk_constructor.py $1 output_1.faa output_1.gbk
echo " Deleting unnessesary files "
rm output_1.faa
rm output_1.gbk
echo " Your file is {$1}.gbk "

