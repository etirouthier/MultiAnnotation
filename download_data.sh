#!/bin/bash

# This script is aimed at downloading the data from ucsc by giving the name of the specie and the name of the assembly. The data are then loaded where they need to be loaded.

mkdir $1

cd $1
echo downloading the data from $2 ...
wget http://hgdownload.soe.ucsc.edu/goldenPath/$2/bigZips/$2.fa.gz
wget http://hgdownload.soe.ucsc.edu/goldenPath/$2/database/refGene.txt.gz

echo converting the .fa file to multi .fa ...
python ../fa_to_multi_fa.py -f $2.fa.gz
rm -f $2.fa.gz

echo adapting refGene.txt to our standards ...
gunzip refGene.txt.gz
python ../analyse_refGene.py -f refGene.txt
rm -f refGene.txt

cd ..

echo moving the data to their directory ...
load_specie.sh $1
rm -rf $1
