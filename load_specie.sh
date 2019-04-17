#!/bin/bash

for specie in `find $1`
do
	if [ -d $specie ]
	then
		mkdir ./Start_data/`basename $specie`
		mkdir ./seq_chr/`basename $specie`
		
		for chrm in `ls $specie | grep chr*`
		do
			cp `dirname $specie`/`basename $specie`/$chrm ./seq_chr/`basename $specie`
		done
		
		for annotation in `ls $specie | grep ref*`
		do
			cp `dirname $specie`/`basename $specie`/$annotation ./Start_data/`basename $specie`
		done
	fi
done
