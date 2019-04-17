#!/bin/bash

for specie in `ls $1`
do
	load_specie.sh `dirname $1`/`basename $1`/$specie
done
