#!/bin/bash

# This script is aimed at downloading the data from USCS or NCBI by giving the name of the specie and the name of the assembly (plus NCBI or USCS). The data are then loaded where they need to be loaded.

download_USCS()
{
	echo downloading the data from $2 ...
	wget http://hgdownload.soe.ucsc.edu/goldenPath/$2/bigZips/$2.fa.gz &> /dev/null
	wget http://hgdownload.soe.ucsc.edu/goldenPath/$2/database/refGene.txt.gz &> /dev/null

	echo converting the .fa file to multi .fa ...
	python ../fa_to_multi_fa.py -f $2.fa.gz
	rm -f $2.fa.gz

	echo adapting refGene.txt to our standards ...
	gunzip refGene.txt.gz
	python ../analyse_refGene.py -f refGene.txt
	rm -f refGene.txt
}

download_NCBI()
{
	echo downloading the annotation for $1 ...
	wget ftp://ftp.ncbi.nih.gov/genomes/$1/GFF/ref_$2_top_level.gff3.gz &> /dev/null
	refGene=ref_$2_top_level.gff3
	gunzip $refGene.gz
	access_num=`grep "GCF_[0-9]\+\.[0-9]\+$" $refGene -o`
	grep "^#" $refGene -v > 'tmpfile'
	mv -f 'tmpfile' $refGene
	
	echo downloading the assembly report ...
	wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${access_num:4:3}/${access_num:7:3}/${access_num:10:3}/$access_num'_'$2/$access_num'_'$2'_assembly_report.txt' &> /dev/null
	assembly_rep=$access_num'_'$2'_assembly_report.txt'
	grep "^#" $assembly_rep -v > 'tmpfile'
	mv -f 'tmpfile' $assembly_rep
	
	echo creating a refGene.csv that match our standards ...
	seq_type=`python ../analyse_refGene.py -f $refGene -a $assembly_rep`

	echo downloading the DNA sequences in fasta ...
	
	if [ $seq_type = chromosome ]
	then 
		for i in `seq 1 $max_chr`
		do
			echo downloading chromosome $i ...
			if [ $i -lt 10 ]
			then
				wget ftp://ftp.ncbi.nih.gov/genomes/$1/CHR_0$i/*.fa.gz &> /dev/null
				mv *chr$i.fa.gz chr$i.fa.gz
			else
				wget ftp://ftp.ncbi.nih.gov/genomes/$1/CHR_$i/*.fa.gz &> /dev/null
				mv *chr$i.fa.gz chr$i.fa.gz
			fi
		done
	else
		wget ftp://ftp.ncbi.nih.gov/genomes/$1/CHR_Un/*.fa.gz &> /dev/null
		mv *chrUn.fa.gz chrUn.fa.gz
		echo converting fa to multi fa ...
		python ../fa_to_multi_fa.py -f chrUn.fa.gz -a $assembly_rep
		rm -f chrUn.fa.gz
	fi
	
	rm -f $refGene
	rm -f $assembly_rep
}

mkdir $1
cd $1

if [ $3 = "NCBI" ]
then
	download_NCBI $1 $2
elif [ $3 = "USCS" ]
then
	download_USCS $1 $2
fi

cd ..
echo moving the data to their directory...
load_specie.sh $1
rm -rf $1

