# MultiAnnotation

This project is aimed at predicting the position of a given annotation along the genome of a specie. 

### What are the data needed ?

For every specie the data need to be the chromosome in .fa.gz format with one file per chromosome. The module fa_to_multi_fa.py can be used to separate a fa file with several chromosome into several uni-chromosome file. We need also a csv file named
refAnnotation.csv that contains columns Chr, Start, Stop, Strand (the model will be trained to predict the begining of the 
annotation in the reading direction). All those data need to be stored in a directory named as the specie.

##### How can I load my data ?

- To load the data from one specie use *~/MultiAnnotation$ load_specie.sh specie* with specie being a directory containing the data 
described in the previous section.
- To load the data from several specie use *~/MultiAnnotation$ load_multi _species.sh directory* with directory containing several 
species directory.
- To download the data needed to predict the TSS from UCSC browser use *~/MultiAnnotation$ download_data.sh specie assembly USCS*. Specie is the name of the specie and assembly is the name of its consensus sequence. The sequence need to be available in chromosome and not scaffold (with chromosome name chr1, ...). The refGene need to contains at least 11 columns (with num, name, chr, strand, start, stop, cdsStart, cdsStop, exonStarts, exonEnds). After using this script the data will be ready for training or prediction.
- To download the data needed to predict the TSS from NCBI browser use *~/MultiAnnotation$ download_data.sh specie assembly NCBI*. The data can be either in scaffold or chromosome (in rare case the name of contig are none of them, the script will not work). When the data are in scaffold the low limit of scaffold length is set to 50000 bp otherwise the number of files is too high (this parameter can nonetheless be tuned in analysis_refGene.py). 

### How to train and predict with a model ?

- **training a model:** *~/MultiAnnotation$ python training_session.py -s specie_name -a Annotation --max_chr num* with specie name 
corresponding to the name of the directory that contains the data of the specie we want to train the model on, Annotation (with a
capital first letter) is the annotation, num is the number of the last chromsome that will be in the training + validation sets.

- **training a model on several species:** *~/MultiAnnotation$ python training_session.py -s specie1 ... specieN groupname -a Annotation --max_chr num1 ... numN*. Be carefull to give the name of the group of species we want to train on after the list of species. It will create the corresponding directory in both *Start_data* and *Results_multi*.

- **predicting with a model:** *~/MultiAnnotation$ python prediction_session.py -t training_specie -p predicted_specie -a Annotation --start num --stop num* with the arguments being the name of the specie on which the model was trained, the name of the specie we want to predict on, the annotation, the first and the last chromosome on which we want to predict.
