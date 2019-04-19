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

### How to train and predict with a model ?

- **training a model:** *~/MultiAnnotation$ python training_session.py -s specie_name -a Annotation --max_chr num* with specie name 
corresponding to the name of the directory that contains the data of the speice we want to train the model on, Annotation (with a
capital first letter) is the annotation, num is the number of the last chromsome that will be in the training + validation sets.

- **predicting with a model:** *~/MultiAnnotation$ python prediction_session.py -t training_specie -p predicted_specie -a Annotation --start num --stop num* with the arguments being the name of the specie on which the model was trained, the name of the specie we want to predict on, the annotation, the first and the last chromosome on which we want to predict.
