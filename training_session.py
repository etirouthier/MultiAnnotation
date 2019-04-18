#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 10:27:48 2019

@author: routhier
"""

import os
import argparse

from MultiAnnotation import training
from MultiAnnotation.DataPipeline import fasta_reader, datagenerator


def _parse_arguments(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s',
                        '--specie',
                        help='''the species on which we want to train a model
                        ''')
    parser.add_argument('-a',
                        '--annotation',
                        help='''the annotation we want to predict''')
    parser.add_argument('--max_chr',
                        help='''the last chromosome to put in the training and
                        validation sets''')

    return parser.parse_args(args)

def main(command_line_arguments=None):
    """
        Make all the manipulation from the row data to a trained CNN model to
        predict a given annotation in a given specie.

        The DNA data should be stored in a directory named as the specie
        in a .fa.gzip format and a .csv file containing the annotation should
        be stored as well in a directory named as the specie with the name
        refAnnotation.csv. This script will convert the .fa data in .hdf5
        and create the training data for the training. The model will then be
        trained to predict weither a sequence contains the annotation or not.
        The data will not be created if they already exist.
        This script is demanding in terms of memory so it should be used in
        local.
    """
    args = _parse_arguments(command_line_arguments)
    dna_directory = os.path.join(os.path.dirname(__file__), 'seq_chr',
                                 args.specie)
    annotation_dir = os.path.join(os.path.dirname(__file__),
                                  'Start_data', args.specie)

    if not os.path.exists(os.path.join(dna_directory, 'chr1.hdf5')):
        fasta_reader.main(['--directory', dna_directory,
                           '--output', dna_directory])

    if not os.path.exists(os.path.join(annotation_dir,
                                       'X0_start_' + args.annotation + '.npy')):
        datagenerator.main(['--directory', dna_directory,
                            '--file', os.path.join(annotation_dir,
                                                   'ref' + args.annotation + '.csv'),
                            '--balance_factor', '100', 
                            '--max_chr', args.max_chr])
    results_dir = os.path.join(os.path.dirname(__file__),
                               'Results_multi', args.specie)

    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    training.main(['--file',
                  os.path.join(results_dir,
                               'weights_CNN_' + args.annotation + '_' + \
                               args.specie + '.hdf5'),
                  '--directory', annotation_dir,
                  '--annotation', args.annotation])

if __name__ == '__main__':
    main()
    