#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 10:27:48 2019

@author: routhier
"""

import os
import argparse
import numpy as np

import training
from DataPipeline import fasta_reader, datagenerator


def _parse_arguments(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-s',
                        '--species',
                        nargs='+',
                        help='''the species on which we want to train a model,
                        if several the last one is the name of the group''')
    parser.add_argument('-a',
                        '--annotation',
                        help='''the annotation we want to predict''')
    parser.add_argument('--max_chr',
                        nargs='+',
                        help='''the last chromosome to put in the training and
                        validation sets for every specie''')

    return parser.parse_args(args)

def _prepare_data(annotation, groupename, **species):
    for specie, max_chr in species.items():
        dna_directory = os.path.join(os.path.dirname(__file__), 'seq_chr',
                                     specie)
        annotation_dir = os.path.join(os.path.dirname(__file__),
                                      'Start_data', specie)
        
        if not os.path.exists(os.path.join(dna_directory, 'chr1.hdf5')):
            fasta_reader.main(['--directory', dna_directory,
                               '--output', dna_directory])
    
        if not os.path.exists(os.path.join(annotation_dir,
                                           'X0_start_' + annotation + '.npy')):
            datagenerator.main(['--directory', dna_directory,
                                '--file', os.path.join(annotation_dir,
                                                       'ref' + annotation + '.csv'),
                                '--balance_factor', '100', 
                                '--max_chr', max_chr])
    
    group_dir = os.path.join(os.path.dirname(__file__),
                             'Start_data', groupename)
    
    if not os.path.exists(group_dir):
        os.mkdir(group_dir)
    
    if not os.path.exists(os.path.join(group_dir,
                                       'X0_start_' + annotation + '.npy')):
        for specie in species.keys():
            annotation_dir = os.path.join(os.path.dirname(__file__),
                                      'Start_data', specie)
            x0_ = np.load(os.path.join(annotation_dir,
                                      'X0_start_' + annotation + '.npy'))
            x1_ = np.load(os.path.join(annotation_dir,
                                      'X1_start_' + annotation + '.npy'))
            
            if specie == species.keys()[0]:
                x0 = x0_
                x1 = x1_
            else:
                x0 = np.append(x0, x0_, axis=0)
                x1 = np.append(x1, x1_, axis=0)
        
        np.save(os.path.join(group_dir, 'X0_start_' + annotation + '.npy'), x0)
        np.save(os.path.join(group_dir, 'X1_start_' + annotation + '.npy'), x1)

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
        Several species can be parsed as arguments, in this case the last
        name of the list must be the name of the group of species.
    """
    args = _parse_arguments(command_line_arguments)
    groupname = args.species[-1]
    species = {}
    
    for i in range(len(args.max_chr)):
        species[args.species[i]] = args.max_chr[i]
    
    _prepare_data(args.annotation, groupname, **species)

    annotation_dir = os.path.join(os.path.dirname(__file__),
                                  'Start_data', groupname)

    results_dir = os.path.join(os.path.dirname(__file__),
                               'Results_multi', groupname)

    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    training.main(['--file',
                  os.path.join(results_dir,
                               'weights_CNN_' + args.annotation + '_' + \
                               groupname + '.hdf5'),
                  '--directory', annotation_dir,
                  '--annotation', args.annotation])

if __name__ == '__main__':
    main()
    