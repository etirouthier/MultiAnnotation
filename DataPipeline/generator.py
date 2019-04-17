#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 14:42:41 2019

@author: routhier
"""

import os
import re
import pandas as pd
import numpy as np
import h5py
import argparse

from MyModuleLibrary.array_modifier import reorganize_random_multi_array


def _parse_arguments(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',
                        '--file',
                        help='''path to the csv file with Start, Stop for
                                every gene on all Chr.
                             ''')
    parser.add_argument('-d',
                        '--directory',
                        help='''path to the directory where the .npy files
                                will be stored.
                             ''')

    return parser.parse_args(args)

def gene_position(refGene_file, directory):
    """
        Create a numpy file for every chromosome containing the inidicative 
        function of the Open reading Frame (without exons).
        
        Args:
        refGene_file: path of the csv file with Start, Stop of every gene for 
                      every Chr.
    """
    refGene_new = pd.read_csv(refGene_file)
    refGene_new = refGene_new.dropna(axis=1)

    size_chr = pd.read_csv(os.path.join(os.path.dirname(refGene_file),
                                        'chr_size.txt'), sep='\t')
    
    for chrname in refGene_new.Chr.unique():
        if re.match('chr\d+$', chrname):

            refGene_parsed = refGene_new[refGene_new.Chr == chrname].copy()

            x = np.zeros((size_chr[size_chr.Chr == chrname].Size.values[0],))

            for start, stop in zip(refGene_parsed.Start.values,
                                   refGene_parsed.Stop.values):
                x[start : stop] = 1.0

            np.save(os.path.join(directory,
                                 'gene_indicator_' + chrname + '.npy'),
                    x)
    
def _nucleotid_arrays(path_to_directory, train, val):
    train_chr = range(train[0], train[1] + 1)
    val_chr = range(val[0], val[1] + 1)

    for i in train_chr:
        path_to_file = os.path.join(path_to_directory, 'chr' + str(i) +
                                    '.hdf5')

        f = h5py.File(path_to_file, 'r')
        nucleotid_ = np.array(f[f.keys()[0]])
        f.close()

        if i == train_chr[0]:
            nucleotid_train = nucleotid_
        else:
            nucleotid_train = np.append(nucleotid_train, nucleotid_)

    for i in val_chr:
        path_to_file = os.path.join(path_to_directory, 'chr' + str(i) +
                                    '.hdf5')

        f = h5py.File(path_to_file, 'r')
        nucleotid_ = np.array(f[f.keys()[0]])
        f.close()

        if i == val_chr[0]:
            nucleotid_val = nucleotid_
        else:
            nucleotid_val = np.append(nucleotid_val, nucleotid_)

    return nucleotid_train, nucleotid_val

def _gene_indicator(path_to_directory, train, val):
    train_chr = range(train[0], train[1] + 1)
    val_chr = range(val[0], val[1] + 1)

    proba_train = np.load(os.path.join(path_to_directory,
                                       'gene_indicator_chr' + \
                                       str(train_chr[0]) + '.npy'))

    for i in train_chr[1:]:
        proba_ = np.load(os.path.join(path_to_directory,
                                       'gene_indicator_chr' + str(i) + '.npy'))
        proba_train = np.append(proba_train, proba_)

    proba_val = np.load(os.path.join(path_to_directory,
                                       'gene_indicator_chr' + \
                                       str(val_chr[0]) + '.npy'))

    for i in val_chr[1:]:
        proba_ = np.load(os.path.join(path_to_directory,
                                       'gene_indicator_chr' + str(i) + '.npy'))
        proba_val = np.append(proba_train, proba_)

    # Preparing the weight (1/frequence of occurrence)
    proba = np.append(proba_train, proba_val)
    weights = [1/float(len(np.where(proba == 0.)[0])), 
               1/float(len(np.where(proba == 1.)[0]))]

    weights_train = np.zeros(proba_train.shape)
    weights_train[proba_train == 0] = weights[0]
    weights_val = np.zeros(proba_val.shape)
    weights_val[proba_val == 1] = weights[1]

    return proba_train, weights_train, proba_val, weights_val

def generator(path_to_seq, path_to_genes, train, val):
    """
        Creates two keras data generator for the train set and the validation
        set.

        :param path_to_seq: the path to the a directory containing the
        DNA sequence of all chromosomes in .hdf5 format.
        :param path_to_genes : the path to the .npy directory cointaining files
        with the gene indicator function.
        :param train: tuple (first chromosome in train set, last one)
        :param val: tuple (first chromosome in val set, last one)
    """
    nucleotid_train, nucleotid_val = _nucleotid_arrays(path_to_seq,
                                                       train, val)
    proba_train, weights_train, proba_val, weights_val = \
    _gene_indicator(path_to_genes, train, val)

    positions_train = np.arange(0, nucleotid_train.shape[0])
    positions_val = np.arange(0, nucleotid_val.shape[0])
    batch_size = 512
    number_of_set_train = positions_train.shape[0] // batch_size
    number_of_set_val = positions_val.shape[0] // batch_size

    positions_train = positions_train[1500 : - 1501]
    positions_val = positions_val[1500 : - 1501]

    def generator_function(positions, nucleotid, proba, weights):
        """
            The generator which will be used by the keras model to train.
        """
        window = 299
        number_of_set = positions.shape[0] // batch_size
        half_wx = int((window-1)/2.)
        length = int(positions.shape[0] // number_of_set)

        while True:

            # reshuffled the train set after an epoch
            position = reorganize_random_multi_array(positions)

            for num in range(0, number_of_set):
                if window % 2 == 0:
                    raise ValueError("window must be an odd number")

                positions_ = position[num*length : (num + 1) * length]
                X_ = np.zeros((positions_.shape[0], window, 4, 1))

                for i in range(0, positions_.shape[0]):
                    nucleotid_ = nucleotid[positions_[i] - half_wx :
                                           positions_[i] + half_wx + 1]
                    nucleotid_ = nucleotid_.reshape(nucleotid_.shape[0], 1)
                    X_one_hot = (np.arange(4) == nucleotid_[..., None]-1).astype(int)
                    _X_ = X_one_hot.reshape(X_one_hot.shape[0],
                                            X_one_hot.shape[1] * X_one_hot.shape[2], 1)
                    X_[i] = _X_

                y = proba[positions_]
                w = weights[positions_]
                y = y.reshape(y.shape[0], 1)

                yield X_, y, w

    return generator_function(positions_train,
                              nucleotid_train,
                              proba_train,
                              weights_train), \
           number_of_set_train, \
           generator_function(positions_val,
                              nucleotid_val,
                              proba_val,
                              weights_val), \
           number_of_set_val

def main(command_line_arguments=None):
    """
        Create a csv file with the gene indicative function for every chr.
    """
    args = _parse_arguments(command_line_arguments)
    gene_position(args.file, args.directory)
    
if __name__ == '__main__':
    main()
    