#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:08:19 2019

@author: routhier
"""
import os
import argparse
import pandas as pd
import h5py

import numpy as np
np.random.seed(1337)  # for reproducibility

from MyModuleLibrary.array_modifier import reorganize_random_multi_array
from MyModuleLibrary.array_modifier import rolling_window

def _convert_array_to_multi(myArray, number_of_lines, number_of_column):
    """
       Convert a numpy array to multi array and if the shape is not correct
       then reshape it (removing starting elements in most of case).
    """
    if len(myArray) > number_of_lines * number_of_column:
        # if the array has not the right shape,
        # then reshape it by removing x starting elements
        resized_array = np.delete(myArray,
                                  range(0, len(myArray) - (number_of_lines * number_of_column)),
                                  0)
        res = np.reshape(resized_array, (number_of_lines, number_of_column))
        
    elif len(myArray) < number_of_lines * number_of_column:
        myArray = rolling_window(myArray, window=number_of_column)
        np.random.shuffle(myArray)
        res = myArray[ : number_of_lines]
    else:
        res = np.reshape(myArray, (number_of_lines, number_of_column))

    return res

def _parse_arguments(args=None):
    """
        Parsing arguments.
        The csv file needs to contains the following columns : Chr, Strand,
        Start and Stop.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',
                        '--directory',
                        help='''Directory containing the DNA sequence
                                  chromosome by chromosome in .hdf5
                               ''')
    parser.add_argument('-f',
                        '--file',
                        help='''csv file containing the annotation positions, 
                        its name should be refAnnotation.csv''')
    parser.add_argument('-n',
                        '--balance_factor',
                        help='''n times more non annotation than annotation
                        data''')
    parser.add_argument('-m',
                        '--max_chr',
                        help='''last chromosome number in training and
                                validation set
                             ''')

    return parser.parse_args(args)

def generating_tss(max_chr, refGene_file, directory):
    """
        Creating an array containing the DNA sequence around every annotation
        (299 bp)

        An csv file containing the start and stop position of every gene for
        every chromosome need to be passed as well as the last chromosome that
        we want to put in the training and validation set.
        In addition the directory containing the DNA sequence needs to be
        passed.
        A numpy file with sequences of 299 bp long will be created in the same
        directory as refAnnotation.csv

        max_chr: number of the last chromosome in the training and validation
                 set.
        refGene_new: csv file with the genes positions on the specie.
        directory: name of the directory containing the DNA sequence 
    """
    path_to_directory = os.path.abspath(directory)
    index = []
    
    # suppressing gaps in chromosome naming (for ex: chr2A will lead to the 
    #lack of chr2 or for scaffold could be difficult to anticipate)
    for i in range(0, int(max_chr) + 1):
        if os.path.exists(os.path.join(path_to_directory,
                                       'chr' + str(i) + '.hdf5')):
            index.append(i)
    
    try:
        refGene_new = pd.read_csv(refGene_file)
    except:
        refGene_new = pd.read_csv(refGene_file, sep='\t')
    refGene_new = refGene_new.dropna(axis=1)

    for i in index:
        print 'generating positively labeled sequences for chr' + str(i)
        f = h5py.File(os.path.join(path_to_directory,'chr' + str(i) + '.hdf5'),
                      'r')
        seq = np.array(f[f.keys()[0]])
        f.close()

        refGene_parsed_ = refGene_new[(refGene_new.Chr == 'chr' + str(i))]
        refGene_parsed_ = refGene_parsed_.drop_duplicates(subset=['Start', 'Stop'],
                                                          keep='last')

        refGene_parsed_start = refGene_parsed_[(refGene_parsed_.Strand == '+')]
        refGene_parsed_stop = refGene_parsed_[(refGene_parsed_.Strand == '-')]

        start = refGene_parsed_start['Start'].values
        stop = refGene_parsed_stop['Stop'].values

        sequ = seq.reshape(seq.shape[0],)
        adn_sq = sequ.astype('int')

        X_slide = np.array([])
        for x in start:
            n = 149
            X_slide = np.append(X_slide, adn_sq[x - n : x + n + 1])
        X_slide_start = X_slide.reshape(X_slide.shape[0] / 299, 299)

        X_slide = np.array([])
        for x in stop:
            n = 149
            X_slide = np.append(X_slide, adn_sq[x - n : x + n + 1])
        X_slide_stop = X_slide.reshape(X_slide.shape[0] / 299, 299)

        np.place(X_slide_stop, X_slide_stop == 1., [5])
        np.place(X_slide_stop, X_slide_stop == 2., [6])
        np.place(X_slide_stop, X_slide_stop == 3., [7])
        np.place(X_slide_stop, X_slide_stop == 4., [8])

        np.place(X_slide_stop, X_slide_stop == 5., [2])
        np.place(X_slide_stop, X_slide_stop == 6., [1])
        np.place(X_slide_stop, X_slide_stop == 7., [4])
        np.place(X_slide_stop, X_slide_stop == 8., [3])

        reverse = np.flip(X_slide_stop, axis=1)
        X1 = np.append(X_slide_start, reverse, axis=0)

        if i == index[0]:
            res = X1
        else:
            res = np.append(res, X1, axis=0)

    return res

def generating_non_tss(max_chr, n, refGene_file, directory):
    """
        Creating an array containing the DNA sequence far from any TSS (299 bp)

        An csv file containing the start and stop position of every gene for
        every chromosome need to be passed as well as the last chromosome that
        we want to put in the training and validation set.
        In addition the directory containing the DNA sequence needs to be
        passed.
        A numpy file with sequences of 299 bp long will be created in the same
        directory as refGene_new.

        max_chr: number of the last chromosome in the training and validation
                 set.
        n: generating n times as much non TSS sequence as TSS sequence.
        refGene_new: csv file with the genes positions on the specie.
        directory: name of the directory containing the DNA sequence
    """
    path_to_directory = os.path.abspath(directory)
    index = []
    
    # suppressing gaps in chromosome naming (for ex: chr2A will lead to the 
    #lack of chr2)
    for i in range(1, int(max_chr) + 1):
        if os.path.exists(os.path.join(path_to_directory,
                                       'chr' + str(i) + '.hdf5')):
            index.append(i)

    try:
        refGene_new = pd.read_csv(refGene_file)
    except:
        refGene_new = pd.read_csv(refGene_file, sep='\t')
    refGene_new = refGene_new.dropna(axis=1)

    for i in index:
        print 'generating negatively labeled sequences for chr' + str(i)
        f = h5py.File(os.path.join(path_to_directory,'chr' + str(i) + '.hdf5'),
                      'r')
        seq = np.array(f[f.keys()[0]])
        f.close()

        refGene_parsed = refGene_new[(refGene_new.Chr == 'chr' + str(i))]
        refGene_parsed_ = refGene_parsed.drop_duplicates(subset=['Start', 'Stop'],
                                                         keep='last')
        refGene_parsed_start = refGene_parsed_[(refGene_parsed_.Strand == '+')]
        refGene_parsed_stop = refGene_parsed_[(refGene_parsed_.Strand == '-')]

        start = refGene_parsed_start['Start'].values
        stop = refGene_parsed_stop['Stop'].values

        positions = np.append(start, stop)

        del_range = 149 * 2

        del_arr_inc = np.array([])
        del_arr_dec = np.array([])

        for num in range(1, del_range + 1):
            del_arr_inc = np.concatenate((del_arr_inc,
                                          [x + num for x in positions]),
                                         axis=0)
            del_arr_dec = np.concatenate((del_arr_dec,
                                          [x - num for x in positions]),
                                         axis=0)

        del_arr = np.concatenate((del_arr_dec, positions), axis=0)
        del_arr = np.concatenate((del_arr, del_arr_inc), axis=0)
        del_arr = del_arr[del_arr >= 0] # Remove Negative

        C = np.delete(seq, del_arr)

        conv_array = _convert_array_to_multi(C, len(positions)*int(n), 299)

        if i == index[0]:
            res = reorganize_random_multi_array(conv_array)
        else:
            res = np.append(res,
                            reorganize_random_multi_array(conv_array),
                            axis=0)

    return res

def main(command_line_arguments=None):
    """
        Generating both annotation and non annotation class.
    """
    args = _parse_arguments(command_line_arguments)
    
    X1_start = generating_tss(args.max_chr,
                              args.file,
                              args.directory)
    X0_start = generating_non_tss(args.max_chr,
                                  args.balance_factor,
                                  args.file,
                                  args.directory)
    annotation = os.path.basename(args.file)[3 : -4]
    
    np.save(os.path.join(os.path.dirname(args.file),
                         'X0_start_' + annotation + '.npy'), X0_start)

    np.save(os.path.join(os.path.dirname(args.file),
                         'X1_start_' + annotation + '.npy'), X1_start)

if __name__ == '__main__':
    main()
