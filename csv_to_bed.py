#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:56:34 2019

@author: routhier
"""

import argparse
import os
import pandas as pd
import numpy as np

def _parse_arguments(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file',
                        help='''Name of the file containing the prediction of 
                        the CNN model''')
    return parser.parse_args(args)

def main(command_line_arguments=None):
    """
        Function used to transform the csv file containing the prediction of 
        the model into a .bed file. 
        
        The genome is divided in sequence of 100 bp and the corresponding 
        annotation is the mean prediction score at this posiition.
    """
    args = _parse_arguments(command_line_arguments)
    
    if command_line_arguments:
        # used in a prediction session
        path_to_file = args.file
    else:
        # used frealy in command line
        path_to_file = os.path.join(os.path.dirname(args.file),
                                    'Results_multi',
                                    args.file)
    prediction = pd.read_csv(path_to_file)
    bedfile = pd.DataFrame()
    
    for chrom in np.unique(prediction.Chr.values):
        bedfile_ = pd.DataFrame()
        y_pred = prediction[prediction.Chr == chrom].prediction.values
        
        start = np.arange(0, len(y_pred) - 100, 100)
        stop = np.arange(100, len(y_pred), 100)
        
        score = np.array([np.mean(y_pred[i : i + 100]) for i in start])
        
        bedfile_['chrom'] = np.repeat(chrom, start.shape[0])
        bedfile_['chromStart'] = start
        bedfile_['chromStop'] = stop
        bedfile_['score'] = score
        
        bedfile = bedfile.append(bedfile_)
        
    bedfile.to_csv(path_or_buf=path_to_file[:-3] + 'bed',
                   sep='\t')
        
if __name__ == '__main__':
    main()