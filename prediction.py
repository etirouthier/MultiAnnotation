#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 15:20:32 2019

@author: routhier
"""
import os
import argparse
import numpy as np
import pandas as pd
import h5py

from keras.models import load_model

from MyModuleLibrary.array_modifier import rolling_window
from MyModuleLibrary.mykeras.losses import MCC


class DataLoader:
    """
        Loading the data for the prediction.

        directory: contains the DNA sequence en .fa or .fa.gz format.
        num_chr: num of the chromosome on which the data will be taken.
    """
    def __init__(self, directory, num_chr):
        self.directory = directory
        self.num_chr = num_chr

    def _one_hot_encoder(self, nucleotid):
        res = (np.arange(nucleotid.max()) == nucleotid[..., None]-1).astype(int)
        res = res.reshape(res.shape[0], res.shape[2])
        return res
    
    def _reverse_dna(self, nucleotid):
        np.place(nucleotid, nucleotid == 1., [5])
        np.place(nucleotid, nucleotid == 2., [6])
        np.place(nucleotid, nucleotid == 3., [7])
        np.place(nucleotid, nucleotid == 4., [8])

        np.place(nucleotid, nucleotid == 5., [2])
        np.place(nucleotid, nucleotid == 6., [1])
        np.place(nucleotid, nucleotid == 7., [4])
        np.place(nucleotid, nucleotid == 8., [3])

        reverse = np.flip(nucleotid, axis=0)
        return reverse

    def process(self):
        """
            Load the data from a .fa file and process it with a rolling window
            so that the model can handle it.
        """
        WX = 299

        path_to_directory = os.path.dirname(os.path.dirname(self.directory))
        path_to_directory = os.path.join(path_to_directory,
                                         'seq_chr',
                                         self.directory + '/')

        f = h5py.File(path_to_directory + 'chr' + str(self.num_chr) + '.hdf5')
        seq = np.array(f['data'])
        f.close()
        
        reverse = self._reverse_dna(seq)
        
        seq = self._one_hot_encoder(seq)
        reverse = self._one_hot_encoder(reverse)
        seq_slide = rolling_window(seq, window=(WX, 4),
                                   asteps=None, wsteps=None,
                                   axes=None, toend=True)
        seq = seq_slide.reshape(seq_slide.shape[0], WX, 4, 1)
        
        reverse_slide = rolling_window(reverse, window=(WX, 4),
                                   asteps=None, wsteps=None,
                                   axes=None, toend=True)
        reverse = reverse_slide.reshape(reverse_slide.shape[0], WX, 4, 1)

        return seq, reverse

def _parse_arguments(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-d',
                        '--directory',
                        help='''Directory containing the DNA sequence
                                  chromosome by chromosome in .hdf5
                               ''')
    parser.add_argument('-w',
                        '--weights',
                        help='''file containing the weights of the model''')
    parser.add_argument('-o',
                        '--output_dir',
                        help='''directory to register the prediction, if
                        different from the directory that contains the trained
                        model''')
    parser.add_argument('--start',
                        help='''first chromosome used for prediction''')
    parser.add_argument('--stop',
                        help='''last chromosome used for prediction''')
    parser.add_argument('-r', '--reverse', action='store_true',
                        help='''predict the annotation in both strand''')
    return parser.parse_args(args)

def main(command_line_arguments=None):
    """
        Predict the annotation positions on the asked chromosomes.
    """
    args = _parse_arguments(command_line_arguments)
    directory = os.path.abspath(args.directory)
    
    if not command_line_arguments:
        # if we are using the prediction outside of a prediction session
        path_to_model = os.path.join(os.path.dirname(__file__),
                                     'Results_multi',
                                     os.path.basename(os.path.dirname(args.weights)),
                                     os.path.basename(args.weights))
    else:
        # if we are using in a prediction session that already treat the file 
        # name.
        path_to_model = args.weights

    model = load_model(path_to_model, custom_objects={'MCC':MCC})
    bedfile = pd.DataFrame()

    for num_chr in range(int(args.start), int(args.stop) + 1):
        data = DataLoader(directory, num_chr)
        x_test, x_test_rev = data.process()
        y_pred = model.predict(x_test)

        bedfile_ = pd.DataFrame()
        bedfile_['prediction'] = y_pred[:, 0]
                                 
        bedfile_['Chr'] = 'chr' + str(num_chr)
        
        if args.reverse:
            y_pred_rev = model.predict(x_test_rev)
            bedfile_['prediction_rev'] = np.flip(y_pred_rev[:, 0], axis=0)
        bedfile = bedfile.append(bedfile_)
    
    if args.output_dir:
        bedfile.to_csv(os.path.join(args.output_dir,
                                    'y_pred_' + \
                                    os.path.basename(args.weights)[8:-5] + \
                                    '_predict_' + os.path.basename(directory) + \
                                    '.csv'))
    else:
        bedfile.to_csv(os.path.join(os.path.dirname(path_to_model),
                                    'y_pred_' + \
                                    os.path.basename(args.weights)[8:-5] + \
                                    '_predict_' + os.path.basename(directory) + \
                                    '.csv'))

if __name__ == '__main__':
    main()
    