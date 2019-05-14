#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 15:17:09 2019

@author: routhier
"""

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.metrics import roc_curve, roc_auc_score

import prediction
import csv_to_bed
from DataPipeline.fasta_reader import faconverter


def _parse_arguments(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-t',
                        '--training_specie',
                        help='''The specie on which the training was made''')
    parser.add_argument('-p',
                        '--predicted_specie',
                        help='''The specie on which we want make prediction''')
    parser.add_argument('-a',
                        '--annotation',
                        help='''The annotation we want to predict''')
    parser.add_argument('--start',
                        help='''first chromosome used for prediction''')
    parser.add_argument('--stop',
                        help='''last chromosome used for prediction''')
    return parser.parse_args(args)

def prediction_analysis(filename, annotation_file, num):
    """
        If the file corresponding to the specie on which the prediction are 
        made contains a refAnnotation.csv file then the analysis consist in 
        plotting the heat_map and the mean prediction in Annotation region.
        
        Args:
            filename: the name of the prediction file that we want to analyse
            annotation_file: a csv file with the position of the annotation
            num: str, the chromosome number on which the analysis will be made
    """
    WX = 299
    HALF_WX = WX // 2

    y_pred = pd.read_csv(filename)
    y_pred = y_pred[y_pred.Chr == 'chr' + num].prediction.values
    z = y_pred - np.mean(y_pred)
    z /= np.std(y_pred)

    refGene = pd.read_csv(annotation_file)
    refGene_parsed_ = refGene[(refGene.Chr == 'chr' + num)]
    refGene_parsed_ = refGene_parsed_.drop_duplicates(subset=['Start', 'Stop'],
                                                      keep='last')

    refGene_parsed_start = refGene_parsed_[(refGene_parsed_.Strand == '+')]
    refGene_parsed_stop = refGene_parsed_[(refGene_parsed_.Strand == '-')]

    start = refGene_parsed_start['Start'].values - HALF_WX
    stop = refGene_parsed_stop['Stop'].values - HALF_WX

    positions = np.append(start, stop)
    matrix = np.array([z[pos - 5000 : pos + 5000] for pos in positions])
    
    y_true = np.zeros(y_pred.shape)
    for pos in positions:
        y_true[pos - HALF_WX : pos + HALF_WX] = 1
        
    roc_points = roc_curve(y_true, y_pred)
    auroc = roc_auc_score(y_true, y_pred)
    
    mean_matrix = np.mean(matrix, axis=0)
    lambda_score = np.sum(mean_matrix[mean_matrix > 2.0]) / \
                   float(len(mean_matrix[mean_matrix > 2]))

    fig = plt.figure(figsize=(7,12))

    ax1 = fig.add_axes([0.1,0.42,0.8,0.1])
    ax2 = fig.add_axes([0.13,0.55,0.92,0.3])

    ax2.tick_params(labelbottom=False,labelleft=True)

    ax2.imshow(matrix, aspect='auto')
    ax2.set_title('Heat_map of the prediction near TSS')
    ax1.plot(range(-5000, 5000), np.mean(matrix, axis=0))
    ax1.text(2000, 5, 'Score : {}'.format(round(100*lambda_score) / 100.),
             fontsize=15)
    ax1.set_title('Averaged zscore of prediction around TSS')

    pl = ax2.pcolormesh(matrix)
    fig.colorbar(pl, ax=ax2)

    ax = fig.add_axes([0.13, 0.1, 0.8, 0.25])

    ax.plot(roc_points[0], roc_points[1])
    ax.set(xlabel='False Positive Rate', ylabel='True Positive Rate',
           title='Receiver Operating Characteristics')
    ax.text(0.4, 0.4,
            'AUROC : {}'.format(round(100*auroc) / 100.),
            fontsize=15)

    plt.savefig(os.path.join(os.path.dirname(filename),
                             'analysis' + os.path.basename(filename)[6 : -4] + \
                             '.png'))
    
def main(command_line_arguments=None):
    """
        Predict the desired annotation on the desired chrosomomes of the specie.
        The model was trained on a specie that is specified.
    """
    args = _parse_arguments(command_line_arguments)
    dna_directory = os.path.join(os.path.dirname(__file__),
                                 'seq_chr', args.predicted_specie)
    annotation_dir = os.path.join(os.path.dirname(__file__),
                                  'Start_data',
                                  args.predicted_specie)
    path_to_model = os.path.join(os.path.dirname(__file__), 'Results_multi',
                                 args.training_specie,
                                 'weights_CNN_' + args.annotation + '_' + \
                                 args.training_specie + '.hdf5' )
    results_dir = os.path.join(os.path.dirname(__file__), 'Results_multi',
                               args.predicted_specie)

    # Verification of the existence of the data on the chromosome we want to
    # predict (i.e hdf5 file for every chromosome). If not we convert the .fa
    for num in range(int(args.start), int(args.stop) + 1):
        if not os.path.exists(os.path.join(dna_directory,
                                      'chr' + str(num) + '.hdf5')):
            faconverter(os.path.join(dna_directory,
                                     'chr' + str(num) + '.fa.gz'),
                        os.path.join(dna_directory,
                                     'chr' + str(num) + '.hdf5'))
    
    # Verification of the existence of the directory to store the results
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
        
    prediction.main(['--directory', dna_directory,
                     '--weights', path_to_model,
                     '--output_dir', results_dir,
                     '--start', args.start,
                     '--stop', args.stop])

    csv_to_bed.main(['--file', os.path.join(results_dir,
                                            'y_pred_CNN_' + args.annotation + \
                                            '_' + args.training_specie + \
                                            '_predict_' + args.predicted_specie + \
                                            '.csv')])

    if os.path.exists(os.path.join(annotation_dir, 'ref' + \
                                   args.annotation + '.csv')):
        prediction_analysis(os.path.join(results_dir,
                                         'y_pred_CNN_' + args.annotation + \
                                         '_' + args.training_specie + \
                                         '_predict_' + args.predicted_specie + \
                                         '.csv'),
                            os.path.join(annotation_dir,
                                         'ref' + args.annotation + '.csv'),
                            args.start)
    
    os.remove(os.path.join(results_dir, 'y_pred_CNN_' + args.annotation + \
                                         '_' + args.training_specie + \
                                         '_predict_' + args.predicted_specie + \
                                         '.csv'))    

if __name__ == '__main__':
    main()
