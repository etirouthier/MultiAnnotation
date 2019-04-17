#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 15:19:55 2019

@author: routhier
"""

import argparse
import os
import numpy as np

from keras.callbacks import ModelCheckpoint, EarlyStopping, TensorBoard

from MyModuleLibrary.mykeras.losses import MCC, correlate
from CustomModel.cnn_model import cnn_model
from DataPipeline.generator import generator

def _parse_arguments(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',
                        '--file',
                        help='''file to register the weights''')
    parser.add_argument('-d',
                        '--directory',
                        help='''directory that contains the inputs of the model
                        ''')
    parser.add_argument('-a',
                        '--annotation',
                        help='''the annotation we want to be able to predict''')
    parser.add_argument('-c', '--continuous',
                        action='store_true',
                        help='''learn the indicative function of gene''')

    return parser.parse_args(args)

def _data_loader(directory, annotation, batch_size):
    path_to_directory = os.path.dirname(os.path.dirname(directory))
    path_to_directory = os.path.join(path_to_directory,
                                     'Start_data',
                                     directory + '/')
    
    x_0 = np.load(path_to_directory + 'X0_start_' + annotation + '.npy')
    x_1 = np.load(path_to_directory + 'X1_start_' + annotation + '.npy')

    x_ = np.append(x_0, x_1, axis=0)
    y_1 = np.ones((x_1.shape[0]))
    y_0 = np.zeros((x_0.shape[0]))

    y = np.append(y_0,  y_1, axis=0)
    y = y.reshape(y.shape[0],1)

    indices = np.arange(x_.shape[0])
    np.random.shuffle(indices)
    position_train = indices[: int(indices.shape[0] * 0.85)]
    position_val = indices[int(indices.shape[0] * 0.85):]

    number_of_set_train = position_train.shape[0] // batch_size
    number_of_set_val = position_val.shape[0] // batch_size

    def _generator_function(position):
        number_of_set = position.shape[0] // batch_size
        length = int(position.shape[0] // number_of_set)
        position_ = position

        while True:
            np.random.shuffle(position_)

            for i in range(number_of_set):
                local_pos = position_[i * length : (i + 1) * length]
                nucleotid = x_[local_pos]
                X_ = (np.arange(nucleotid.max()) == nucleotid[...,None]-1).astype(int)
                images = X_.reshape(X_.shape[0], X_.shape[1], X_.shape[2], 1)
                targets = y[local_pos]

                yield images, targets
    
    return _generator_function(position_train), number_of_set_train, \
           _generator_function(position_val), number_of_set_val

def main(command_line_arguments=None):
    """
        Train a CNN to predict weither a sequence contains the annotation
        or not.
    """
    args = _parse_arguments(args=command_line_arguments)
    num_epochs = 200
    batch_size = 512
    
    if not command_line_arguments:
        # if we use directly this script without using a training session
        path_to_file = os.path.dirname(os.path.dirname(args.file))
        path_to_file = os.path.join(path_to_file, 'Results_multi', args.file)
    else:
        # using this script in a training session (which give the correct path)
        path_to_file = args.file
    
    model = cnn_model()

    checkpointer = ModelCheckpoint(filepath=path_to_file,
                                   monitor='val_loss',
                                   verbose=0,
                                   save_best_only=True,
                                   save_weights_only=False,
                                   mode='min',
                                   period=1)
    early = EarlyStopping(monitor='val_loss',
                          min_delta=0,
                          patience=10,
                          verbose=0,
                          mode='auto')
    tensorboard = TensorBoard(log_dir=os.path.join(os.getcwd(), 'Tensorboard'),
                              update_freq=200)

    if args.continuous:
        model.compile(loss='mae',
              optimizer='adam',
              metrics=[correlate, MCC])

        directory = input('The directory containing the DNA seq in .hdf5')
        filename = input('''The file that contains the indicator function
                         of the gene in .csv ''')
        train = input('A tuple (first chromosome, last chr) in training set')
        val = input('A tuple (first chromosome, last chr) in validation set')
        
        generator_train, number_of_set_train, \
        generator_val, number_of_set_val = generator(directory, filename,
                                                     train, val)
        model.fit_generator(generator=generator_train,
                            steps_per_epoch=500,
                            epochs=num_epochs,
                            validation_data=generator_val,
                            validation_steps=200,
                            callbacks=[checkpointer, early])
    else:
        model.compile(loss='binary_crossentropy',
              optimizer='adam',
              metrics=['accuracy',MCC])
        
        directory = args.directory
        generator_train, number_of_set_train, \
        generator_val, number_of_set_val = _data_loader(directory,
                                                        args.annotation,
                                                        batch_size)
    
        model.fit_generator(generator=generator_train,
                            steps_per_epoch=number_of_set_train,
                            epochs=num_epochs,
                            validation_data=generator_val,
                            validation_steps=number_of_set_val,
                            callbacks=[checkpointer, early, tensorboard])
        
        # remove the Tensorboard summary after the training
        for element in os.listdir(os.path.join(os.getcwd(), 'Tensorboard')):
            os.remove(os.path.join(os.getcwd(), 'Tensorboard', element))

if __name__ == '__main__':
    main()