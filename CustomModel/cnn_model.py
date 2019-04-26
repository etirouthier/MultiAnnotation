#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 16:38:56 2019

@author: routhier
"""

from keras.models import Sequential, Model
from keras.layers import Dropout, Flatten, BatchNormalization, Add
from keras.layers import Dense, Conv2D, MaxPooling2D
from keras.layers.advanced_activations import LeakyReLU
from keras import Input


def cnn_model() :
    """
        Create a convolutional model with 3 convolutional layers before a final 
        dense a layer with one node used to make the final prediction.
        
        ..notes: the precision of the prediction does not depend strongly 
        with the architecture.
    """
    WINDOW = 299
    num_classes = 1

    fashion_model = Sequential()
    fashion_model.add(Conv2D(32, 
                             kernel_size=(3,3),
                             activation='relu',
                             input_shape=(WINDOW,4,1),
                             padding='same'))
    fashion_model.add(LeakyReLU(alpha=0.1))
    fashion_model.add(MaxPooling2D((2, 2), padding='same'))
    fashion_model.add(Dropout(0.2))

    fashion_model.add(Conv2D(64,
                             (3, 3),
                             activation='relu',
                             padding='same'))
    fashion_model.add(LeakyReLU(alpha=0.1))
    fashion_model.add(MaxPooling2D(pool_size=(2, 2),padding='same'))
    fashion_model.add(Dropout(0.2))

    fashion_model.add(Conv2D(128,
                             (3, 3),
                             activation='relu',
                             padding='same'))
    fashion_model.add(LeakyReLU(alpha=0.1))                  
    fashion_model.add(MaxPooling2D(pool_size=(2, 2),padding='same'))
    fashion_model.add(Dropout(0.2))

    fashion_model.add(Flatten())

    fashion_model.add(Dropout(0.2))
    fashion_model.add(LeakyReLU(alpha=0.1))                  

    fashion_model.add(Dense(num_classes, activation='sigmoid'))

    return fashion_model 

def resnet_model() :
    """
        Create a ResNet to predict the nucleosome density along the genome.
    """
    WINDOW = 299
    num_classes = 1

    dna_input = Input(shape=(WINDOW, 4, 1), name='dna_sequence')

    x = Conv2D(48, kernel_size=(3,4),
                   activation='relu',
                   padding='valid')(dna_input)
    x = Dropout(0.2)(x)
    x = BatchNormalization()(x)
    
    x = Conv2D(64, kernel_size=(3,1),
                   activation='relu',
                   padding='same')(x)
    x = Dropout(0.2)(x)
    x = MaxPooling2D((2,1),padding='same')(x)
    x = BatchNormalization()(x)
    
    fx = Conv2D(64, kernel_size=(3,1),
                   activation='relu',
                   padding='same')(x)
    fx = Dropout(0.2)(fx)
    fx = BatchNormalization()(fx)
    fx = Conv2D(64, kernel_size=(3,1),
                   activation='relu',
                   padding='same')(fx)
    fx = Dropout(0.2)(fx)
    x = Add()([fx, x])
    x = MaxPooling2D((2,1),padding='same')(fx)
    x = BatchNormalization()(x)
    
    fx = Conv2D(64, kernel_size=(7,1),
                   activation='relu',
                   padding='same')(x)
    fx = BatchNormalization()(fx)
    fx = Dropout(0.2)(fx)
    fx = Conv2D(64, kernel_size=(7,1),
                   activation='relu',
                   padding='same')(fx)
    fx = Dropout(0.2)(fx)
    x = Add()([fx, x])
    x = MaxPooling2D((2,1),padding='same')(fx)
    x = BatchNormalization()(x)
    
    x = Flatten()(x)
    x = Dense(16, activation = 'relu')(x)

    output = Dense(num_classes, activation='relu')(x)
    fashion_model = Model([dna_input], output)

    return fashion_model
