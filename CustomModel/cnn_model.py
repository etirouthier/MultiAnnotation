#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 16:38:56 2019

@author: routhier
"""

from keras.models import Sequential
from keras.layers import Dropout, Flatten
from keras.layers import Dense, Conv2D, MaxPooling2D
from keras.layers.advanced_activations import LeakyReLU

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