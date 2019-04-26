#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 13:45:01 2019

@author: routhier
"""

import pandas as pd
import argparse
import os

def txt_to_csv(filename):
    df = pd.read_csv(filename, sep='\t')
    df = df.drop(df.columns[11:], axis=1)
    
    df.columns = ['num', 'name', 'Chr', 'Strand', 'Start', 'Stop',
                  'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']
    
    df.to_csv(os.path.join(os.path.dirname(filename),
                           os.path.basename(filename)[:-3] + 'csv'),
              index=False)
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file',
                        help='''the refGene.txt file to be converted in csv''')
    args = parser.parse_args()

    txt_to_csv(args.file)

if __name__ == '__main__':
    main()