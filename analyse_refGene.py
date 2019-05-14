#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 13:45:01 2019

@author: routhier
"""

import pandas as pd
import argparse
import os
import re
import numpy as np

def txt_to_csv(filename):
    """
        This function is aimed at converting the txt annotation file that we 
        obtain on USCS to a refGene.csv file with the shape used in our 
        program. ('Chr', 'Strand', 'Start', 'Stop' at least)
    """
    df = pd.read_csv(filename, sep='\t')
    df = df.drop(df.columns[11:], axis=1)

    df.columns = ['num', 'name', 'Chr', 'Strand', 'Start', 'Stop',
                  'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds']

    df.to_csv(os.path.join(os.path.dirname(filename),
                           os.path.basename(filename)[:-3] + 'csv'),
              index=False)

def _analysis_report(assembly_report):
    df = pd.read_csv(assembly_report, sep='\t',
                     names=['Sequence-Name', 'Sequence-Role',
                            'Assigned-Molecule',
                            'Assigned-Molecule-Location/Type',
                            'GenBank-Accn', 'Relationship',
                            'RefSeq-Accn', 'Assembly-Unit', 'Sequence-Length',
                            'USCS-style-name'])

    dico = dict()
    for name in df['Sequence-Name'].values:
        if re.search('([Ss]caffold)?([Cc]hr)?(omosome)?.?\d+$', name) \
        and df[df['Sequence-Name'] == name]['Sequence-Length'].values[0] > 120000000:
            num = str(int(re.search('\d+$', name).group()))
            dico[df[df['Sequence-Name'] == name]['RefSeq-Accn'].values[0]] = \
            'chr' + num
    
    if re.match('[Cc]hr(omosome)?',
                df['Assigned-Molecule-Location/Type'].values[0]):
        seq_type = 'chromosome'
    else:
        seq_type = 'scaffold'

    return dico, seq_type

def gff3_to_csv(filename, assembly_report):
    """
        This function is aimed at converting the gff3 annotation file that we 
        obtain on NCBI to a refGene.csv file with the shape used in our 
        program. ('Chr', 'Strand', 'Start', 'Stop' at least)
    """
    dico, seq_type = _analysis_report(assembly_report)

    df = pd.read_csv(filename,
                     sep='\t',
                     names=['name', 'RefSeq', 'type', 'Start',
                             'Stop', '1', 'Strand', '2', '3'])
    df = df.drop(['1', '2', '3'], axis=1)
    df = df[df.type == 'gene']

    chrom = list()
    for name in df.name.values:
        try:
            chrom.append(dico[name])
        except KeyError:
            chrom.append('non_chr')

    chrom = np.array(chrom)
    df['Chr'] = chrom

    df.to_csv(os.path.join(os.path.dirname(filename),
                           'refGene.csv'),
              index=False)
    return seq_type

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file',
                        help='''the refGene.txt file to be converted in csv''')
    parser.add_argument('-a', '--assembly_report',
                        help='''the assembly_report.txt given by NCBI''')
    args = parser.parse_args()

    if not args.assembly_report:
        txt_to_csv(args.file)
    else:
        seq_type = gff3_to_csv(args.file, args.assembly_report)
        print seq_type

if __name__ == '__main__':
    main()
