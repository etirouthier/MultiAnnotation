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
        if re.match('([Ss]caffold_?\d+$)', name) or \
        re.match('[Cc]hr(omosome)?_?\d+$', name) or re.match('^\d+$', name) \
        and df[df['Sequence-Name'] == name]['Sequence-Length'].values[0] > 50000:
            num = str(int(re.search('\d+$', name).group()))
            dico[df[df['Sequence-Name'] == name]['RefSeq-Accn'].values[0]] = \
            'chr' + num
    
    if re.match('[Cc]hr(omosome)?',
                df['Assigned-Molecule-Location/Type'].values[0]):
        seq_type = 'chromosome'
    else:
        seq_type = 'scaffold'

    return dico, seq_type

def _invert_strand(dataframe):
    df = dataframe.copy()
    
    array = df.Strand.values
    array = (array == '+').astype(int)
    array = array.astype('|S1')
    array[array == '0'] = '+'
    array[array == '1'] = '-'
    df.Strand = array
    return df

def _change_chr_name(dataframe, dico):
    df = dataframe.copy()
    chrom = list()
    for name in df.name.values:
        try:
            chrom.append(dico[name])
        except KeyError:
            chrom.append('non_chr')

    chrom = np.array(chrom)
    df['Chr'] = chrom
    return df

def gff3_to_csv(filename, assembly_report):
    """
        This function is aimed at converting the gff3 annotation file that we 
        obtain on NCBI to a refTSS.csv and a refTTS.csv files with the shape 
        used in our program. ('Chr', 'Strand', 'Start', 'Stop' at least)
    """
    dico, seq_type = _analysis_report(assembly_report)

    df = pd.read_csv(filename,
                     sep='\t',
                     names=['name', 'RefSeq', 'type', 'Start',
                             'Stop', '1', 'Strand', '2', '3'])
    df = df.drop(['1', '2', '3'], axis=1)
    dfg = df[df.type == 'gene']
    dfe = df[df.type == 'exon']

    dfg = _change_chr_name(dfg, dico)
    dfe = _change_chr_name(dfe, dico)

    dfg.to_csv(os.path.join(os.path.dirname(filename),
                            'refTSS.csv'),
              index=False)
    dfe.to_csv(os.path.join(os.path.dirname(filename),
                            'refESS.csv'),
              index=False)

    dfg = _invert_strand(dfg)
    dfg.to_csv(os.path.join(os.path.dirname(filename),
                           'refTTS.csv'),
              index=False)
    dfe = _invert_strand(dfe)
    dfe.to_csv(os.path.join(os.path.dirname(filename),
                           'refETS.csv'),
              index=False)
    return seq_type, max([int(dico.values()[i][3:]) for i in range(len(dico))])

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
        seq_type, max_chr = gff3_to_csv(args.file, args.assembly_report)
        print seq_type, max_chr

if __name__ == '__main__':
    main()
