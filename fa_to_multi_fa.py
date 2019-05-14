#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:27:07 2019

@author: routhier
"""

import gzip
import re
import argparse
import os

from Bio import SeqIO

from analyse_refGene import _analysis_report

def file_to_multi_files(filenamein, assembly_report=None):
    """
       Takes a .fa file containing all the DNA sequence and creates one fasta
       per chromosome (compressed in .gz format)

       filenamein: the .fa file that need to be converted (or .gz file)
       assembly_report: for data taken from NCBI an assembly report is added.
       If this assembly is passed then it means the data came from NCBI.
    """
    if re.match(r'.*\.fa$', filenamein):
        file_in = open(filenamein, 'r')
    elif re.match(r'.*\.fa\.gz$', filenamein):
        file_in = gzip.open(filenamein, 'r')
    else:
        raise ValueError("file must be a fasta file (or .fa.gz)")

    directory = os.path.dirname(filenamein)
    
    if not assembly_report:
    
        for seq_record in SeqIO.parse(file_in, 'fasta'):
            if re.match(r'chr\d+$', seq_record.id):
                file_out = os.path.join(directory, seq_record.id + '.fa.gz')
    
                with gzip.open(file_out, 'w') as f_out:
                    SeqIO.write(seq_record, f_out, 'fasta')
    else:
        dico, _ = _analysis_report(assembly_report)
        for seq_record in SeqIO.parse(file_in, 'fasta'):
            try:
                name_out = dico[re.search('N._\d{9}\.\d+',
                                          seq_record.id).group()]
                file_out = os.path.join(directory, name_out + '.fa.gz')
                with gzip.open(file_out, 'w') as f_out:
                    SeqIO.write(seq_record, f_out, 'fasta')
            except KeyError:
                pass

def main():
    """
       Parsing arguments and applying file_to_multi_files()
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file',
                        help='''fasta file to be separated in several files''')
    parser.add_argument('-a', '--assembly_report',
                        help='''assembly_report to convert the NCBI name in 
                        our format of chromosome name''')
    args = parser.parse_args()

    file_to_multi_files(args.file, args.assembly_report)

if __name__ == '__main__':
    main()
    