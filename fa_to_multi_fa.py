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

def file_to_multi_files(filenamein):
    """
       Takes a .fa file containing all the DNA sequence and creates one fasta
       per chromosome (compressed in .gz format)

       filenamein: the .fa file that need to be converted (or .gz file)
    """
    if re.match(r'.*\.fa$', filenamein):
        file_in = open(filenamein, 'r')
    elif re.match(r'.*\.fa\.gz$', filenamein):
        file_in = gzip.open(filenamein, 'r')
    else:
        raise ValueError("file must be a fasta file (or .fa.gz)")

    directory = os.path.dirname(filenamein)

    for seq_record in SeqIO.parse(file_in, 'fasta'):
        if re.match(r'chr\d+$', seq_record.id):
            file_out = os.path.join(directory, seq_record.id + '.fa.gz')

            with gzip.open(file_out, 'w') as f_out:
                SeqIO.write(seq_record, f_out, 'fasta')

def main():
    """
       Parsing arguments and applying file_to_multi_files()
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file',
                        help='''fasta file to be separated in several files''')
    args = parser.parse_args()

    file_to_multi_files(args.file)

if __name__ == '__main__':
    main()
    