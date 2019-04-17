#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 15:14:24 2019

@author: routhier
"""

from Bio import SeqIO
import h5py
import numpy as np

import fasta_reader as script


class CustomDictionnary:
    
    def __init__(self, identity, sequence):
        self.id = identity
        self.seq = sequence


def test_faconverter(tmpdir, monkeypatch):

    def mock_seqIO_parse(filename, extension):
        seq_record = (CustomDictionnary('chr1', 'ATCGGCTA') for i in range(0, 1))
        return seq_record

    monkeypatch.setattr(SeqIO, 'parse', mock_seqIO_parse)

    p = tmpdir.mkdir("tmpStart_data")
    d = '/users/invites/routhier/Documents/MultiAnnotation/seq_chr/human'
    print(p)
    script.main(["--directory", d, "--output", str(p)])

    for i in range(1,22):
        fh5 = h5py.File(p.join('chr' + str(i) + '.hdf5'))
        seq = np.array(fh5[fh5.keys()[0]])
        fh5.close()
        assert (seq[:,0] == np.array([1., 2., 4., 3., 3., 4., 2., 1.])).all()
