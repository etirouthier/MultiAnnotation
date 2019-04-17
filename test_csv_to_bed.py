#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 16:57:52 2019

@author: routhier
"""

import pandas as pd
import numpy as np

import csv_to_bed as script

def test_main(tmpdir, monkeypatch):

    def mock_read_csv(filename):
        mock_pred = pd.DataFrame()
        mock_pred['prediction'] = np.ones((1023,))
        mock_pred['Chr'] = 'chr12'
        return mock_pred
    
    monkeypatch.setattr(pd, 'read_csv', mock_read_csv)

    f = tmpdir.mkdir("tmpResults_multi").join('tmpfile.csv')
    
    script.main(["--file", f])
    
    bedfile = pd.read_table(f, sep='\t')
    
    assert bedfile.chrom.values[0] == 'chr12'
    assert bedfile.score.values[0] == 1
    assert len(bedfile) == 10