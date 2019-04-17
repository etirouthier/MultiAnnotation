#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 11:02:21 2019

@author: routhier
"""

import os
import pandas as pd
import numpy as np

import datagenerator as script

def test_generator(tmpdir, monkeypatch):
    
    def mock_read_csv(file_in, sep):
        table = pd.DataFrame()
        table['Start'] = np.array([500, 8000])
        table['Stop'] = np.array([1500, 7600])
        table['Strand'] = np.array(['+', '-'])
        table['Chr'] = 'chr1'
        return table

    monkeypatch.setattr(pd, 'read_csv', mock_read_csv)

    p = tmpdir.mkdir("tmpStart_data").join("tmpfile.csv")
    d ='/users/invites/routhier/Documents/' + \
                        'Projet_nucleosomes/' + \
                        'Programme/seq_chr_sacCer3/sacCer3'

    # run script
    script.main(["--directory", d, "--file", str(p),
                 "--balance_factor", "4", "--max_chr", "1"])

    local_x0 = np.load(os.path.dirname(str(p)) + '/X0_start.npy')
    local_x1 = np.load(os.path.dirname(str(p)) + '/X1_start.npy')
    assert local_x1.shape == (2, 299) and local_x0.shape == (8, 299)