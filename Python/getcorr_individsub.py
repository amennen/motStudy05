#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 16:28:24 2018
Made to take matrix of subjects and print out each
From last one, changing the way it prints out allr to be (nsub,nwin)
@author: amennen
"""

import numpy as np
import scipy
import numpy
from scipy import stats
from scipy.stats import norm

def getcorr_matrix(subjectMatrix,TRmatrix):
    nwin = np.shape(TRmatrix)[1]
    nsub = np.shape(subjectMatrix)[1]
    allr = np.zeros((nsub,nwin))
    for s in np.arange(nsub):
        thissim = subjectMatrix[:, s]
        for w in np.arange(nwin):
            thisTR = TRmatrix[:, w, s]
            nas = np.logical_or(np.isnan(thissim), np.isnan(thisTR))
            if len(np.argwhere(~nas)[:,0]) > 2:
                allr[s,w], p = scipy.stats.pearsonr(thissim[~nas],thisTR[~nas])
            else: # not enough samples
                allr[s,w] = np.nan
    return allr
# taken from: http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html


