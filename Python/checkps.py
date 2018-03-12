#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 14:40:20 2018

@author: amennen
Purpose: look at why the BOLD PS and BETA pattern similarity is so different
- see if there's any correlation between things that have high/low PS before and after MOT
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import scipy
import scipy.io
import os,glob
from matplotlib.pyplot import cm
import matplotlib
matplotlib.rcParams.update({'font.size': 22})
import pandas as pd
import itertools
import seaborn as sns
from scipy import stats
from scipy.stats import norm
from math import exp, sqrt
import pickle
from sklearn.neighbors import KernelDensity
from sklearn import linear_model
# reset random seed just in case
import random
from datetime import datetime
random.seed(datetime.now())
flatui = ["#DB5461", "#593C8F"]


def mean_confidence_interval(data, confidence=0.95):
    d_sorted = np.sort(data)
    n = len(data)
    n_in_middle = confidence*n
    first_index = (n-n_in_middle)/2 - 1
    m = np.mean(data)
    low_val = d_sorted[np.int(first_index)]
    high_val = d_sorted[np.int(first_index+n_in_middle)]
    return m, low_val, high_val

def nanzscore(inputdata):
    if len(np.shape(inputdata)) == 2:
        nsub = np.shape(inputdata)[1]
        nstim = np.shape(inputdata)[0]
        zdata = np.zeros((nstim,nsub))
        for s in np.arange(nsub):
            zdata[:,s] = (inputdata[:,s] - np.nanmean(inputdata[:,s]))/np.nanstd(inputdata[:,s])
    return zdata

#%% specify subject numbers
all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
nsub= np.int(len(all_sub))
npairs = np.int(nsub/2)
nSub= np.int(len(all_sub))
npairs = np.int(nSub/2)
RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)
subarray = np.zeros(nSub)
subarray[YC_ind] = 1
nstim = 10
# variables
detailthreshold = 2
zscoreIV = 0
bw = 0.1
#%% load BOLD PS
BOLD_RTsim = np.zeros((nstim, nsub))
BOLD_OMITsim = np.zeros((nstim, nsub))
rdata_dir = '/Volumes/norman/amennen/PythonMot5/RecallPat/'
for s in np.arange(nsub):
    subj = all_sub[s]
    # 1/31 after sfn: adding z to zscore differently before taking out timepoints
    filename = 'recallPATz%i.mat' % subj
    filepath = os.path.join(rdata_dir, filename)
    print(filepath)
    d = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
    OMIT = d['OMITPAT']
    RT = d['RTPAT']
    nstim = RT.shape[0]
    nvox = RT.shape[1]

    for stim in np.arange(nstim):
        r1 = RT[stim, :, 0]
        r2 = RT[stim, :, 1]
        BOLD_RTsim[stim, s] = np.corrcoef(r1, r2)[1, 0]
        r1 = OMIT[stim, :, 0]
        r2 = OMIT[stim, :, 1]
        BOLD_OMITsim[stim, s] = np.corrcoef(r1, r2)[1, 0] 
#%%load beta PS
with open("/Volumes/norman/amennen/PythonMot5/betas_recall_orderedstim.pickle", "rb") as f:  # Python 3: open(..., 'rb')
    betasbystim_RT, betasbystim_OM = pickle.load(f)
BETA_RTsim = np.zeros((nstim,nsub))
BETA_OMITsim = np.zeros((nstim,nsub))
for s in np.arange(nsub):
    subj = all_sub[s]
    sub = "Subject%01d" % subj
    OMIT = betasbystim_OM[sub]
    RT = betasbystim_RT[sub]
    R1_RT = np.mean(RT[0:4, :, :], axis=0)
    R1_OM = np.mean(OMIT[0:4, :, :], axis=0)
    R2_RT = np.mean(RT[4:8, :, :], axis=0)
    R2_OM = np.mean(OMIT[4:8, :, :], axis=0)
    nstim = RT.shape[2]
    nvox = RT.shape[1]

    for stim in np.arange(nstim):
        r1 = R1_RT[:, stim]
        r2 = R2_RT[:, stim]
        BETA_RTsim[stim, s] = np.corrcoef(r1, r2)[1, 0]
        r1 = R1_OM[:, stim]
        r2 = R2_OM[:, stim]
        BETA_OMITsim[stim, s] = np.corrcoef(r1, r2)[1, 0]
#%% now compare with BETAS/BOLD
        
RT_BOLD = BOLD_RTsim.flatten()
RT_BETA = BETA_RTsim.flatten()
plt.figure(figsize=(10,7))
plt.plot(RT_BOLD,RT_BETA, '.')
plt.title('RT Stimuli')
plt.xlabel('BOLD')
plt.ylabel('BETA')
plt.xlim(-.6,.6)
OM_BOLD = BOLD_OMITsim.flatten()
OM_BETA = BETA_OMITsim.flatten()
plt.figure(figsize=(10,7))
plt.plot(OM_BOLD,OM_BETA, '.')
plt.title('OMIT Stimuli')
plt.xlabel('BOLD')
plt.ylabel('BETA')
plt.xlim(-.6,.6)
scipy.stats.pearsonr(RT_BOLD,RT_BETA)
scipy.stats.pearsonr(OM_BOLD,OM_BETA)
