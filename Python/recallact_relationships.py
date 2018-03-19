#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 12:08:53 2018

Purpose: compare recall activations and memory results
Measurements of memory:
    Neural:
        - Recall activation (post - pre) (or just post)
        - Pattern similarity (post-pre) (raw/GLM)
    Behavioral:
        - Lure RT post only-->(regular/LOG) choosing log
        - Ratings (post - pre)
        - Word vector cosine similarity (cosine/softmax)

Things this script does:
    1. loads all the data
    2. plots histograms between RT and YC evidence distribution
    3. calculates the maximum bin for that stimulus and replots distirbution
    4. calculates correlations & plots relationship between DV and maxbins
    5. NEW: will plot and caculate correlations between different DV's (maybe average by subject first?)
    6. Will make a large matrix to show all the correlations
    

@author: amennen
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
import statsmodels.api as sm # import statsmodels 
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
# variablesll
detailthreshold = 1
zscoreIV = 0
bw = 0.1
#%% LOAD ALL DATA
pickle_in = open("/Volumes/norman/amennen/PythonMot5/evidencebystim_glmclassifier_alpha100_intercept_motion.pickle","rb")
evbystim = pickle.load(pickle_in)

data_dir = '/Volumes/norman/amennen/PythonMot5/'
filename = 'compareExp5.mat'
filepath = os.path.join(data_dir,filename)
d = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
allSep = d['sepbystimD']

# now specify path for betas ps
motpath = '/Volumes/norman/amennen/motStudy05_transferred/'
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'
savepath = '/Volumes/norman/amennen/PythonMot5/'
targRT = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
lureRT = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
targAcc = np.load('/Volumes/norman/amennen/PythonMot5/targAcc.npy')
lureAcc = np.load('/Volumes/norman/amennen/PythonMot5/lureAcc.npy')
lureAcc_bool = lureAcc==1
targAcc_bool = targAcc==1
lureRT_correctOnly = np.log(np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy'))
lureRT_correctOnly[~lureAcc_bool] = np.nan
lureRT_incorrectOnly = np.log(np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy'))
lureRT_incorrectOnly[lureAcc_bool] = np.nan
targRT_correctOnly = np.log(np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy'))
targRT_correctOnly[~targAcc_bool] = np.nan
targRT_incorrectOnly = np.log(np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy'))
targRT_incorrectOnly[targAcc_bool] = np.nan
hard_sm = np.load('/Volumes/norman/amennen/wordVec/hard_smF.npy')
simHard = np.load('/Volumes/norman/amennen/wordVec/simHardF.npy')
corDetHard = np.load('/Volumes/norman/amennen/wordVec/corDetHard.npy')
INcorDetHard = np.load('/Volumes/norman/amennen/wordVec/INcorDetHard.npy')
corDetHard = corDetHard.T
INcorDetHard = INcorDetHard.T
hard_sm = hard_sm.T # this is the soft max
simHard = simHard.T # this is just the version of it regualr cosines

# load PS
with open("/Volumes/norman/amennen/PythonMot5/betas_recall_orderedstim_PHG2.pickle", "rb") as f:  # Python 3: open(..., 'rb')
    betasbystim_RT, betasbystim_OM = pickle.load(f)
BETA_RTsim = np.zeros((nstim, nsub))
BETA_avgRTsim = np.zeros(nsub)
#if usebetas_ps:
for s in np.arange(nsub):
    subj = all_sub[s]
    sub = "Subject%01d" % subj
    RT = betasbystim_RT[sub]
    R1_RT = np.mean(RT[0:4, :, :], axis=0)
    R2_RT = np.mean(RT[4:8, :, :], axis=0)
    nstim = RT.shape[2]
    nvox = RT.shape[1]
    for stim in np.arange(nstim):
        r1 = R1_RT[:, stim]
        r2 = R2_RT[:, stim]
        BETA_RTsim[stim, s] = np.corrcoef(r1, r2)[1, 0]

rdata_dir = '/Volumes/norman/amennen/PythonMot5/RecallPat/'
RTsim = np.zeros((nstim, nsub))
OMITsim = np.zeros((nstim, nsub))
avgRTsim = np.zeros(nsub)
avgOMITsim = np.zeros(nsub)
for s in np.arange(nsub):
    subj = all_sub[s]
    # 1/31 after sfn: adding z to zscore differently before taking out timepoints
    filename = 'recallPATz_PHG2_%i.mat' % subj
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
        RTsim[stim, s] = np.corrcoef(r1, r2)[1, 0]


# load detail ratings
diffEasy = np.zeros((10,npairs*2))
diffHard = np.zeros((10,npairs*2))
postHard = np.zeros((10,npairs*2))
goodRTstim = {}
easyR = {}
hardR = {}
for s in np.arange(npairs*2):
    subjPath = behavioral_data + str(all_sub[s]) + '/'
    subjName = 'subj' + str(s)
    sub = "Subject%01d" % all_sub[s]
    fn = glob.glob(subjPath + 'ratings'+ '.mat')
    d = scipy.io.loadmat(fn[0])
    easyR[subjName] = d['easyScores']
    # find any nans
    nas = np.isnan(easyR[subjName])
    easyR[subjName] = easyR[subjName].astype(np.float16)
    easyR[subjName][nas] = np.nan
    hardR[subjName] = d['hardScores']
    nas = np.isnan(hardR[subjName])
    hardR[subjName] = hardR[subjName].astype(np.float16)
    hardR[subjName][nas] = np.nan
    # int 16 converts the scores to integers, but the problem is that nan's become zeros (when they should stay as nan's)
    diffEasy[:,s] = np.diff(easyR[subjName],axis=0)
    diffHard[:,s] = np.diff(hardR[subjName],axis=0)
    postHard[:,s] = hardR[subjName][1,:]
    goodRTstim[sub] = hardR[subjName][0,:] > detailthreshold

# load activity for recall    
hardAct = np.zeros((nstim,4,2,nsub))
inteldir = '/Volumes/norman/amennen/motStudy05_transferred/datafromintelrt/data/'
for s in np.arange(npairs*2):
    subjPath = inteldir + str(all_sub[s]) + '/'
    subjName = 'subj' + str(s)
    sub = "Subject%01d" % all_sub[s]
    fn = glob.glob(subjPath + 'recallactivations'+ '.mat')
    d = scipy.io.loadmat(fn[0])
    hardAct[:,:,:,s] = d['hard_activation']
avg_hardAct = np.mean(hardAct,axis=0)
hardAct_diff = np.diff(hardAct,axis=2).squeeze()
avg_preAct = np.mean(hardAct[:,:,0,:],axis=1)
avg_postAct = np.mean(hardAct[:,:,1,:],axis=1)
avg_hardAct_diff = np.mean(hardAct_diff,axis=1)

#%% filter by first rating
FILTERED_BOLD_RTsim = RTsim
FILTERED_diffhard = diffHard # behavioral ratings difference
FILTERED_lureRT_CO = lureRT_correctOnly #lureAcc + targAcc
FILTERED_lureRT = lureRT
FILTERED_targRT = targRT
FILTERED_recallact = avg_hardAct_diff
FILTERED_postrecallact = avg_postAct
FILTERED_BETA_RTsim = BETA_RTsim
## CHANGE HERE WHAT KIND OF SIM YOU WANT
FILTERED_wv_sm = hard_sm
FILTERED_wv_cos = simHard
FILTERED_lureRT_CO_Z  = nanzscore(lureRT_correctOnly)
for s in np.arange(nsub):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    stimkeep = goodRTstim[sub]
    FILTERED_BOLD_RTsim[~stimkeep,s] = np.nan
    FILTERED_diffhard[~stimkeep,s] = np.nan
    FILTERED_lureRT_CO[~stimkeep, s] = np.nan
    FILTERED_lureRT[~stimkeep, s] = np.nan
    FILTERED_targRT[~stimkeep, s] = np.nan
    FILTERED_wv_cos[~stimkeep, s] = np.nan
    FILTERED_wv_sm[~stimkeep, s] = np.nan
    FILTERED_recallact[~stimkeep,s] = np.nan
    FILTERED_postrecallact[~stimkeep,s] = np.nan
    FILTERED_BETA_RTsim[~stimkeep,s] = np.nan
#%% calculate maxbins per stimulus per subject
windowsize = 0.05

min=-1.5
max = -1*min + windowsize # to go one over
catrange = np.arange(min,max+windowsize,windowsize)
cr2 = np.reshape(catrange,(len(catrange),1))
nwin = catrange.shape[0]
TRmatrix_kde = np.zeros((nstim,nwin,npairs*2))
maxbins = np.zeros((nstim,nsub))
for s in np.arange(nsub):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    subvals_z = stats.zscore(allSep[:,:,s])
    for st in np.arange(nstim):
        thissep = allSep[st,:,s]
        if zscoreIV:
            thissep = subvals_z[st,:]
        x2 = np.reshape(thissep, (len(thissep), 1))
        kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(x2)
        allvals = np.exp(kde.score_samples(cr2))
        maxbins[st,s] = cr2[np.argmax(allvals)]
# %% general distribution first
        
nTRs = 15
nTR_total = nTRs*3*nstim
allSepVec_RT = allSep[:,:,RT_ind].flatten()
allSepVec_YC = allSep[:,:,YC_ind].flatten()
data = np.concatenate((allSepVec_RT[:,np.newaxis],allSepVec_YC[:,np.newaxis]),axis=0)
AB = np.concatenate((np.zeros((npairs*nTR_total,1)),np.ones((npairs*nTR_total,1))),axis=0)
data2b = np.concatenate((data,AB),axis=1)
df = pd.DataFrame(data2b,columns = ['Evidence','AB'])
fig, ax = plt.subplots(figsize=(10,7))
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "RT"
labels[1] = "YC"
min=-1.5
max=1.5
binw = .1
bins = np.arange(min,max+binw,binw)
# plot mean instead of ideal range
sns.distplot(allSepVec_YC, bins=bins,color=flatui[1], label='YC',norm_hist=False,kde=False,hist_kws={'alpha':0.6})
sns.distplot(allSepVec_RT, bins=bins,color=flatui[0], label='RT',norm_hist=False,kde=False,hist_kws={'alpha':0.6})
plt.plot([np.median(allSepVec_YC), np.median(allSepVec_YC)],[0,2000], color=flatui[1], linestyle='--', linewidth=3)
plt.plot([np.median(allSepVec_RT), np.median(allSepVec_RT)],[0,2000], color=flatui[0], linestyle='--', linewidth=3)
plt.title('Distribution of Evidence During MOT')
plt.xlabel('Retrieval evidence bin')
range = np.array([0,.17])
scale = range*len(allSepVec_RT)
plt.ylabel('Fraction of TRs in range')
plt.ylim(scale)
plt.xlim(-.8,.8)
labels2 = np.arange(0,0.2,0.05)
scaled_labels = len(allSepVec_RT)*labels2
result2 = ["{:.2f}".format(x) for x in labels2]
plt.yticks( scaled_labels, result2 )
sns.despine()
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(30)
plt.legend()


#%% compare between groups
# instead maybe go through each trial and say each trial's most common point
maxbins_RT = maxbins[:,RT_ind].flatten()
maxbins_YC = maxbins[:,YC_ind].flatten()
data = np.concatenate((maxbins_RT[:,np.newaxis],maxbins_YC[:,np.newaxis]),axis=0)
AB = np.concatenate((np.zeros((npairs*nstim,1)),np.ones((npairs*nstim,1))),axis=0)
data2b = np.concatenate((data,AB),axis=1)
df = pd.DataFrame(data2b,columns = ['Evidence','AB'])
fig, ax = plt.subplots(figsize=(10,7))
#plt.hist(maxbins_YC)
#plt.hist(maxbins_RT)
sns.distplot(maxbins_YC, bins=bins, color=flatui[1], label='YC',norm_hist=False,kde=False,hist_kws={'alpha':0.6})
sns.distplot(maxbins_RT, bins=bins, color=flatui[0], label='RT',norm_hist=False,kde=False,hist_kws={'alpha':0.6})
plt.plot([np.mean(maxbins_YC), np.mean(maxbins_YC)],[0,2000], color=flatui[1], linestyle='--', linewidth=3)
plt.plot([np.mean(maxbins_RT), np.mean(maxbins_RT)],[0,2000], color=flatui[0], linestyle='--', linewidth=3)
plt.legend()
range = np.array([0,.5])
scale = range*len(maxbins_YC)
plt.ylabel('Fraction of TRs in range')
plt.ylim(scale)
labels2 = np.arange(0,0.5,0.05)
scaled_labels = len(maxbins_YC)*labels2
# gets rid of probelm too many floating points
result2 = ["{:.2f}".format(x) for x in labels2]
plt.yticks( scaled_labels, result2 )
plt.xlim(-.65,.65)
res = stats.ttest_rel(maxbins_YC,maxbins_RT)
# for the question is YC greater than RT we take the pvalue over 2
#real pvalue
res.pvalue/2 # almost significantly greater!

#%% plot: for stimulus: max bin & (1) recall activation difference
plt.figure()
plt.plot(maxbins,avg_hardAct_diff, '.')
plt.xlabel('Most often classifier bin during RT')
plt.ylabel('Post-Pre Recall Activation')
scipy.stats.pearsonr(maxbins.flatten(),avg_hardAct_diff.flatten())

#%% plot: for stimulus: max bin & (1) recall activation PRE ONLY
plt.figure()
plt.plot(maxbins,avg_preAct, '.')
plt.xlabel('Most often classifier bin during RT')
plt.ylabel('Pre Recall Activation')
scipy.stats.pearsonr(maxbins.flatten(),avg_preAct.flatten())
#%% plot: for stimulus: max bin & (1) recall activation POST ONLY
plt.figure()
plt.plot(maxbins,avg_postAct, '.')
plt.xlabel('Most often classifier bin during RT')
plt.ylabel('Pre Recall Activation')
scipy.stats.pearsonr(maxbins.flatten(),avg_postAct.flatten())
#%% plot: maxbin and pattern similarty
# similar relationship with pattern similarity?
plt.figure()
x = maxbins.flatten()
y = RTsim.flatten()
nas = np.logical_or(np.isnan(x),np.isnan(y))
plt.plot(x,y, '.')
scipy.stats.pearsonr(x[~nas],y[~nas])

#%% lure RT correct only?
plt.figure()
x = maxbins.flatten()
y = lureRT_correctOnly.flatten()
plt.plot(x,y, '.')
nas = np.logical_or(np.isnan(x),np.isnan(y))
plt.show()
scipy.stats.pearsonr(x[~nas],y[~nas])
#%% lure RT correct only?-zscored
plt.figure()
x = maxbins.flatten()
y = FILTERED_lureRT_CO_Z.flatten()
plt.plot(x,y, '.')
nas = np.logical_or(np.isnan(x),np.isnan(y))
plt.show()
scipy.stats.pearsonr(x[~nas],y[~nas])
#%% plot: maxbin and rating difference
plt.figure()
plt.plot(maxbins,diffHard, '.')
x = maxbins.flatten()
y = diffHard.flatten()
nas = np.logical_or(np.isnan(x),np.isnan(y))
scipy.stats.pearsonr(x[~nas],y[~nas])

#%% now: relate recall to different DV's:
# recall act (post-pre) vs. RT sim(post:pre)
x = FILTERED_recallact.flatten()
y = FILTERED_BOLD_RTsim.flatten()
plt.figure()
nas = np.logical_or(np.isnan(x),np.isnan(y))
plt.plot(x,y, '.')
plt.show()
z = scipy.stats.pearsonr(x[~nas],y[~nas])
# take correlation by subject
corrbysubj = np.zeros((nsub))
for s in np.arange(nsub):
    # compute average correlation
    x = FILTERED_BOLD_RTsim[:,s]
    y = FILTERED_recallact[:,s]
    nas = np.logical_or(np.isnan(x),np.isnan(y))
    corrbysubj[s] = scipy.stats.pearsonr(x[~nas],y[~nas])[0]
print('total correlation is %.2f, p val: %.2f' % (z[0],z[1]))
print('average correlation is %.2f +/- %.2f' % (np.mean(corrbysubj),scipy.stats.sem(corrbysubj)))
    
#%% now: recall vs. reaction time  
x = FILTERED_recallact.flatten()
y = FILTERED_lureRT_CO.flatten()
plt.figure()
nas = np.logical_or(np.isnan(x),np.isnan(y))
plt.plot(x,y, '.')
plt.show()
z = scipy.stats.pearsonr(x[~nas],y[~nas])
# take correlation by subject
corrbysubj = np.zeros((nsub))
for s in np.arange(nsub):
    # compute average correlation
    y = FILTERED_lureRT_CO[:,s]
    x = FILTERED_recallact[:,s]
    nas = np.logical_or(np.isnan(x),np.isnan(y))
    corrbysubj[s] = scipy.stats.pearsonr(x[~nas],y[~nas])[0]
print('total correlation is %.2f, p val: %.2f' % (z[0],z[1]))
print('average correlation is %.2f +/- %.2f' % (np.mean(corrbysubj),scipy.stats.sem(corrbysubj)))
#%% now: recall vs. word vector 
x = FILTERED_postrecallact.flatten()
y = FILTERED_simHard.flatten()
plt.figure()
nas = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),y==0)
plt.plot(x[~nas],y[~nas], '.')
plt.show()
z = scipy.stats.pearsonr(x[~nas],y[~nas])
# take correlation by subject
corrbysubj = np.zeros((nsub))
for s in np.arange(nsub):
    # compute average correlation
    y = FILTERED_simHard[:,s]
    x = FILTERED_recallact[:,s]
    nas = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),y==0)
    corrbysubj[s] = scipy.stats.pearsonr(x[~nas],y[~nas])[0]
print('total correlation is %.2f, p val: %.2f' % (z[0],z[1]))
print('average correlation is %.2f +/- %.2f' % (np.mean(corrbysubj),scipy.stats.sem(corrbysubj)))
#%% now: reaction time vs. word vector
x = FILTERED_BOLD_RTsim.flatten()
y = FILTERED_simHard.flatten()
plt.figure()
nas = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),y==0)
plt.plot(x[~nas],y[~nas], '.')
plt.show()
z = scipy.stats.pearsonr(x[~nas],y[~nas])
# take correlation by subject
corrbysubj = np.zeros((nsub))
for s in np.arange(nsub):
    # compute average correlation
    y = FILTERED_simHard[:,s]
    x = FILTERED_BOLD_RTsim[:,s]
    nas = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),y==0)
    corrbysubj[s] = scipy.stats.pearsonr(x[~nas],y[~nas])[0]
print('total correlation is %.2f, p val: %.2f' % (z[0],z[1]))
print('average correlation is %.2f +/- %.2f' % (np.mean(corrbysubj),scipy.stats.sem(corrbysubj)))
#%% now: reaction time vs. word vector
x = FILTERED_lureRT_CO.flatten()
y = FILTERED_BOLD_RTsim.flatten()
plt.figure()
nas = np.logical_or(np.isnan(x),np.isnan(y))
plt.plot(x[~nas],y[~nas], '.')
plt.show()
z = scipy.stats.pearsonr(x[~nas],y[~nas])
# take correlation by subject
corrbysubj = np.zeros((nsub))
for s in np.arange(nsub):
    # compute average correlation
    y = FILTERED_BOLD_RTsim[:,s]
    x = FILTERED_lureRT_CO[:,s]
    nas = np.logical_or(np.isnan(x),np.isnan(y))
    corrbysubj[s] = scipy.stats.pearsonr(x[~nas],y[~nas])[0]
print('total correlation is %.2f, p val: %.2f' % (z[0],z[1]))
print('average correlation is %.2f +/- %.2f' % (np.mean(corrbysubj),scipy.stats.sem(corrbysubj)))
      
#%% now: reaction time vs. word vector
x = FILTERED_lureRT_CO.flatten()
y = FILTERED_simHard.flatten()
plt.figure()
nas = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),y==0)
plt.plot(x[~nas],y[~nas], '.')
plt.show()
z = scipy.stats.pearsonr(x[~nas],y[~nas])
# take correlation by subject
corrbysubj = np.zeros((nsub))
for s in np.arange(nsub):
    # compute average correlation
    y = FILTERED_simHard[:,s]
    x = FILTERED_lureRT_CO[:,s]
    nas = np.logical_or(np.logical_or(np.isnan(x),np.isnan(y)),y==0)
    corrbysubj[s] = scipy.stats.pearsonr(x[~nas],y[~nas])[0]
print('total correlation is %.2f, p val: %.2f' % (z[0],z[1]))
print('average correlation is %.2f +/- %.2f' % (np.mean(corrbysubj),scipy.stats.sem(corrbysubj)))

#%% now: recall activity vs. ratings
x = FILTERED_recallact.flatten()
y = FILTERED_diffhard.flatten()
plt.figure()
nas = np.logical_or(np.isnan(x),np.isnan(y))
plt.plot(x[~nas],y[~nas], '.')
plt.show()
z = scipy.stats.pearsonr(x[~nas],y[~nas])
# take correlation by subject
corrbysubj = np.zeros((nsub))
for s in np.arange(nsub):
    # compute average correlation
    y = FILTERED_recallact[:,s]
    x = FILTERED_diffhard[:,s]
    nas = np.logical_or(np.isnan(x),np.isnan(y))
    corrbysubj[s] = scipy.stats.pearsonr(x[~nas],y[~nas])[0]
print('total correlation is %.2f, p val: %.2f' % (z[0],z[1]))
print('average correlation is %.2f +/- %.2f' % (np.nanmean(corrbysubj),scipy.stats.sem(corrbysubj)))

#%% maybe do multiple linear regression?      
# build data frame--first nothing zscored
nas = np.isnan(FILTERED_lureRT_CO.flatten())
datamatrix = np.concatenate((FILTERED_lureRT_CO.flatten()[~nas,np.newaxis],FILTERED_BOLD_RTsim.flatten()[~nas,np.newaxis],FILTERED_recallact.flatten()[~nas,np.newaxis],FILTERED_simHard.flatten()[~nas,np.newaxis]),axis=1)
maxbins_vector = maxbins.flatten()
df = pd.DataFrame(datamatrix,columns = ['RT','PS', 'recallact', 'wv'])
X = df
X = sm.add_constant(X)
# Note the difference in argument order
model = sm.OLS(maxbins_vector[~nas], X).fit() ## sm.OLS(output, input)
predictions = model.predict(X)

# Print out the statistics
model.summary()

#%% now calculate all correlations
ncat = 7
FIL