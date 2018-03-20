#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 18:39:20 2018

@author: amennen
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 16:56:11 2018

@author: amennen
"""

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

    6. To gauge how strong relationships are: --will build TR matrix and ask the correlations between -0.5+.5 where is it the most pos (RT) or neg(all others)
    

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

#%% build TR matrix
windowsize = 0.05
min=-1.5
max = -1*min + windowsize # to go one over
catrange = np.arange(min,max+windowsize,windowsize)
cr2 = np.reshape(catrange,(len(catrange),1))
nwin = catrange.shape[0]
TRmatrix_BETA = np.zeros((nstim,nwin,npairs*2))
TRmatrix_BOLD = np.zeros((nstim,nwin,npairs*2))

maxbinsBETA = np.zeros((nstim,nsub))
maxbinsBOLD = np.zeros((nstim,nsub))

for s in np.arange(nsub):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    for st in np.arange(nstim):
        thissepBETA = evbystim[sub][:,st]
        x2 = np.reshape(thissepBETA, (len(thissepBETA), 1))
        kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(x2)
        allvals = np.exp(kde.score_samples(cr2))
        maxbinsBETA[st,s] = cr2[np.argmax(allvals)]
        TRmatrix_BETA[st, :, s] = allvals
        
        thissepBOLD = allSep[st,:,s]
        x2 = np.reshape(thissepBOLD, (len(thissepBOLD), 1))
        kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(x2)
        allvals = np.exp(kde.score_samples(cr2))
        maxbinsBOLD[st,s] = cr2[np.argmax(allvals)]
        TRmatrix_BOLD[st, :, s] = allvals
#%% filter by first rating
FILTERED_BOLD_RTsim = RTsim
FILTERED_diffhard = diffHard # behavioral ratings difference
FILTERED_lureRT_CO = lureRT_correctOnly #lureAcc + targAcc
FILTERED_recallact = avg_hardAct_diff
FILTERED_postrecallact = avg_postAct
FILTERED_BETA_RTsim = BETA_RTsim
FILTERED_wv_sm = hard_sm
FILTERED_wv_cos = simHard
FILTERED_TRmatrix_BETA = TRmatrix_BETA
FILTERED_TRmatrix_BOLD = TRmatrix_BOLD
FILTERED_maxbins_BOLD = maxbinsBOLD
FILTERED_maxbins_BETA = maxbinsBETA
for s in np.arange(nsub):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    stimkeep = goodRTstim[sub]
    FILTERED_BOLD_RTsim[~stimkeep,s] = np.nan
    FILTERED_diffhard[~stimkeep,s] = np.nan
    FILTERED_lureRT_CO[~stimkeep, s] = np.nan
    FILTERED_wv_cos[~stimkeep, s] = np.nan
    FILTERED_wv_sm[~stimkeep, s] = np.nan
    FILTERED_recallact[~stimkeep,s] = np.nan
    FILTERED_postrecallact[~stimkeep,s] = np.nan
    FILTERED_BETA_RTsim[~stimkeep,s] = np.nan
    FILTERED_TRmatrix_BETA[~stimkeep,:,s]  = np.nan
    FILTERED_TRmatrix_BOLD[~stimkeep,:,s]  = np.nan
    FILTERED_maxbins_BOLD[~stimkeep,s] = np.nan
    FILTERED_maxbins_BETA[~stimkeep,s] = np.nan
#%% correlate maxbins with everything

matplotlib.rcParams.update({'font.size': 18})
ncategories = 9
subjcorrelations = np.zeros((ncategories,ncategories,nsub))
for s in np.arange(nsub):
    subjmat = np.zeros((ncategories,nstim))
    subjmat[0,:] = FILTERED_BOLD_RTsim[:,s]
    subjmat[1,:] = FILTERED_BETA_RTsim[:,s]
    subjmat[2,:] = FILTERED_lureRT_CO[:,s]
    subjmat[3,:] = FILTERED_postrecallact[:,s]
    subjmat[4,:] = FILTERED_recallact[:,s]
    subjmat[5,:] = FILTERED_wv_sm[:,s]
    subjmat[6,:] = FILTERED_wv_cos[:,s]
    subjmat[7,:] = FILTERED_maxbins_BOLD[:,s]
    subjmat[8,:] = FILTERED_maxbins_BETA[:,s]
    df = pd.DataFrame(subjmat.T)
    subjcorrelations[:,:,s] = df.corr()
    
mean_corr = np.mean(subjcorrelations,axis=2)
std_corr = np.std(subjcorrelations, axis=2)
plt.figure()
plt.imshow(mean_corr,vmin=-1,vmax=1,cmap='coolwarm')
plt.colorbar()
plt.xticks(np.arange(ncategories),('boldPS', 'betaPS', 'RTCO', 'pra', 'dra', 'wvsm', 'wvcos','MB_BOLD', 'MB_BETA'))
plt.yticks(np.arange(ncategories),('boldPS', 'betaPS', 'RTCO', 'post_rec_act', 'diff_rec_act', 'wvsm', 'wvcos', 'MBBOLD', 'MBBETA'))
plt.title('Mean Correlations over subjects')
plt.show()

#plt.figure()
#plt.imshow(std_corr,vmin=0,vmax=1,cmap='coolwarm')
#plt.colorbar()
#plt.xticks(np.arange(ncategories),('boldPS', 'betaPS', 'RTCO', 'pra', 'dra', 'wvsm', 'wvcos'))
#plt.yticks(np.arange(ncategories),('boldPS', 'betaPS', 'RTCO', 'pra', 'dra', 'wvsm', 'wvcos'))
#plt.show()
z= scipy.stats.ttest_1samp(subjcorrelations,0,axis=2)
plt.figure()
plt.imshow(z.pvalue,cmap='coolwarm',vmin=0,vmax=0.3)
plt.colorbar()
plt.xticks(np.arange(ncategories),('boldPS', 'betaPS', 'RTCO', 'pra', 'dra', 'wvsm', 'wvcos','MB_BOLD', 'MB_BETA'))
plt.yticks(np.arange(ncategories),('boldPS', 'betaPS', 'RTCO', 'post_rec_act', 'diff_rec_act', 'wvsm', 'wvcos', 'MBBOLD', 'MBBETA'))
plt.title('Pvals--all correlations diff from 0')
plt.show()

# not reliable results
plt.figure()
plt.plot(FILTERED_maxbins_BOLD.flatten(),FILTERED_recallact.flatten(), '.')
nas = np.logical_or(np.isnan(FILTERED_recallact.flatten()),np.isnan(maxbinsBOLD.flatten()))
plt.xlabel('Most often classifier bin during RT')
plt.ylabel('Post-Pre Recall Activation')
scipy.stats.pearsonr(FILTERED_maxbins_BOLD.flatten()[~nas],FILTERED_recallact.flatten()[~nas])
# larger correlation when not averaged by individual subjects
allc = np.zeros((nsub))
for s in np.arange(nsub):
    nas = np.logical_or(np.isnan(FILTERED_recallact[:,s]),np.isnan(maxbinsBOLD[:,s]))
    allc[s] = scipy.stats.pearsonr(FILTERED_maxbins_BOLD[~nas,s],FILTERED_recallact[~nas,s])[0]
    

#%% get at consistency of TRmatrix points
from getcorr_individsub import getcorr_matrix
from getcorr_allsub import getcorr_vector

nboot = 1000
nboot=1
if nboot > 1:
    correlations_BOLDps = np.zeros((nboot,nwin))
    correlations_BETAps = np.zeros((nboot,nwin))
    correlations_lureRTCO = np.zeros((nboot,nwin))
    correlations_wvcos = np.zeros((nboot,nwin))
    correlations_WVsm = np.zeros((nboot,nwin))
    correlations_recallact = np.zeros((nboot,nwin)) 
    correlations_postrecallact = np.zeros((nboot,nwin))
    
    for b in np.arange(nboot):
    
        bootsubjects = np.random.randint(0,high=nsub,size=nsub)
        correlations_BOLDps[b, :] = getcorr_vector(FILTERED_BOLD_RTsim[:,bootsubjects],FILTERED_TRmatrix_BOLD[:,:,bootsubjects])
        correlations_BETAps[b, :] = getcorr_vector(FILTERED_BETA_RTsim[:,bootsubjects],FILTERED_TRmatrix_BOLD[:,:,bootsubjects])
        correlations_lureRTCO[b, :] = getcorr_vector(FILTERED_lureRT_CO[:,bootsubjects],FILTERED_TRmatrix_BOLD[:,:,bootsubjects])
        correlations_wvcos[b, :] = getcorr_vector(FILTERED_wv_cos[:,bootsubjects],FILTERED_TRmatrix_BOLD[:,:,bootsubjects])
        correlations_WVsm[b, :] = getcorr_vector(FILTERED_wv_sm[:,bootsubjects],FILTERED_TRmatrix_BOLD[:,:,bootsubjects])
        correlations_recallact[b, :] = getcorr_vector(FILTERED_recallact[:,bootsubjects],FILTERED_TRmatrix_BOLD[:,:,bootsubjects])
        correlations_postrecallact[b, :] = getcorr_vector(FILTERED_postrecallact[:, bootsubjects],FILTERED_TRmatrix_BOLD[:, :, bootsubjects])
    
    b1 = 20
    b2 = 40
    BOLDPSMIN_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_BOLDps,axis=0)[b1:b2])][0]
    BETAPSMIN_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_BETAps,axis=0)[b1:b2])][0]
    LURERTMIN_BOLD = cr2[b1+np.argmax(np.nanmean(correlations_lureRTCO,axis=0)[b1:b2])][0]
    WVCOS_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_wvcos,axis=0)[b1:b2])][0]
    WMSM_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_WVsm,axis=0)[b1:b2])][0]
    RECALLACT_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_recallact,axis=0)[b1:b2])][0]
    POSTRECALL_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_postrecallact,axis=0)[b1:b2])][0]
    plt.figure()
    y = np.array((BOLDPSMIN_BOLD, BETAPSMIN_BOLD, LURERTMIN_BOLD,WVCOS_BOLD,WMSM_BOLD,RECALLACT_BOLD,POSTRECALL_BOLD))
    plt.plot(np.arange(7),y, '.',label='bold data',markersize=20)
    
    for b in np.arange(nboot):
        bootsubjects = np.random.randint(0,high=nsub,size=nsub)
        correlations_BOLDps[b, :] = getcorr_vector(FILTERED_BOLD_RTsim[:,bootsubjects],FILTERED_TRmatrix_BETA[:,:,bootsubjects])
        correlations_BETAps[b, :] = getcorr_vector(FILTERED_BETA_RTsim[:,bootsubjects],FILTERED_TRmatrix_BETA[:,:,bootsubjects])
        correlations_lureRTCO[b, :] = getcorr_vector(FILTERED_lureRT_CO[:,bootsubjects],FILTERED_TRmatrix_BETA[:,:,bootsubjects])
        correlations_wvcos[b, :] = getcorr_vector(FILTERED_wv_cos[:,bootsubjects],FILTERED_TRmatrix_BETA[:,:,bootsubjects])
        correlations_WVsm[b, :] = getcorr_vector(FILTERED_wv_sm[:,bootsubjects],FILTERED_TRmatrix_BETA[:,:,bootsubjects])
        correlations_recallact[b, :] = getcorr_vector(FILTERED_recallact[:,bootsubjects],FILTERED_TRmatrix_BETA[:,:,bootsubjects])
        correlations_postrecallact[b, :] = getcorr_vector(FILTERED_postrecallact[:, bootsubjects],FILTERED_TRmatrix_BETA[:, :, bootsubjects])  
   
    BOLDPSMIN_BETA = cr2[b1+np.argmin(np.nanmean(correlations_BOLDps,axis=0)[b1:b2])][0]
    BETAPSMIN_BETA = cr2[b1+np.argmin(np.nanmean(correlations_BETAps,axis=0)[b1:b2])][0]
    LURERTMIN_BETA = cr2[b1+np.argmax(np.nanmean(correlations_lureRTCO,axis=0)[b1:b2])][0]
    WVCOS_BETA = cr2[b1+np.argmin(np.nanmean(correlations_wvcos,axis=0)[b1:b2])][0]
    WMSM_BETA = cr2[b1+np.argmin(np.nanmean(correlations_WVsm,axis=0)[b1:b2])][0]
    RECALLACT_BETA = cr2[b1+np.argmin(np.nanmean(correlations_recallact,axis=0)[b1:b2])][0]
    POSTRECALL_BETA = cr2[b1+np.argmin(np.nanmean(correlations_postrecallact,axis=0)[b1:b2])][0]
    ybeta = np.array((BOLDPSMIN_BETA, BETAPSMIN_BETA, LURERTMIN_BETA,WVCOS_BETA,WMSM_BETA,RECALLACT_BETA,POSTRECALL_BETA))
    plt.plot(np.arange(7),ybeta, 'r.',label='glm data',markersize=20)
    plt.title('Correlation  minimums')
    plt.ylim(-.6,.6)
    plt.legend()
    
else:
    correlations_BOLDps = getcorr_matrix(FILTERED_BOLD_RTsim,FILTERED_TRmatrix_BOLD)
    correlations_BETAps = getcorr_matrix(FILTERED_BETA_RTsim,FILTERED_TRmatrix_BOLD)
    correlations_lureRTCO = getcorr_matrix(FILTERED_lureRT_CO,FILTERED_TRmatrix_BOLD)
    correlations_wvcos = getcorr_matrix(FILTERED_wv_cos,FILTERED_TRmatrix_BOLD)
    correlations_WVsm = getcorr_matrix(FILTERED_wv_sm,FILTERED_TRmatrix_BOLD)
    correlations_recallact = getcorr_matrix(FILTERED_recallact,FILTERED_TRmatrix_BOLD)  
    correlations_postrecallact = getcorr_matrix(FILTERED_postrecallact,FILTERED_TRmatrix_BOLD)  
    
    # within certain zones, where is the minimum point
    b1 = 20
    b2 = 40
    BOLDPSMIN_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_BOLDps,axis=0)[b1:b2])][0]
    BETAPSMIN_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_BETAps,axis=0)[b1:b2])][0]
    LURERTMIN_BOLD = cr2[b1+np.argmax(np.nanmean(correlations_lureRTCO,axis=0)[b1:b2])][0]
    WVCOS_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_wvcos,axis=0)[b1:b2])][0]
    WMSM_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_WVsm,axis=0)[b1:b2])][0]
    RECALLACT_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_recallact,axis=0)[b1:b2])][0]
    POSTRECALL_BOLD = cr2[b1+np.argmin(np.nanmean(correlations_postrecallact,axis=0)[b1:b2])][0]
    plt.figure()
    y = np.array((BOLDPSMIN_BOLD, BETAPSMIN_BOLD, LURERTMIN_BOLD,WVCOS_BOLD,WMSM_BOLD,RECALLACT_BOLD,POSTRECALL_BOLD))
    plt.plot(np.arange(7),y, '.',label='bold data',markersize=20)
    
    correlations_BOLDps = getcorr_matrix(FILTERED_BOLD_RTsim,FILTERED_TRmatrix_BETA)
    correlations_BETAps = getcorr_matrix(FILTERED_BETA_RTsim,FILTERED_TRmatrix_BETA)
    correlations_lureRTCO = getcorr_matrix(FILTERED_lureRT_CO,FILTERED_TRmatrix_BETA)
    correlations_wvcos = getcorr_matrix(FILTERED_wv_cos,FILTERED_TRmatrix_BETA)
    correlations_WVsm = getcorr_matrix(FILTERED_wv_sm,FILTERED_TRmatrix_BETA)
    correlations_recallact = getcorr_matrix(FILTERED_recallact,FILTERED_TRmatrix_BETA)  
    correlations_postrecallact = getcorr_matrix(FILTERED_postrecallact,FILTERED_TRmatrix_BETA)  
    
    # within certain zones, where is the minimum point
    b1 = 20
    b2 = 40
    BOLDPSMIN_BETA = cr2[b1+np.argmin(np.nanmean(correlations_BOLDps,axis=0)[b1:b2])][0]
    BETAPSMIN_BETA = cr2[b1+np.argmin(np.nanmean(correlations_BETAps,axis=0)[b1:b2])][0]
    LURERTMIN_BETA = cr2[b1+np.argmax(np.nanmean(correlations_lureRTCO,axis=0)[b1:b2])][0]
    WVCOS_BETA = cr2[b1+np.argmin(np.nanmean(correlations_wvcos,axis=0)[b1:b2])][0]
    WMSM_BETA = cr2[b1+np.argmin(np.nanmean(correlations_WVsm,axis=0)[b1:b2])][0]
    RECALLACT_BETA = cr2[b1+np.argmin(np.nanmean(correlations_recallact,axis=0)[b1:b2])][0]
    POSTRECALL_BETA = cr2[b1+np.argmin(np.nanmean(correlations_postrecallact,axis=0)[b1:b2])][0]
    ybeta = np.array((BOLDPSMIN_BETA, BETAPSMIN_BETA, LURERTMIN_BETA,WVCOS_BETA,WMSM_BETA,RECALLACT_BETA,POSTRECALL_BETA))
    plt.plot(np.arange(7),ybeta, 'r.',label='glm data',markersize=20)
    plt.title('Correlation  minimums')
    plt.ylim(-.6,.6)
    plt.legend()