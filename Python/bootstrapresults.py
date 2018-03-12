# purpose: bootstrap/Ghootstrap w/ replacement subject population to get estimate of uncertainty of correlations
import numpy as np
import matplotlib.pyplot as plt
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.metrics import r2_score
import pickle
import pdb
import scipy
import scipy.io
from sklearn import metrics
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
from usefulfns import getcorr
import pickle
from sklearn.neighbors import KernelDensity
from sklearn import linear_model
from getcorr_allsub import getcorr_vector
# reset random seed just in case
import random
from datetime import datetime
random.seed(datetime.now())

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


# using the original patterns the effect is worse than using the betas
nboot = 1000
bw = 0.1 # set it here for everyone!!
detailthreshold = 2
usebetas_ps = 0 # whether or not to use the betas for the "Y" PS
usebetas_mot = 0 # whether or not to use betas for classifier evidence for MOT
zscoreDV = 0 # if true zscore all DV
zscoreIV = 0 # if you should zscore all classifier values


# LOAD ALL DATA
pickle_in = open("/Volumes/norman/amennen/PythonMot5/evidencebystim_glmclassifier_alpha100_intercept_motion.pickle","rb")
evbystim = pickle.load(pickle_in)

data_dir = '/Volumes/norman/amennen/PythonMot5/'
filename = 'compareExp5.mat'
filepath = os.path.join(data_dir,filename)
d = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
allSep = d['sepbystimD']

# now specify path for betas ps
with open("/Volumes/norman/amennen/PythonMot5/betas_recall_orderedstim.pickle", "rb") as f:  # Python 3: open(..., 'rb')
    betasbystim_RT, betasbystim_OM = pickle.load(f)
motpath = '/Volumes/norman/amennen/motStudy05_transferred/'
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'
savepath = '/Volumes/norman/amennen/PythonMot5/'
targRT = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
lureRT = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
targAcc = np.load('/Volumes/norman/amennen/PythonMot5/targAcc.npy')
lureAcc = np.load('/Volumes/norman/amennen/PythonMot5/lureAcc.npy')
lureAcc_bool = lureAcc==1
targAcc_bool = targAcc==1
lureRT_correctOnly = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
lureRT_correctOnly[~lureAcc_bool] = np.nan
lureRT_incorrectOnly = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
lureRT_incorrectOnly[lureAcc_bool] = np.nan
targRT_correctOnly = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
targRT_correctOnly[~targAcc_bool] = np.nan
targRT_incorrectOnly = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
targRT_incorrectOnly[targAcc_bool] = np.nan
hard_sm = np.load('/Volumes/norman/amennen/wordVec/hard_smF.npy')
simHard = np.load('/Volumes/norman/amennen/wordVec/simHardF.npy')
corDetHard = np.load('/Volumes/norman/amennen/wordVec/corDetHard.npy')
INcorDetHard = np.load('/Volumes/norman/amennen/wordVec/INcorDetHard.npy')
corDetHard = corDetHard.T
INcorDetHard = INcorDetHard.T
hard_sm = hard_sm.T # this is the soft max
simHard = simHard.T # this is just the version of it regualr cosines
all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
nsub= np.int(len(all_sub))
npairs = np.int(nsub/2)
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
nstim = 10
data_dir = '/Volumes/norman/amennen/PythonMot5/RecallPat/'
RTsim = np.zeros((nstim, nsub))
OMITsim = np.zeros((nstim, nsub))
OMITsim_diff = np.zeros((nstim, nsub))
RTsim_diff = np.zeros((nstim, nsub))
avgRTsim = np.zeros(nsub)
avgOMITsim = np.zeros(nsub)
avgRTsim_diff = np.zeros(nsub)
avgOMITsim_diff = np.zeros(nsub)
if usebetas_ps:
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
            RTsim[stim, s] = np.corrcoef(r1, r2)[1, 0]
            r1 = R1_OM[:, stim]
            r2 = R2_OM[:, stim]
            OMITsim[stim, s] = np.corrcoef(r1, r2)[1, 0]

        # now look at how different they are to others
        for stim in np.arange(nstim):
            r1 = R1_RT[:, stim]
            cor_dif = np.zeros(nstim - 1)
            otherstim = np.delete(np.arange(nstim), stim)
            for j in np.arange(nstim - 1):
                other = otherstim[j]
                r2 = R2_RT[:, other]
                cor_dif[j] = np.corrcoef(r1, r2)[1, 0]
            RTsim_diff[stim, s] = np.mean(cor_dif)

        for stim in np.arange(nstim):
            r1 = R1_OM[:, stim]
            cor_dif = np.zeros(nstim - 1)
            otherstim = np.delete(np.arange(nstim), stim)
            for j in np.arange(nstim - 1):
                other = otherstim[j]
                r2 = R2_OM[:, other]
                cor_dif[j] = np.corrcoef(r1, r2)[1, 0]
            OMITsim_diff[stim, s] = np.mean(cor_dif)
else:
    for s in np.arange(nsub):
        subj = all_sub[s]
        # 1/31 after sfn: adding z to zscore differently before taking out timepoints
        filename = 'recallPATz%i.mat' % subj
        filepath = os.path.join(data_dir, filename)
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
            r1 = OMIT[stim, :, 0]
            r2 = OMIT[stim, :, 1]
            OMITsim[stim, s] = np.corrcoef(r1, r2)[1, 0]
        # now look at how different they are to others
        for stim in np.arange(nstim):
            r1 = RT[stim, :, 0]
            cor_dif = np.zeros(nstim - 1)
            otherstim = np.delete(np.arange(nstim), stim)
            for j in np.arange(nstim - 1):
                other = otherstim[j]
                r2 = RT[other, :, 1]
                cor_dif[j] = np.corrcoef(r1, r2)[1, 0]
            RTsim_diff[stim, s] = np.mean(cor_dif)

        for stim in np.arange(nstim):
            r1 = OMIT[stim, :, 0]
            cor_dif = np.zeros(nstim - 1)
            otherstim = np.delete(np.arange(nstim), stim)
            for j in np.arange(nstim - 1):
                other = otherstim[j]
                r2 = OMIT[other, :, 1]
                cor_dif[j] = np.corrcoef(r1, r2)[1, 0]
            OMITsim_diff[stim, s] = np.mean(cor_dif)
# specify analyses we'll do and number of bootstrap iterations

# now load activation score differences: recall retrieval
easyAct = np.zeros((nstim,4,2,nsub))
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
# now average over the four TR's to get one number
avg_hardAct_diff = np.mean(hardAct_diff,axis=1)
postAct = hardAct[:,:,1,:]
avg_hardAct_post = np.mean(postAct,axis=1)
windowsize = 0.05
#min = -1
#max = -1*min + windowsize # to go one over
#min=-0.8
#max=0.9
min=-1.5
max = -1*min + windowsize # to go one over
if zscoreIV:
    windowsize = 0.1
    min=-1.5
    max=-1*min
catrange = np.arange(min,max+windowsize,windowsize)
cr2 = np.reshape(catrange,(len(catrange),1))
nwin = catrange.shape[0]
TRmatrix_kde = np.zeros((nstim,nwin,npairs*2))
maxbins = np.zeros((nstim,nsub))
for s in np.arange(nsub):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    if usebetas_mot:
        subvals_z = stats.zscore(evbystim[sub])
    else:
        subvals_z = stats.zscore(allSep[:,:,s])
    for st in np.arange(nstim):
        if usebetas_mot:
            thissep = evbystim[sub][:,st]
        else:
            thissep = allSep[st,:,s]
        if zscoreIV:
            if usebetas_mot:
                thissep = subvals_z[:,st]
            else:
                thissep = subvals_z[st,:]
        x2 = np.reshape(thissep, (len(thissep), 1))
        kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(x2)
        allvals = np.exp(kde.score_samples(cr2))
        maxbins[st,s] = cr2[np.argmax(allvals)]
        TRmatrix_kde[st, :, s] = allvals

FILTERED_RTsim = RTsim
FILTERED_diffhard = diffHard # behavioral ratings difference
FILTERED_TRmatrix_kde = TRmatrix_kde
FILTERED_lureRT_CO = lureRT_correctOnly #lureAcc + targAcc
FILTERED_lureRT = lureRT
FILTERED_targRT = targRT
FILTERED_recallact = avg_hardAct_diff
#FILTERED_lureRT_CO = lureRT_correctOnly
#FILTERED_lureRT = lureRT
FILTERED_simHard = simHard

for s in np.arange(nsub):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    stimkeep = goodRTstim[sub]
    FILTERED_TRmatrix_kde[~stimkeep,:,s] = np.nan
    FILTERED_RTsim[~stimkeep,s] = np.nan
    FILTERED_diffhard[~stimkeep,s] = np.nan
    FILTERED_lureRT_CO[~stimkeep, s] = np.nan
    FILTERED_lureRT[~stimkeep, s] = np.nan
    FILTERED_targRT[~stimkeep, s] = np.nan
    FILTERED_simHard[~stimkeep, s] = np.nan
    FILTERED_recallact[~stimkeep,s] = np.nan

if zscoreDV:
    FILTERED_diffhard = nanzscore(FILTERED_diffhard) # zscore over subjects
    #FILTERED_TRmatrix_kde[~stimkeep,:,s] = np.nan
    FILTERED_RTsim = nanzscore(FILTERED_RTsim)
    FILTERED_lureRT_CO = nanzscore(FILTERED_lureRT_CO)
    FILTERED_lureRT = nanzscore(FILTERED_lureRT)
    FILTERED_simHard = nanzscore(FILTERED_simHard)
    FILTERED_recallact = nanzscore(FILTERED_recallact)

# analysis: Evidence w/ (1) pattern similarity (2) lure RT (3) detail difference (4) word vector similarity
correlations_ps = np.zeros((nboot,nwin)) # each result will be for a specific window the mean for that bootstrap run
correlations_lureRT = np.zeros((nboot,nwin))
correlations_lureRTCO = np.zeros((nboot,nwin))
correlations_detaildiff = np.zeros((nboot,nwin))
correlations_WVsoftmax = np.zeros((nboot,nwin))
correlations_recallact = np.zeros((nboot,nwin))

for b in np.arange(nboot):
    bootsubjects = np.random.randint(0,high=nsub,size=nsub)
    correlations_ps[b, :] = getcorr_vector(FILTERED_RTsim[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
    #correlations_ps[b,:] = np.nanmean(allr,axis=1)
    correlations_lureRT[b, :] = getcorr_vector(FILTERED_lureRT[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
    #correlations_lureRT[b,:] = np.nanmean(allr,axis=1)
    correlations_lureRTCO[b, :] = getcorr_vector(FILTERED_lureRT_CO[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
    #correlations_lureRTCO[b,:] = np.nanmean(allr,axis=1)
    correlations_detaildiff[b, :] = getcorr_vector(FILTERED_diffhard[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
    #correlations_detaildiff[b,:] = np.nanmean(allr,axis=1)
    correlations_WVsoftmax[b, :] = getcorr_vector(FILTERED_simHard[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
    correlations_recallact[b, :] = getcorr_vector(FILTERED_recallact[:, bootsubjects],FILTERED_TRmatrix_kde[:, :, bootsubjects])
    #correlations_WVsoftmax[b,:] = np.nanmean(allr,axis=1)

# now analyze results of bootstrap!!
# find the 95% confidence interval and plot that
print("done!")

# get confidence intervals
ps_mean = np.zeros((nwin))
ps_errL = np.zeros((nwin))
ps_errH = np.zeros((nwin))

lureRT_mean = np.zeros((nwin))
lureRT_errL = np.zeros((nwin))
lureRT_errH = np.zeros((nwin))

lureRTCO_mean = np.zeros((nwin))
lureRTCO_errL = np.zeros((nwin))
lureRTCO_errH = np.zeros((nwin))

detaildiff_mean = np.zeros((nwin))
detaildiff_errL = np.zeros((nwin))
detaildiff_errH = np.zeros((nwin))

WVsoftmax_mean = np.zeros((nwin))
WVsoftmax_errL = np.zeros((nwin))
WVsoftmax_errH = np.zeros((nwin))

recallact_mean = np.zeros((nwin))
recallact_errL = np.zeros((nwin))
recallact_errH = np.zeros((nwin))

for w in np.arange(nwin):
    ps_mean[w],ps_errL[w],ps_errH[w] = mean_confidence_interval(correlations_ps[:,w])
    lureRT_mean[w], lureRT_errL[w], lureRT_errH[w] = mean_confidence_interval(correlations_lureRT[:, w])
    lureRTCO_mean[w], lureRTCO_errL[w],lureRTCO_errH[w],  = mean_confidence_interval(correlations_lureRTCO[:, w])
    detaildiff_mean[w], detaildiff_errL[w], detaildiff_errH[w] = mean_confidence_interval(correlations_detaildiff[:, w])
    WVsoftmax_mean[w], WVsoftmax_errL[w],WVsoftmax_errH[w] = mean_confidence_interval(correlations_WVsoftmax[:, w])
    recallact_mean[w], recallact_errL[w],recallact_errH[w] = mean_confidence_interval(correlations_recallact[:, w])
# now calculate pvalues!


# now plot results!!

fig, ax = plt.subplots(figsize=(7,5))
plt.title('PATTERN SIMILARITY')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, ps_errL, ps_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,ps_mean, color='r')
plt.ylim(-.25,.25)
ax.set_yticks([-.2,-.1,0,.1,.2])


fig, ax = plt.subplots(figsize=(7,5))
plt.title('LURE RT CO')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, lureRTCO_errL, lureRTCO_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,lureRTCO_mean, color='r')
plt.ylim(-.25,.25)
ax.set_yticks([-.2,-.1,0,.1,.2])

fig, ax = plt.subplots(figsize=(7,5))
plt.title('LURE RT ALL')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, lureRT_errL, lureRT_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,lureRT_mean, color='r')
plt.ylim(-.25,.25)
ax.set_yticks([-.2,-.1,0,.1,.2])

fig, ax = plt.subplots(figsize=(7,5))
plt.title('WV softmax')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, WVsoftmax_errL, WVsoftmax_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,WVsoftmax_mean, color='r')
plt.ylim(-.25,.25)
ax.set_yticks([-.2,-.1,0,.1,.2])

fig, ax = plt.subplots(figsize=(7,5))
plt.title('Recall Act')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, recallact_errL, recallact_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,recallact_mean, color='r')
plt.ylim(-.25,.25)
ax.set_yticks([-.2,-.1,0,.1,.2])

fig, ax = plt.subplots(figsize=(7,5))
plt.title('Detail diff')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, detaildiff_errL, detaildiff_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,detaildiff_mean, color='r')
plt.ylim(-.25,.25)
ax.set_yticks([-.2,-.1,0,.1,.2])

# as a check, look at the kde of all subjects evidence
megadatamatrix= np.empty((0,12))
for s in np.arange(len(all_sub)):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    # calculate individual bw for that subject
    allvals = evbystim[sub]
    stimkeep = goodRTstim[sub]
    allvalsbysubj = allvals[:,stimkeep].T
    megadatamatrix = np.concatenate((megadatamatrix,allvalsbysubj),axis=0)
vectorevidence = megadatamatrix.flatten()
kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(vectorevidence[:,np.newaxis])
allvals = np.exp(kde.score_samples(cr2))

plt.figure()
plt.plot(catrange,allvals, 'r')
bw2=0.05
bins = np.arange(-1.5,1.6,bw2)
x = plt.hist(vectorevidence[:,np.newaxis],bins)
nx = x[0]/np.sum(x[0])
plt.plot(bins[0:-1]+(bw2/2),nx, 'b.')
plt.show()
ev_mean = np.mean(vectorevidence)
ev_std = np.std(vectorevidence)
ev_mean + 3*ev_std
plt.xlabel('Evidence')
plt.ylabel('Proportion in that range')
plt.title('Normalized counts: all classifier evidence')

# PLOT DOR RANGE OF VALUES--TAKE THE SUM AND THEN PLOT CORRELATION
b1 = 20
b2 = 25

sum_filt = np.sum(FILTERED_TRmatrix_kde[:,b1:b2,:],axis=1)
vec_ev = sum_filt.flatten()
vec_rt = FILTERED_RTsim.flatten()
nas = np.logical_or(np.isnan(vec_ev), np.isnan(vec_rt))
m, b = np.polyfit(vec_ev[~nas], vec_rt[~nas], 1)
plt.figure()
plt.plot(vec_ev[~nas],vec_rt[~nas], 'k.')
x = np.arange(9)
plt.plot(x,x*m+b, 'r',linewidth=3, alpha=0.5)
ms = np.zeros(nsub)
for s in np.arange(nsub):
    thisev = sum_filt[:,s]
    thisrt = FILTERED_RTsim[:,s]
    nas = np.logical_or(np.isnan(thisev), np.isnan(thisrt))
    x2 = np.arange(9)
    ms[s],bs = np.polyfit(thisev[~nas],thisrt[~nas],1)
    #plt.plot(x2,x2*ms[s]+bs, 'r', alpha=0.25)
plt.show()
m_c,p = scipy.stats.pearsonr(vec_ev[~nas],vec_rt[~nas])
plt.xlabel('Probability density in range')
plt.ylabel('Pattern similarity pre:post')
plt.title('Pattern similarity vs. time spent in range')

plt.plot()
#xh = plt.hist(ms)
#bww=0.01
#bins = np.arange(-.11,.1,bww)
#xvals = plt.hist(ms,bins)
#nx = 100*x[0]/np.sum(x[0])
#plt.plot(bins[0:-1]+(bw2/2),nx, '.')
sns.stripplot(ms,orient="v",color='k')
sns.barplot(ms,orient="v",ci=68)
plt.title('Slopes for individual subjects')

# now what if we do that with word vector similarity? might be too hard to explain
vec_dt = FILTERED_diffhard.flatten()
nas = np.logical_or(np.isnan(vec_ev), np.isnan(vec_dt))
m, b = np.polyfit(vec_ev[~nas], vec_dt[~nas], 1)
plt.figure()
plt.plot(vec_ev[~nas],vec_dt[~nas], 'k.')
x = np.arange(9)
plt.plot(x,x*m+b, 'r',linewidth=3, alpha=0.5)
ms = np.zeros(nsub)
for s in np.arange(nsub):
    thisev = sum_filt[:,s]
    thisrt = FILTERED_diffhard[:,s]
    nas = np.logical_or(np.isnan(thisev), np.isnan(thisrt))
    x2 = np.arange(9)
    ms[s],bs = np.polyfit(thisev[~nas],thisrt[~nas],1)
plt.show()
m_c,p = scipy.stats.pearsonr(vec_ev[~nas],vec_dt[~nas])
plt.xlabel('Probability density in range')
plt.ylabel('Pattern similarity pre:post')
plt.title('Pattern similarity vs. time spent in range')

# what can subtract to get this value?
# can subtract from subject's mean reaction time for everything?
# calculate subject reaction times
tmeans = np.nanmean(FILTERED_targRT,axis=0)
lmeans = np.nanmean(FILTERED_lureRT,axis=0)
allmeans = (tmeans + lmeans)/2
sub_rt = FILTERED_lureRT_CO - allmeans
vec_RT = np.log(FILTERED_lureRT_CO.flatten())
vec_subrt = sub_rt.flatten()

nas = np.logical_or(np.isnan(vec_ev), np.isnan(vec_subrt))
m, b = np.polyfit(vec_ev[~nas], vec_subrt[~nas], 1)
plt.figure()
plt.plot(vec_ev[~nas],vec_subrt[~nas], 'k.')
x = np.arange(9)
plt.plot(x,x*m+b, 'r',linewidth=3, alpha=0.5)
ms = np.zeros(nsub)
for s in np.arange(nsub):
    thisev = sum_filt[:,s]
    thisrt = sub_rt[:,s]
    nas = np.logical_or(np.isnan(thisev), np.isnan(thisrt))
    x2 = np.arange(9)
    ms[s],bs = np.polyfit(thisev[~nas],thisrt[~nas],1)
plt.show()
m_c,p = scipy.stats.pearsonr(vec_ev[~nas],vec_subrt[~nas])
plt.xlabel('Probability density in range')
plt.ylabel('Reaction time for lure image')
plt.title('RT vs. time spent in range')

plt.figure()
sns.stripplot(ms,orient="v",color='k')
sns.barplot(ms,orient="v",ci=68)
plt.title('Slopes for individual subjects')


vec_wv= FILTERED_simHard.flatten()
nas = np.logical_or(np.isnan(vec_ev), np.isnan(vec_wv))
m, b = np.polyfit(vec_ev[~nas], vec_wv[~nas], 1)
plt.figure()
plt.plot(vec_ev[~nas],vec_wv[~nas], '.')
x = np.arange(9)
plt.plot(x,x*m+b, 'k')
plt.show()
m_c,p = scipy.stats.pearsonr(vec_ev[~nas],vec_wv[~nas])
plt.xlabel('Probability density in range')
plt.ylabel('Pattern similarity pre:post')
plt.title('Pattern similarity vs. time spent in range')

# look for relationship with activations
avgsep = np.median(allSep,axis=1)
avg_hardAct_diff = np.median(hardAct_diff,axis=1)
plt.figure()
plt.plot(avgsep,avg_hardAct_diff, '.')
scipy.stats.pearsonr(avgsep.flatten(),avg_hardAct_diff.flatten())
# you get positive correlation: more higher their region is during MOT, the more likely positive increase afterwards?
# could also see instead of a difference just post scores- relationship a lot oworse if jus tlooking at post only

# look at the max bins
avg_hardAct_diff = np.mean(hardAct_diff,axis=1)
plt.figure()
plt.plot(maxbins,avg_hardAct_diff, '.')
scipy.stats.pearsonr(maxbins.flatten(),avg_hardAct_diff.flatten())
# stronger effect with larger bandwith

# similar relationship with pattern similarity?
plt.figure()
plt.plot(maxbins,RTsim, '.')
scipy.stats.pearsonr(maxbins.flatten(),RTsim.flatten())


plt.figure()
x = maxbins.flatten()
y = lureRT_correctOnly.flatten()
plt.plot(x,y, '.')
nas = np.logical_or(np.isnan(x),np.isnan(y))
scipy.stats.pearsonr(x[~nas],y[~nas])

# is it because YC have higher classifier evidence mean?
nTRs = 15
nTR_total = nTRs*3*nstim
allSepVec_RT = allSep[:,:,RT_ind].flatten()
allSepVec_YC = allSep[:,:,YC_ind].flatten()
data = np.concatenate((allSepVec_RT[:,np.newaxis],allSepVec_YC[:,np.newaxis]),axis=0)
AB = np.concatenate((np.zeros((npairs*nTR_total,1)),np.ones((npairs*nTR_total,1))),axis=0)
data2b = np.concatenate((data,AB),axis=1)
df = pd.DataFrame(data2b,columns = ['Evidence','AB'])
fig, ax = sns.plt.subplots(figsize=(10,7))
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "RT"
labels[1] = "YC"
min=-1.5
max=1.5
binw = .1
bins = np.arange(min,max+binw,binw)
#sns.distplot(allSepVec_YC, bins=bins,color=flatui[1], label='YC',norm_hist=False,kde=False,hist_kws={'alpha':0.6})
#sns.distplot(allSepVec_RT, bins=bins,color=flatui[0], label='RT',norm_hist=False,kde=False,hist_kws={'alpha':0.6})
sns.pointplot(all)
# plot mean instead of ideal range
sns.plt.plot([np.median(allSepVec_YC), np.median(allSepVec_YC)],[0,2000], color=flatui[1], linestyle='--', linewidth=3)
sns.plt.plot([np.median(allSepVec_RT), np.median(allSepVec_RT)],[0,2000], color=flatui[0], linestyle='--', linewidth=3)
sns.plt.title('Distribution of Evidence During MOT')
sns.plt.xlabel('Retrieval evidence bin')
range = np.array([0,.17])
scale = range*len(allSepVec_RT)
sns.plt.ylabel('Fraction of TRs in range')
sns.plt.ylim(scale)
sns.plt.xlim(-.8,.8)
labels2 = np.arange(0,0.2,0.05)
scaled_labels = len(allSepVec_RT)*labels2
result2 = [str(x) for x in labels2]
sns.plt.yticks( scaled_labels, result2 )
sns.despine()
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(30)
sns.plt.legend()

# instead maybe go through each trial and say each trial's most common point
maxbins_RT = maxbins[:,RT_ind].flatten()
maxbins_YC = maxbins[:,YC_ind].flatten()
data = np.concatenate((maxbins_RT[:,np.newaxis],maxbins_YC[:,np.newaxis]),axis=0)
AB = np.concatenate((np.zeros((npairs*nstim,1)),np.ones((npairs*nstim,1))),axis=0)
data2b = np.concatenate((data,AB),axis=1)
df = pd.DataFrame(data2b,columns = ['Evidence','AB'])
fig, ax = sns.plt.subplots(figsize=(10,7))
#plt.hist(maxbins_YC)
#plt.hist(maxbins_RT)
sns.distplot(maxbins_YC, bins=bins, color=flatui[1], label='YC',norm_hist=False,kde=False,hist_kws={'alpha':0.6})
sns.distplot(maxbins_RT, bins=bins, color=flatui[0], label='RT',norm_hist=False,kde=False,hist_kws={'alpha':0.6})
sns.plt.plot([np.mean(maxbins_YC), np.mean(maxbins_YC)],[0,2000], color=flatui[1], linestyle='--', linewidth=3)
sns.plt.plot([np.mean(maxbins_RT), np.mean(maxbins_RT)],[0,2000], color=flatui[0], linestyle='--', linewidth=3)
sns.plt.legend()
range = np.array([0,.5])
scale = range*len(maxbins_YC)
sns.plt.ylabel('Fraction of TRs in range')
sns.plt.ylim(scale)
labels2 = np.arange(0,0.5,0.05)
scaled_labels = len(maxbins_YC)*labels2
result2 = [str(x) for x in labels2]
sns.plt.yticks( scaled_labels, result2 )
plt.xlim(-.65,.65)
res = stats.ttest_rel(maxbins_YC,maxbins_RT)
# for the question is YC greater than RT we take the pvalue over 2
#real pvalue
res.pvalue/2 # almost significantly greater!
# in