# purpose: bootstrap/Ghootstrap w/ replacement subject population to get estimate of uncertainty of correlations

# let's just make this script ONLY bootstrapping results

# calculates bootstrap correlations of mot5 based on input parameters (zscore/betas)
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
from usefulfns import getcorr
import pickle
from sklearn.neighbors import KernelDensity
from sklearn import linear_model
from getcorr_allsub import getcorr_vector
from getcorr_individsub import getcorr_matrix
# reset random seed just in case
import random
from datetime import datetime
random.seed(datetime.now())
flatui = ["#DB5461", "#593C8F"]

# using the original patterns the effect is worse than using the betas
nboot = 1000 # put nboot as == 1 to say only run once
bw = 0.1 # set it here for everyone!!
detailthreshold = 1
usebetas_ps = 0 # whether or not to use the betas for the "Y" PS
usebetas_mot = 0# whether or not to use betas for classifier evidence for MOT
zscoreDV = 0 # if true zscore all DV
zscoreIV = 0 # if you should zscore all classifier values
# specify type of recall evidence to do
DIFF = 1
POST = 2
recalltype = DIFF

def mean_confidence_interval(data, confidence=0.95):
    if confidence == 0.95:
        d_sorted = np.sort(data)
        d_nonnan = d_sorted[np.argwhere(~np.isnan(d_sorted))]
        n = len(d_nonnan)
        n_in_middle = np.int(np.round(confidence*n))
        first_index = (n-n_in_middle)/2 - 1
        m = np.nanmean(data)
        low_val = d_sorted[np.int(first_index)]
        high_val = d_sorted[np.int(first_index+n_in_middle)]
    elif confidence == 0.68:
        # use scipy.stats.sem
        m = np.nanmean(data)
        d_sem = scipy.stats.sem(data,nan_policy='omit')
        low_val = m - d_sem
        high_val = m + d_sem
    return m, low_val, high_val

def nanzscore(inputdata):
    if len(np.shape(inputdata)) == 2:
        nsub = np.shape(inputdata)[1]
        nstim = np.shape(inputdata)[0]
        zdata = np.zeros((nstim,nsub))
        for s in np.arange(nsub):
            zdata[:,s] = (inputdata[:,s] - np.nanmean(inputdata[:,s]))/np.nanstd(inputdata[:,s])
    return zdata

#%%
# LOAD ALL DATA
pickle_in = open("/Volumes/norman/amennen/PythonMot5/evidencebystim_glmclassifier_alpha100_intercept_motion.pickle","rb")
evbystim = pickle.load(pickle_in)

data_dir = '/Volumes/norman/amennen/PythonMot5/'
filename = 'compareExp5.mat'
filepath = os.path.join(data_dir,filename)
d = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
allSep = d['sepbystimD']

# now specify path for betas ps
with open("/Volumes/norman/amennen/PythonMot5/betas_recall_orderedstim_PHG2.pickle", "rb") as f:  # Python 3: open(..., 'rb')
    betasbystim_RT, betasbystim_OM = pickle.load(f)
motpath = '/Volumes/norman/amennen/motStudy05_transferred/'
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'
savepath = '/Volumes/norman/amennen/PythonMot5/'
# new: make logs so that RT is normally distributed**
targRT = np.log(np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy'))
lureRT = np.log(np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')) 
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
nSub= np.int(len(all_sub))
npairs = np.int(nSub/2)
RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)
subarray = np.zeros(nSub)
subarray[YC_ind] = 1
diffEasy = np.zeros((10,npairs*2))
diffHard = np.zeros((10,npairs*2))
postHard = np.zeros((10,npairs*2))
#%%
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
BOLD_RTsim = np.zeros((nstim, nsub))
BOLD_avgRTsim = np.zeros(nsub)
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

#else:
for s in np.arange(nsub):
    subj = all_sub[s]
    # 1/31 after sfn: adding z to zscore differently before taking out timepoints
    # 3/13: using more strict mask
    filename = 'recallPATz_PHG2_%i.mat' % subj
    filepath = os.path.join(data_dir, filename)
    #print(filepath)
    d = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
    RT = d['RTPAT']
    nstim = RT.shape[0]
    nvox = RT.shape[1]
    for stim in np.arange(nstim):
        r1 = RT[stim, :, 0]
        r2 = RT[stim, :, 1]
        BOLD_RTsim[stim, s] = np.corrcoef(r1, r2)[1, 0]

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
#%%
FILTERED_BOLD_RTsim = BOLD_RTsim
FILTERED_BETA_RTsim = BETA_RTsim
FILTERED_diffhard = diffHard # behavioral ratings difference
FILTERED_TRmatrix_kde = TRmatrix_kde
FILTERED_lureRT_CO = lureRT_correctOnly #lureAcc + targAcc
FILTERED_lureRT = lureRT
FILTERED_targRT = targRT

FILTERED_recallact = avg_hardAct_diff
FILTERED_simHard = simHard

for s in np.arange(nsub):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    stimkeep = goodRTstim[sub]
    FILTERED_TRmatrix_kde[~stimkeep,:,s] = np.nan
    FILTERED_BOLD_RTsim[~stimkeep,s] = np.nan
    FILTERED_BETA_RTsim[~stimkeep,s] = np.nan
    FILTERED_diffhard[~stimkeep,s] = np.nan
    FILTERED_lureRT_CO[~stimkeep, s] = np.nan
    FILTERED_lureRT[~stimkeep, s] = np.nan
    FILTERED_targRT[~stimkeep, s] = np.nan
    FILTERED_simHard[~stimkeep, s] = np.nan
    FILTERED_recallact[~stimkeep,s] = np.nan

if zscoreDV:
    FILTERED_diffhard = nanzscore(FILTERED_diffhard) # zscore over subjects
    #FILTERED_TRmatrix_kde[~stimkeep,:,s] = np.nan
    FILTERED_BOLD_RTsim = nanzscore(FILTERED_BOLD_RTsim)
    FILTERED_BETA_RTsim = nanzscore(FILTERED_BETA_RTsim)
    FILTERED_lureRT_CO = nanzscore(FILTERED_lureRT_CO)
    FILTERED_lureRT = nanzscore(FILTERED_lureRT)
    FILTERED_simHard = nanzscore(FILTERED_simHard)
    FILTERED_recallact = nanzscore(FILTERED_recallact)
#%% check that BOLD and BETA are related
#allcorr = np.zeros((nsub))
#allp = np.zeros((nsub))
#plt.figure(figsize=(10,7))    
#for s in np.arange(nsub):
#    x = FILTERED_BETA_RTsim[:,s]
#    y = FILTERED_BOLD_RTsim[:,s]
#    allcorr[s],allp[s] = scipy.stats.pearsonr(x,y)
#    plt.plot(x,y, '.')
#plt.xlabel('beta PS')
#plt.ylabel('BOLD PS')
#plt.title('BOLD vs. beta PS')
#plt.figure(figsize=(10,7))    
##plt.hist(allcorr)
##plt.xlabel('Correlation')
##plt.yticks([0,8,16], ['0', '.25', '.5'])
##plt.ylabel('Frequency')
##plt.title('Histogram of PS correlations')
### find subjects where relationship is less than .5
#badsubj = (allcorr<.5)
#%% now run bootstrap
# analysis: Evidence w/ (1) pattern similarity (2) lure RT (3) detail difference (4) word vector similarity
if nboot > 1:
    correlations_BOLDps = np.zeros((nboot,nwin)) # each result will be for a specific window the mean for that bootstrap run
    correlations_BETAps = np.zeros((nboot,nwin)) 
    correlations_lureRT = np.zeros((nboot,nwin))
    correlations_lureRTCO = np.zeros((nboot,nwin))
    correlations_detaildiff = np.zeros((nboot,nwin))
    correlations_WVsoftmax = np.zeros((nboot,nwin))
    correlations_recallact = np.zeros((nboot,nwin))
    
    for b in np.arange(nboot):
    
        bootsubjects = np.random.randint(0,high=nsub,size=nsub)
        correlations_BOLDps[b, :] = getcorr_vector(FILTERED_BOLD_RTsim[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
        correlations_BETAps[b, :] = getcorr_vector(FILTERED_BETA_RTsim[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
        correlations_lureRT[b, :] = getcorr_vector(FILTERED_lureRT[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
        correlations_lureRTCO[b, :] = getcorr_vector(FILTERED_lureRT_CO[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
        correlations_detaildiff[b, :] = getcorr_vector(FILTERED_diffhard[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
        correlations_WVsoftmax[b, :] = getcorr_vector(FILTERED_simHard[:,bootsubjects],FILTERED_TRmatrix_kde[:,:,bootsubjects])
        correlations_recallact[b, :] = getcorr_vector(FILTERED_recallact[:, bootsubjects],FILTERED_TRmatrix_kde[:, :, bootsubjects])
else:
    correlations_BOLDps = np.zeros((nsub,nwin)) # each result will be for a specific window the mean for that bootstrap run
    correlations_BETAps = np.zeros((nsub,nwin)) # each result will be for a specific window the mean for that bootstrap run
    correlations_lureRT = np.zeros((nsub,nwin))
    correlations_lureRTCO = np.zeros((nsub,nwin))
    correlations_detaildiff = np.zeros((nsub,nwin))
    correlations_WVsoftmax = np.zeros((nsub,nwin))
    correlations_recallact = np.zeros((nsub,nwin))
    correlations_BOLDps = getcorr_matrix(FILTERED_BOLD_RTsim,FILTERED_TRmatrix_kde)
    correlations_BETAps = getcorr_matrix(FILTERED_BETA_RTsim,FILTERED_TRmatrix_kde)
    correlations_lureRT = getcorr_matrix(FILTERED_lureRT,FILTERED_TRmatrix_kde)
    correlations_lureRTCO = getcorr_matrix(FILTERED_lureRT_CO,FILTERED_TRmatrix_kde)
    correlations_detaildiff = getcorr_matrix(FILTERED_diffhard,FILTERED_TRmatrix_kde)
    correlations_WVsoftmax= getcorr_matrix(FILTERED_simHard,FILTERED_TRmatrix_kde)
    correlations_recallact = getcorr_matrix(FILTERED_recallact,FILTERED_TRmatrix_kde)  
#%% check that the correlations are correlated
#x = correlations_BETAps.flatten()
#y = correlations_BOLDps.flatten()
#plt.figure()
#plt.plot(x,y, '.')
#scipy.stats.pearsonr(x,y)
#allcorr = np.zeros((nsub))
## these are also strongy correlated can check per subj
#plt.figure(figsize=(10,7))
#for s in np.arange(nsub):
#    x = correlations_BETAps[s,:]
#    y = correlations_BOLDps[s,:]
#    plt.plot(x,y,'.')
#    allcorr[s],p = scipy.stats.pearsonr(x,y,)
#plt.xlabel('beta correlations')
#plt.ylabel('BOLD correlations')
#plt.title('BOLD vs. beta PS correlations')
#
#plt.figure()
#plt.hist(allcorr)
#
#x = np.mean(correlations_BETAps,axis=0)
#y = np.mean(correlations_BOLDps,axis=0)
#alldifferences = x-y
#plt.figure()
#plt.plot(catrange,alldifferences, '.')
#plt.figure()
#plt.plot(x,y, '.')
## now define bad subj on this this relationship
#badsubj = allcorr< 0.5
#%% get confidence intervals
BOLDps_mean = np.zeros((nwin))
BOLDps_errL = np.zeros((nwin))
BOLDps_errH = np.zeros((nwin))
BETAps_mean = np.zeros((nwin))
BETAps_errL = np.zeros((nwin))
BETAps_errH = np.zeros((nwin))
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
if nboot == 1: #just get 68% confidence interval
    CI = 0.68
else:
    CI = 0.95

for w in np.arange(nwin):
    BOLDps_mean[w],BOLDps_errL[w],BOLDps_errH[w] = mean_confidence_interval(correlations_BOLDps[:,w],CI)
    BETAps_mean[w],BETAps_errL[w],BETAps_errH[w] = mean_confidence_interval(correlations_BETAps[:,w],CI)
    lureRT_mean[w], lureRT_errL[w], lureRT_errH[w] = mean_confidence_interval(correlations_lureRT[:, w],CI)
    lureRTCO_mean[w], lureRTCO_errL[w],lureRTCO_errH[w],  = mean_confidence_interval(correlations_lureRTCO[:, w],CI)
    detaildiff_mean[w], detaildiff_errL[w], detaildiff_errH[w] = mean_confidence_interval(correlations_detaildiff[:, w],CI)
    WVsoftmax_mean[w], WVsoftmax_errL[w],WVsoftmax_errH[w] = mean_confidence_interval(correlations_WVsoftmax[:, w],CI)
    recallact_mean[w], recallact_errL[w],recallact_errH[w] = mean_confidence_interval(correlations_recallact[:, w],CI)
# now calculate pvalues!


#%% now plot results!! pattern similarity first!

fig, ax = plt.subplots(figsize=(12,7))
plt.title('BOLD PATTERN SIMILARITY')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, BOLDps_errL, BOLDps_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,BOLDps_mean, color='r')
#plt.ylim(-.25,.25)
#ax.set_yticks([-.2,-.1,0,.1,.2])


fig, ax = plt.subplots(figsize=(12,7))
plt.title('BETA PATTERN SIMILARITY')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, BETAps_errL, BETAps_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,BETAps_mean, color='r')
#plt.ylim(-.25,.25)
#ax.set_yticks([-.2,-.1,0,.1,.2])

#%%
fig, ax = plt.subplots(figsize=(12,7))
plt.title('LURE RT CO')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, lureRTCO_errL, lureRTCO_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,lureRTCO_mean, color='r')
plt.ylim(-.25,.25)
#ax.set_yticks([-.2,-.1,0,.1,.2])

fig, ax = plt.subplots(figsize=(12,7))
plt.title('WV correlations')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, WVsoftmax_errL, WVsoftmax_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,WVsoftmax_mean, color='r')
plt.ylim(-.25,.25)
#ax.set_yticks([-.2,-.1,0,.1,.2])

fig, ax = plt.subplots(figsize=(12,7))
plt.title('Recall Act')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
sns.despine()
plt.fill_between(catrange, recallact_errL, recallact_errH,facecolor='r',alpha=0.3)
plt.plot(catrange,recallact_mean, color='r')
#plt.ylim(-.25,.25)
#ax.set_yticks([-.2,-.1,0,.1,.2])

#fig, ax = plt.subplots(figsize=(7,5))
#plt.title('Detail diff')
#plt.ylabel('Correlation')
#plt.xlabel('Retrieval evidence bin-kde')
#palette = itertools.cycle(sns.color_palette("husl",8))
#sns.despine()
#plt.fill_between(catrange, detaildiff_errL, detaildiff_errH,facecolor='r',alpha=0.3)
#plt.plot(catrange,detaildiff_mean, color='r')
#plt.ylim(-.25,.25)
#ax.set_yticks([-.2,-.1,0,.1,.2])
