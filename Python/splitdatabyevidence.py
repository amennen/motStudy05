# purpose: split data by evidence then organize by other categories

import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn import svm
from sklearn.metrics import r2_score
import pickle
import pdb
import sklearn
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
from usefulfns import *
import pickle
from sklearn.neighbors import KernelDensity
from sklearn import linear_model

RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)

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
usebetas = 0
zscoreDV = 0 # if true zscore all DV
zscoreIV = 0 # if you should zscore all classifier values


# LOAD ALL DATA
pickle_in = open("/Volumes/norman/amennen/PythonMot5/evidencebystim_glmclassifier_alpha100_intercept_motion.pickle","rb")
evbystim = pickle.load(pickle_in)
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
    hardR[subjName] = d['hardScores']
    diffEasy[:,s] = np.diff(easyR[subjName].astype(np.int16),axis=0)
    diffHard[:,s] = np.diff(hardR[subjName].astype(np.int16),axis=0)
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
if usebetas:
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



windowsize = 0.05
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
for s in np.arange(nsub):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    subvals_z = stats.zscore(evbystim[sub])
    for st in np.arange(nstim):
        thissep = evbystim[sub][:,st]
        if zscoreIV:
            thissep = subvals_z[:,st]
        x2 = np.reshape(thissep, (len(thissep), 1))
        kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(x2)
        allvals = np.exp(kde.score_samples(cr2))
        TRmatrix_kde[st, :, s] = allvals

sumforkde = np.sum(TRmatrix_kde[:,20:24,:],axis=1).flatten()
stim_ordered = np.argsort(sumforkde) # not sure if it's okay to combine separate kde's like this--otherwise could zscore first
# goes from low to high most time spent in that range
low_stim = stim_ordered[0:len(stim_ordered)/2]
high_stim = stim_ordered[len(stim_ordered)/2:]

FILTERED_RTsim = RTsim
FILTERED_diffhard = diffHard # behavioral ratings difference
FILTERED_TRmatrix_kde = TRmatrix_kde
FILTERED_lureRT_CO = lureRT_correctOnly #lureAcc + targAcc
FILTERED_lureRT = lureRT
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
    FILTERED_simHard[~stimkeep, s] = np.nan

wv_vector = FILTERED_simHard.flatten()
ps_vector = FILTERED_RTsim.flatten()
lr_vector = FILTERED_lureRT_CO.flatten()

allstim = np.zeros((len(stim_ordered)))
allstim[high_stim] = 1
alldata = np.concatenate((ps_vector[:,np.newaxis],wv_vector[:,np.newaxis],lr_vector[:,np.newaxis],allstim[:,np.newaxis]),axis=1)
df = pd.DataFrame(alldata,columns=['ps','wv','lr','ord'])

plt.figure()
#sns.stripplot(data=df,x='ord',y='ps', color='black',alpha=0.5)
sns.barplot(data=df,x='ord',y='ps',ci=68)

plt.figure()
#sns.stripplot(data=df,x='ord',y='ps', color='black',alpha=0.5)
sns.barplot(data=df,x='ord',y='wv',ci=68)

plt.figure()
#sns.stripplot(data=df,x='ord',y='ps', color='black',alpha=0.5)
sns.barplot(data=df,x='ord',y='lr',ci=68)