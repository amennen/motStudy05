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
from usefulfns import getcorr

# specify now which computer you're using!
motpath = '/Volumes/norman/amennen/motStudy05_transferred/'
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'

all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
nSub= np.int(len(all_sub))
npairs = np.int(nSub/2)
RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)
subarray = np.zeros(nSub)
subarray[YC_ind] = 1

# get the recall data from RECALL patterns coming from nifti files
# now look at the difference in PS between RT and OM
nsub = 32
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
for s in np.arange(nsub):
    subj = all_sub[s]
    filename = 'recallPAT%i.mat' % subj
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

avgRTsim = np.mean(RTsim, axis=0)
avgOMITsim = np.mean(OMITsim, axis=0)
avgRTsim_diff = np.mean(RTsim_diff, axis=0)
avgOMITsim_diff = np.mean(OMITsim_diff, axis=0)


######################################################### LOAD EVIDENCE DATA ############################################################################
data_dir = '/Volumes/norman/amennen/PythonMot5/'
filename = 'compareExp5.mat'
filepath = os.path.join(data_dir,filename)
d = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
secondDiff = d['allSecondDiffD']
goodRange = d['nGoodRangeD']
nLow = d['nLowD']
nHigh = d['nHighD']
vectorSep = d['vectorSepD']
lowAverage = d['avg_lowD']
highAverage = d['avg_highD']
nConsec = d['nConsecD']
allSep = d['sepbystimD']
allSpeed = d['speedbystimD']

# separate by RT/YC
sepRT = vectorSep[RT_ind,:]
sepYC = vectorSep[YC_ind,:]
goodRange_RT = goodRange[RT_ind]
goodRange_YC = goodRange[YC_ind]
lowA_RT = lowAverage[RT_ind,:]
lowA_YC = lowAverage[YC_ind,:]
highA_RT = highAverage[RT_ind,:]
highA_YC = highAverage[YC_ind,:]
nConsec_RT = nConsec[RT_ind]
nConsec_YC = nConsec[YC_ind]
secondDiff_RT = secondDiff[RT_ind,:]
secondDiff_YC = secondDiff[YC_ind,:]
allSep_RT = allSep[:,:,RT_ind]
allSep_YC = allSep[:,:,YC_ind]
nTRs = 450
allSepVec_RT = np.reshape(allSep_RT,(npairs*nTRs,1))
allSepVec_YC = np.reshape(allSep_YC,(npairs*nTRs,1))


# BUILD TR MATRIX FOR EVIDENCE
windowsize = 0.1
min = -0.4
max = -1*min + windowsize # to go one over
catrange = np.arange(min,max,windowsize)
nwin = catrange.shape[0] - 1
TRmatrix = np.zeros((nstim,nwin,npairs*2))
for s in np.arange(npairs*2):
    for st in np.arange(nstim):
        thissep = allSep[st,:,s]
        for w in np.arange(nwin):
            TRmatrix[st,w,s] = np.where((thissep >= catrange[w]) & (thissep < catrange[w+1]))[0].shape[0]

allr = getcorr(RTsim,TRmatrix)

fig, ax = plt.subplots()
plt.title('Pattern Similarity Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
sem = stats.sem(allr.T)
yerr = sem
y = np.mean(allr.T,axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
sns.tsplot(data=allr.T,color=next(palette), err_style=None)
#sns.barplot(data=allr.T,ci=0)
#sns.stripplot(data=allr.T)
#plt.plot(pvals)

l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)

######################################################### LOAD BEHAVIORAL DATA ############################################################################
diffEasy = np.zeros((10,npairs*2))
diffHard = np.zeros((10,npairs*2))
easyR = {}
hardR = {}
for s in np.arange(npairs*2):
    subjPath = behavioral_data + str(all_sub[s]) + '/'
    subjName = 'subj' + str(s)
    fn = glob.glob(subjPath + 'ratings'+ '.mat')
    d = scipy.io.loadmat(fn[0])
    easyR[subjName] = d['easyScores']
    hardR[subjName] = d['hardScores']
    diffEasy[:,s] = np.diff(easyR[subjName].astype(np.int16),axis=0)
    diffHard[:,s] = np.diff(hardR[subjName].astype(np.int16),axis=0)
allr = getcorr(diffHard,TRmatrix)
# because in this case some subjects are bad!!!
bad = np.unique(np.where(np.isnan(allr))[1])
mask = np.ones(32, dtype=bool)
mask[bad] = False
allr = allr[:,mask]
fig, ax = plt.subplots()
plt.title('Detail Difference Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
sem = stats.sem(allr.T)
yerr = sem
y = np.mean(allr.T,axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
sns.tsplot(data=allr.T,color=next(palette), err_style=None)
#sns.barplot(data=allr.T,ci=0)
#sns.stripplot(data=allr.T)
#plt.plot(pvals)

l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)

######################################################### LOAD WORD VECTOR DATA ############################################################################
hard_sm = np.load('/Volumes/norman/amennen/wordVec/hard_smF.npy')
simHard = np.load('/Volumes/norman/amennen/wordVec/simHardF.npy')
corDetHard = np.load('/Volumes/norman/amennen/wordVec/corDetHard.npy')
INcorDetHard = np.load('/Volumes/norman/amennen/wordVec/INcorDetHard.npy')
corDetHard = corDetHard.T
INcorDetHard = INcorDetHard.T
hard_sm = hard_sm.T
simHard = simHard.T

allr = getcorr(hard_sm,TRmatrix)
fig, ax = plt.subplots()
plt.title('Softmax Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
sem = stats.sem(allr.T)
yerr = sem
y = np.mean(allr.T,axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
sns.tsplot(data=allr.T,color=next(palette), err_style=None)
l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)


allr = getcorr(simHard,TRmatrix)
fig, ax = plt.subplots()
plt.title('Mturk Sim Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
sem = stats.sem(allr.T)
yerr = sem
y = np.mean(allr.T,axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
sns.tsplot(data=allr.T,color=next(palette), err_style=None)
l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)


fig, ax = plt.subplots()
plt.plot(simHard,RTsim, '.')





plt.show()