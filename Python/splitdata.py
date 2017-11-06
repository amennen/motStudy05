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
nstim = 10
# first calculate TRmatrix_consec
# sort the subjects by the highest number of bins consecutive
# see how the distribution varies by RT/YC subjects

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
TRmatrix_consec = np.zeros((nstim,nwin,npairs*2))

for s in np.arange(npairs*2):
    for st in np.arange(nstim):
        thissep = allSep[st,:,s]
        for w in np.arange(nwin):
            TRmatrix[st,w,s] = np.where((thissep >= catrange[w]) & (thissep < catrange[w+1]))[0].shape[0]
            z = np.where(np.diff(np.where((thissep >= catrange[w]) & (thissep < catrange[w + 1])))[0] < 4)
            TRmatrix_consec[st, w, s] = z[0].size

win = 5 # this is the range we care about
n_consec_low = TRmatrix_consec[:,win,:]


plt.legend()
sub_consec_low = np.median(n_consec_low,axis=0)
#sub_consec_low = np.sum(n_consec_low,axis=0)
ord_low_high = np.argsort(sub_consec_low)
ordered_sub = sub_consec_low[ord_low_high]
subject_ordered = all_sub[ord_low_high]
low_sub = subject_ordered[0:np.int(nSub/2)]
high_sub = subject_ordered[np.int(nSub/2):]

# WOWOWOW data is split evenly! Reshuffle!
low_index = ord_low_high[0:np.int(nSub/2)]
high_index = ord_low_high[np.int(nSub/2):]
subarray_ord = np.zeros(nSub)
subarray_ord[high_index] = 1

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

MOT_sim_minus_diff = RTsim - RTsim_diff
OMIT_sim_minus_diff = OMITsim - OMITsim_diff
sub_MOT_sim = np.mean(RTsim,axis=0)
sub_OMIT_sim = np.mean(OMITsim,axis=0)
sub_MOT_sim_minus_diff = np.mean(MOT_sim_minus_diff,axis=0)
sub_OMIT_sim_minus_diff = np.mean(OMIT_sim_minus_diff,axis=0)

master_diff = MOT_sim_minus_diff - OMIT_sim_minus_diff
sub_master_diff = np.mean(master_diff,axis=0)

n_consec_vector = np.reshape(n_consec_low,(10*32,1))
largeord = np.argsort(n_consec_vector,axis=0)[:,0]
ind_1 = largeord[0:80]
ind_2 = largeord[80:160]
ind_3 = largeord[160:240]
ind_4 = largeord[240:]

# CHANGE TO TWO GROUPS
ind_1 = largeord[0:160]
ind_2 = largeord[160:]

RTsim_vector = np.reshape(RTsim-RTsim_diff,(10*32,1))
allInd = np.zeros((len(largeord),1))
allInd[ind_2] = 1
#allInd[ind_3] = 2
#allInd[ind_4] = 3

data = RTsim_vector
group = allInd
data2b = np.concatenate((data,group),axis=1)
df = pd.DataFrame(data2b, columns=['data','group'])
fig, ax = plt.subplots()
#RTsim_vector[ind_1]
sns.barplot(data=df,x='group', y='data',ci=68)
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "IND1"
labels[1] = "IND2"

lureRT = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
lureRT_vector = np.reshape(lureRT,(10*32,1))
data = lureRT_vector
group = allInd
data2b = np.concatenate((data,group),axis=1)
df = pd.DataFrame(data2b, columns=['data','group'])
fig, ax = plt.subplots()
#RTsim_vector[ind_1]
sns.barplot(data=df,x='group', y='data',ci=68)
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "IND1"
labels[1] = "IND2"
plt.ylim(.5,.8)

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
diffHard_vector = np.reshape(diffHard,(10*32,1))
data = diffHard_vector
group = allInd
data2b = np.concatenate((data,group),axis=1)
df = pd.DataFrame(data2b, columns=['data','group'])
fig, ax = plt.subplots()
#RTsim_vector[ind_1]
sns.barplot(data=df,x='group', y='data',ci=68)
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "IND1"
labels[1] = "IND2"






plt.plot(n_consec_vector,RTsim_vector, '.')
for s in np.arange(nsub):
    #plt.plot(n_consec_low[:,s], label='Sub %i' %s)
    thissub = n_consec_low[:,s]
    stimOrd = np.argsort(thissub)


# plot subtracted values for each category by group
#data = np.concatenate((np.reshape(sub_MOT_sim_minus_diff,(nsub,1)),np.reshape(sub_OMIT_sim_minus_diff,(nsub,1))),axis=0)
data = np.concatenate((np.reshape(sub_MOT_sim,(nsub,1)),np.reshape(sub_OMIT_sim,(nsub,1))),axis=0)
group = np.concatenate((np.reshape(subarray_ord,(nsub,1)),np.reshape(subarray_ord,(nsub,1))),axis=0)
type = np.concatenate((0*np.ones((nsub,1)),np.ones((nsub,1))),axis=0)
data2b = np.concatenate((data,group,type),axis=1)
df = pd.DataFrame(data2b, columns=['data','group','type'])
fig, ax = plt.subplots()
sns.boxplot(data=df,x='type', y='data', hue="group")
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "MOT: SAME-ALL"
labels[1] = "OMIT: SAME-ALL"


data = np.reshape(sub_master_diff,(nsub,1))
group = np.reshape(subarray_ord,(nsub,1))
data2b = np.concatenate((data,group),axis=1)
df = pd.DataFrame(data2b, columns=['data','group'])
fig, ax = plt.subplots()
sns.boxplot(data=df,x='group', y='data')
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "IND1"
labels[1] = "IND2"

data = np.reshape(sub_MOT_sim,(nsub,1))
group = np.reshape(subarray_ord,(nsub,1))
data2b = np.concatenate((data,group),axis=1)
df = pd.DataFrame(data2b, columns=['data','group'])
fig, ax = plt.subplots()
sns.violinplot(data=df,x='group', y='data')
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "IND1"
labels[1] = "IND2"

# want RT same vs. RT
#sns.barplot(data = df,hue='type', y='data', x='group',split=True,color='k',  alpha=0.7)
#plt.title('Pattern Similarity Post vs. Pre MOT')
#ax.set_xticklabels(labels)
#plt.ylabel('Pattern Similarity')
#plt.ylim(-.1,.25)
#plt.show()
data = np.concatenate((np.reshape(sub_consec_low,(32,1)),np.reshape(subarray_ord,(32,1))),axis=1)
df = pd.DataFrame(data,columns=['data', 'group'])
plt.hist(sub_consec_low[low_index], alpha=0.5)
plt.hist(sub_consec_low[high_index], alpha=0.5)

