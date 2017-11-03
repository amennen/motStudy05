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


# first plot for RT subjects
avgRTsimG = avgRTsim[RT_ind]
avgOMITsimG = avgOMITsim[RT_ind]
avgRTsim_diffG = avgRTsim_diff[RT_ind]
avgOMITsim_diffG = avgOMITsim_diff[RT_ind]

# now make plots--different subtraction look back on meeting notes!!
# now plot conditions: realtime similar, realtime dissimilar, omit similar, omit dissimilar
sns.set(style="whitegrid", color_codes=True)
data = np.concatenate((np.reshape(avgRTsimG,(npairs,1)),np.reshape(avgRTsim_diffG,(npairs,1)),np.reshape(avgOMITsimG,(npairs,1)),np.reshape(avgOMITsim_diffG,(npairs,1))),axis=0)
type = np.concatenate((0*np.ones((npairs,1)),0*np.ones((npairs,1)),np.ones((npairs,1)),np.ones((npairs,1))),axis=0) # MOT or OMIT
simdiff = np.concatenate((0*np.ones((npairs,1)),np.ones((npairs,1)),np.zeros((npairs,1)),np.ones((npairs,1))),axis=0) # SAME vs DIFF
data2b = np.concatenate((data,type,simdiff),axis=1)

df = pd.DataFrame(data2b, columns=['data','type', 'simdiff'])
fig, ax = plt.subplots()
sns.violinplot(data=df,x='type', y='data', hue="simdiff")
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "RT"
labels[1] = "OMIT"
# want RT same vs. RT
#sns.swarmplot(data = df,hue='simdiff', y='data', x='type',split=True,color='k',  alpha=0.7)
x1 = np.array([-.2,0.2])
x2 = np.array([.8,1.2])
for s in np.arange(npairs):
    y1 = np.array([avgRTsimG[s],avgRTsim_diffG[s]])
    y2= np.array([avgOMITsimG[s],avgOMITsim_diffG[s]])
    plt.plot(x1,y1)
    plt.plot(x2,y2)

plt.title('Pattern Similarity Post vs. Pre MOT')
ax.set_xticklabels(labels)
plt.ylabel('Pattern Similarity')
#plt.ylim(-.1,.25)
plt.show()

#########################################################
# now make plots--different subtraction look back on meeting notes!!
# now plot conditions: realtime similar, realtime dissimilar, omit similar, omit dissimilar
sns.set(style="whitegrid", color_codes=True)
data = np.concatenate((np.reshape(avgRTsim,(nsub,1)),np.reshape(avgRTsim_diff,(nsub,1)),np.reshape(avgOMITsim,(nsub,1)),np.reshape(avgOMITsim_diff,(nsub,1))),axis=0)
type = np.concatenate((0*np.ones((nsub,1)),np.ones((nsub,1)),2*np.ones((nsub,1)),3*np.ones((nsub,1))),axis=0)
group = np.concatenate((np.reshape(subarray,(nsub,1)),np.reshape(subarray,(nsub,1)),np.reshape(subarray,(nsub,1)),np.reshape(subarray,(nsub,1))),axis=0)
data2b = np.concatenate((data,group,type),axis=1)

df = pd.DataFrame(data2b, columns=['data','group','type'])
fig, ax = plt.subplots()
sns.violinplot(data=df,x='type', y='data', hue="group")
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "RT: Same Item"
labels[1] = "RT: Different Item"
labels[2] = "OMIT: Same Item"
labels[3] = "OMIT: Different Item"
# want RT same vs. RT
sns.swarmplot(data = df,hue='type', y='data', x='group',split=True,color='k',  alpha=0.7)
plt.title('Pattern Similarity Post vs. Pre MOT')
ax.set_xticklabels(labels)
plt.ylabel('Pattern Similarity')
#plt.ylim(-.1,.25)
plt.show()