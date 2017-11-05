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
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'
save_dir = '/Volumes/norman/amennen/PythonMot5/'
RECALL = np.array([20,24])
RECOGNITION = 26
Z = norm.ppf

all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
nSub= np.int(len(all_sub))
npairs = np.int(nSub/2)
RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)
subarray = np.zeros(nSub)
subarray[YC_ind] = 1
responses = {}
nTrials = 43
# stimID
# cond
# RT
# accuracy
# resp
# cresp (1 = target, 2 = lure)
for s in np.arange(nSub):
    subjPath = behavioral_data + str(all_sub[s]) + '/'
    subjName = 'subj' + str(s)
    print(all_sub[s])
    fn = glob.glob(subjPath + 'recog'+ '*.mat')
    d = scipy.io.loadmat(fn[0])

    responses[subjName] = np.zeros((nTrials,5))
    responses[subjName][:,0] = d['stimID'][:,0]
    responses[subjName][:,1] = d['cond'][:,0]
    responses[subjName][:,2] = d['rt'][:,0]
    responses[subjName][:,3] = d['acc'][:,0]
    responses[subjName][:,4] = d['cresp'][:,0]
lureRT = np.zeros((10,nSub))
targRT = np.zeros((10,nSub))
for s in np.arange(nSub):
    data = responses['subj'+str(s)]
    RT = data[:,2]
    acc = data[:,3]
    cond = data[:,1]
    cresp = data[:,4]
    id = data[:,0]
    hard = np.argwhere(cond == 1)
    target = np.argwhere(cresp == 1)
    lure = np.argwhere(cresp == 2)
    HL = np.intersect1d(hard,lure)
    HT = np.intersect1d(hard,target)
    id_HL = id[HL]
    x_HL = np.argsort(id_HL)
    RT_HL = RT[HL]
    lureRT[:,s] = RT_HL[x_HL]
    id_HT = id[HT]
    x_HT = np.argsort(id_HT)
    RT_HT = RT[HT]
    targRT[:,s] = RT_HT[x_HT]


np.save(save_dir + 'targRT',targRT)
np.save(save_dir + 'lureRT',lureRT)

## LOOK AT RATINGS
# LOOK AT RATINGS NOW
easyR = {}
hardR = {}
diffEasy = np.zeros((10,nSub))
diffHard = np.zeros((10,nSub))
nTrials = 43

# compute for RT and YC separately
for s in np.arange(nSub):
    subjPath = behavioral_data + str(all_sub[s]) + '/'
    subjName = 'subj' + str(s)
    fn = glob.glob(subjPath + 'ratings'+ '.mat')
    d = scipy.io.loadmat(fn[0])

    easyR[subjName] = d['easyScores']
    hardR[subjName] = d['hardScores']
    diffEasy[:,s] = np.diff(easyR[subjName].astype(np.int16),axis=0)
    diffHard[:,s] = np.diff(hardR[subjName].astype(np.int16),axis=0)

easyAvg = np.mean(diffEasy,axis=0)
hardAvg = np.mean(diffHard,axis=0)

# first do for RT
easyAvgG = easyAvg[RT_ind]
hardAvgG = hardAvg[RT_ind]
eA = np.mean(easyAvgG)
hA = np.mean(hardAvgG)
fig, ax = plt.subplots(figsize=(10,7))
palette = itertools.cycle(sns.color_palette("husl",8))
# first show all of them
data = np.concatenate((np.reshape(hardAvgG,(npairs,1)),np.reshape(easyAvgG,(npairs,1))),axis=0)
subj = np.concatenate((np.reshape(np.arange(npairs),(npairs,1)),np.reshape(np.arange(npairs),(npairs,1))),axis=0)
ABC = np.concatenate((np.zeros((npairs,1)),np.ones((npairs,1))),axis=0)
data2b = np.concatenate((data,subj,ABC),axis=1)
df = pd.DataFrame(data2b,columns = ['Similarity','Subject','RTOM'])

labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "Real-time"
labels[1] = "Omit"
plt.plot([0,1],[hA,eA],color='k', linewidth=5.0)

sns.pointplot(data=df,x="RTOM", y="Similarity", hue="Subject", linestyles="--")
plt.title('Ratings Differences')
ax.set_xticklabels(labels)
plt.ylabel('Ratings')
plt.xlabel('Condition')
plt.ylim([-1.5,1])
plt.show()

# repeat for YC gropu
easyAvgG = easyAvg[YC_ind]
hardAvgG = hardAvg[YC_ind]
eA = np.mean(easyAvgG)
hA = np.mean(hardAvgG)
fig, ax = plt.subplots(figsize=(10,7))
palette = itertools.cycle(sns.color_palette("husl",8))
# first show all of them
data = np.concatenate((np.reshape(hardAvgG,(npairs,1)),np.reshape(easyAvgG,(npairs,1))),axis=0)
subj = np.concatenate((np.reshape(np.arange(npairs),(npairs,1)),np.reshape(np.arange(npairs),(npairs,1))),axis=0)
ABC = np.concatenate((np.zeros((npairs,1)),np.ones((npairs,1))),axis=0)
data2b = np.concatenate((data,subj,ABC),axis=1)
df = pd.DataFrame(data2b,columns = ['Similarity','Subject','RTOM'])

labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "Real-time"
labels[1] = "Omit"
plt.plot([0,1],[hA,eA],color='k', linewidth=5.0)

sns.pointplot(data=df,x="RTOM", y="Similarity", hue="Subject", linestyles="--")
plt.title('Ratings Differences')
ax.set_xticklabels(labels)
plt.ylabel('Ratings')
plt.xlabel('Condition')
plt.ylim([-1.5,1])
plt.show()


# taken from: http://lindeloev.net/calculating-d-in-python-and-php/
def dPrime(hits, misses, fas, crs):
    # Floors an ceilings are replaced by half hits and half FA's
    halfHit = 0.5 / (hits + misses)
    halfFa = 0.5 / (fas + crs)

    # Calculate hitrate and avoid d' infinity
    hitRate = hits / (hits + misses)
    if hitRate == 1: hitRate = 1 - halfHit
    if hitRate == 0: hitRate = halfHit

    # Calculate false alarm rate and avoid d' infinity
    faRate = fas / (fas + crs)
    if faRate == 1: faRate = 1 - halfFa
    if faRate == 0: faRate = halfFa

    # Return d', beta, c and Ad'
    out = {}
    out['d'] = Z(hitRate) - Z(faRate)
    out['beta'] = exp((Z(faRate) ** 2 - Z(hitRate) ** 2) / 2)
    out['c'] = -(Z(hitRate) + Z(faRate)) / 2
    out['Ad'] = norm.cdf(out['d'] / sqrt(2))
    return out