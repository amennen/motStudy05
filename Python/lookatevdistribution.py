# purpose: check distribution of all evidence
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

nboot = 1000
bw = 0.07 # set it here for everyone!!
detailthreshold = 2
usebetas = 0
zscoreDV = 1 # if true zscore all DV
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



windowsize = 0.05
#min = -1
#max = -1*min + windowsize # to go one over
min=-0.7
max=0.9
if zscoreIV:
    windowsize = 0.1
    min=-1.5
    max=-1*min
catrange = np.arange(min,max+windowsize,windowsize)
cr2 = np.reshape(catrange,(len(catrange),1))
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
bw2=0.05
bins = np.arange(-1,1.1,bw2)
x = plt.hist(vectorevidence[:,np.newaxis],bins)



plt.figure()
plt.plot(catrange,allvals, 'r')
nx = x[0]/np.sum(x[0])
plt.plot(bins[0:-1]+(bw2/2),nx*10, 'b.')
plt.show()
ev_mean = np.mean(vectorevidence)
ev_std = np.std(vectorevidence)
ev_mean + 3*ev_std
plt.xlabel('Evidence')
plt.ylabel('Proportion in that range')
plt.title('Normalized counts: all classifier evidence')

nstim=10
bw = 0.05
# plot for individual person
bins = np.arange(min,max,windowsize)
s = np.random.randint(0,high=nsub,size=1)
st = np.random.randint(0,nstim,size=1)
s_ind = all_sub[s]
sub = "Subject%01d" % s_ind
bw = 0.1
thissep = evbystim[sub][:,st]
x2 = np.reshape(thissep, (len(thissep), 1))
x = plt.hist(x2,bins, normed=True)
plt.figure()
nx = x[0]/np.sum(x[0])
plt.plot(catrange[0:-2],nx, 'b.')
kde = KernelDensity(kernel='gaussian', bandwidth=bw).fit(x2)
allvals = np.exp(kde.score_samples(cr2))
plt.plot(catrange,allvals, 'r')
plt.legend(['normalized counts', 'kde'])
plt.xlabel('Classifier evidence')
plt.ylabel('Frequency')
plt.title('KDE vs. HISTOGRAM')