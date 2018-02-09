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
import pickle
from sklearn.neighbors import KernelDensity

# goal: plot correlation results-maybe each subject's average as line plot
# (do the same for plotting averages)

sns.set(font_scale = 1.5)
custom = {'axes.linewidth':5,'font.family':'sans-serif','font.sans-serif':['STHeiti']}
sns.set_style('white',custom)
bw = 0.15
# this is the one where we're going to take GLM classifier
pickle_in = open("/Volumes/norman/amennen/PythonMot5/evidencebystim_glmclassifier_alpha100_intercept_motion.pickle","rb")
evbystim = pickle.load(pickle_in)
# specify now which computer you're using!
motpath = '/Volumes/norman/amennen/motStudy05_transferred/'
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'
savepath = '/Volumes/norman/amennen/PythonMot5/'
flatui = ["#DB5461", "#593C8F"]
all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
subtouse = all_sub
nSub= np.int(len(subtouse))
npairs = np.int(nSub/2)
RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)
subarray = np.zeros(nSub)
subarray[YC_ind] = 1
allpallet = ["#DB5461", "#FFD9CE", "593C8F", "#8EF9F3", "#171738"]
# we want matrix of good stimuli for each subject
diffEasy = np.zeros((10,npairs*2))
diffHard = np.zeros((10,npairs*2))
postHard = np.zeros((10,npairs*2))
goodRTstim = {}
easyR = {}
hardR = {}
for s in np.arange(npairs*2):
    subjPath = behavioral_data + str(subtouse[s]) + '/'
    subjName = 'subj' + str(s)
    sub = "Subject%01d" % subtouse[s]
    fn = glob.glob(subjPath + 'ratings'+ '.mat')
    d = scipy.io.loadmat(fn[0])
    easyR[subjName] = d['easyScores']
    hardR[subjName] = d['hardScores']
    diffEasy[:,s] = np.diff(easyR[subjName].astype(np.int16),axis=0)
    diffHard[:,s] = np.diff(hardR[subjName].astype(np.int16),axis=0)
    postHard[:,s] = hardR[subjName][1,:]
    goodRTstim[sub] = hardR[subjName][0,:] > 1

nsub = nSub
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
    subj = subtouse[s]
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

windowsize = 0.5
min = -1
max = -1*min + windowsize # to go one over
catrange = np.arange(min,max,windowsize)
cr2 = np.reshape(catrange,(len(catrange),1))
nwin = catrange.shape[0] - 1 # add the -1 when NOT doing KDE
TRmatrix_kde = np.zeros((nstim,nwin,npairs*2))

for s in np.arange(len(subtouse)):
    s_ind = subtouse[s]
    sub = "Subject%01d" % s_ind
    # calculate individual bw for that subject
    allvals_z = evbystim[sub]
    #allvals_z = scipy.stats.zscore(allvals)

    for st in np.arange(nstim):
        #thissep = allSep[st,:,s]
        #thissep = evbystim[sub][:,st]
        thissep = allvals_z[:,st]
        bw = 1.06 * np.std(allvals_z) * (12 ** -.2)
        x2 = np.reshape(thissep, (len(thissep), 1))
        kde = KernelDensity(kernel='gaussian',bandwidth=bw).fit(x2)
        allvals = np.exp(kde.score_samples(cr2))
        for w in np.arange(nwin):
            TRmatrix_kde[st,w,s] = np.where((thissep >= catrange[w]) & (thissep < catrange[w+1]))[0].shape[0]
            #TRmatrix_kde[st, w, s] = allvals[w]
            #z = np.where(np.diff(np.where((thissep >= catrange[w]) & (thissep < catrange[w + 1])))[0] < 3)
            #TRmatrix_consec[st, w, s] = z[0].size

targRT = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
lureRT = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
targAcc = np.load('/Volumes/norman/amennen/PythonMot5/targAcc.npy')
lureAcc = np.load('/Volumes/norman/amennen/PythonMot5/lureAcc.npy')
lureAcc_bool = lureAcc==1
targAcc_bool = targAcc==1
#allr_lureRT = getcorr(lureRT,TRmatrix, 'lureRT', 'TRmatrix')
#allr_lureRT_consec = getcorr(lureRT,TRmatrix_consec, 'lureRT', 'TRmatrix_consec')

# now want to treat correct/incorrect lure trials separately

lureRT_correctOnly = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
lureRT_correctOnly[~lureAcc_bool] = np.nan
lureRT_incorrectOnly = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
lureRT_incorrectOnly[lureAcc_bool] = np.nan

targRT_correctOnly = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
targRT_correctOnly[~targAcc_bool] = np.nan
targRT_incorrectOnly = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
targRT_incorrectOnly[targAcc_bool] = np.nan

FILTERED_lureRT_CO = lureRT_correctOnly
FILTERED_TRmatrix_kde = TRmatrix_kde
FILTERED_RTsim = RTsim
FILTERED_diffhard = diffHard # behavioral ratings difference
FILTERED_posthard = postHard
for s in np.arange(len(subtouse)):
    s_ind = subtouse[s]
    sub = "Subject%01d" % s_ind
    stimkeep = goodRTstim[sub]
    FILTERED_TRmatrix_kde[~stimkeep,:,s] = np.nan
    FILTERED_RTsim[~stimkeep,s] = np.nan
    FILTERED_diffhard[~stimkeep,s] = np.nan
    FILTERED_posthard[~stimkeep,s] = np.nan
    FILTERED_lureRT_CO[~stimkeep, s] = np.nan
allr = getcorr(FILTERED_RTsim,FILTERED_TRmatrix_kde, 'RTsim', 'TRmatrix')

corrplot = allr.T
plt.figure()
for s in np.arange(nSub):
    plt.plot(catrange[0:-1],corrplot[s,:], '-.')
plt.show()

allr_lureRT = getcorr(FILTERED_lureRT_CO,FILTERED_TRmatrix_kde, 'diffHard', 'TRmatrix')

corrplot = allr_lureRT.T
plt.figure()
for s in np.arange(nSub):
    plt.plot(catrange[0:-1],corrplot[s,:], '-.')
plt.ylim([-2,2])
plt.show()
