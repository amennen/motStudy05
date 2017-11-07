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
import usefulfns
from usefulfns import getcorr
from astropy.modeling.models import Lorentz1D
from astropy.convolution import convolve, Gaussian1DKernel, Box1DKernel
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
TRmatrix_smooth = np.zeros((nstim,nwin,npairs*2))
for s in np.arange(npairs*2):
    for st in np.arange(nstim):
        thissep = allSep[st,:,s]
        box_kernel = Box1DKernel(2)
        smoothed_data_box = convolve(thissep, box_kernel)
        for w in np.arange(nwin):
            #TRmatrix[st,w,s] = np.where((thissep >= catrange[w]) & (thissep < catrange[w+1]))[0].shape[0]
            TRmatrix[st, w, s] = np.where((smoothed_data_box >= catrange[w]) & (smoothed_data_box < catrange[w+1]))[0].shape[0]
            z = np.where(np.diff(np.where((smoothed_data_box >= catrange[w]) & (smoothed_data_box < catrange[w + 1])))[0] < 4)
            TRmatrix_smooth[st, w, s] = z[0].size
# # plot evidence and see how spikey it is//how smoothing helps
s = 27
fig, axarr = plt.subplots(5, 2)
for st in np.arange(nstim):
    if st < 5:
        col=0
        row = st
    else:
        col=1
        row = st - 5
    y = allSep[st,:,s]
    axarr[row,col].plot(y, '-o', color='r', label='orig')
    #axarr[st].ylim(-1,1)
    n=len(allSep[st,:,s])
    box_kernel = Box1DKernel(2)
    smoothed_data_box = convolve(y, box_kernel)
    axarr[row,col].plot(smoothed_data_box, '-o', color='g', label='boxwidth3')
    gauss_kernel = Gaussian1DKernel(1)
    smoothed_data_gauss = convolve(y, gauss_kernel)
    axarr[row,col].plot(smoothed_data_gauss, '-o', color='b', label='gauss1')
    axarr[row,col].set_ylim(-1,1)
plt.legend()


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




allr = getcorr(RTsim,TRmatrix, 'RTsim', 'TRmatrix')
allr_consec = getcorr(RTsim,TRmatrix_smooth, 'RTsim', 'TRmatrix_smooth')

fig, ax = plt.subplots()
plotting_data = allr.T
plot2 = allr_consec.T
plt.title('Pattern Similarity Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
ye2 = stats.sem(plot2, nan_policy='omit')
y2 = np.nanmean(plot2, axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(y, color='r')
plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor='b',alpha=0.3)
plt.plot(y2, color='b')
plt.xlim(0,7)
plt.ylim(-.2,.2)
l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
plt.legend(['Time', 'Consecutive Time'])
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
allr_ratings = getcorr(diffHard,TRmatrix, 'diffHard', 'TRmatrix')
allr_ratings_consec = getcorr(diffHard,TRmatrix_smooth, 'diffHard', 'TRmatrix_smooth')

fig, ax = plt.subplots()
plotting_data = allr_ratings.T
plot2 = allr_ratings_consec.T
plt.title('Detail Difference Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
ye2 = stats.sem(plot2, nan_policy='omit')
y2 = np.nanmean(plot2, axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(y, color='r')
plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor='b',alpha=0.3)
plt.plot(y2, color='b')
plt.xlim(0,7)
plt.ylim(-.2,.2)

l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
plt.legend(['Time', 'Consecutive Time'])
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

allr_sim_sm = getcorr(hard_sm,TRmatrix, 'sm', 'TRmatrix')
allr_sim_sm_consec = getcorr(hard_sm,TRmatrix_smooth, 'sm', 'TRmatrix_smooth')
fig, ax = plt.subplots()
plotting_data = allr_sim_sm.T
plot2 = allr_sim_sm_consec.T
plt.title('Mturk Similarity Softmax Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
ye2 = stats.sem(plot2, nan_policy='omit')
y2 = np.nanmean(plot2, axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(y, color='r')
plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor='b',alpha=0.3)
plt.plot(y2, color='b')
plt.xlim(0,7)
plt.ylim(-.2,.2)
l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
plt.legend(['Time', 'Consecutive Time'])
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)



allr_sim = getcorr(simHard,TRmatrix, 'simHard', 'TRmatrix')
allr_sim_consec = getcorr(simHard,TRmatrix_smooth, 'simHard', 'TRmatrix_smooth')
fig, ax = plt.subplots()
plotting_data = allr_sim.T
plot2 = allr_sim_consec.T
plt.title('Mturk Similarity Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
ye2 = stats.sem(plot2, nan_policy='omit')
y2 = np.nanmean(plot2, axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(y, color='r')
plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor='b',alpha=0.3)
plt.plot(y2, color='b')
plt.xlim(0,7)
plt.ylim(-.2,.2)
l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
plt.legend(['Time', 'Consecutive Time'])
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)


allr_cordet = getcorr(corDetHard,TRmatrix, 'cordet', 'TRmatrix')
allr_cordet_consec = getcorr(corDetHard,TRmatrix_smooth, 'cordet', 'TRmatrix_smooth')
fig, ax = plt.subplots()
plotting_data = allr_cordet.T
plot2 = allr_cordet_consec.T
plt.title('Correct # Details Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
ye2 = stats.sem(plot2, nan_policy='omit')
y2 = np.nanmean(plot2, axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(y, color='r')
plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor='b',alpha=0.3)
plt.plot(y2, color='b')
plt.xlim(0,7)
plt.ylim(-.2,.2)

l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
plt.legend(['Time', 'Consecutive Time'])
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)

allr_INcordet = getcorr(INcorDetHard,TRmatrix, 'incordet', 'TRmatrix')
allr_INcordet_consec = getcorr(INcorDetHard,TRmatrix_smooth, 'incordet', 'TRmatrix_smooth')
fig, ax = plt.subplots()
plotting_data = allr_INcordet.T
plot2 = allr_INcordet_consec.T
plt.title('INCorrect # Details Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
ye2 = stats.sem(plot2, nan_policy='omit')
y2 = np.nanmean(plot2, axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(y, color='r')
plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor='b',alpha=0.3)
plt.plot(y2, color='b')
plt.xlim(0,7)
plt.ylim(-.2,.2)

l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
plt.legend(['Time', 'Consecutive Time'])
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)

fig, ax = plt.subplots()
plt.plot(simHard,RTsim, '.')

######################################################### LOAD RECOG DATA ############################################################################
targRT = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
lureRT = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')

allr_lureRT = getcorr(lureRT,TRmatrix, 'lureRT', 'TRmatrix')
allr_lureRT_consec = getcorr(lureRT,TRmatrix_smooth, 'lureRT', 'TRmatrix_smooth')
fig, ax = plt.subplots()
plotting_data = allr_lureRT.T
plot2 = allr_lureRT_consec.T
plt.title('Lure RT Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
ye2 = stats.sem(plot2, nan_policy='omit')
y2 = np.nanmean(plot2, axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(y, color='r')
plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor='b',alpha=0.3)
plt.plot(y2, color='b')
plt.xlim(0,7)
plt.ylim(-.2,.2)

l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
plt.legend(['Time', 'Consecutive Time'])
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)

allr_targRT = getcorr(targRT,TRmatrix, 'targRT', 'TRmatrix')
allr_targRT_consec = getcorr(targRT,TRmatrix_smooth, 'targRT', 'TRmatrix_smooth')
fig, ax = plt.subplots()
plotting_data = allr_targRT.T
plot2 = allr_targRT_consec.T
plt.title('Target RT Correlations')
plt.ylabel('Correlation')
plt.xlabel('MOT Evidence Bin')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
ye2 = stats.sem(plot2, nan_policy='omit')
y2 = np.nanmean(plot2, axis=0)
plt.fill_between(np.arange(nwin), y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(y, color='r')
plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor='b',alpha=0.3)
plt.plot(y2, color='b')
plt.xlim(0,7)
plt.ylim(-.2,.2)

l2 = [item.get_text() for item in ax.get_xticklabels()]
l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
ax.set_xticklabels(l2)
plt.legend(['Time', 'Consecutive Time'])
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
    item.set_fontsize(20)
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(15)


plt.show()



plt.show()