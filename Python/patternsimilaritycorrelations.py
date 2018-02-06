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
sns.set(font_scale = 3)
custom = {'axes.linewidth':5,'font.family':'sans-serif','font.sans-serif':['STHeiti']}
sns.set_style('white',custom)

# this is the one where we're going to take GLM classifier
pickle_in = open("/Volumes/norman/amennen/PythonMot5/evidencebystim_glmclassifieralpha3000.pickle","rb")
evbystim = pickle.load(pickle_in)
# specify now which computer you're using!
motpath = '/Volumes/norman/amennen/motStudy05_transferred/'
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'
savepath = '/Volumes/norman/amennen/PythonMot5/'
flatui = ["#DB5461", "#593C8F"]
all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
subtouse = all_sub
#subtouse =np.array([1,3,4,5,6,8,10,12,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
nSub= np.int(len(subtouse))
npairs = np.int(nSub/2)
RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)


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
    goodRTstim[sub] = hardR[subjName][0,:] > 0

# get the recall data from RECALL patterns coming from nifti files
# now look at the difference in PS between RT and OM
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
# changing 2/5/18: going from original classifier evidence to loaded classifier GLM evidence
windowsize = 0.05
min = -1.5
max = -1*min + windowsize # to go one over
catrange = np.arange(min,max,windowsize)
cr2 = np.reshape(catrange,(len(catrange),1))
nwin = catrange.shape[0]
TRmatrix = np.zeros((nstim,nwin,npairs*2))
TRmatrix_consec = np.zeros((nstim,nwin,npairs*2))
TRmatrix_kde = np.zeros((nstim,nwin,npairs*2))

for s in np.arange(len(subtouse)):
    s_ind = subtouse[s]
    sub = "Subject%01d" % s_ind
    for st in np.arange(nstim):
        thissep = allSep[st,:,s]
        #thissep = evbystim[sub][:,st]
        x2 = np.reshape(thissep, (len(thissep), 1))
        kde = KernelDensity(kernel='gaussian').fit(x2)
        allvals = np.exp(kde.score_samples(cr2))
        for w in np.arange(nwin):
            TRmatrix_kde[st, w, s] = allvals[w]
            #TRmatrix[st,w,s] = np.where((thissep >= catrange[w]) & (thissep < catrange[w+1]))[0].shape[0]
            #z = np.where(np.diff(np.where((thissep >= catrange[w]) & (thissep < catrange[w + 1])))[0] < 3)
            #TRmatrix_consec[st, w, s] = z[0].size
    print(TRmatrix_consec[:,5,s])
allr = getcorr(RTsim,TRmatrix_kde, 'RTsim', 'TRmatrix')



######################################################### LOAD BEHAVIORAL DATA ############################################################################


######################################################### LOAD WORD VECTOR DATA ############################################################################
hard_sm = np.load('/Volumes/norman/amennen/wordVec/hard_smF.npy')
simHard = np.load('/Volumes/norman/amennen/wordVec/simHardF.npy')
corDetHard = np.load('/Volumes/norman/amennen/wordVec/corDetHard.npy')
INcorDetHard = np.load('/Volumes/norman/amennen/wordVec/INcorDetHard.npy')
corDetHard = corDetHard.T
INcorDetHard = INcorDetHard.T
hard_sm = hard_sm.T
simHard = simHard.T

######################################################### LOAD RECOG DATA ############################################################################
targRT = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
lureRT = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
targAcc = np.load('/Volumes/norman/amennen/PythonMot5/targAcc.npy')
lureAcc = np.load('/Volumes/norman/amennen/PythonMot5/lureAcc.npy')
lureAcc_bool = lureAcc==1
targAcc_bool = targAcc==1
allr_lureRT = getcorr(lureRT,TRmatrix, 'lureRT', 'TRmatrix')
allr_lureRT_consec = getcorr(lureRT,TRmatrix_consec, 'lureRT', 'TRmatrix_consec')

# now want to treat correct/incorrect lure trials separately

lureRT_correctOnly = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
lureRT_correctOnly[~lureAcc_bool] = np.nan
lureRT_incorrectOnly = np.load('/Volumes/norman/amennen/PythonMot5/lureRT.npy')
lureRT_incorrectOnly[lureAcc_bool] = np.nan

targRT_correctOnly = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
targRT_correctOnly[~targAcc_bool] = np.nan
targRT_incorrectOnly = np.load('/Volumes/norman/amennen/PythonMot5/targRT.npy')
targRT_incorrectOnly[targAcc_bool] = np.nan

######################################################### Filter at the end ############################################################################
FILTERED_TRmatrix_kde = TRmatrix_kde
FILTERED_RTsim = RTsim
FILTERED_diffhard = diffHard # behavioral ratings difference
FILTERED_lureRT_CO = lureRT_correctOnly
for s in np.arange(len(subtouse)):
    s_ind = subtouse[s]
    sub = "Subject%01d" % s_ind
    stimkeep = goodRTstim[sub]
    FILTERED_TRmatrix_kde[~stimkeep,:,s] = np.nan
    FILTERED_RTsim[~stimkeep,s] = np.nan
    FILTERED_diffhard[~stimkeep,s] = np.nan
    FILTERED_lureRT_CO[~stimkeep,s] = np.nan

######################################################### PLOT TR MAT AND RT SIM ############################################################################
allr = getcorr(FILTERED_RTsim,FILTERED_TRmatrix_kde, 'RTsim', 'TRmatrix')
plotting_data = allr.T
fig, ax = plt.subplots(figsize=(7,5))
plt.title('Pattern Similarity Correlations')
plt.ylabel('Correlation')
plt.xlabel('Retrieval evidence bin-kde')
palette = itertools.cycle(sns.color_palette("husl",8))
yerr = stats.sem(plotting_data, nan_policy='omit')
y = np.nanmean(plotting_data,axis=0)
#ye2 = stats.sem(plot2, nan_policy='omit')
#y2 = np.nanmean(plot2, axis=0)
sns.despine()
plt.fill_between(catrange, y-yerr, y+yerr,facecolor='r',alpha=0.3)
plt.plot(catrange,y, color='r')
#plt.fill_between(np.arange(nwin), y2-ye2, y2+ye2,facecolor=allpallet[3],alpha=0.3)
#plt.plot(y2, color=allpallet[4], linewidth=6)
#l2 = [item.get_text() for item in ax.get_xticklabels()]
#l2 = ["-.5,-.4","-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4",".4,.5"]
#l2 = ["-.4,-.3","-.3,-.2" , "-.2,-.1", "-.1,0",".0,.1",".1,.2",".2,.3" , ".3,.4"]
#ax.set_xticklabels(l2)
#plt.xlim(0,7)
plt.ylim(-.25,.3)
ax.set_yticks([-.2,-.1,0,.1,.2])



plt.show()