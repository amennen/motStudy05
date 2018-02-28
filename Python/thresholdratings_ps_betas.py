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

sns.set(font_scale = 1.5)
custom = {'axes.linewidth':5,'font.family':'sans-serif','font.sans-serif':['STHeiti']}
sns.set_style('white',custom)
bw = 0.1
# this is the one where we're going to take GLM classifier
pickle_in = open("/Volumes/norman/amennen/PythonMot5/evidencebystim_glmclassifier_alpha100_intercept_motion.pickle","rb")
evbystim = pickle.load(pickle_in)

# now specify path for betas ps
with open("/Volumes/norman/amennen/PythonMot5/betas_recall_orderedstim.pickle", "rb") as f:  # Python 3: open(..., 'rb')
    betasbystim_RT, betasbystim_OM = pickle.load(f)


def nanzscore(inputdata):
    zdata = (inputdata - np.nanmean(inputdata))/np.nanstd(inputdata)
    return zdata

# specify now which computer you're using!
motpath = '/Volumes/norman/amennen/motStudy05_transferred/'
behavioral_data = '/Volumes/norman/amennen/motStudy05_transferred/BehavioralData/'
savepath = '/Volumes/norman/amennen/PythonMot5/'
flatui = ["#DB5461", "#593C8F"]
all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
nSub= np.int(len(all_sub))
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
    goodRTstim[sub] = hardR[subjName][0,:] > 2
# now repeat other things only considering those stimuli
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

laptop = 0
if laptop:
    motpath = '/Users/amennen/motStudy05/'
else:
    motpath = '/Users/amennen/Documents/Norman/MOT/motStudy05/'
bd = motpath + 'BehavioralData/'
# get responses during MOT-detail ratings
allsubresponse = {}
for s in np.arange(len(all_sub)):
    subj = all_sub[s]
    sub = "Subject%01d" % subj
    thisev = evbystim[sub]
    subjPath = bd + str(subj) + '/'
    fn = subjPath + 'RTresponses'+ '.mat'
    d = scipy.io.loadmat(fn)
    responses = d['allresp']
    responses = responses.T
    allsubresponse[sub] = responses
# first we have to make TR matrix and filter out ones we don't want
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

hard_sm = np.load('/Volumes/norman/amennen/wordVec/hard_smF.npy')
simHard = np.load('/Volumes/norman/amennen/wordVec/simHardF.npy')
corDetHard = np.load('/Volumes/norman/amennen/wordVec/corDetHard.npy')
INcorDetHard = np.load('/Volumes/norman/amennen/wordVec/INcorDetHard.npy')
corDetHard = corDetHard.T
INcorDetHard = INcorDetHard.T
hard_sm = hard_sm.T # this is the soft max
simHard = simHard.T # this is just the version of it regualr cosines

# so maybe take all points together?
# Filter all points together- first take kde
# filter ev by stim and filter RTsim
# start with not doing the kde because would have to pool all the data together
megadatamatrix = np.empty((0,12))
megapsmatrix = np.empty((0))
megartmatrix = np.empty((0))
megaTrtmatrix = np.empty((0))
megadiffmatrix = np.empty((0))
megasumrt = np.empty((0))
megasm = np.empty((0))
megadet = np.empty((0,12))

# ZSCORED BY SUBJECT
megadatamatrixZ = np.empty((0,12))
megapsmatrixZ = np.empty((0))
megartmatrixZ = np.empty((0))
megaTrtmatrixZ = np.empty((0))
megadiffmatrixZ = np.empty((0))
megasumrtZ = np.empty((0))
megasmZ = np.empty((0))
megadetZ = np.empty((0,12))
for s in np.arange(len(all_sub)):
    s_ind = all_sub[s]
    sub = "Subject%01d" % s_ind
    # calculate individual bw for that subject
    allvals = evbystim[sub]
    stimkeep = goodRTstim[sub]
    allvalsbysubj = allvals[:,stimkeep].T
    subdet = allsubresponse[sub]
    alldetbysub = subdet[:,stimkeep].T
    allpsbysubj = RTsim[stimkeep,s]
    allrtbysubj = lureRT[stimkeep,s]
    allrtTbysubj = targRT[stimkeep, s]
    alldiffbysubj = diffHard[stimkeep,s]
    allwv = hard_sm[stimkeep,s]
    allrt = lureRT[stimkeep,s] + targRT[stimkeep,s]
    # now it's nstim x 12
    megadatamatrix = np.concatenate((megadatamatrix,allvalsbysubj),axis=0)
    megapsmatrix = np.concatenate((megapsmatrix,allpsbysubj))
    megartmatrix = np.concatenate((megartmatrix,allrtbysubj))
    megadiffmatrix = np.concatenate((megadiffmatrix,alldiffbysubj))
    megaTrtmatrix = np.concatenate((megaTrtmatrix,allrtTbysubj))
    megasumrt = np.concatenate((megasumrt,allrt))
    megasm = np.concatenate((megasm,allwv))
    megadet = np.concatenate((megadet,alldetbysub),axis=0)

    megadatamatrixZ = np.concatenate((megadatamatrixZ,nanzscore(allvalsbysubj)),axis=0)
    megapsmatrixZ = np.concatenate((megapsmatrixZ,nanzscore(allpsbysubj)))
    megartmatrixZ = np.concatenate((megartmatrixZ,nanzscore(allrtbysubj)))
    megadiffmatrixZ = np.concatenate((megadiffmatrixZ,nanzscore(alldiffbysubj)))
    megaTrtmatrixZ = np.concatenate((megaTrtmatrixZ,nanzscore(allrtTbysubj)))
    megasumrtZ = np.concatenate((megasumrtZ,nanzscore(allrt)))
    megasmZ = np.concatenate((megasmZ,nanzscore(allwv)))
    megadetZ = np.concatenate((megadetZ,nanzscore(alldetbysub)),axis=0)
# now can just look for correlations!
nstimT = np.shape(megapsmatrix)[0]
windowsize = 0.2
min = -1
max = -1*min + windowsize # to go one over
catrange = np.arange(min,max,windowsize)
cr2 = np.reshape(catrange,(len(catrange),1))
nwin = catrange.shape[0]-1
TRmatrix = np.zeros((nstimT,nwin))

for st in np.arange(nstim):
    thissep = megadatamatrix[st,:]
    for w in np.arange(nwin):
        TRmatrix[st,w] = np.where((thissep >= catrange[w]) & (thissep < catrange[w+1]))[0].shape[0]

# now get total correlation
def getcorr_T(memoryeffect,datamatrix):
    nwin = np.shape(datamatrix)[1]
    allr2 = np.zeros(nwin)
    for w in np.arange(nwin):
        thisTR = datamatrix[:,w]
        nas = np.logical_or(np.isnan(memoryeffect),np.isnan(thisTR))
        allr2[w] = scipy.stats.pearsonr(memoryeffect[~nas],thisTR[~nas])[0]
    return allr2

allr = getcorr_T(megapsmatrix,TRmatrix)
corrplot = allr.T
plt.figure()
for s in np.arange(nSub):
    plt.plot(catrange[0:-1],corrplot, '-.')
plt.show()

avgev = np.median(megadatamatrix,axis=1)
sum_ev = np.sum(megadatamatrix,axis=1)/np.shape(megadatamatrix)[1]
plt.figure()
plt.plot(avgev,megapsmatrix, '.')

#order by least to most similar
ps_ord = np.argsort(megapsmatrix)
data_ord = megadatamatrix[ps_ord,:]
data_all = np.reshape(data_ord,(nstimT*12,1))
x = np.arange(nstimT)
x2 = np.repeat(x,12)
x2 = np.reshape(x2,(len(x2),1))
df = np.concatenate((data_all,x2),axis=1)
plotdata = pd.DataFrame(df,columns=['ps','order'])
plt.figure()
sns.boxplot(data=plotdata,x="order",y="ps")
# could split by half?
# mulitply variability times mean
#allstd = np.reshape(np.std(megadatamatrix,axis=1),(294,1))
allstd = np.std(data_ord,axis=1)
allmn = np.mean(data_ord,axis=1)
weighted = allmn/allstd
plt.figure()
plt.plot(np.arange(nstimT),weighted, '.')


# do lure and ps relate?
plt.figure()
plt.plot(megartmatrix,megapsmatrix, '.')
nas = np.logical_or(np.isnan(megartmatrix), np.isnan(megapsmatrix))
scipy.stats.pearsonr(megartmatrix[~nas],megapsmatrix[~nas])

plt.figure()
plt.plot(megaTrtmatrix,megapsmatrix, '.')
nas = np.logical_or(np.isnan(megaTrtmatrix), np.isnan(megapsmatrix))
scipy.stats.pearsonr(megaTrtmatrix[~nas],megapsmatrix[~nas])


# QUESTION: is pattern similarity related to other behavioral measures?
# seems to be related to reaction time--stronger when combine both reaction times
# describe: problem with averaging individual subjects is that not a lot of sample points get treated as equally as other points
plt.figure()
plt.plot(megartmatrixZ,megapsmatrixZ, '.')
plt.xlabel('Lure RT + Targ RT (all trials)')
plt.ylabel('Pattern Similarity')
plt.title('RT vs. PS')

nas = np.logical_or(np.isnan(megartmatrixZ), np.isnan(megapsmatrixZ))
scipy.stats.pearsonr(megartmatrixZ[~nas],megapsmatrixZ[~nas])

# so what do those stimuli have in common? do they have similar distributions of mot scores

# check with word vector descriptions
plt.figure()
plt.plot(megasmZ,megapsmatrixZ, '.')
plt.xlabel('Softmax word vector correlation')
plt.ylabel('Pattern Similarity')
plt.title('Word vector vs. PS')

nas = np.logical_or(np.isnan(megasmZ), np.isnan(megapsmatrixZ))
scipy.stats.pearsonr(megasmZ[~nas],megapsmatrixZ[~nas])



# sum everything (some will be nan)
nas = np.logical_or(np.logical_or(np.isnan(megasm), np.isnan(megapsmatrix)),np.isnan(megasumrt))
MEMORYSCORETOTAL = megapsmatrixZ[~nas] - megartmatrixZ[~nas] + megasmZ[~nas]
#MEMORYSCORETOTAL = scipy.stats.zscore(megapsmatrix[~nas]) - scipy.stats.zscore(megasumrt[~nas]) + scipy.stats.zscore(megasm[~nas])

#MEMORYSCORETOTAL = -1* scipy.stats.zscore(megasumrt[~nas])
MEM_SORTED = np.argsort(MEMORYSCORETOTAL)
SCORE_SORTED = MEMORYSCORETOTAL[MEM_SORTED]

TOTALEV1 = megadatamatrixZ[~nas,:]
TOTALEV2 =  megadetZ[~nas, :]
EV_SORTED = TOTALEV1[MEM_SORTED,:]
EV_SORTED2 = TOTALEV1[MEM_SORTED,:]

VAR_EV_SORTED = np.nanmean(EV_SORTED,axis=1) #+ np.nanmean(EV_SORTED2,axis=1)
nex = len(SCORE_SORTED)
# compare evidence distributions
EV1 = EV_SORTED[0:np.floor(nex/3),:].flatten()
EV2 = EV_SORTED[-np.floor(nex/3):,:].flatten()

EV1 = VAR_EV_SORTED[0:np.floor(nex/3)].flatten()
EV2 = VAR_EV_SORTED[-np.floor(nex/3):].flatten()

gr = np.concatenate((np.zeros(len(EV1)),np.ones(len(EV2))),axis=0)
data = np.concatenate((EV1,EV2),axis=0)
data2b = np.concatenate((np.reshape(data,(len(data),1)),np.reshape(gr,(len(gr),1))),axis=1)
df = pd.DataFrame(data2b, columns=['data','group'])
fig, ax = plt.subplots()
sns.stripplot(data=df,y='data',jitter=True,alpha=0.25,x="group",split="True" , color="k")
sns.barplot(data=df,y='data',x="group")

#plt.xlim([-1,1])
plt.show()

scipy.stats.ttest_ind(EV1,EV2)


# NOW TRY WITH DETAIL RATINGS
# zscore within subjects first????
nas = np.logical_or(np.logical_or(np.isnan(megasm), np.isnan(megapsmatrix)),np.isnan(megasumrt))
MEMORYSCORETOTAL = megapsmatrixZ[~nas] - megartmatrixZ[~nas] + megasmZ[~nas]

#MEMORYSCORETOTAL = scipy.stats.zscore(megapsmatrix[~nas]) - scipy.stats.zscore(megasumrt[~nas]) + scipy.stats.zscore(megasm[~nas])
IVTOTAL = scipy.stats.zscore(megadet[~nas,:])
#MEMORYSCORETOTAL = -1* scipy.stats.zscore(megasumrt[~nas])
MEM_SORTED = np.argsort(MEMORYSCORETOTAL)
SCORE_SORTED = MEMORYSCORETOTAL[MEM_SORTED]

TOTALEV = megadetZ[~nas,:]
EV_SORTED = TOTALEV[MEM_SORTED,:]
VAR_EV_SORTED = np.max(EV_SORTED,axis=1)
nex = len(SCORE_SORTED)
# compare evidence distributions
EV1 = EV_SORTED[0:np.floor(nex/3),:].flatten()
EV2 = EV_SORTED[-np.floor(nex/3):,:].flatten()

#EV1 = VAR_EV_SORTED[0:np.floor(nex/3)].flatten()
#EV2 = VAR_EV_SORTED[-np.floor(nex/3):].flatten()

gr = np.concatenate((np.zeros(len(EV1)),np.ones(len(EV2))),axis=0)
data = np.concatenate((EV1,EV2),axis=0)
data2b = np.concatenate((np.reshape(data,(len(data),1)),np.reshape(gr,(len(gr),1))),axis=1)
df = pd.DataFrame(data2b, columns=['data','group'])
fig, ax = plt.subplots()
sns.stripplot(data=df,y='data',jitter=True,alpha=0.25,x="group",split="True" , color="k")
sns.barplot(data=df,y='data',x="group")

#plt.xlim([-1,1])
plt.show()

scipy.stats.ttest_ind(EV1,EV2)

xmatrix = np.concatenate(megapsmatrixZ[~nas],megasumrtZ)
# can you make glm to predict outcome based on x's?
regModel = linear_model.LinearRegression()
regModel.fit

# make histogram of evidence
evRT = np.empty((0,12))
evYC = np.empty((0,12))
for s in np.arange(npairs):
    s_ind = RT_sub[s]
    sub = "Subject%01d" % s_ind
    allvals = evbystim[sub]
    sv = allvals[:, stimkeep].T
    evRT = np.concatenate((evRT, sv), axis=0)
    s_ind = YC_sub[s]
    sub = "Subject%01d" % s_ind
    allvals = evbystim[sub]
    sv = allvals[:, stimkeep].T
    evYC = np.concatenate((evYC, sv), axis=0)

eRT = evRT.flatten()
eYC = evYC.flatten()
allev = np.concatenate((eRT,eYC),axis=0)
groups = np.arange(2)
groups = np.repeat(groups,120*npairs)
data = np.concatenate((allev[:,np.newaxis],groups[:,np.newaxis]),axis=1)
df = pd.DataFrame(data,columns=['y', 's'])
plt.figure()
sns.distplot(eRT, color='r')
sns.distplot(eYC)