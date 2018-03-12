# what this code does:

# - plots factor plot of recall activation by group (RT/YC,H/E)


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
from sklearn.neighbors import KernelDensity
from usefulfns import getcorr
import pickle
sns.set(font_scale = 1.5)
custom = {'axes.linewidth':5,'font.family':'sans-serif','font.sans-serif':['STHeiti']}
sns.set_style('white',custom)

# this is the one where we're going to take GLM classifier
data_dir = '/Volumes/norman/amennen/PythonMot5/'
filename = 'compareExp5.mat'
filepath = os.path.join(data_dir,filename)
d = scipy.io.loadmat(filepath, squeeze_me=True, struct_as_record=False)
allSep = d['sepbystimD']
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
subarray = np.zeros(nSub)
subarray[YC_ind] = 1
allpallet = ["#DB5461", "#FFD9CE", "593C8F", "#8EF9F3", "#171738"]

flatui = ["#DB5461", "#593C8F"]

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

# divide everything into thirds--do they get more related over time?
diffHard = np.zeros((10,npairs*2))
diffEasy = np.zeros((10,npairs*2))
postHard = np.zeros((10,npairs*2))
preHard = np.zeros((10,npairs*2))
easyR = {}
hardR = {}
for s in np.arange(npairs*2):
    subjPath = behavioral_data + str(subtouse[s]) + '/'
    subjName = 'subj' + str(s)
    sub = "Subject%01d" % subtouse[s]
    fn = glob.glob(subjPath + 'ratings'+ '.mat')
    d = scipy.io.loadmat(fn[0])
    easyR[subjName] = d['easyScores']
    # find any nans
    nas = np.isnan(easyR[subjName])
    easyR[subjName] = easyR[subjName].astype(np.float16)
    easyR[subjName][nas] = np.nan
    hardR[subjName] = d['hardScores']
    nas = np.isnan(hardR[subjName])
    hardR[subjName] = hardR[subjName].astype(np.float16)
    hardR[subjName][nas] = np.nan
    # int 16 converts the scores to integers, but the problem is that nan's become zeros (when they should stay as nan's)
    diffEasy[:, s] = np.diff(easyR[subjName], axis=0)
    diffHard[:, s] = np.diff(hardR[subjName], axis=0)
    postHard[:, s] = hardR[subjName][1, :]
    preHard[:, s] = hardR[subjName][0, :]
nruns=3
nstim = 10
all_correlations = np.zeros((nstim,nruns,npairs))
# make vector to normalize correlations
all_other_correlations = np.zeros((nstim,nruns,npairs))
for p in np.arange(npairs):
    s1 = RT_sub[p]
    s2 = YC_sub[p]
    for st in np.arange(nstim):
        timecourse1 = allSep[st,:,RT_ind[p]]
        timecourse2 = allSep[st,:,YC_ind[p]]
        all_correlations[st,0,p] = scipy.stats.pearsonr(timecourse1[0:15],timecourse2[0:15])[0]
        all_correlations[st,1,p] = scipy.stats.pearsonr(timecourse1[15:30],timecourse2[15:30])[0]
        all_correlations[st,2,p] = scipy.stats.pearsonr(timecourse1[30:45],timecourse2[30:45])[0]
        # do for each stimuli
        othermatches = np.delete(np.arange(npairs), p)
        othercor = np.zeros((nruns,npairs-1))
        for j in np.arange(npairs-1):
            # wiat shouldn't be done by stimulus because not other see same stimulus there
            # now match the yoked to another RT person
            timecourse1 = allSep[st, :, RT_ind[othermatches[j]]]
            othercor[0,j] = scipy.stats.pearsonr(timecourse1[0:15],timecourse2[0:15])[0]
            othercor[1,j] = scipy.stats.pearsonr(timecourse1[15:30],timecourse2[15:30])[0]
            othercor[2,j] = scipy.stats.pearsonr(timecourse1[30:45],timecourse2[30:45])[0]
        all_other_correlations[st,0,p] = np.mean(othercor[0,:])
        all_other_correlations[st, 1, p] = np.mean(othercor[1 :])
        all_other_correlations[st, 2, p] = np.mean(othercor[2, :])
# now each pair gets average over stimulus per time block
# average all_other_correlations becuase stimuli don't matter
normalized_correlations = all_correlations/np.mean(all_other_correlations,axis=0)
#pair_avg = np.mean(normalized_correlations,axis=0)
pair_avg = np.mean(all_correlations,axis=0)
pair_total_score = np.mean(pair_avg,axis=0)
# plot data
scores = pair_avg.flatten()
runs = np.concatenate((np.zeros((npairs)),np.ones((npairs)),2*np.ones((npairs))),axis=0)
pairs = np.concatenate((np.arange((npairs)),np.arange((npairs)),np.arange((npairs))),axis=0)
data = np.concatenate((scores[:,np.newaxis],runs[:,np.newaxis],pairs[:,np.newaxis]),axis=1)
df = pd.DataFrame(data=data, columns=['scores', 'runs', 'pairs'])
fig, ax = plt.subplots(figsize=(7,5))
sns.pointplot(data=df,x='runs', y='scores', hue='pairs')


# see stimulus average
# maybe just look at are correlations higher/lower than baseline?
stim_avg = np.mean(all_correlations,axis=1)
RTsim_RT = RTsim[:,RT_ind]
RTsim_YC = RTsim[:,YC_ind]
sim_differences = RTsim_RT - RTsim_YC
plt.figure()
plt.plot(stim_avg, sim_differences, '.')
scipy.stats.pearsonr(stim_avg.flatten(),sim_differences.flatten())

diffhard_RT = diffHard[:,RT_ind]
diffhard_YC = diffHard[:,YC_ind]
ratingdiff = np.abs(diffhard_RT - diffhard_YC)
prediff = np.abs(preHard[:,RT_ind] - preHard[:,YC_ind])
postdiff = np.abs(postHard[:,RT_ind] - postHard[:,YC_ind])
avgscores = (preHard + postHard)/2
avgdiff = np.abs(avgscores[:,RT_ind] - avgscores[:,YC_ind])

plt.figure()
plt.plot(stim_avg, postdiff, '.')
plt.ylim([-1 ,6])
scipy.stats.pearsonr(stim_avg.flatten(),prediff.flatten())
rv= np.zeros((npairs))
pv = np.zeros((npairs))
for p in np.arange(npairs):
    rv[p], pv[p] = scipy.stats.pearsonr(stim_avg[:,p],prediff[:,p])


plt.figure()
y = prediff.flatten()
x = stim_avg.flatten()
keep = np.argwhere(y<10)
plt.plot(x[keep], y[keep], '.')
plt.xlabel('Average correlation during RT MOT')
plt.ylabel('Difference in pre MOT ratings abs(RT-YC)')
plt.title('Difference in behavioral ratings pre MOT vs. similarity during RT')
scipy.stats.pearsonr(stim_avg.flatten(),postdiff.flatten())

plt.figure()
plt.plot(stim_avg, avgdiff, '.')
nas = np.logical_or(np.isnan(stim_avg.flatten()), np.isnan(avgdiff.flatten()))
scipy.stats.pearsonr(stim_avg.flatten()[~nas],avgdiff.flatten()[~nas])

diff_combined = (prediff + postdiff)/2

plt.figure()
plt.plot(stim_avg, avgdiff, '.')
nas = np.logical_or(np.isnan(stim_avg.flatten()), np.isnan(diff_combined.flatten()))
scipy.stats.pearsonr(stim_avg.flatten()[~nas],diff_combined.flatten()[~nas])

plt.figure()
plt.plot(RTsim,preHard, '.')
plt.ylim([0, 6])


# load in into matrices the different actiavtions?
# maybe first just plot time course
easyAct = np.zeros((nstim,4,2,nsub))
hardAct = np.zeros((nstim,4,2,nsub))
inteldir = '/Volumes/norman/amennen/motStudy05_transferred/datafromintelrt/data/'
for s in np.arange(npairs*2):
    subjPath = inteldir + str(subtouse[s]) + '/'
    subjName = 'subj' + str(s)
    sub = "Subject%01d" % subtouse[s]
    fn = glob.glob(subjPath + 'recallactivations'+ '.mat')
    d = scipy.io.loadmat(fn[0])
    easyAct[:,:,:,s] = d['easy_activation']
    hardAct[:,:,:,s] = d['hard_activation']

avg_easyAct = np.mean(easyAct,axis=0)
avg_hardAct = np.mean(hardAct,axis=0)
easyAct_diff = np.diff(easyAct,axis=2).squeeze()
hardAct_diff = np.diff(hardAct,axis=2).squeeze()
avg_hardAct_diff = np.mean(hardAct_diff,axis=0)
avg_easyAct_diff = np.mean(easyAct_diff,axis=0)

# now make labels to plot
easy1 = avg_easyAct[:,0,:].flatten() # goes through each person's first
sublabels = np.tile(np.arange(nsub),4)
TR = np.arange(4)
easy2 = avg_easyAct[:,1,:].flatten()
hard1 = avg_hardAct[:,0,:].flatten()
hard2 = avg_hardAct[:,1,:].flatten()
data = np.concatenate((easy1[:,np.newaxis],easy2[:,np.newaxis],hard1[:,np.newaxis],hard2[:,np.newaxis]),axis=0)
pair = np.concatenate((sublabels[:,np.newaxis],sublabels[:,np.newaxis],sublabels[:,np.newaxis],sublabels[:,np.newaxis]),axis=0)
#allTR = np.tile(TR,(nsub*4))
allTR = np.concatenate((np.repeat(TR,nsub),np.repeat(TR,nsub),np.repeat(TR,nsub),np.repeat(TR,nsub)),axis=0)
easyhard = np.concatenate((np.ones((nsub*4,1)),np.ones((nsub*4,1)),np.zeros((nsub*4,1)),np.zeros((nsub*4,1))),axis=0)
run1run2 = np.concatenate((np.zeros((nsub*4,1)),np.ones((nsub*4,1)),np.zeros((nsub*4,1)),np.ones((nsub*4,1))),axis=0)
subarray_t = np.tile(subarray,4)
rtyc = np.concatenate((subarray_t[:,np.newaxis],subarray_t[:,np.newaxis],subarray_t[:,np.newaxis],subarray_t[:,np.newaxis]),axis=0)
largedata = np.concatenate((data,pair,allTR[:,np.newaxis],easyhard,run1run2,rtyc),axis=1)
df = pd.DataFrame(data=largedata,columns=['data','pair','tr', 'EH', 'R1R2', 'RTYC'])
sns.factorplot(data=df,x='tr',y='data', hue='R1R2', col='RTYC',row='EH')

easy = avg_easyAct_diff.flatten() # cycles over subjects first
hard = avg_hardAct_diff.flatten()
data = np.concatenate((easy[:,np.newaxis],hard[:,np.newaxis]),axis=0) # cycles through all subjects 4 x  EACH
pair = np.concatenate((sublabels[:,np.newaxis],sublabels[:,np.newaxis]),axis=0)
easyhard = np.concatenate((np.ones((nsub*4,1)),np.zeros((nsub*4,1))),axis=0)
rtyc = np.concatenate((subarray_t[:,np.newaxis],subarray_t[:,np.newaxis]),axis=0)
# TR isn't right!! it goes all 0's all 1's all 2's all 3's
#allTR = np.tile(TR,(nsub*2))
allTR = np.concatenate((np.repeat(TR,nsub),np.repeat(TR,nsub)),axis=0)
largedata2 = np.concatenate((data,pair,allTR[:,np.newaxis],easyhard,rtyc),axis=1)
df = pd.DataFrame(data=largedata2,columns=['data','pair','tr', 'EH', 'RTYC'])
sns.factorplot(data=df,x='tr',y='data', col='RTYC',row='EH', ci=68)

# is there a significant differnece?
subjavg_easy_diff = np.mean(avg_easyAct_diff,axis=0)
subjavg_hard_diff = np.mean(avg_hardAct_diff,axis=0)

stats.ttest_rel(subjavg_hard_diff[RT_ind], subjavg_hard_diff[YC_ind])
stats.ttest_rel(subjavg_easy_diff[RT_ind], subjavg_easy_diff[YC_ind])

plt.figure()
sns.pointplot(data=df,x='tr', y='data', hue="R1R2")

# now can correlate for each pair--average MOT evidence and average Recall evidence afterwards?
# or it can be just difference before/afterwards like with recall ratings

# can also correlate recall ratings and average recall evidence
avgbystim_easyAct = np.mean(easyAct,axis=1)
diff_easyAct = np.diff(avgbystim_easyAct,axis=1).squeeze() # it's two minus one
avgbystim_hardAct = np.mean(hardAct,axis=1)
diff_hardAct = np.diff(avgbystim_hardAct,axis=1).squeeze() # it's two minus one

diffhardAct_RT = diff_hardAct[:,RT_ind]
diffhardAct_YC = diff_hardAct[:,YC_ind]
recallactdiff = np.abs(diffhardAct_RT - diffhardAct_YC)

plt.figure()
plt.plot(ratingdiff,recallactdiff, '.')

plt.figure()
plt.plot(stim_avg,recallactdiff, '.')
plt.ylim([-1 ,6])
scipy.stats.pearsonr(stim_avg.flatten(),recallactdiff.flatten())

# in general is there a correlation between stimuli activation and ratings?
allhard_act = diff_hardAct.flatten()
allhard_ratings = diffHard.flatten()

plt.figure()
plt.plot(allhard_act,allhard_ratings, '.')
scipy.stats.pearsonr(allhard_act,allhard_ratings)
