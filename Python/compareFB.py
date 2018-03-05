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
import os
from matplotlib.pyplot import cm
import matplotlib
matplotlib.rcParams.update({'font.size': 22})
import seaborn as sns
from scipy import stats
import pandas as pd
import itertools

sns.set(font_scale = 3)
custom = {'axes.linewidth':5,'font.family':'sans-serif','font.sans-serif':['STHeiti']}
sns.set_style('white',custom)
savepath = '/Users/amennen/Dropbox/sfn2017/'
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

# now divide data into RT and YC groups
all_sub = np.array([1,3,4,5,6,8,10,11,12,13,14,16,17,19,20,21,23,25,26,27,29,30,31,32,33,34,35,36,37,38,39,40])
RT_sub = np.array([1, 3, 4,5,6,8,10,12,13,14,19,21,23,26,29,32])
YC_sub = np.array([20,40,36,34,11,27,17,16,35,30,39,25,37,31,38,33])
RT_ind = np.searchsorted(all_sub,RT_sub)
YC_ind = np.searchsorted(all_sub,YC_sub)

nsub = len(all_sub)
npairs = len(RT_sub)
nPairs = npairs
nTrials = 10
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
# now repeat analysis above but with the new data that includes all 15 TR's!
ALLEVIDENCE = np.array([])
ALLSUBJECT = np.array([])
ALLGROUP = np.array([])
ALLSESSION = np.array([])
target = 0.15
res = 0

nTRs = 15
nCorrect_RT = np.zeros((nPairs, 3))
nCorrect_YC = np.zeros((nPairs, 3))
errorScale_RT = np.empty((0, 2), int)
errorScale_YC = np.empty((0, 2), int)
errorScale_Pos_RT = np.empty((0, 2), int)
errorScale_Pos_YC = np.empty((0, 2), int)
errorScale_Neg_RT = np.empty((0, 2), int)
errorScale_Neg_YC = np.empty((0, 2), int)
nCor_Pos_RT = np.zeros((nPairs, 3))
nCor_Neg_RT = np.zeros((nPairs, 3))
nCor_Pos_YC = np.zeros((nPairs, 3))
nCor_Neg_YC = np.zeros((nPairs, 3))

for s in np.arange(nPairs):

    RT_ev = allSep_RT[:, :, s]
    YC_ev = allSep_YC[:, :, s]

    # RTmat = np.concatenate((RT_ev[:,:,0],RT_ev[:,:,1],RT_ev[:,:,2]))
    # YCmat = np.concatenate((YC_ev[:,:,0],YC_ev[:,:,1],YC_ev[:,:,2]))

    RT_diff = np.diff(RT_ev, axis=1)
    YC_diff = np.diff(YC_ev, axis=1)

    for session in np.arange(3):
        session_ind = np.arange(session * nTRs, (session + 1) * nTRs)
        thisSession_RT = RT_ev[:, session_ind]
        RT_diff = np.diff(thisSession_RT, axis=1)
        thisSession_YC = YC_ev[:, session_ind]
        YC_diff = np.diff(thisSession_YC, axis=1)

        pos_error_RT = 0
        neg_error_RT = 0
        pos_error_YC = 0
        neg_error_YC = 0

        for row in np.arange(nTrials): # iterate over all stimuli
            for col in np.arange(nTRs - 2):
                # RT_thisdiff = RT_diff[row,col+1]
                RT_thisdiff = thisSession_RT[row, col + 2] - thisSession_RT[row, col]
                RT_1 = thisSession_RT[row, col]
                if RT_1 > target - res:  # want to decrease
                    pos_error_RT += 1
                    if RT_thisdiff < 0:
                        nCorrect_RT[s, session] += 1
                        nCor_Pos_RT[s, session] += 1
                    # print(np.array([RT_1-target,RT_thisdiff]))
                    errorScale_Pos_RT = np.concatenate(
                        (errorScale_Pos_RT, np.reshape(np.array([RT_1 - target, RT_thisdiff]), (1, 2))), axis=0)
                elif RT_1 < target + res:  # if RT is less than target
                    neg_error_RT += 1
                    if RT_thisdiff > 0:
                        nCorrect_RT[s, session] += 1
                        nCor_Neg_RT[s, session] += 1
                    errorScale_Neg_RT = np.concatenate(
                        (errorScale_Neg_RT, np.reshape(np.array([RT_1 - target, RT_thisdiff]), (1, 2))), axis=0)
                errorScale_RT = np.concatenate(
                    (errorScale_RT, np.reshape(np.array([RT_1 - target, RT_thisdiff]), (1, 2))), axis=0)
                # YC_thisdiff = YC_diff[row,col+1]
                YC_thisdiff = thisSession_YC[row, col + 2] - thisSession_YC[row, col]
                YC_1 = thisSession_YC[row, col]
                if YC_1 > target - res:  # want to decrease
                    pos_error_YC += 1
                    if YC_thisdiff < 0:
                        nCorrect_YC[s, session] += 1
                        nCor_Pos_YC[s, session] += 1
                    errorScale_Pos_YC = np.concatenate(
                        (errorScale_Pos_YC, np.reshape(np.array([YC_1 - target, YC_thisdiff]), (1, 2))), axis=0)
                elif YC_1 < target + res:  # if RT is less than target
                    neg_error_YC += 1
                    if YC_thisdiff > 0:
                        nCorrect_YC[s, session] += 1
                        nCor_Neg_YC[s, session] += 1
                    errorScale_Neg_YC = np.concatenate(
                        (errorScale_Neg_YC, np.reshape(np.array([YC_1 - target, YC_thisdiff]), (1, 2))), axis=0)
                errorScale_YC = np.concatenate(
                    (errorScale_YC, np.reshape(np.array([YC_1 - target, YC_thisdiff]), (1, 2))), axis=0)
        nCorrect_RT[s, session] = nCorrect_RT[s, session] / (pos_error_RT + neg_error_RT)
        nCorrect_YC[s, session] = nCorrect_YC[s, session] / (pos_error_YC + neg_error_YC)
        nCor_Pos_RT[s, session] = nCor_Pos_RT[s, session] / pos_error_RT
        nCor_Neg_RT[s, session] = nCor_Neg_RT[s, session] / neg_error_RT
        nCor_Pos_YC[s, session] = nCor_Pos_YC[s, session] / pos_error_YC
        nCor_Neg_YC[s, session] = nCor_Neg_YC[s, session] / neg_error_YC

# average over subjects
pos_rt_sessavg = np.mean(nCor_Pos_RT, axis=1)
pos_yc_sessavg = np.mean(nCor_Pos_YC, axis=1)
neg_rt_sessavg = np.mean(nCor_Neg_RT, axis=1)
neg_yc_sessavg = np.mean(nCor_Neg_YC, axis=1)
data = np.concatenate((neg_rt_sessavg,neg_yc_sessavg,pos_rt_sessavg,pos_yc_sessavg),axis=0)
sign = np.concatenate((np.zeros((16*2)),np.ones((16*2))),axis=0)
group = np.concatenate((np.zeros((16)),np.ones((16)),np.zeros((16)),np.ones((16))),axis=0)
data2b = np.concatenate((data[:,np.newaxis],group[:,np.newaxis],sign[:,np.newaxis]),axis=1)
df = pd.DataFrame(data2b,columns = ['Correct','Group','Sign'])
flatui = ["#DB5461", "#593C8F"]

fig, ax = plt.subplots(figsize=(11,9))
p = sns.swarmplot(data = df,x = "Sign",y="Correct",hue="Group",split=True,palette=(['k', 'k']),size=6)
p.legend_.remove()
g = sns.barplot(data = df,x = "Sign",y="Correct",hue="Group",ci=68,palette =flatui,errwidth=6 )
plt.title('Proportion Correct Changes')
sns.despine()
leg = g.axes.get_legend()
handles,labels = g.axes.get_legend_handles_labels()
handles = [handles[2], handles[3]]
new_labels = ['RT', 'YC']
g.legend(handles,new_labels)
new_title = 'Group'
leg.set_title(new_title)
labels = [item.get_text() for item in g.get_xticklabels()]
labels[0] = "retrieval < target"
labels[1] = "retrieval > target"
g.set_xticklabels(labels)
plt.ylabel('Fraction correct change')
plt.ylim(.45,1.05)
plt.xlabel('Type of error')
ax.set_yticks([0.5,0.75,1])
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(27)
fn = savepath + 'correctchanges_avgsess.pdf'
plt.savefig(fn)
# MAKE HISTOGRAM OF TOTAL EVIDENCE!!

sns.palplot(sns.color_palette(flatui))


nTRs = 15
nTR_total = nTRs*3*nTrials
data = np.concatenate((allSepVec_RT,allSepVec_YC,),axis=0)
AB = np.concatenate((np.zeros((npairs*nTR_total,1)),np.ones((npairs*nTR_total,1))),axis=0)
data2b = np.concatenate((data,AB),axis=1)
df = pd.DataFrame(data2b,columns = ['Evidence','AB'])

fig, ax = sns.plt.subplots(figsize=(15,10))
labels = [item.get_text() for item in ax.get_xticklabels()]
labels[0] = "RT"
labels[1] = "YC"
min=-1.5
max=1.5
binw = .1
bins = np.arange(min,max+binw,binw)
sns.distplot(allSepVec_RT, bins=bins,color=flatui[0], label='RT',norm_hist=False,kde=False,hist_kws={'alpha':0.8})
sns.distplot(allSepVec_YC, bins=bins,color=flatui[1], label='YC',norm_hist=False,kde=False,hist_kws={'alpha':0.8})
sns.plt.plot([0.1, .1],[0,2000], color='k', linestyle='--', linewidth=6)
sns.plt.plot([.2, .2],[0,2000], color='k', linestyle='--', linewidth=6)
sns.plt.title('Distribution of Evidence During MOT')
sns.plt.xlabel('Retrieval evidence bin')
range = np.array([0,.17])
scale = range*len(allSepVec_RT)
sns.plt.ylabel('Fraction of TRs in range')
sns.plt.ylim(scale)
sns.plt.xlim(-1,1)
labels2 = np.arange(0,0.2,0.05)
scaled_labels = len(allSepVec_RT)*labels2
result2 = [str(x) for x in labels2]
sns.plt.yticks( scaled_labels, result2 )
sns.despine()
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(30)
sns.plt.legend()
fn = savepath + 'histevidence.pdf'
plt.savefig(fn)

pair = np.array([])
session = np.array([])
group = np.array([])

data = np.concatenate((np.reshape(nCorrect_RT.T,(16*3,1)),np.reshape(nCorrect_YC.T,(16*3,1))))
for i in np.arange(2):
    thispair = np.concatenate((np.reshape(np.arange(nPairs),(nPairs,1)),np.reshape(np.arange(nPairs),(nPairs,1)),np.reshape(np.arange(nPairs),(nPairs,1))),axis=0)
    pair = np.append(pair,thispair)
    thissession = np.concatenate((np.zeros((nPairs,1)),np.ones((nPairs,1)),2*np.ones((nPairs,1))),axis=0)
    session = np.append(session,thissession)
    thisgroup = i*np.ones((nPairs*3,1))
    group = np.append(group,thisgroup)
pair = np.reshape(pair,(len(pair),1))
group = np.reshape(group,(len(pair),1))
session = np.reshape(session,(len(pair),1))
data2b = np.concatenate((data,pair,group,session),axis=1)
df = pd.DataFrame(data2b,columns = ['Correct','Pair','Group','Run'])

fig, ax = plt.subplots(figsize=(10,7))
#sns.barplot(data=data, palette=sns.color_palette("Set2", 10),ci=68)
#labels = [item.get_text() for item in ax.get_xticklabels()]
#labels[0] = "RT"
#labels[1] = "YC"
sns.stripplot(data = df,x = "Run",y="Correct",hue="Group",split=True,palette = itertools.cycle(sns.color_palette("husl",3)))

sns.boxplot(data = df,x = "Run",y="Correct",hue="Group",palette = itertools.cycle(sns.color_palette("husl",9)))
plt.title('Proportion Correct Changes')
#ax.set_xticklabels(labels)
plt.ylabel('Average Correct Change')
plt.ylim(.5,.8)
plt.xlabel('Run Number')
print(stats.ttest_rel(nCorrect_RT,nCorrect_YC))

# plot both on the same graph, collapse across runs
pair = np.array([])
session = np.array([])
group = np.array([])
data = np.concatenate((np.reshape(nCor_Neg_RT.T,(16*3,1)),np.reshape(nCor_Neg_YC.T,(16*3,1)),np.reshape(nCor_Pos_RT.T,(16*3,1)),np.reshape(nCor_Pos_YC.T,(16*3,1))))
sign = np.concatenate((np.zeros((16*3*2,1)),np.ones((16*3*2,1))),axis=0)
for i in np.arange(2):
    thispair = np.concatenate((np.reshape(np.arange(nPairs),(nPairs,1)),np.reshape(np.arange(nPairs),(nPairs,1)),np.reshape(np.arange(nPairs),(nPairs,1))),axis=0)
    pair = np.append(pair,thispair)
    thisgroup = i*np.ones((nPairs*3,1))
    group = np.append(group,thisgroup)
pair = np.append(pair,pair)
group = np.append(group,group)
pair = np.reshape(pair,(len(pair),1))
group = np.reshape(group,(len(pair),1))
data2b = np.concatenate((data,pair,group,sign),axis=1)
df = pd.DataFrame(data2b,columns = ['Correct','Pair','Group','Sign'])
fig, ax = plt.subplots(figsize=(15,10))
p = sns.swarmplot(data = df,x = "Sign",y="Correct",hue="Group",split=True,palette=(['k', 'k']),size=6)
p.legend_.remove()
g = sns.barplot(data = df,x = "Sign",y="Correct",hue="Group",ci=68,palette =flatui,errwidth=6 )
plt.title('Proportion Correct Changes')
sns.despine()
leg = g.axes.get_legend()
handles,labels = g.axes.get_legend_handles_labels()
handles = [handles[2], handles[3]]
new_labels = ['RT', 'YC']
g.legend(handles,new_labels)
new_title = 'Group'
leg.set_title(new_title)
labels = [item.get_text() for item in g.get_xticklabels()]
labels[0] = "retrieval < target"
labels[1] = "retrieval > target"
g.set_xticklabels(labels)
plt.ylabel('Fraction correct change')
plt.ylim(.45,1.05)
plt.xlabel('Type of error')
ax.set_yticks([0.5,0.75,1])
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(30)
fn = savepath + 'correctchanges.eps'
plt.savefig(fn)

fakeres = np.reshape(np.array([0.6,.8,.7,.7]),(4,1))
group = np.reshape(np.array([0,1,0,1]),(4,1))
type = np.reshape(np.array([0,0,1,1]),(4,1))
data2b = np.concatenate((fakeres,group,type),axis=1)
df = pd.DataFrame(data2b,columns = ['data', 'group', 'type'])
fig, ax = plt.subplots(figsize=(15,10))
g = sns.barplot(data = df,x = "type",y="data",hue="group",palette =flatui )
plt.title('Memory results by group')
sns.despine()
leg = g.axes.get_legend()
handles,labels = g.axes.get_legend_handles_labels()
new_labels = ['RT', 'YC']
g.legend(handles,new_labels)
new_title = 'Group'
leg.set_title(new_title)
labels = [item.get_text() for item in g.get_xticklabels()]
labels[0] = "MOT"
labels[1] = "Omit"
g.set_xticklabels(labels)
plt.ylabel('Memory score')
plt.ylim(.45,1.05)
ax.set_yticks([0.5,0.75,1])
plt.xlabel('Stimulus group')
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(30)
fn = savepath + 'exdata.eps'
plt.savefig(fn)


print(stats.ttest_rel(nCor_Neg_RT,nCor_Neg_YC))

# now separate when evidence is too low
pair = np.array([])
session = np.array([])
group = np.array([])
data = np.concatenate((np.reshape(nCor_Neg_RT.T,(16*3,1)),np.reshape(nCor_Neg_YC.T,(16*3,1))))
for i in np.arange(2):
    thispair = np.concatenate((np.reshape(np.arange(nPairs),(nPairs,1)),np.reshape(np.arange(nPairs),(nPairs,1)),np.reshape(np.arange(nPairs),(nPairs,1))),axis=0)
    pair = np.append(pair,thispair)
    thissession = np.concatenate((np.zeros((nPairs,1)),np.ones((nPairs,1)),2*np.ones((nPairs,1))),axis=0)
    session = np.append(session,thissession)
    thisgroup = i*np.ones((nPairs*3,1))
    group = np.append(group,thisgroup)
pair = np.reshape(pair,(len(pair),1))
group = np.reshape(group,(len(pair),1))
session = np.reshape(session,(len(pair),1))
data2b = np.concatenate((data,pair,group,session),axis=1)
df = pd.DataFrame(data2b,columns = ['Correct','Pair','Group','Run'])

fig, ax = plt.subplots(figsize=(10,7))
#sns.barplot(data=data, palette=sns.color_palette("Set2", 10),ci=68)
#labels = [item.get_text() for item in ax.get_xticklabels()]
#labels[0] = "RT"
#labels[1] = "YC"
sns.stripplot(data = df,x = "Run",y="Correct",hue="Group",split=True,palette = itertools.cycle(sns.color_palette("husl",3)))

sns.boxplot(data = df,x = "Run",y="Correct",hue="Group",palette = itertools.cycle(sns.color_palette("husl",9)))
plt.title('Proportion Correct Changes when too low')
#ax.set_xticklabels(labels)
plt.ylabel('Average Correct Change')
plt.ylim(.45,1)
plt.xlabel('Run Number')
print(stats.ttest_rel(nCor_Neg_RT,nCor_Neg_YC))


# POSITIVE CHANGES ONLY
pair = np.array([])
session = np.array([])
group = np.array([])

data = np.concatenate((np.reshape(nCor_Pos_RT.T,(16*3,1)),np.reshape(nCor_Pos_YC.T,(16*3,1))))
for i in np.arange(2):
    thispair = np.concatenate((np.reshape(np.arange(nPairs),(nPairs,1)),np.reshape(np.arange(nPairs),(nPairs,1)),np.reshape(np.arange(nPairs),(nPairs,1))),axis=0)
    pair = np.append(pair,thispair)
    thissession = np.concatenate((np.zeros((nPairs,1)),np.ones((nPairs,1)),2*np.ones((nPairs,1))),axis=0)
    session = np.append(session,thissession)
    thisgroup = i*np.ones((nPairs*3,1))
    group = np.append(group,thisgroup)
pair = np.reshape(pair,(len(pair),1))
group = np.reshape(group,(len(pair),1))
session = np.reshape(session,(len(pair),1))
data2b = np.concatenate((data,pair,group,session),axis=1)
df = pd.DataFrame(data2b,columns = ['Correct','Pair','Group','Run'])

fig, ax = plt.subplots(figsize=(10,7))
#sns.barplot(data=data, palette=sns.color_palette("Set2", 10),ci=68)
#labels = [item.get_text() for item in ax.get_xticklabels()]
#labels[0] = "RT"
#labels[1] = "YC"
sns.stripplot(data = df,x = "Run",y="Correct",hue="Group",split=True,palette = itertools.cycle(sns.color_palette("husl",3)))

sns.boxplot(data = df,x = "Run",y="Correct",hue="Group",palette = itertools.cycle(sns.color_palette("husl",9)))
plt.title('Proportion Correct Changes when too high')
#ax.set_xticklabels(labels)
plt.ylabel('Average Correct Change')
plt.ylim(.45,1)
plt.xlabel('Run Number')
print(stats.ttest_rel(nCor_Pos_RT,nCor_Pos_YC))

# now look at how the magnitude of changes in evidence compare (FOR INREASES AND DECREASES IN EVIDENCE)
# first for increases in evidence
bins=np.arange(-1.5,1.6,.1)
y = errorScale_Pos_RT[:,1]
y2 = errorScale_Pos_YC[:,1]
fig, ax = plt.subplots(figsize=(10,7))
plt.title('Magnitude of Evidence Change when E>target')
plt.hist(y,alpha=0.2,bins=bins,color='b', weights=np.zeros_like(y) + 1/ y.size, label='RT')
plt.hist(y2,alpha=0.2,bins=bins, color='r',weights=np.zeros_like(y2) + 1/ y2.size, label='YC')
plt.xlim(bins[0],bins[-1])
plt.ylim(0,.2)
plt.legend()
print(stats.ttest_ind(y,y2,equal_var=False))

# first for decreases in evidence
bins=np.arange(-1.5,1.6,.1)
y = errorScale_Neg_RT[:,1]
y2 = errorScale_Neg_YC[:,1]
fig, ax = plt.subplots(figsize=(10,7))
plt.title('Magnitude of Evidence Change when E<target')
plt.hist(y,alpha=0.2,bins=bins,color='b', weights=np.zeros_like(y) + 1/ y.size, label='RT')
plt.hist(y2,alpha=0.2,bins=bins, color='r',weights=np.zeros_like(y2) + 1/ y2.size, label='YC')
plt.xlim(bins[0],bins[-1])
plt.ylim(0,.2)
plt.legend()
print(stats.ttest_ind(y,y2,equal_var=False))

# NOW : try to see how the ratio of good evidence/feedback changes by different time points either behind or ahead
# store for RT/YC: (E, DE at whatever timepoint)
# nposDE, nnegDE
# nposE, nnegE

########################################################################################################################
# NOW COMPARE WITH DIFFERENT SHIFTS IN TIME

target = 0.15
nTRs = 15

test = np.arange(-5, 6, 1)
npoints = len(test)
n_negDE = np.zeros((npairs, npoints, 2))
n_posDE = np.zeros((npairs, npoints, 2))
n_posE = np.zeros((npairs, npoints, 2))
n_negE = np.zeros((npairs, npoints, 2))
n_negDE_posE = np.zeros((npairs, npoints, 2))
n_negDE_negE = np.zeros((npairs, npoints, 2))
n_posDE_posE = np.zeros((npairs, npoints, 2))
n_posDE_negE = np.zeros((npairs, npoints, 2))
corr = np.zeros((npairs, npoints, 2))
# first iterate over all subjects
for s in np.arange(nPairs):
    RT_ev = allSep_RT[:, :, s]
    YC_ev = allSep_YC[:, :, s]
    # then iterate over all possible shifting values
    for ti in np.arange(npoints):
        errorScale_RT = np.empty((0, 2), int)
        errorScale_YC = np.empty((0, 2), int)
        t = test[ti]
        # now look at the data and get the way it changed comparing
        for session in np.arange(3):
            session_ind = np.arange(session * nTRs, (session + 1) * nTRs)
            thisSession_RT = RT_ev[:, session_ind]
            thisSession_YC = YC_ev[:, session_ind]
            for row in np.arange(nTrials):
                if t > 0:
                    scan = np.arange(nTRs - t - 1)
                elif t == 0:
                    scan = np.arange(1, nTRs - 1)
                elif t < 0:
                    scan = np.arange(np.abs(t) + 1, nTRs, 1)
                for col in scan:
                    # RT_thisdiff = thisSession_RT[row,col+t] - thisSession_RT[row,col]
                    RT_thisdiff = (thisSession_RT[row, col + t + 1] - thisSession_RT[row, col + t - 1]) / 2
                    RT_1 = thisSession_RT[row, col]
                    errorScale_RT = np.concatenate(
                        (errorScale_RT, np.reshape(np.array([RT_1 - target, RT_thisdiff]), (1, 2))), axis=0)
                    # YC_thisdiff = thisSession_YC[row,col+t] - thisSession_YC[row,col]
                    YC_thisdiff = (thisSession_YC[row, col + t + 1] - thisSession_YC[row, col + t - 1]) / 2
                    YC_1 = thisSession_YC[row, col]
                    errorScale_YC = np.concatenate(
                        (errorScale_YC, np.reshape(np.array([YC_1 - target, YC_thisdiff]), (1, 2))), axis=0)
        x = errorScale_RT[:, 0]
        y = errorScale_RT[:, 1]
        corr[s, ti, 0] = stats.pearsonr(x, y)[0]
        x = errorScale_YC[:, 0]
        y = errorScale_YC[:, 1]
        corr[s, ti, 1] = stats.pearsonr(x, y)[0]

        negDE = np.argwhere(errorScale_RT[:, 1] < 0)
        n_negDE[s, ti, 0] = len(negDE)
        x1 = errorScale_RT[negDE, 0]
        y1 = errorScale_RT[negDE, 1]

        n_negDE_posE[s, ti, 0] = np.shape(np.argwhere(x1[:, 0] > (target)))[0]
        n_negDE_negE[s, ti, 0] = np.shape(np.argwhere(x1[:, 0] < (target)))[0]

        posDE = np.argwhere(errorScale_RT[:, 1] > 0)
        n_posDE[s, ti, 0] = len(posDE)
        x1 = errorScale_RT[posDE, 0]
        y1 = errorScale_RT[posDE, 1]
        n_posDE_posE[s, ti, 0] = np.shape(np.argwhere(x1[:, 0] > (target)))[0]
        n_posDE_negE[s, ti, 0] = np.shape(np.argwhere(x1[:, 0] < (target)))[0]
        negE = np.argwhere(errorScale_RT[:, 0] < target)
        n_negE[s, ti, 0] = len(negE)
        posE = np.argwhere(errorScale_RT[:, 0] > target)
        n_posE[s, ti, 0] = len(posE)

        negDE = np.argwhere(errorScale_YC[:, 1] < 0)
        n_negDE[s, ti, 1] = len(negDE)
        x1 = errorScale_YC[negDE, 0]
        y1 = errorScale_YC[negDE, 1]
        n_negDE_posE[s, ti, 1] = np.shape(np.argwhere(x1[:, 0] > (target)))[0]
        n_negDE_negE[s, ti, 1] = np.shape(np.argwhere(x1[:, 0] < (target)))[0]

        posDE = np.argwhere(errorScale_YC[:, 1] > 0)
        n_posDE[s, ti, 1] = len(posDE)
        x1 = errorScale_YC[posDE, 0]
        y1 = errorScale_YC[posDE, 1]
        n_posDE_posE[s, ti, 1] = np.shape(np.argwhere(x1[:, 0] > (target)))[0]
        n_posDE_negE[s, ti, 1] = np.shape(np.argwhere(x1[:, 0] < (target)))[0]
        negE = np.argwhere(errorScale_YC[:, 0] < target)
        n_negE[s, ti, 1] = len(negE)
        posE = np.argwhere(errorScale_YC[:, 0] > target)
        n_posE[s, ti, 1] = len(posE)
fig, ax = plt.subplots(figsize=(10,7))
palette = itertools.cycle(sns.color_palette("husl",8))
label = (['RT', 'YC'])
colors = ['b', 'r']
for gr in np.arange(2):
    y = corr[:,:,gr]
    desc = label[gr]
    sns.tsplot(y,condition=desc,color=colors[gr],ci=68, linewidth=3.0)
    for s in np.arange(nPairs):
        plt.plot(y[s,:], '--', color=colors[gr], alpha=0.5)
plt.title('Correlation (E(t) and dE/dt at t + shift)')
plt.xlabel('Shift (TRs)')
plt.ylabel('Correlation')
labels2 = test
result2 = [str(x) for x in labels2]
scaled_labels = np.arange(0,11)
plt.xticks( scaled_labels, result2 )

# now test significance for all time points
for t in np.arange(npoints):
    gr1 = corr[:,t,0]
    gr2 = corr[:,t,1]

    print('for %i, p value is %2.2f' % (test[t],stats.ttest_rel(gr1,gr2).pvalue))

plt.show()
